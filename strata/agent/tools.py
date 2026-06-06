"""Claude tool-use schemas and dispatcher for the methylation analysis engine.

``TOOLS`` is a list of Anthropic tool-schema dicts that can be passed directly
to the ``tools=`` parameter of ``client.messages.create``.

``dispatch_tool(name, args)`` routes an incoming tool call to the right engine
function and returns a JSON-serializable dict.

Caching
-------
The methylation matrix (~1055 x 14608) takes ~15 s to load from disk.  The two
loaders are wrapped with ``lru_cache`` so the first call pays the I/O cost;
subsequent calls are instant.  Tests can monkeypatch ``_meth`` / ``_drug`` at
the *module level* because ``dispatch_tool`` and the private helpers always look
them up by name (``import strata.agent.tools as _mod; _mod._meth()``) rather
than capturing them at import time.  The ``lru_cache`` wrapper is replaced along
with the function object when a test does ``monkeypatch.setattr(tools, "_meth",
_fake_meth)``, so the real data files are never touched.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Any

import pandas as pd

# ---------------------------------------------------------------------------
# Cached loaders
# ---------------------------------------------------------------------------


@lru_cache(maxsize=1)
def _meth() -> pd.DataFrame:
    from strata.data import gdsc  # deferred so tests can monkeypatch before load

    return gdsc.load_gdsc_methylation()


@lru_cache(maxsize=1)
def _drug() -> pd.DataFrame:
    from strata.data import gdsc

    return gdsc.load_gdsc_drug_response()


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

import strata.agent.tools as _self  # noqa: E402  (import self for monkeyable lookup)


def _auc_for(drug_name: str) -> pd.Series:
    """Return a per-COSMIC_ID AUC Series for *drug_name* (case-insensitive).

    If multiple rows exist for the same COSMIC_ID, the first is kept (they
    should be identical after GDSC dedup, but this guard keeps the contract).
    Raises ``KeyError`` if the drug is not found.
    """
    drug = _self._drug()
    mask = drug["drug_name"].str.lower() == drug_name.lower()
    sub = drug.loc[mask, ["COSMIC_ID", "auc"]].drop_duplicates("COSMIC_ID").set_index("COSMIC_ID")["auc"]
    if sub.empty:
        raise KeyError(drug_name)
    return sub


# ---------------------------------------------------------------------------
# Anthropic tool schemas
# ---------------------------------------------------------------------------

TOOLS: list[dict[str, Any]] = [
    {
        "name": "list_drugs",
        "description": (
            "List drug names available in the GDSC2 drug-response dataset. "
            "Optionally filter by a substring. Returns up to 40 matches. "
            "Use this tool first to discover valid drug names before calling "
            "rank_markers or evaluate_subgroup."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "contains": {
                    "type": "string",
                    "description": (
                        "Case-insensitive substring to filter drug names. "
                        "Omit or pass an empty string to list the first 40 drugs."
                    ),
                }
            },
            "required": [],
        },
    },
    {
        "name": "rank_markers",
        "description": (
            "Rank methylation marker genes by their Spearman correlation with "
            "drug AUC (lower AUC = more sensitive). Returns the top_n genes with "
            "rho and q-value. A negative rho means high methylation marks "
            "sensitivity; a positive rho means high methylation marks resistance. "
            "Only genes passing FDR q < 0.1 are returned."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "drug_name": {
                    "type": "string",
                    "description": "Exact drug name as returned by list_drugs.",
                },
                "top_n": {
                    "type": "integer",
                    "description": "Number of top markers to return (default 15).",
                    "default": 15,
                },
            },
            "required": ["drug_name"],
        },
    },
    {
        "name": "evaluate_subgroup",
        "description": (
            "Median-split cell lines on a marker gene's methylation level to "
            "define a high-methylation vs low-methylation subgroup, then test "
            "whether the two groups differ in drug AUC (Kruskal-Wallis). "
            "Returns grade (++ / + / neutral / - / --), p-value, effect size, "
            "responder cluster label, and sample count."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "drug_name": {
                    "type": "string",
                    "description": "Exact drug name as returned by list_drugs.",
                },
                "marker": {
                    "type": "string",
                    "description": "Gene symbol to use as the methylation marker (e.g. 'BRCA1').",
                },
            },
            "required": ["drug_name", "marker"],
        },
    },
    {
        "name": "write_opportunity_brief",
        "description": (
            "Record a commercial-development opportunity brief for a drug+marker "
            "pairing. The agent composes the rationale text; this tool echoes it "
            "back as a confirmation. Call this after evaluate_subgroup confirms a "
            "meaningful responder subgroup."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "drug_name": {
                    "type": "string",
                    "description": "Drug name for the opportunity.",
                },
                "marker": {
                    "type": "string",
                    "description": "Methylation marker gene that defines the subgroup.",
                },
                "rationale": {
                    "type": "string",
                    "description": (
                        "1-3 sentence explanation of why this drug+marker pairing "
                        "represents a commercial opportunity."
                    ),
                },
            },
            "required": ["drug_name", "marker", "rationale"],
        },
    },
]

# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------


def dispatch_tool(name: str, args: dict[str, Any]) -> dict[str, Any]:
    """Route a tool call to the appropriate engine function.

    Returns a JSON-serializable dict.  Never raises; errors are returned as
    ``{"error": "..."}``.
    """
    try:
        if name == "list_drugs":
            return _dispatch_list_drugs(args)
        if name == "rank_markers":
            return _dispatch_rank_markers(args)
        if name == "evaluate_subgroup":
            return _dispatch_evaluate_subgroup(args)
        if name == "write_opportunity_brief":
            return _dispatch_write_opportunity_brief(args)
        return {"error": f"unknown tool: {name}"}
    except Exception as exc:  # pragma: no cover – belt-and-suspenders
        return {"error": str(exc)}


# ---------------------------------------------------------------------------
# Per-tool dispatch helpers
# ---------------------------------------------------------------------------


def _dispatch_list_drugs(args: dict[str, Any]) -> dict[str, Any]:
    drug = _self._drug()
    contains: str = args.get("contains") or ""
    names = drug["drug_name"].unique().tolist()
    names = sorted(set(n for n in names if contains.lower() in n.lower()))
    return {"drugs": names[:40]}


def _dispatch_rank_markers(args: dict[str, Any]) -> dict[str, Any]:
    from strata.engine.associate import rank_response_markers

    drug_name: str = args.get("drug_name", "")
    top_n: int = int(args.get("top_n", 15))

    try:
        auc = _auc_for(drug_name)
    except KeyError:
        return {"error": f"drug not found: {drug_name!r}"}

    meth = _self._meth()
    markers_df = rank_response_markers(meth, auc)

    rows = []
    for _, row in markers_df.head(top_n).iterrows():
        rows.append(
            {
                "gene": str(row["gene"]),
                "rho": float(row["rho"]),
                "qvalue": float(row["qvalue"]),
            }
        )
    n_significant = int(len(markers_df))
    n_cell_lines = int(markers_df["n"].iloc[0]) if len(markers_df) else 0
    return {
        "drug": drug_name,
        "markers": rows,
        "n_significant": n_significant,
        "n": n_cell_lines,
    }


def _dispatch_evaluate_subgroup(args: dict[str, Any]) -> dict[str, Any]:
    from strata.engine.associate import discover_responder_subgroup

    drug_name: str = args.get("drug_name", "")
    marker: str = args.get("marker", "")

    try:
        auc = _auc_for(drug_name)
    except KeyError:
        return {"error": f"drug not found: {drug_name!r}"}

    meth = _self._meth()
    if marker not in meth.columns:
        return {"error": f"marker not found in methylation matrix: {marker!r}"}

    clusters, assoc = discover_responder_subgroup(meth, auc, marker)
    n = int(clusters.notna().sum())

    return {
        "drug": drug_name,
        "marker": marker,
        "grade": assoc.grade,
        "pvalue": float(assoc.pvalue),
        "effect_size": float(assoc.effect_size),
        "responder_label": int(assoc.responder_label),
        "n": n,
    }


def _dispatch_write_opportunity_brief(args: dict[str, Any]) -> dict[str, Any]:
    drug_name: str = args.get("drug_name", "")
    marker: str = args.get("marker", "")
    # rationale is recorded by the agent; we simply echo confirmation
    return {"ok": True, "drug": drug_name, "marker": marker}
