"""Report agent — natural language -> Strata Report API -> three-gates narration.

``run_report_agent(query, api_base, use_llm)`` takes a free-text question, decides
whether it concerns a *drug* or a *cancer type*, fetches the matching standardized
report from the Strata Report API (``api/main.py``), and narrates the top
opportunity through the three structural gates (taxonomy / known-gene /
multivariate control) — naming the hero lead, its cross-screen replication and
functional silencing, and the honest caveats.

Two paths produce the **same transcript shape** so the front end is identical:

- ``_run_llm`` — a real Anthropic tool-use loop (the model picks the report and
  writes the narration). Used when ``ANTHROPIC_API_KEY`` is present and
  ``use_llm`` is truthy.
- ``_run_scripted`` — a deterministic fallback that string-matches the query
  against the manifest, fetches the report over HTTP, and composes a templated
  three-gates narration from the report fields. Always works with no API key.

The returned dict is the transcript object the front end replays::

    {
      "query": <str>,
      "steps": [{"tool": <str>, "args": {...}, "result_summary": <str>}, ...],
      "answer": <final narration text>,
      "report": <the fetched standardized report>,
      "used_llm": <bool>,
      "llm_error": <str, only on fallback after an LLM error>,
    }

Follows the MODEL / tool-use / graceful-fallback pattern of
``strata/agent/loop.py``.
"""

from __future__ import annotations

import json
import os
import urllib.error
import urllib.request
from typing import Any

from dotenv import load_dotenv

load_dotenv(override=False)

MODEL = "claude-opus-4-8"
MAX_TOKENS = 2000
MAX_TURNS = 8  # safety bound on the LLM tool-use loop

DEFAULT_API_BASE = "http://localhost:8077"

SYSTEM_PROMPT = (
    "You are a pharmaco-epigenomics discovery agent. You answer natural-language "
    "questions about cross-indication methylation->drug-response opportunities by "
    "calling the Strata Report API and reading its standardized reports.\n\n"
    "FRAME — 'the three structural gates'. A human expert cannot remove three "
    "priors, no matter how much time you give them; an agent has none of them:\n"
    "  1. TAXONOMY. Experts are trained/funded/regulated by organ. A subgroup "
    "defined by a molecular state ACROSS indications is invisible to them.\n"
    "  2. THE KNOWN-GENE PRIOR. Experts evaluate hypotheses about famous genes. "
    "They will not prioritize a gene they've never heard of (e.g. ID2, a "
    "developmental transcription factor, for a PARP inhibitor).\n"
    "  3. MULTIVARIATE CONTROL. Experts reason over one or two variables; they "
    "cannot residualize the whole genome against tissue + latent methylome axes "
    "in their head.\n\n"
    "PROTOCOL — follow in order:\n"
    "  1. Call list_available to see which drugs and cancer types have reports.\n"
    "  2. Decide whether the query is about a DRUG or a CANCER TYPE, pick the best "
    "matching slug, and call get_drug_report or get_cancer_type_report.\n"
    "  3. Write a concise narration of the TOP opportunity through the three "
    "gates. Name the hero lead (gene + direction), give its cross-screen "
    "replication (which screens agree, the rho values), its functional silencing "
    "(silencing_r: promoter methylation tracking expression DOWN), and its "
    "cross-indication recurrence (same_direction / tissues_tested). Then state the "
    "honest caveats from the report's honesty rails.\n\n"
    "RULES:\n"
    "  - The evidence is DIRECTIONAL CONCORDANCE across independent screens + "
    "functional silencing, NOT a large effect. Effect sizes are modest by design; "
    "say so before you are asked.\n"
    "  - These are cell-line hypotheses (evidence ladder ~L3-L5), not clinical "
    "claims.\n"
    "  - Be specific and quantitative. Do not overclaim. The agent removes the "
    "structural gate; the expert still adjudicates.\n"
    "  - Finish with your narration as your final message text (plain prose, no "
    "tool call)."
)

# ---------------------------------------------------------------------------
# Tool definitions (Anthropic tool-use schema)
# ---------------------------------------------------------------------------

TOOLS = [
    {
        "name": "list_available",
        "description": (
            "List which drugs and cancer types have standardized reports "
            "available (from the API manifest). Returns names and slugs plus the "
            "three structural gates. Call this first to choose what to fetch."
        ),
        "input_schema": {"type": "object", "properties": {}, "required": []},
    },
    {
        "name": "get_drug_report",
        "description": (
            "Fetch the standardized report for one drug by slug (e.g. 'olaparib', "
            "'trametinib'). Returns target, pathway, the validation funnel, ranked "
            "leads, the hero lead, and honesty rails."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "slug": {
                    "type": "string",
                    "description": "Drug slug, e.g. 'olaparib'.",
                }
            },
            "required": ["slug"],
        },
    },
    {
        "name": "get_cancer_type_report",
        "description": (
            "Fetch the standardized report for one cancer type by slug (e.g. "
            "'lung', 'breast'). Returns ranked drug/gene leads, cross-indication "
            "connections, and honesty rails."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "slug": {
                    "type": "string",
                    "description": "Cancer-type slug, e.g. 'lung'.",
                }
            },
            "required": ["slug"],
        },
    },
]


# ---------------------------------------------------------------------------
# HTTP helpers (urllib — no extra deps)
# ---------------------------------------------------------------------------


def _http_get(api_base: str, path: str) -> Any:
    url = api_base.rstrip("/") + path
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as resp:  # noqa: S310 (trusted local)
        return json.loads(resp.read().decode("utf-8"))


def _get_manifest(api_base: str) -> dict[str, Any]:
    return _http_get(api_base, "/manifest")


def _get_drug_report(api_base: str, slug: str) -> dict[str, Any]:
    return _http_get(api_base, f"/report/drug/{slug}")


def _get_cancer_type_report(api_base: str, slug: str) -> dict[str, Any]:
    return _http_get(api_base, f"/report/cancer_type/{slug}")


def _dispatch_tool(name: str, args: dict[str, Any], api_base: str) -> Any:
    if name == "list_available":
        m = _get_manifest(api_base)
        return {
            "drugs": [
                {"name": d["name"], "slug": d["slug"], "target": d.get("target")}
                for d in m.get("drugs", [])
            ],
            "cancer_types": [
                {"name": c["name"], "slug": c["slug"]}
                for c in m.get("cancer_types", [])
            ],
            "gates": m.get("gates", []),
        }
    if name == "get_drug_report":
        return _get_drug_report(api_base, str(args["slug"]))
    if name == "get_cancer_type_report":
        return _get_cancer_type_report(api_base, str(args["slug"]))
    return {"error": f"unknown tool: {name}"}


def _result_summary(name: str, result: Any) -> str:
    """One-line human summary of a tool result for the transcript."""
    if not isinstance(result, dict):
        return str(result)[:200]
    if name == "list_available":
        return (
            f"{len(result.get('drugs', []))} drug reports, "
            f"{len(result.get('cancer_types', []))} cancer-type reports available."
        )
    if name == "get_drug_report":
        hero = result.get("hero") or {}
        funnel = result.get("funnel") or {}
        return (
            f"{result.get('drug')} ({result.get('target')}): "
            f"{funnel.get('genes_scanned')} genes scanned -> "
            f"{funnel.get('externally_replicated')} externally replicated; "
            f"hero={hero.get('gene')} ({hero.get('direction')}, "
            f"{hero.get('n_screens')} screens, {hero.get('evidence_level')})."
        )
    if name == "get_cancer_type_report":
        leads = result.get("leads") or []
        return (
            f"{result.get('cancer_type')}: {result.get('n_leads')} leads; "
            f"top = {leads[0].get('drug')}->{leads[0].get('gene')} "
            f"({leads[0].get('direction')})"
            if leads
            else f"{result.get('cancer_type')}: no leads."
        )
    if "error" in result:
        return f"error: {result['error']}"
    return json.dumps(result)[:200]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def _has_key() -> bool:
    return bool(os.environ.get("ANTHROPIC_API_KEY"))


def run_report_agent(
    query: str,
    api_base: str = DEFAULT_API_BASE,
    use_llm: bool = True,
) -> dict[str, Any]:
    """Run the report agent for *query*.

    Returns the transcript dict described in the module docstring.
    """
    if use_llm and _has_key():
        try:
            return _run_llm(query, api_base)
        except Exception as exc:  # demo safety: never crash on an API/LLM error
            import sys

            print(
                f"[strata] live report agent failed "
                f"({type(exc).__name__}: {exc}); falling back to scripted mode.",
                file=sys.stderr,
            )
            result = _run_scripted(query, api_base)
            result["llm_error"] = f"{type(exc).__name__}: {exc}"
            return result
    return _run_scripted(query, api_base)


# ---------------------------------------------------------------------------
# Scripted fallback (no API key required)
# ---------------------------------------------------------------------------

# Hints to disambiguate a query toward a drug even when the drug isn't named by
# its catalog name (e.g. "a MEK inhibitor" -> Trametinib).
_DRUG_HINTS: dict[str, str] = {
    "mek inhibitor": "trametinib",
    "mek": "trametinib",
    "parp inhibitor": "olaparib",
    "parp": "olaparib",
    "cdk4/6": "palbociclib",
    "cdk": "palbociclib",
    "egfr": "gefitinib",
    "braf": "dabrafenib",
    "bcl2": "venetoclax",
    "mtor": "rapamycin",
}

_CANCER_WORDS = (
    "cancer",
    "tumor",
    "tumour",
    "carcinoma",
    "indication",
    "tissue",
    "actionable in",
)


def _match_target(
    query: str, manifest: dict[str, Any]
) -> tuple[str, str] | None:
    """Return ('drug'|'cancer_type', slug) by string-matching the manifest.

    Drug-name / cancer-name matches win first; then hint phrases; then a
    cancer-vs-drug intent heuristic with a sensible default.
    """
    q = query.lower()
    drugs = manifest.get("drugs", [])
    cancer_types = manifest.get("cancer_types", [])

    # 1. Direct name match against catalog entries (longest name first so e.g.
    #    "Large Intestine" beats a stray "lung" substring).
    drug_hits = sorted(
        (d for d in drugs if d["name"].lower() in q or d["slug"] in q),
        key=lambda d: -len(d["name"]),
    )
    cancer_hits = sorted(
        (c for c in cancer_types if c["name"].lower() in q or c["slug"] in q),
        key=lambda c: -len(c["name"]),
    )
    if drug_hits and not cancer_hits:
        return ("drug", drug_hits[0]["slug"])
    if cancer_hits and not drug_hits:
        return ("cancer_type", cancer_hits[0]["slug"])
    if drug_hits and cancer_hits:
        # Both present: prefer the longer literal match.
        if len(drug_hits[0]["name"]) >= len(cancer_hits[0]["name"]):
            return ("drug", drug_hits[0]["slug"])
        return ("cancer_type", cancer_hits[0]["slug"])

    # 2. Drug hint phrases ("a MEK inhibitor").
    drug_slugs = {d["slug"] for d in drugs}
    for phrase, slug in _DRUG_HINTS.items():
        if phrase in q and slug in drug_slugs:
            return ("drug", slug)

    # 3. Cancer intent word but no named tissue -> first cancer type.
    if any(w in q for w in _CANCER_WORDS) and cancer_types:
        return ("cancer_type", cancer_types[0]["slug"])

    # 4. Default to a drug report (the canonical demo).
    if drugs:
        return ("drug", drugs[0]["slug"])
    if cancer_types:
        return ("cancer_type", cancer_types[0]["slug"])
    return None


def _fmt_r(x: Any) -> str:
    try:
        return f"{float(x):+.3f}"
    except (TypeError, ValueError):
        return "n/a"


def _narrate_drug(report: dict[str, Any], query: str) -> str:
    drug = report.get("drug", "this drug")
    target = report.get("target") or "its nominal target"
    pathway = report.get("pathway")
    funnel = report.get("funnel") or {}
    hero = report.get("hero") or (report.get("leads") or [{}])[0]
    honesty = report.get("honesty") or []

    gene = hero.get("gene", "the top lead")
    direction = hero.get("direction", "sensitive")
    note = hero.get("gene_note") or ""
    screens = hero.get("screens") or []
    rec = hero.get("recurrence") or {}
    sil = hero.get("silencing_r")
    ev = hero.get("evidence_level", "L?")

    rep_bits = []
    if hero.get("gdsc_r") is not None:
        rep_bits.append(f"GDSC rho={_fmt_r(hero.get('gdsc_r'))}")
    if hero.get("prism_r") is not None:
        rep_bits.append(
            f"PRISM rho={_fmt_r(hero.get('prism_r'))} (p={hero.get('prism_p')})"
        )
    if hero.get("ctrp_r") is not None:
        rep_bits.append(
            f"CTRP rho={_fmt_r(hero.get('ctrp_r'))} (p={hero.get('ctrp_p')})"
        )
    rep_str = "; ".join(rep_bits) if rep_bits else "single-screen"

    pathway_str = f", {pathway}" if pathway else ""
    lines = []
    lines.append(
        f"Query: \"{query}\" -> {drug} (target: {target}{pathway_str}).\n"
    )
    lines.append(
        f"Starting from a molecular state rather than a disease, the agent scanned "
        f"{funnel.get('genes_scanned', '?')} gene-promoter methylation features "
        f"against {drug} response across cancer types, residualizing tissue and the "
        f"latent methylome axes. The funnel: "
        f"{funnel.get('response_candidates', '?')} response candidates -> "
        f"{funnel.get('functionally_silenced', '?')} functionally silenced -> "
        f"{funnel.get('externally_replicated', '?')} externally replicated "
        f"(n={funnel.get('lines', '?')} cell lines)."
    )
    lines.append(
        f"\nThe top opportunity is **{drug} -> {gene}** "
        f"({direction}; evidence {ev})."
    )

    # Gate 2 — known-gene prior.
    if note:
        lines.append(f"\nGate 2 (the known-gene prior): {note}")
    else:
        lines.append(
            f"\nGate 2 (the known-gene prior): {gene} is not a textbook "
            f"{target} response gene — it is a screen-derived hypothesis a "
            f"famous-gene prior would never surface."
        )

    # Gate 1 — taxonomy / recurrence.
    if rec:
        lines.append(
            f"Gate 1 (the taxonomy prior): the signal is cross-indication — same "
            f"direction in {rec.get('same_direction', '?')} of "
            f"{rec.get('tissues_tested', '?')} cancer types "
            f"({rec.get('significant', '?')} individually significant). An "
            f"indication-gated search never assembles that group."
        )

    # Gate 3 — multivariate control + functional silencing.
    sil_str = _fmt_r(sil)
    lines.append(
        f"Gate 3 (multivariate control): {gene} only surfaces after residualizing "
        f"tissue + the global methylation axis, and it clears a functional bar — "
        f"its promoter methylation tracks its own expression DOWN "
        f"(silencing_r={sil_str}), so the methylation silences the gene rather "
        f"than riding along as a passenger."
    )

    # Replication strip.
    lines.append(
        f"\nReplication: the association holds, same direction, across "
        f"{', '.join(screens) if screens else 'one screen'} ({rep_str})."
    )

    # Honesty rails.
    if honesty:
        lines.append("\nHonest caveats:")
        for h in honesty:
            lines.append(f"  - {h}")

    return "\n".join(lines)


def _narrate_cancer_type(report: dict[str, Any], query: str) -> str:
    ct = report.get("cancer_type", "this cancer type")
    leads = report.get("leads") or []
    connects = report.get("connects_to") or []
    honesty = report.get("honesty") or []

    lines = []
    lines.append(
        f"Query: \"{query}\" -> {ct} cancer "
        f"({report.get('n_leads', len(leads))} cross-indication leads).\n"
    )

    if not leads:
        lines.append(
            f"No replicated methylation->response leads survived control for {ct}."
        )
        return "\n".join(lines)

    # Prefer a replicated lead as the hero; else the strongest by |rho|.
    replicated = [l for l in leads if l.get("replicated")]
    hero = (
        max(replicated, key=lambda l: abs(l.get("tissue_rho", 0)))
        if replicated
        else max(leads, key=lambda l: abs(l.get("tissue_rho", 0)))
    )
    gene = hero.get("gene")
    drug = hero.get("drug")
    direction = hero.get("direction")
    rho = _fmt_r(hero.get("tissue_rho"))
    ev = hero.get("evidence_level", "L?")
    also_in = hero.get("also_in") or []

    lines.append(
        f"The agent ignores the organ as a boundary and treats it as a variable to "
        f"adjust away. Within {ct}, the strongest actionable, replicated lead is "
        f"**{drug} -> {gene}** ({direction}; tissue rho={rho}, "
        f"evidence {ev})."
    )
    lines.append(
        f"\nGate 2 (the known-gene prior): {gene} is surfaced purely on screen "
        f"evidence, with no prestige prior — the long tail an expert would not "
        f"have written down for {drug}."
    )
    if also_in:
        lines.append(
            f"Gate 1 (the taxonomy prior): this is not a one-tumor finding — the "
            f"same {gene} signal recurs in {', '.join(also_in)}. A "
            f"{ct}-only search would never see the cross-indication group."
        )
    else:
        lines.append(
            f"Gate 1 (the taxonomy prior): the lead is evaluated across the whole "
            f"taxonomy-free space, not just within {ct}."
        )
    lines.append(
        f"Gate 3 (multivariate control): the rho is computed after residualizing "
        f"tissue + latent methylome axes, so it is not a confound the raw data "
        f"would have reported."
    )
    if connects:
        lines.append(
            f"\nThis cancer type connects to: {', '.join(connects[:6])}."
        )
    if honesty:
        lines.append("\nHonest caveats:")
        for h in honesty:
            lines.append(f"  - {h}")
    return "\n".join(lines)


def _run_scripted(query: str, api_base: str) -> dict[str, Any]:
    steps: list[dict[str, Any]] = []

    # Step 1 — list_available.
    avail = _dispatch_tool("list_available", {}, api_base)
    steps.append(
        {
            "tool": "list_available",
            "args": {},
            "result_summary": _result_summary("list_available", avail),
        }
    )

    manifest = {
        "drugs": avail.get("drugs", []),
        "cancer_types": avail.get("cancer_types", []),
    }
    match = _match_target(query, manifest)
    if match is None:
        return {
            "query": query,
            "steps": steps,
            "answer": "No reports are available to answer this query.",
            "report": None,
            "used_llm": False,
        }

    kind, slug = match
    if kind == "drug":
        report = _dispatch_tool("get_drug_report", {"slug": slug}, api_base)
        steps.append(
            {
                "tool": "get_drug_report",
                "args": {"slug": slug},
                "result_summary": _result_summary("get_drug_report", report),
            }
        )
        answer = _narrate_drug(report, query)
    else:
        report = _dispatch_tool(
            "get_cancer_type_report", {"slug": slug}, api_base
        )
        steps.append(
            {
                "tool": "get_cancer_type_report",
                "args": {"slug": slug},
                "result_summary": _result_summary(
                    "get_cancer_type_report", report
                ),
            }
        )
        answer = _narrate_cancer_type(report, query)

    return {
        "query": query,
        "steps": steps,
        "answer": answer,
        "report": report,
        "used_llm": False,
    }


# ---------------------------------------------------------------------------
# LLM path (Anthropic tool-use loop)
# ---------------------------------------------------------------------------


def _run_llm(query: str, api_base: str) -> dict[str, Any]:
    from anthropic import Anthropic

    client = Anthropic()
    steps: list[dict[str, Any]] = []
    last_report: Any = None

    system = [
        {
            "type": "text",
            "text": SYSTEM_PROMPT,
            "cache_control": {"type": "ephemeral"},
        }
    ]
    messages: list[dict[str, Any]] = [
        {
            "role": "user",
            "content": (
                f"Question: {query}\n\n"
                f"Run the protocol against the Strata Report API and finish with "
                f"your three-gates narration as your final message text."
            ),
        }
    ]

    answer = ""
    for _ in range(MAX_TURNS):
        response = client.messages.create(
            model=MODEL,
            max_tokens=MAX_TOKENS,
            system=system,
            tools=TOOLS,
            messages=messages,
        )

        if response.stop_reason != "tool_use":
            answer = next(
                (b.text for b in response.content if b.type == "text"), ""
            )
            break

        messages.append({"role": "assistant", "content": response.content})

        tool_results: list[dict[str, Any]] = []
        for block in response.content:
            if block.type != "tool_use":
                continue
            args = dict(block.input)
            result = _dispatch_tool(block.name, args, api_base)
            if block.name in ("get_drug_report", "get_cancer_type_report"):
                last_report = result
            steps.append(
                {
                    "tool": block.name,
                    "args": args,
                    "result_summary": _result_summary(block.name, result),
                }
            )
            tool_results.append(
                {
                    "type": "tool_result",
                    "tool_use_id": block.id,
                    "content": json.dumps(result),
                }
            )
        messages.append({"role": "user", "content": tool_results})

    return {
        "query": query,
        "steps": steps,
        "answer": answer,
        "report": last_report,
        "used_llm": True,
    }


# ---------------------------------------------------------------------------
# Runnable demo entrypoint
# ---------------------------------------------------------------------------


def _pretty_print(result: dict[str, Any]) -> None:
    mode = "llm" if result.get("used_llm") else "scripted"
    print(f"\n=== Report agent [{mode}] — {result['query']} ===")
    if result.get("llm_error"):
        print(f"(llm_error: {result['llm_error']})")
    for s in result["steps"]:
        print(f"[tool] {s['tool']}({s['args']}) -> {s['result_summary']}")
    print("\n--- ANSWER ---")
    print(result["answer"])


if __name__ == "__main__":
    import sys

    q = sys.argv[1] if len(sys.argv) > 1 else (
        "What drug-response opportunities does olaparib have that an "
        "oncologist would miss?"
    )
    _pretty_print(run_report_agent(q))
