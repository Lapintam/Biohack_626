"""Live discovery agent — self-contained for the API image.

``run_ask(query)`` takes a free-text question about a drug, a cancer type, or a
patient profile, and either:

- ``_run_llm`` — a real Anthropic tool-use loop (the model picks the tool and
  writes the narration), used when ``ANTHROPIC_API_KEY`` is set; or
- ``_run_scripted`` — a deterministic fallback that string-matches the query,
  fetches/scores in-process, and composes a templated three-gates narration.
  Always works with no API key (or after any LLM/API error, incl. low-credit 400).

All tools call LOCAL functions directly (read report JSON, run the annotate
scorer) — no HTTP. This module imports only ``anthropic`` (optional) + stdlib +
the sibling ``annotate`` module, so it ships in the ``api/`` image without the
``strata`` package.

Returns the transcript dict::

    {"query", "steps":[{"tool","args","result_summary"}], "answer",
     "used_llm": bool, "llm_error"?: str}
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Any

from annotate import annotate_profile, load_annotation_model

REPORTS_DIR = Path(__file__).parent / "reports"
DEMO_DIR = REPORTS_DIR / "demo_profiles"

MODEL = "claude-opus-4-8"
MAX_TOKENS = 2000
MAX_TURNS = 8  # safety bound on the LLM tool-use loop

SYSTEM_PROMPT = (
    "You are a pharmaco-epigenomics discovery agent. You answer natural-language "
    "questions about cross-indication methylation->drug-response opportunities and "
    "patient methylation profiles by calling local tools and reading standardized "
    "reports.\n\n"
    "FRAME — 'the three structural gates'. A human expert cannot remove three "
    "priors; an agent has none of them:\n"
    "  1. TAXONOMY. Experts are trained/funded/regulated by organ. A subgroup "
    "defined by a molecular state ACROSS indications is invisible to them.\n"
    "  2. THE KNOWN-GENE PRIOR. Experts evaluate hypotheses about famous genes; "
    "they will not prioritize a gene they've never heard of (e.g. ID2 for a PARP "
    "inhibitor).\n"
    "  3. MULTIVARIATE CONTROL. Experts reason over one or two variables; they "
    "cannot residualize the whole genome against tissue + latent methylome axes.\n\n"
    "PROTOCOL — follow in order:\n"
    "  1. Call list_available to see which drugs and cancer types have reports.\n"
    "  2. Decide whether the query is about a DRUG, a CANCER TYPE, or a PATIENT "
    "PROFILE. For a drug/cancer type, pick the best slug and call get_drug_report "
    "or get_cancer_type_report. For a patient (words like 'patient', 'annotate', "
    "'profile', 'this sample', or a demo id), call annotate_patient.\n"
    "  3. Narrate the TOP opportunity through the three gates. Name the hero lead "
    "or the top predicted-sensitive drug + its driver genes, give the cross-screen "
    "replication / cross-indication recurrence, and state the honest caveats.\n\n"
    "RULES:\n"
    "  - The evidence is DIRECTIONAL CONCORDANCE across independent screens + "
    "functional silencing, NOT a large effect. Held-out validation rho is modest "
    "(~ -0.25 to -0.37); say so before you are asked.\n"
    "  - These are cell-line hypotheses, research-use only, not clinical claims; "
    "validation is on held-out lines.\n"
    "  - Be specific and quantitative. Do not overclaim. The agent removes the "
    "structural gate; the expert still adjudicates.\n"
    "  - Finish with your narration as your final message text (plain prose)."
)

TOOLS = [
    {
        "name": "list_available",
        "description": (
            "List which drugs and cancer types have standardized reports, plus the "
            "available demo patient profiles and the three structural gates. Call "
            "this first to choose what to fetch."
        ),
        "input_schema": {"type": "object", "properties": {}, "required": []},
    },
    {
        "name": "get_drug_report",
        "description": (
            "Fetch the standardized report for one drug by slug (e.g. 'olaparib'). "
            "Returns target, pathway, the validation funnel, ranked leads, the hero "
            "lead, and honesty rails."
        ),
        "input_schema": {
            "type": "object",
            "properties": {"slug": {"type": "string", "description": "Drug slug, e.g. 'olaparib'."}},
            "required": ["slug"],
        },
    },
    {
        "name": "get_cancer_type_report",
        "description": (
            "Fetch the standardized report for one cancer type by slug (e.g. "
            "'lung'). Returns ranked drug/gene leads, cross-indication connections, "
            "and honesty rails."
        ),
        "input_schema": {
            "type": "object",
            "properties": {"slug": {"type": "string", "description": "Cancer-type slug, e.g. 'lung'."}},
            "required": ["slug"],
        },
    },
    {
        "name": "annotate_patient",
        "description": (
            "Score a demo patient methylation profile against the validated leads. "
            "Returns ranked per-drug predictions (sensitive/resistant/neutral) with "
            "driver genes, plus the sample label and ground-truth AUCs. Use for "
            "patient / profile / annotate queries."
        ),
        "input_schema": {
            "type": "object",
            "properties": {"demo_id": {"type": "string", "description": "Demo profile id, e.g. 'nci-h209'."}},
            "required": ["demo_id"],
        },
    },
]


# ---------------------------------------------------------------------------
# Local data access (no HTTP)
# ---------------------------------------------------------------------------


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text())


def load_manifest() -> dict[str, Any]:
    path = REPORTS_DIR / "manifest.json"
    if not path.exists():
        return {"drugs": [], "cancer_types": [], "gates": []}
    return _read_json(path)


def list_demo_profiles() -> list[dict[str, Any]]:
    path = DEMO_DIR / "index.json"
    if not path.exists():
        return []
    return _read_json(path).get("profiles", [])


def load_demo_profile(demo_id: str) -> dict[str, Any] | None:
    path = DEMO_DIR / f"{demo_id}.json"
    if not path.exists():
        return None
    return _read_json(path)


def list_available() -> dict[str, Any]:
    m = load_manifest()
    return {
        "drugs": [{"name": d["name"], "slug": d["slug"], "target": d.get("target")} for d in m.get("drugs", [])],
        "cancer_types": [{"name": c["name"], "slug": c["slug"]} for c in m.get("cancer_types", [])],
        "demo_profiles": [{"id": p["id"], "label": p.get("label"), "tissue": p.get("tissue")} for p in list_demo_profiles()],
        "gates": m.get("gates", []),
    }


def get_drug_report(slug: str) -> dict[str, Any]:
    path = REPORTS_DIR / "drug" / f"{slug}.json"
    if not path.exists():
        return {"error": f"drug report not found: {slug}"}
    return _read_json(path)


def get_cancer_type_report(slug: str) -> dict[str, Any]:
    path = REPORTS_DIR / "cancer_type" / f"{slug}.json"
    if not path.exists():
        return {"error": f"cancer_type report not found: {slug}"}
    return _read_json(path)


def annotate_patient(demo_id: str) -> dict[str, Any]:
    prof = load_demo_profile(demo_id)
    if prof is None:
        return {"error": f"demo profile not found: {demo_id}"}
    out = annotate_profile(prof["profile"], load_annotation_model())
    out["source_label"] = prof.get("label")
    out["tissue"] = prof.get("tissue")
    out["ground_truth_auc"] = prof.get("ground_truth_auc")
    return out


def _dispatch_tool(name: str, args: dict[str, Any]) -> Any:
    if name == "list_available":
        return list_available()
    if name == "get_drug_report":
        return get_drug_report(str(args["slug"]))
    if name == "get_cancer_type_report":
        return get_cancer_type_report(str(args["slug"]))
    if name == "annotate_patient":
        return annotate_patient(str(args["demo_id"]))
    return {"error": f"unknown tool: {name}"}


def _result_summary(name: str, result: Any) -> str:
    if not isinstance(result, dict):
        return str(result)[:200]
    if "error" in result:
        return f"error: {result['error']}"
    if name == "list_available":
        return (
            f"{len(result.get('drugs', []))} drug reports, "
            f"{len(result.get('cancer_types', []))} cancer-type reports, "
            f"{len(result.get('demo_profiles', []))} demo profiles available."
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
        if leads:
            return (
                f"{result.get('cancer_type')}: {result.get('n_leads')} leads; "
                f"top = {leads[0].get('drug')}->{leads[0].get('gene')} "
                f"({leads[0].get('direction')})"
            )
        return f"{result.get('cancer_type')}: no leads."
    if name == "annotate_patient":
        preds = result.get("predictions") or []
        top = preds[0] if preds else {}
        return (
            f"{result.get('source_label')}: {result.get('n_drugs_scored')} drugs "
            f"scored; top = {top.get('drug')} ({top.get('call')}, "
            f"score={top.get('score')})."
        )
    return json.dumps(result)[:200]


# ---------------------------------------------------------------------------
# Public entrypoint
# ---------------------------------------------------------------------------


def _has_key() -> bool:
    return bool(os.environ.get("ANTHROPIC_API_KEY"))


def run_ask(query: str) -> dict[str, Any]:
    """Run the discovery agent for *query*. Returns the transcript dict."""
    if _has_key():
        try:
            return _run_llm(query)
        except Exception as exc:  # demo safety: never crash on an API/LLM error
            print(
                f"[strata] live agent failed ({type(exc).__name__}: {exc}); "
                f"falling back to scripted mode.",
                file=sys.stderr,
            )
            result = _run_scripted(query)
            result["llm_error"] = f"{type(exc).__name__}: {exc}"
            return result
    return _run_scripted(query)


# ---------------------------------------------------------------------------
# Scripted fallback (no API key required)
# ---------------------------------------------------------------------------

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

_CANCER_WORDS = ("cancer", "tumor", "tumour", "carcinoma", "indication", "tissue", "actionable in")
_PATIENT_WORDS = ("patient", "annotate", "profile", "this sample", "my sample", "methylation profile")


def _fmt_r(x: Any) -> str:
    try:
        return f"{float(x):+.3f}"
    except (TypeError, ValueError):
        return "n/a"


def _match_patient(query: str) -> str | None:
    """Return a demo_id if the query is a patient/profile query (or names a demo)."""
    q = query.lower()
    profiles = list_demo_profiles()
    # Direct id / name match wins.
    for p in profiles:
        pid = p["id"].lower()
        name = (p.get("label") or "").lower().split(" ")[0]
        if pid in q or (name and name in q):
            return p["id"]
    if any(w in q for w in _PATIENT_WORDS) and profiles:
        return profiles[0]["id"]
    return None


def _match_target(query: str, manifest: dict[str, Any]) -> tuple[str, str] | None:
    q = query.lower()
    drugs = manifest.get("drugs", [])
    cancer_types = manifest.get("cancer_types", [])

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
        if len(drug_hits[0]["name"]) >= len(cancer_hits[0]["name"]):
            return ("drug", drug_hits[0]["slug"])
        return ("cancer_type", cancer_hits[0]["slug"])

    drug_slugs = {d["slug"] for d in drugs}
    for phrase, slug in _DRUG_HINTS.items():
        if phrase in q and slug in drug_slugs:
            return ("drug", slug)

    if any(w in q for w in _CANCER_WORDS) and cancer_types:
        return ("cancer_type", cancer_types[0]["slug"])
    if drugs:
        return ("drug", drugs[0]["slug"])
    if cancer_types:
        return ("cancer_type", cancer_types[0]["slug"])
    return None


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
        rep_bits.append(f"PRISM rho={_fmt_r(hero.get('prism_r'))} (p={hero.get('prism_p')})")
    if hero.get("ctrp_r") is not None:
        rep_bits.append(f"CTRP rho={_fmt_r(hero.get('ctrp_r'))} (p={hero.get('ctrp_p')})")
    rep_str = "; ".join(rep_bits) if rep_bits else "single-screen"
    pathway_str = f", {pathway}" if pathway else ""

    lines = [
        f"Query: \"{query}\" -> {drug} (target: {target}{pathway_str}).\n",
        (
            f"Starting from a molecular state rather than a disease, the agent scanned "
            f"{funnel.get('genes_scanned', '?')} gene-promoter methylation features "
            f"against {drug} response across cancer types, residualizing tissue and the "
            f"latent methylome axes. The funnel: "
            f"{funnel.get('response_candidates', '?')} response candidates -> "
            f"{funnel.get('functionally_silenced', '?')} functionally silenced -> "
            f"{funnel.get('externally_replicated', '?')} externally replicated "
            f"(n={funnel.get('lines', '?')} cell lines)."
        ),
        f"\nThe top opportunity is **{drug} -> {gene}** ({direction}; evidence {ev}).",
    ]
    if note:
        lines.append(f"\nGate 2 (the known-gene prior): {note}")
    else:
        lines.append(
            f"\nGate 2 (the known-gene prior): {gene} is not a textbook {target} "
            f"response gene — it is a screen-derived hypothesis a famous-gene prior "
            f"would never surface."
        )
    if rec:
        lines.append(
            f"Gate 1 (the taxonomy prior): the signal is cross-indication — same "
            f"direction in {rec.get('same_direction', '?')} of "
            f"{rec.get('tissues_tested', '?')} cancer types "
            f"({rec.get('significant', '?')} individually significant). An "
            f"indication-gated search never assembles that group."
        )
    lines.append(
        f"Gate 3 (multivariate control): {gene} only surfaces after residualizing "
        f"tissue + the global methylation axis, and it clears a functional bar — its "
        f"promoter methylation tracks its own expression DOWN "
        f"(silencing_r={_fmt_r(sil)}), so the methylation silences the gene rather "
        f"than riding along as a passenger."
    )
    lines.append(
        f"\nReplication: the association holds, same direction, across "
        f"{', '.join(screens) if screens else 'one screen'} ({rep_str}). Held-out "
        f"validation rho is modest by design (~ -0.25 to -0.37) — directional "
        f"concordance, not a large effect."
    )
    if honesty:
        lines.append("\nHonest caveats:")
        lines.extend(f"  - {h}" for h in honesty)
    return "\n".join(lines)


def _narrate_cancer_type(report: dict[str, Any], query: str) -> str:
    ct = report.get("cancer_type", "this cancer type")
    leads = report.get("leads") or []
    connects = report.get("connects_to") or []
    honesty = report.get("honesty") or []

    lines = [
        f"Query: \"{query}\" -> {ct} cancer "
        f"({report.get('n_leads', len(leads))} cross-indication leads).\n"
    ]
    if not leads:
        lines.append(f"No replicated methylation->response leads survived control for {ct}.")
        return "\n".join(lines)

    replicated = [l for l in leads if l.get("replicated")]
    hero = (
        max(replicated, key=lambda l: abs(l.get("tissue_rho", 0)))
        if replicated
        else max(leads, key=lambda l: abs(l.get("tissue_rho", 0)))
    )
    gene, drug = hero.get("gene"), hero.get("drug")
    direction = hero.get("direction")
    rho = _fmt_r(hero.get("tissue_rho"))
    ev = hero.get("evidence_level", "L?")
    also_in = hero.get("also_in") or []

    lines.append(
        f"The agent ignores the organ as a boundary and treats it as a variable to "
        f"adjust away. Within {ct}, the strongest actionable, replicated lead is "
        f"**{drug} -> {gene}** ({direction}; tissue rho={rho}, evidence {ev})."
    )
    lines.append(
        f"\nGate 2 (the known-gene prior): {gene} is surfaced purely on screen "
        f"evidence, with no prestige prior — the long tail an expert would not have "
        f"written down for {drug}."
    )
    if also_in:
        lines.append(
            f"Gate 1 (the taxonomy prior): not a one-tumor finding — the same {gene} "
            f"signal recurs in {', '.join(also_in)}. A {ct}-only search never sees the "
            f"cross-indication group."
        )
    else:
        lines.append(
            f"Gate 1 (the taxonomy prior): the lead is evaluated across the whole "
            f"taxonomy-free space, not just within {ct}."
        )
    lines.append(
        f"Gate 3 (multivariate control): the rho is computed after residualizing "
        f"tissue + latent methylome axes (held-out rho modest by design, "
        f"~ -0.25 to -0.37), so it is not a confound the raw data would have reported."
    )
    if connects:
        lines.append(f"\nThis cancer type connects to: {', '.join(connects[:6])}.")
    if honesty:
        lines.append("\nHonest caveats:")
        lines.extend(f"  - {h}" for h in honesty)
    return "\n".join(lines)


def _narrate_patient(result: dict[str, Any], query: str) -> str:
    label = result.get("source_label", "this sample")
    tissue = result.get("tissue")
    preds = result.get("predictions") or []
    gt = result.get("ground_truth_auc") or {}
    honesty = result.get("honesty") or []

    tissue_str = f" ({tissue})" if tissue else ""
    lines = [
        f"Query: \"{query}\" -> patient methylation profile {label}{tissue_str}.\n",
        (
            f"The agent scored this promoter-methylation profile against the validated "
            f"leads, percentiling each driver gene within the cell-line cohort and "
            f"aggregating per drug. {result.get('n_drugs_scored', 0)} drugs were "
            f"scorable."
        ),
    ]
    sens = [p for p in preds if p.get("call") == "sensitive"]
    top = (sens[0] if sens else (preds[0] if preds else None))
    if top is None:
        lines.append("No drug reached the call threshold for this profile.")
        return "\n".join(lines)

    drv = ", ".join(
        f"{d['gene']} (pct={d['patient_pct']}, contrib={d['contribution']:+})"
        for d in (top.get("drivers") or [])
    )
    lines.append(
        f"\nTop predicted-{top.get('call')} drug: **{top.get('drug')}** "
        f"(score={top.get('score'):+}, {top.get('n_leads_scored')} leads scored). "
        f"Drivers: {drv}."
    )
    lines.append(
        f"\nGate 2 (the known-gene prior): the drivers are screen-derived methylation "
        f"leads, not textbook response genes — the agent reads the long tail an expert "
        f"would not have flagged."
    )
    lines.append(
        f"Gate 1 (the taxonomy prior): the call is made from the molecular state "
        f"itself, not the tissue of origin{(' (' + tissue + ')') if tissue else ''} — a "
        f"cross-indication subgroup an organ-gated workup would miss."
    )
    lines.append(
        f"Gate 3 (multivariate control): each driver is percentiled against the cohort "
        f"distribution and weighted by cross-screen evidence, aggregating many weak "
        f"signals rather than reading one gene."
    )
    if gt:
        gt_str = ", ".join(f"{d}={a}" for d, a in gt.items())
        lines.append(f"\nGround-truth AUCs for this illustrative line: {gt_str} (lower AUC = more sensitive).")
    if honesty:
        lines.append("\nHonest caveats:")
        lines.extend(f"  - {h}" for h in honesty)
    return "\n".join(lines)


def _run_scripted(query: str) -> dict[str, Any]:
    steps: list[dict[str, Any]] = []
    avail = _dispatch_tool("list_available", {})
    steps.append({"tool": "list_available", "args": {},
                  "result_summary": _result_summary("list_available", avail)})

    # Patient route first (most specific intent).
    demo_id = _match_patient(query)
    if demo_id is not None:
        result = _dispatch_tool("annotate_patient", {"demo_id": demo_id})
        steps.append({"tool": "annotate_patient", "args": {"demo_id": demo_id},
                      "result_summary": _result_summary("annotate_patient", result)})
        return {"query": query, "steps": steps,
                "answer": _narrate_patient(result, query), "used_llm": False}

    manifest = {"drugs": avail.get("drugs", []), "cancer_types": avail.get("cancer_types", [])}
    match = _match_target(query, manifest)
    if match is None:
        return {"query": query, "steps": steps,
                "answer": "No reports are available to answer this query.", "used_llm": False}

    kind, slug = match
    if kind == "drug":
        report = _dispatch_tool("get_drug_report", {"slug": slug})
        steps.append({"tool": "get_drug_report", "args": {"slug": slug},
                      "result_summary": _result_summary("get_drug_report", report)})
        answer = _narrate_drug(report, query)
    else:
        report = _dispatch_tool("get_cancer_type_report", {"slug": slug})
        steps.append({"tool": "get_cancer_type_report", "args": {"slug": slug},
                      "result_summary": _result_summary("get_cancer_type_report", report)})
        answer = _narrate_cancer_type(report, query)
    return {"query": query, "steps": steps, "answer": answer, "used_llm": False}


# ---------------------------------------------------------------------------
# LLM path (Anthropic tool-use loop, in-process tools)
# ---------------------------------------------------------------------------


def _run_llm(query: str) -> dict[str, Any]:
    from anthropic import Anthropic

    client = Anthropic()
    steps: list[dict[str, Any]] = []

    system = [{"type": "text", "text": SYSTEM_PROMPT, "cache_control": {"type": "ephemeral"}}]
    messages: list[dict[str, Any]] = [
        {
            "role": "user",
            "content": (
                f"Question: {query}\n\n"
                f"Run the protocol against the local tools and finish with your "
                f"three-gates narration as your final message text."
            ),
        }
    ]

    answer = ""
    for _ in range(MAX_TURNS):
        response = client.messages.create(
            model=MODEL, max_tokens=MAX_TOKENS, system=system, tools=TOOLS, messages=messages
        )
        if response.stop_reason != "tool_use":
            answer = next((b.text for b in response.content if b.type == "text"), "")
            break

        messages.append({"role": "assistant", "content": response.content})
        tool_results: list[dict[str, Any]] = []
        for block in response.content:
            if block.type != "tool_use":
                continue
            args = dict(block.input)
            result = _dispatch_tool(block.name, args)
            steps.append({"tool": block.name, "args": args,
                          "result_summary": _result_summary(block.name, result)})
            tool_results.append({"type": "tool_result", "tool_use_id": block.id,
                                 "content": json.dumps(result)})
        messages.append({"role": "user", "content": tool_results})

    return {"query": query, "steps": steps, "answer": answer, "used_llm": True}


if __name__ == "__main__":
    q = sys.argv[1] if len(sys.argv) > 1 else "what's actionable for olaparib?"
    out = run_ask(q)
    print(f"=== ask [{'llm' if out['used_llm'] else 'scripted'}] — {out['query']} ===")
    if out.get("llm_error"):
        print(f"(llm_error: {out['llm_error']})")
    for s in out["steps"]:
        print(f"[tool] {s['tool']}({s['args']}) -> {s['result_summary']}")
    print("\n--- ANSWER ---")
    print(out["answer"])
