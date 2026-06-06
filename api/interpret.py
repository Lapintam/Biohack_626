"""LLM-mediated "interpret this result" — applies the fixed interpretation key
to an ACTUAL result already shown on screen.

``run_interpret(kind, ref)`` fetches the real result via the in-process
accessors in ``agent.py`` (drug report, cancer-type report, or the patient
scorer), then either:

- ``_interpret_llm`` — one Anthropic call (same SDK/MODEL as ``agent.run_ask``)
  that APPLIES the key to the real numbers; used when ``ANTHROPIC_API_KEY`` is
  set; or
- ``_interpret_scripted`` — a deterministic templated interpretation built from
  the result fields applying the key. ALWAYS works with no API key (or after any
  LLM/API error).

Self-contained for the API image: imports only ``anthropic`` (optional) +
stdlib + the sibling ``agent`` module (which itself only needs ``annotate``).

Returns::

    {"kind", "ref", "interpretation": str, "used_llm": bool, "llm_error"?: str}
"""

from __future__ import annotations

import os
import sys
from typing import Any

import agent

MODEL = agent.MODEL
MAX_TOKENS = 1200

# The fixed interpretation key — embedded verbatim in the LLM system prompt.
INTERPRETATION_KEY = (
    "THE INTERPRETATION KEY (fixed; apply it to the numbers):\n"
    "- rho is corr(promoter methylation, drug AUC). Lower AUC = more sensitive. "
    "So rho NEGATIVE = hypermethylation marks SENSITIVITY (response); rho "
    "POSITIVE = marks RESISTANCE.\n"
    "- The `direction` label (\"sensitive\"/\"resistant\") is authoritative — "
    "trust it over re-deriving the sign.\n"
    "- Patient annotation `score`: POSITIVE = predicted SENSITIVE (sign-flipped "
    "vs rho).\n"
    "- Evidence level: L3 = GDSC only; L4 = GDSC + 1 external screen; L5 = all 3 "
    "(GDSC+PRISM+CTRP). Every lead is already response-associated AND "
    "functionally silenced (promoter methylation lowers its own expression).\n"
    "- Effects are modest; everything is research-use hypotheses (not clinical), "
    "cell-line-derived; in-cohort cell-line annotations carry training leakage, a "
    "real patient sample does not."
)

SYSTEM_PROMPT = (
    "You are a pharmaco-epigenomics result interpreter. You are given the fixed "
    "interpretation key below, then an ACTUAL result (JSON) already shown to the "
    "user. Apply the key to THESE numbers.\n\n"
    + INTERPRETATION_KEY
)


# ---------------------------------------------------------------------------
# Fetch the actual result (in-process, no HTTP)
# ---------------------------------------------------------------------------


def _fetch_result(kind: str, ref: str) -> dict[str, Any]:
    if kind == "drug":
        return agent.get_drug_report(ref)
    if kind == "cancer_type":
        return agent.get_cancer_type_report(ref)
    if kind == "annotation":
        return agent.annotate_patient(ref)
    return {"error": f"unknown kind: {kind!r} (expected drug|cancer_type|annotation)"}


def _trim_result(kind: str, result: dict[str, Any], top: int = 6) -> dict[str, Any]:
    """Trim to the top ~6 leads/predictions to keep the prompt tight."""
    r = dict(result)
    if kind == "drug":
        if isinstance(r.get("leads"), list):
            r["leads"] = r["leads"][:top]
    elif kind == "cancer_type":
        if isinstance(r.get("leads"), list):
            r["leads"] = r["leads"][:top]
    elif kind == "annotation":
        if isinstance(r.get("predictions"), list):
            r["predictions"] = r["predictions"][:top]
    return r


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------


def _fmt_r(x: Any) -> str:
    try:
        return f"{float(x):+.3f}"
    except (TypeError, ValueError):
        return "n/a"


_LEVEL_MEANS = {
    "L3": "GDSC only (single screen)",
    "L4": "GDSC + one external screen",
    "L5": "replicated across all three screens (GDSC + PRISM + CTRP)",
}


def _level_means(level: str | None) -> str:
    return _LEVEL_MEANS.get(str(level or "").upper(), "an unspecified evidence tier")


def _dir_means(direction: str | None) -> str:
    d = str(direction or "").lower()
    if d == "sensitive":
        return "hypermethylation marks SENSITIVITY (predicted response)"
    if d == "resistant":
        return "hypermethylation marks RESISTANCE"
    return "an unspecified direction"


def _rho_means(rho: Any) -> str:
    try:
        v = float(rho)
    except (TypeError, ValueError):
        return ""
    if v < 0:
        return "negative -> hypermethylation marks sensitivity"
    if v > 0:
        return "positive -> hypermethylation marks resistance"
    return "near zero"


# ---------------------------------------------------------------------------
# Scripted fallback (no API key required) — applies the key deterministically
# ---------------------------------------------------------------------------


def _interpret_scripted(kind: str, ref: str, result: dict[str, Any]) -> str:
    if "error" in result:
        return f"Could not interpret: {result['error']}"
    if kind == "drug":
        return _scripted_drug(result)
    if kind == "cancer_type":
        return _scripted_cancer_type(result)
    if kind == "annotation":
        return _scripted_annotation(result)
    return f"Could not interpret an unknown result kind: {kind!r}."


def _scripted_drug(r: dict[str, Any]) -> str:
    drug = r.get("drug", "this drug")
    target = r.get("target") or "its target"
    leads = r.get("leads") or []
    hero = r.get("hero") or (leads[0] if leads else {})
    if not hero:
        return f"{drug} has no replicated methylation->response leads to interpret."

    gene = hero.get("gene", "the top lead")
    rho = hero.get("gdsc_r")
    direction = hero.get("direction")
    ev = hero.get("evidence_level", "L?")
    sil = hero.get("silencing_r")

    sentences = [
        f"For {drug} ({target}), the top lead is **{gene}**: GDSC rho="
        f"{_fmt_r(rho)} ({_rho_means(rho)}), labelled '{direction}', evidence "
        f"{ev}.",
        f"Per the key, {ev} means {_level_means(ev)}, and the direction is "
        f"authoritative: {_dir_means(direction)}.",
        f"Its promoter methylation also tracks its own expression DOWN "
        f"(silencing_r={_fmt_r(sil)}), so the methylation functionally silences "
        f"{gene} rather than riding along as a passenger.",
    ]

    others = [l for l in leads[1:4] if l.get("gene")]
    if others:
        bits = ", ".join(
            f"{l['gene']} (rho={_fmt_r(l.get('gdsc_r'))}, {l.get('evidence_level','L?')}, "
            f"{l.get('direction')})"
            for l in others
        )
        sentences.append(f"Supporting leads point the same way: {bits}.")

    sensitive_n = sum(1 for l in leads if str(l.get("direction")).lower() == "sensitive")
    sentences.append(
        f"As a research-use hypothesis: across {len(leads)} replicated leads "
        f"({sensitive_n} sensitivity-associated), promoter hypermethylation of "
        f"these genes nominates cell lines likely to respond to {drug}, "
        f"independent of tissue of origin."
    )
    sentences.append(
        "Honest caveat: these are modest, cell-line-derived effects (directional "
        "concordance, not a large effect) and research-use hypotheses, not clinical "
        "claims."
    )
    return " ".join(sentences)


def _scripted_cancer_type(r: dict[str, Any]) -> str:
    ct = r.get("cancer_type", "this cancer type")
    leads = r.get("leads") or []
    if not leads:
        return f"No replicated methylation->response leads survived control for {ct}."

    replicated = [l for l in leads if l.get("replicated")]
    pool = replicated or leads
    hero = max(pool, key=lambda l: abs(l.get("tissue_rho") or 0))
    gene, drug = hero.get("gene"), hero.get("drug")
    rho = hero.get("tissue_rho")
    direction = hero.get("direction")
    ev = hero.get("evidence_level", "L?")
    also_in = hero.get("also_in") or []

    sentences = [
        f"In {ct}, the strongest lead is **{drug} -> {gene}**: tissue rho="
        f"{_fmt_r(rho)} ({_rho_means(rho)}), labelled '{direction}', evidence "
        f"{ev}.",
        f"Per the key, the direction is authoritative ({_dir_means(direction)}) "
        f"and {ev} means {_level_means(ev)}.",
    ]
    if also_in:
        sentences.append(
            f"The same {gene} signal recurs across indications "
            f"({', '.join(also_in[:5])}), so it is a cross-indication state, not a "
            f"one-tumor finding."
        )

    others = [l for l in pool if l is not hero][:3]
    if others:
        bits = ", ".join(
            f"{l.get('drug')}->{l.get('gene')} (rho={_fmt_r(l.get('tissue_rho'))}, "
            f"{l.get('evidence_level','L?')}, {l.get('direction')})"
            for l in others
        )
        sentences.append(f"Other leads: {bits}.")

    sentences.append(
        f"As a research-use hypothesis: {ct} carries {r.get('n_leads', len(leads))} "
        f"methylation->response leads, surfaced by molecular state rather than organ, "
        f"that an indication-gated search would not assemble."
    )
    sentences.append(
        "Honest caveat: effects are modest and cell-line-derived — research-use "
        "hypotheses, not clinical claims."
    )
    return " ".join(sentences)


def _scripted_annotation(r: dict[str, Any]) -> str:
    label = r.get("source_label", "this sample")
    tissue = r.get("tissue")
    preds = r.get("predictions") or []
    if not preds:
        return f"No drug reached the call threshold for {label}."

    sens = [p for p in preds if str(p.get("call")).lower() == "sensitive"]
    top = sens[0] if sens else preds[0]
    drug = top.get("drug")
    score = top.get("score")
    call = top.get("call")
    drivers = top.get("drivers") or []

    score_means = (
        "positive -> predicted SENSITIVE"
        if isinstance(score, (int, float)) and score > 0
        else ("negative -> predicted RESISTANT" if isinstance(score, (int, float)) and score < 0 else "near zero")
    )

    sentences = [
        f"For {label}{f' ({tissue})' if tissue else ''}, the top prediction is "
        f"**{drug}**: score={_fmt_r(score)} ({score_means}), call '{call}'.",
        "Per the key, a positive annotation score is sign-flipped vs rho and means "
        "predicted sensitivity.",
    ]
    if drivers:
        lvls = sorted({str(d.get("evidence_level")) for d in drivers if d.get("evidence_level")})
        bits = ", ".join(
            f"{d.get('gene')} ({d.get('evidence_level','L?')}, contrib={_fmt_r(d.get('contribution'))})"
            for d in drivers[:4]
        )
        sentences.append(
            f"Driver genes ({', '.join(lvls) if lvls else 'mixed evidence'}): {bits} — "
            f"each is a functionally silenced, response-associated lead, not a textbook gene."
        )

    other = [p for p in preds if p is not top][:3]
    if other:
        bits = ", ".join(
            f"{p.get('drug')} (score={_fmt_r(p.get('score'))}, {p.get('call')})"
            for p in other
        )
        sentences.append(f"Other calls: {bits}.")

    sentences.append(
        f"As a research-use hypothesis: this profile's methylation pattern nominates "
        f"{drug} as the leading candidate across {r.get('n_drugs_scored', len(preds))} "
        f"scored drugs."
    )
    if str(r.get("tissue")) and "cell" in str(label).lower() or "illustrative" in str(label).lower():
        sentences.append(
            "Honest caveat: this is an in-cohort cell line (training leakage) — a real "
            "patient sample would be leakage-free; effects are modest, research-use only."
        )
    else:
        sentences.append(
            "Honest caveat: effects are modest and cell-line-derived; research-use only, "
            "not a clinical recommendation."
        )
    return " ".join(sentences)


# ---------------------------------------------------------------------------
# LLM path
# ---------------------------------------------------------------------------


def _interpret_llm(kind: str, ref: str, result: dict[str, Any]) -> str:
    import json

    from anthropic import Anthropic

    client = Anthropic()
    trimmed = _trim_result(kind, result)
    system = [{"type": "text", "text": SYSTEM_PROMPT, "cache_control": {"type": "ephemeral"}}]
    user = (
        f"Result kind: {kind}; reference: {ref}.\n\n"
        f"ACTUAL RESULT (JSON):\n{json.dumps(trimmed)}\n\n"
        "Apply the key to THESE numbers. In 4-7 sentences, plain language: name the "
        "top lead(s)/prediction(s) with their actual gene + rho (or score) + evidence "
        "level, state what each sign/level MEANS per the key (sensitive vs resistant, "
        "how many screens replicated), what the overall result implies as a research-use "
        "hypothesis, and the key honest caveat. Be concrete; do not restate the key "
        "generically. You may use **bold** for the lead gene/drug."
    )
    response = client.messages.create(
        model=MODEL,
        max_tokens=MAX_TOKENS,
        system=system,
        messages=[{"role": "user", "content": user}],
    )
    text = next((b.text for b in response.content if b.type == "text"), "")
    if not text.strip():
        raise RuntimeError("empty LLM response")
    return text.strip()


# ---------------------------------------------------------------------------
# Public entrypoint
# ---------------------------------------------------------------------------


def _has_key() -> bool:
    return bool(os.environ.get("ANTHROPIC_API_KEY"))


def run_interpret(kind: str, ref: str) -> dict[str, Any]:
    """Interpret the actual result for (*kind*, *ref*) by applying the key."""
    result = _fetch_result(kind, ref)

    if _has_key():
        try:
            text = _interpret_llm(kind, ref, result)
            return {"kind": kind, "ref": ref, "interpretation": text, "used_llm": True}
        except Exception as exc:  # demo safety: never crash on an API/LLM error
            print(
                f"[strata] interpret LLM failed ({type(exc).__name__}: {exc}); "
                f"falling back to scripted.",
                file=sys.stderr,
            )
            text = _interpret_scripted(kind, ref, result)
            return {
                "kind": kind,
                "ref": ref,
                "interpretation": text,
                "used_llm": False,
                "llm_error": f"{type(exc).__name__}: {exc}",
            }

    return {
        "kind": kind,
        "ref": ref,
        "interpretation": _interpret_scripted(kind, ref, result),
        "used_llm": False,
    }


if __name__ == "__main__":
    k = sys.argv[1] if len(sys.argv) > 1 else "drug"
    rf = sys.argv[2] if len(sys.argv) > 2 else "olaparib"
    out = run_interpret(k, rf)
    print(f"=== interpret [{'llm' if out['used_llm'] else 'scripted'}] — {k}/{rf} ===")
    if out.get("llm_error"):
        print(f"(llm_error: {out['llm_error']})")
    print(out["interpretation"])
