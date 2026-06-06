"""Discovery agent loop — the demo centerpiece.

``run_discovery_agent(drug_name, goal, use_llm)`` runs a small agentic loop that
scans methylation markers against a drug's response, picks the most biologically
interpretable marker, evaluates a median-split responder subgroup, and writes a
commercial-development opportunity brief.

Two execution paths produce the **same output shape** so the UI is identical:

- ``_run_llm`` — a real Anthropic tool-use loop (model decides the steps). Used
  when ``ANTHROPIC_API_KEY`` is present.
- ``_run_scripted`` — a deterministic emulation that calls ``dispatch_tool``
  directly and records natural-language reasoning between calls. Used when no
  key is set, so the demo always works offline.

Framing rules (encoded in both paths):

- The real signal is the **continuous Spearman rho / q-value**, not the
  median-split subgroup effect (which is small). Lead with rho/q.
- The median split is the *operational companion-marker cut*, not the headline.
- Novel hits are screen-derived **hypotheses**; their credibility comes from the
  same screen recovering the known CDKN2A/9p21 -> palbociclib biology as a
  positive control.
"""

from __future__ import annotations

import os
from typing import Any

from dotenv import load_dotenv

from strata.agent.brief import write_brief
from strata.agent.tools import TOOLS, dispatch_tool
from strata.engine.associate import AssociationResult

load_dotenv(override=False)

MODEL = "claude-opus-4-8"
MAX_TOKENS = 2000
MAX_TURNS = 10  # safety bound on the LLM tool-use loop

DEFAULT_GOAL = (
    "Find a methylation-defined responder subgroup for this drug and write a "
    "commercial-development opportunity brief."
)

# Genes that, when surfaced as top markers, recover known CDK4/6-inhibitor
# biology (CDKN2A/p16 and its 9p21 co-deleted neighbor MTAP). Recovering this
# axis is what licenses trust in the novel hits.
_KNOWN_CONTROL_GENES = {"CDKN2A", "MTAP"}

SYSTEM_PROMPT = (
    "You are a pharmaco-epigenomics discovery agent working for a pharmaceutical "
    "portfolio team (think Pfizer business development). Your job: given a drug, "
    "find a DNA-methylation-defined responder subgroup and write a commercial-"
    "development opportunity brief.\n\n"
    "PROTOCOL — follow these steps in order:\n"
    "  1. Call rank_markers to scan all ~14,600 gene-methylation features against "
    "drug response (Spearman correlation with AUC; lower AUC = more sensitive).\n"
    "  2. Examine the top markers and pick the single most biologically "
    "interpretable one. Recovering KNOWN biology builds trust — e.g. for a "
    "CDK4/6 inhibitor, CDKN2A (p16) or its 9p21 co-deleted neighbor MTAP is the "
    "textbook axis, so prefer it when it appears near the top.\n"
    "  3. Call evaluate_subgroup on that marker to define the median-split "
    "high- vs low-methylation companion-marker cut.\n"
    "  4. Write the opportunity brief in a portfolio-ROI framing: reduced "
    "opportunity-identification time and cost, trial-enrichment and peak-revenue "
    "upside.\n\n"
    "FRAMING RULES (critical):\n"
    "  - LEAD with the continuous Spearman rho and q-value as the evidence. That "
    "is the real signal.\n"
    "  - Describe the median-split subgroup as the OPERATIONAL companion-marker "
    "cut, not the headline effect — the split effect size is deliberately small.\n"
    "  - Present any novel marker as a SCREEN-DERIVED HYPOTHESIS. Its credibility "
    "comes from the same unbiased screen also recovering the known "
    "CDKN2A/9p21 -> palbociclib positive control.\n"
)


def _has_key() -> bool:
    return bool(os.environ.get("ANTHROPIC_API_KEY"))


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def run_discovery_agent(
    drug_name: str,
    goal: str | None = None,
    use_llm: bool | None = None,
) -> dict[str, Any]:
    """Run the discovery agent for *drug_name*.

    Returns a dict with keys: ``drug``, ``goal``, ``transcript`` (ordered list of
    ``{"type": "thought"|"tool_call"|"tool_result", ...}`` steps), ``brief``
    (markdown), and ``mode`` (``"llm"`` | ``"scripted"``).
    """
    goal = goal or DEFAULT_GOAL
    if use_llm is None:
        use_llm = _has_key()

    if use_llm:
        try:
            return _run_llm(drug_name, goal)
        except Exception as exc:  # demo safety: never crash the stage on an API error
            import sys

            print(
                f"[strata] live agent failed ({type(exc).__name__}: {exc}); "
                f"falling back to scripted mode.",
                file=sys.stderr,
            )
            result = _run_scripted(drug_name, goal)
            result["llm_error"] = f"{type(exc).__name__}: {exc}"
            return result
    return _run_scripted(drug_name, goal)


# ---------------------------------------------------------------------------
# Scripted fallback (no API key required)
# ---------------------------------------------------------------------------


def _run_scripted(drug_name: str, goal: str) -> dict[str, Any]:
    transcript: list[dict[str, Any]] = []

    # Step 1 — thought
    transcript.append(
        {
            "type": "thought",
            "text": (
                f"Goal: find a methylation-defined responder subgroup for "
                f"{drug_name}. I'll scan all ~14.6k gene-methylation features "
                f"against drug response (Spearman rho vs AUC; lower AUC = more "
                f"sensitive)."
            ),
        }
    )

    # Step 2 — rank_markers
    rank_args = {"drug_name": drug_name, "top_n": 15}
    transcript.append({"type": "tool_call", "name": "rank_markers", "args": rank_args})
    rank_res = dispatch_tool("rank_markers", rank_args)
    transcript.append({"type": "tool_result", "name": "rank_markers", "result": rank_res})

    markers = rank_res.get("markers", [])
    if not markers:
        transcript.append(
            {
                "type": "thought",
                "text": (
                    f"No FDR-significant methylation markers for {drug_name}. "
                    f"No defensible responder subgroup to report."
                ),
            }
        )
        return {
            "drug": drug_name,
            "goal": goal,
            "transcript": transcript,
            "brief": (
                f"# Commercial-Development Opportunity — {drug_name}\n\n"
                f"No FDR-significant methylation markers were recovered for "
                f"{drug_name}; no responder subgroup is supported by the screen."
            ),
            "mode": "scripted",
        }

    # Step 3 — pick the marker (top by |rho|) and interpret
    top = markers[0]
    marker = str(top["gene"])
    is_known = marker.upper() in _KNOWN_CONTROL_GENES
    if is_known:
        interpretation = (
            f"Top marker is {marker} (rho={top['rho']:.3f}, q={top['qvalue']:.2g}). "
            f"This RECOVERS the known 9p21 / CDK4/6-inhibitor biology "
            f"(CDKN2A-p16 and its co-deleted neighbor MTAP) — a built-in positive "
            f"control that licenses trust in the rest of the screen. I'll carry "
            f"{marker} forward as the companion marker."
        )
    else:
        known_in_top = [
            str(m["gene"])
            for m in markers
            if str(m["gene"]).upper() in _KNOWN_CONTROL_GENES
        ]
        ctrl = (
            f" The screen also surfaces {known_in_top[0]}, recovering the known "
            f"CDKN2A/9p21 positive control, which licenses trust in this hit."
            if known_in_top
            else ""
        )
        interpretation = (
            f"Top marker is {marker} (rho={top['rho']:.3f}, q={top['qvalue']:.2g}). "
            f"This is a SCREEN-DERIVED HYPOTHESIS rather than a textbook axis.{ctrl} "
            f"I'll carry {marker} forward as the candidate companion marker."
        )
    transcript.append({"type": "thought", "text": interpretation})

    # Step 4 — evaluate_subgroup
    eval_args = {"drug_name": drug_name, "marker": marker}
    transcript.append({"type": "tool_call", "name": "evaluate_subgroup", "args": eval_args})
    eval_res = dispatch_tool("evaluate_subgroup", eval_args)
    transcript.append(
        {"type": "tool_result", "name": "evaluate_subgroup", "result": eval_res}
    )

    # Step 5 — interpret (lead with rho/q; subgroup is the operational cut)
    transcript.append(
        {
            "type": "thought",
            "text": (
                f"The real signal is the continuous association: rho={top['rho']:.3f}, "
                f"q={top['qvalue']:.2g} across the full panel. The median-split "
                f"{marker} subgroup (grade={eval_res.get('grade')}, "
                f"p={eval_res.get('pvalue'):.2g}, n={eval_res.get('n')}) is the "
                f"OPERATIONAL companion-marker cut a trial would actually use — its "
                f"effect size is small by design; the continuous rho/q is the case."
            ),
        }
    )

    # Step 6 — build the brief
    assoc = AssociationResult(
        effect_size=float(eval_res.get("effect_size", 0.0)),
        pvalue=float(eval_res.get("pvalue", 1.0)),
        responder_label=int(eval_res.get("responder_label", 0)),
        grade=str(eval_res.get("grade", "neutral")),
    )
    driver_genes = [str(m["gene"]) for m in markers]
    n_samples = int(eval_res.get("n", 0))
    brief = write_brief(drug_name, assoc, driver_genes, n_samples)

    return {
        "drug": drug_name,
        "goal": goal,
        "transcript": transcript,
        "brief": brief,
        "mode": "scripted",
    }


# ---------------------------------------------------------------------------
# LLM path (Anthropic tool-use loop) — structurally correct; needs a key
# ---------------------------------------------------------------------------


def _run_llm(drug_name: str, goal: str) -> dict[str, Any]:
    from anthropic import Anthropic

    client = Anthropic()
    transcript: list[dict[str, Any]] = []

    system = [
        {
            "type": "text",
            "text": SYSTEM_PROMPT,
            "cache_control": {"type": "ephemeral"},
        }
    ]

    user_text = (
        f"Drug: {drug_name}\n"
        f"Goal: {goal}\n\n"
        f"Run the protocol and finish with the opportunity brief as your final "
        f"message text (markdown)."
    )
    messages: list[dict[str, Any]] = [
        {"role": "user", "content": user_text}
    ]

    brief = ""
    for _ in range(MAX_TURNS):
        response = client.messages.create(
            model=MODEL,
            max_tokens=MAX_TOKENS,
            system=system,
            tools=TOOLS,
            messages=messages,
        )

        # Record any thought text the model emitted alongside tool calls.
        for block in response.content:
            if block.type == "text" and block.text.strip():
                transcript.append({"type": "thought", "text": block.text})

        if response.stop_reason != "tool_use":
            # Final answer — the last text block is the brief.
            brief = next(
                (b.text for b in response.content if b.type == "text"), ""
            )
            break

        # Append the assistant turn (including tool_use blocks), then execute.
        messages.append({"role": "assistant", "content": response.content})

        tool_results: list[dict[str, Any]] = []
        for block in response.content:
            if block.type != "tool_use":
                continue
            args = dict(block.input)
            transcript.append(
                {"type": "tool_call", "name": block.name, "args": args}
            )
            result = dispatch_tool(block.name, args)
            transcript.append(
                {"type": "tool_result", "name": block.name, "result": result}
            )
            tool_results.append(
                {
                    "type": "tool_result",
                    "tool_use_id": block.id,
                    "content": _json(result),
                }
            )

        messages.append({"role": "user", "content": tool_results})

    return {
        "drug": drug_name,
        "goal": goal,
        "transcript": transcript,
        "brief": brief,
        "mode": "llm",
    }


def _json(obj: Any) -> str:
    import json

    return json.dumps(obj)


# ---------------------------------------------------------------------------
# Runnable demo entrypoint
# ---------------------------------------------------------------------------


def _pretty_print(result: dict[str, Any]) -> None:
    print(f"\n=== Discovery agent [{result['mode']}] — {result['drug']} ===")
    print(f"Goal: {result['goal']}\n")
    for step in result["transcript"]:
        kind = step["type"]
        if kind == "thought":
            print(f"[thought] {step['text']}\n")
        elif kind == "tool_call":
            print(f"[tool_call] {step['name']}({step['args']})")
        elif kind == "tool_result":
            res = step["result"]
            if step["name"] == "rank_markers":
                markers = res.get("markers", [])[:6]
                rendered = ", ".join(
                    f"{m['gene']}(rho={m['rho']:.3f}, q={m['qvalue']:.2g})"
                    for m in markers
                )
                print(f"[tool_result] {step['name']}: {rendered}")
            else:
                print(f"[tool_result] {step['name']}: {res}")
            print()
    print("--- BRIEF ---")
    print(result["brief"])


if __name__ == "__main__":
    _pretty_print(run_discovery_agent("Palbociclib"))
