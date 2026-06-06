from strata.agent import loop


def test_scripted_run_returns_brief(monkeypatch):
    # force scripted mode regardless of env; mock dispatch_tool to avoid heavy data load
    def fake_dispatch(name, args):
        if name == "rank_markers":
            return {
                "drug": args["drug_name"],
                "markers": [
                    {"gene": "CDKN2A", "rho": -0.336, "qvalue": 7.8e-25},
                    {"gene": "MTAP", "rho": -0.44, "qvalue": 6e-43},
                ],
            }
        if name == "evaluate_subgroup":
            return {
                "drug": args["drug_name"],
                "marker": args["marker"],
                "grade": "neutral",
                "pvalue": 7.3e-18,
                "effect_size": 0.033,
                "responder_label": 1,
                "n": 968,
            }
        return {"error": "x"}

    monkeypatch.setattr(loop, "dispatch_tool", fake_dispatch)
    out = loop.run_discovery_agent("Palbociclib", use_llm=False)
    assert out["mode"] == "scripted"
    assert out["drug"] == "Palbociclib"
    assert out["brief"] and "Palbociclib" in out["brief"]
    assert any(s["type"] == "tool_call" for s in out["transcript"])


def test_scripted_transcript_shape(monkeypatch):
    def fake_dispatch(name, args):
        if name == "rank_markers":
            return {
                "drug": args["drug_name"],
                "markers": [{"gene": "CDKN2A", "rho": -0.336, "qvalue": 7.8e-25}],
            }
        if name == "evaluate_subgroup":
            return {
                "drug": args["drug_name"],
                "marker": args["marker"],
                "grade": "neutral",
                "pvalue": 7.3e-18,
                "effect_size": 0.033,
                "responder_label": 1,
                "n": 968,
            }
        return {"error": "x"}

    monkeypatch.setattr(loop, "dispatch_tool", fake_dispatch)
    out = loop.run_discovery_agent("Palbociclib", use_llm=False)
    types = [s["type"] for s in out["transcript"]]
    # at least one thought, a tool_call, and a tool_result
    assert "thought" in types
    assert "tool_call" in types
    assert "tool_result" in types
    # CDKN2A should be recognized as the known positive control
    text = " ".join(s.get("text", "") for s in out["transcript"] if s["type"] == "thought")
    assert "CDKN2A" in text


def test_llm_failure_falls_back_to_scripted(monkeypatch):
    """A live-agent API error must degrade to the scripted demo, never crash."""

    def boom(*_a, **_k):
        raise RuntimeError("api down")

    def fake_dispatch(name, args):
        if name == "rank_markers":
            return {
                "drug": args["drug_name"],
                "markers": [{"gene": "CDKN2A", "rho": -0.336, "qvalue": 7.8e-25}],
            }
        if name == "evaluate_subgroup":
            return {
                "drug": args["drug_name"],
                "marker": args["marker"],
                "grade": "neutral",
                "pvalue": 7.3e-18,
                "effect_size": 0.033,
                "responder_label": 1,
                "n": 968,
            }
        return {"error": "x"}

    monkeypatch.setattr(loop, "_run_llm", boom)
    monkeypatch.setattr(loop, "dispatch_tool", fake_dispatch)
    # force the LLM path; it raises -> must fall back to scripted with the error recorded
    out = loop.run_discovery_agent("Palbociclib", use_llm=True)
    assert out["mode"] == "scripted"
    assert "api down" in out.get("llm_error", "")
    assert out["brief"] and "Palbociclib" in out["brief"]
