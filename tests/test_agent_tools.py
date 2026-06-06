"""Unit tests for strata.agent.tools — all heavy loaders mocked out."""

import numpy as np
import pandas as pd
import pytest

from strata.agent import tools


# ---------------------------------------------------------------------------
# Fake data factories
# ---------------------------------------------------------------------------


def _fake_meth():
    rng = np.random.default_rng(0)
    n = 60
    idx = [str(1000 + i) for i in range(n)]
    df = pd.DataFrame(rng.uniform(0, 1, (n, 20)), index=idx, columns=[f"G{j}" for j in range(20)])
    # "GOOD" is a gene that should correlate strongly with AUC (set below)
    df["GOOD"] = np.clip(1 - rng.uniform(0, 1, n), 0, 1)
    return df


def _fake_drug():
    rng = np.random.default_rng(1)
    n = 60
    idx = [str(1000 + i) for i in range(n)]
    return pd.DataFrame(
        {
            "COSMIC_ID": idx,
            "drug_name": ["DrugX"] * n,
            "ln_ic50": rng.normal(0, 1, n),
            "auc": rng.uniform(0, 1, n),
        }
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_dispatch_rank_markers(monkeypatch):
    monkeypatch.setattr(tools, "_meth", _fake_meth)
    monkeypatch.setattr(tools, "_drug", _fake_drug)
    out = tools.dispatch_tool("rank_markers", {"drug_name": "DrugX", "top_n": 5})
    assert "markers" in out, f"unexpected output: {out}"
    assert isinstance(out["markers"], list)
    # Each marker dict must have the three required keys
    assert all({"gene", "rho", "qvalue"} <= set(m) for m in out["markers"])
    # Values must be Python scalars (JSON-serialisable, not numpy types)
    for m in out["markers"]:
        assert isinstance(m["rho"], float)
        assert isinstance(m["qvalue"], float)
        assert isinstance(m["gene"], str)


def test_dispatch_rank_markers_drug_not_found(monkeypatch):
    monkeypatch.setattr(tools, "_meth", _fake_meth)
    monkeypatch.setattr(tools, "_drug", _fake_drug)
    out = tools.dispatch_tool("rank_markers", {"drug_name": "NoSuchDrug", "top_n": 5})
    assert "error" in out
    assert "NoSuchDrug" in out["error"]


def test_dispatch_list_drugs(monkeypatch):
    monkeypatch.setattr(tools, "_drug", _fake_drug)
    out = tools.dispatch_tool("list_drugs", {})
    assert "drugs" in out
    assert "DrugX" in out["drugs"]


def test_dispatch_list_drugs_filter(monkeypatch):
    monkeypatch.setattr(tools, "_drug", _fake_drug)
    out = tools.dispatch_tool("list_drugs", {"contains": "Drug"})
    assert "drugs" in out and "DrugX" in out["drugs"]

    out_empty = tools.dispatch_tool("list_drugs", {"contains": "ZZZNOMATCH"})
    assert out_empty["drugs"] == []


def test_dispatch_evaluate_subgroup(monkeypatch):
    monkeypatch.setattr(tools, "_meth", _fake_meth)
    monkeypatch.setattr(tools, "_drug", _fake_drug)
    out = tools.dispatch_tool("evaluate_subgroup", {"drug_name": "DrugX", "marker": "G0"})
    assert "grade" in out
    assert "pvalue" in out
    assert isinstance(out["pvalue"], float)
    assert isinstance(out["effect_size"], float)
    assert isinstance(out["responder_label"], int)
    assert isinstance(out["n"], int)


def test_dispatch_evaluate_subgroup_bad_marker(monkeypatch):
    monkeypatch.setattr(tools, "_meth", _fake_meth)
    monkeypatch.setattr(tools, "_drug", _fake_drug)
    out = tools.dispatch_tool("evaluate_subgroup", {"drug_name": "DrugX", "marker": "NOTAREAL"})
    assert "error" in out
    assert "NOTAREAL" in out["error"]


def test_dispatch_write_opportunity_brief(monkeypatch):
    out = tools.dispatch_tool(
        "write_opportunity_brief",
        {"drug_name": "DrugX", "marker": "G0", "rationale": "It looks promising."},
    )
    assert out == {"ok": True, "drug": "DrugX", "marker": "G0"}


def test_dispatch_unknown(monkeypatch):
    out = tools.dispatch_tool("nope", {})
    assert "error" in out
    assert "nope" in out["error"]


def test_tools_schema_structure():
    """TOOLS must be a non-empty list of dicts with the required Anthropic keys."""
    assert isinstance(tools.TOOLS, list) and len(tools.TOOLS) > 0
    names = {t["name"] for t in tools.TOOLS}
    for required in ("list_drugs", "rank_markers", "evaluate_subgroup", "write_opportunity_brief"):
        assert required in names, f"missing tool schema: {required}"
    for t in tools.TOOLS:
        assert "description" in t
        assert "input_schema" in t
        assert t["input_schema"]["type"] == "object"
