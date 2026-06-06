"""Tests for the report agent (scripted fallback path).

These exercise the deterministic ``use_llm=False`` path, which fetches from a
locally running Strata Report API (``api/main.py`` on port 8077). If the API is
not reachable, the tests skip rather than fail.
"""

from __future__ import annotations

import urllib.error
import urllib.request

import pytest

from strata.agent.report_agent import run_report_agent

API_BASE = "http://localhost:8077"


def _api_up() -> bool:
    try:
        with urllib.request.urlopen(f"{API_BASE}/health", timeout=3) as resp:
            return resp.status == 200
    except (urllib.error.URLError, OSError):
        return False


pytestmark = pytest.mark.skipif(
    not _api_up(), reason="Strata Report API not running on localhost:8077"
)


def test_scripted_olaparib_returns_full_transcript() -> None:
    out = run_report_agent(
        "What drug-response opportunities does olaparib have that an "
        "oncologist would miss?",
        api_base=API_BASE,
        use_llm=False,
    )
    # Shape.
    assert out["used_llm"] is False
    assert out["answer"] and isinstance(out["answer"], str)
    assert out["report"] and isinstance(out["report"], dict)
    assert out["steps"] and isinstance(out["steps"], list)

    # It fetched the olaparib drug report.
    assert out["report"].get("type") == "drug"
    assert out["report"].get("drug") == "Olaparib"

    # The narration names the hero gene and the three-gates frame.
    assert "ID2" in out["answer"]
    assert "Gate 1" in out["answer"]
    assert "Gate 2" in out["answer"]
    assert "Gate 3" in out["answer"]

    # Steps include both tool calls.
    tools = [s["tool"] for s in out["steps"]]
    assert "list_available" in tools
    assert "get_drug_report" in tools
    for s in out["steps"]:
        assert s["result_summary"]


def test_scripted_cancer_type_query() -> None:
    out = run_report_agent(
        "What's actionable in lung cancer?", api_base=API_BASE, use_llm=False
    )
    assert out["used_llm"] is False
    assert out["report"] and out["report"].get("type") == "cancer_type"
    assert out["report"].get("cancer_type") == "Lung"
    assert out["answer"]
    assert "get_cancer_type_report" in [s["tool"] for s in out["steps"]]


def test_scripted_mek_inhibitor_hint() -> None:
    # "a MEK inhibitor" (no catalog name) should resolve to Trametinib.
    out = run_report_agent(
        "Show me a cross-indication methylation marker for a MEK inhibitor.",
        api_base=API_BASE,
        use_llm=False,
    )
    assert out["used_llm"] is False
    assert out["report"] and out["report"].get("drug") == "Trametinib"
    assert out["answer"]
