"""Strata Report API — thin FastAPI server over precomputed standardized JSON reports.

Serves the reports written by ``scripts/build_reports.py`` into ``api/reports/``:
  - ``api/reports/drug/<slug>.json``        — one per drug
  - ``api/reports/cancer_type/<slug>.json``  — one per cancer type
  - ``api/reports/manifest.json``            — index

The Pydantic models below mirror the exact JSON shape emitted by
``build_reports.py`` so the auto-generated OpenAPI schema is the standardized
contract an agent can call against.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

from agent import run_ask
from annotate import annotate_profile, load_annotation_model
from interpret import run_interpret

REPORTS_DIR = Path(__file__).parent / "reports"
DEMO_DIR = REPORTS_DIR / "demo_profiles"

# ---------------------------------------------------------------------------
# Models — mirror scripts/build_reports.py exactly
# ---------------------------------------------------------------------------


class Recurrence(BaseModel):
    tissues_tested: int
    same_direction: int
    significant: int


class Lead(BaseModel):
    gene: str
    direction: str  # "sensitive" | "resistant"
    gdsc_r: float
    prism_r: Optional[float] = None
    prism_p: Optional[float] = None
    ctrp_r: Optional[float] = None
    ctrp_p: Optional[float] = None
    silencing_r: float
    screens: List[str]
    n_screens: int
    evidence_level: str
    replicated: bool
    recurrence: Recurrence
    gene_note: str


class Funnel(BaseModel):
    genes_scanned: int
    response_candidates: int
    functionally_silenced: int
    externally_replicated: int
    lines: int


class DrugReport(BaseModel):
    type: str = "drug"
    drug: str
    target: Optional[str] = None
    pathway: Optional[str] = None
    funnel: Funnel
    leads: List[Lead]
    hero: Optional[Lead] = None
    honesty: List[str]


class CancerTypeLead(BaseModel):
    drug: str
    gene: str
    direction: str
    tissue_rho: float
    tissue_p: float
    evidence_level: str
    replicated: bool
    also_in: List[str]


class CancerTypeReport(BaseModel):
    type: str = "cancer_type"
    cancer_type: str
    n_leads: int
    leads: List[CancerTypeLead]
    connects_to: List[str]
    honesty: List[str]


class DrugManifestEntry(BaseModel):
    name: str
    slug: str
    target: Optional[str] = None
    n_leads: int
    replicated: int


class CancerTypeManifestEntry(BaseModel):
    name: str
    slug: str
    n_leads: int


class Manifest(BaseModel):
    drugs: List[DrugManifestEntry]
    cancer_types: List[CancerTypeManifestEntry]
    gates: List[Dict[str, Any]]


# --- Annotation (patient scorer) models ---


class Driver(BaseModel):
    gene: str
    patient_pct: float
    contribution: float
    evidence_level: str


class Prediction(BaseModel):
    drug: str
    score: float
    call: str  # "sensitive" | "resistant" | "neutral"
    n_leads_scored: int
    drivers: List[Driver]


class AnnotationResult(BaseModel):
    n_drugs_scored: int
    predictions: List[Prediction]
    honesty: List[str]
    source_label: Optional[str] = None
    tissue: Optional[str] = None
    ground_truth_auc: Optional[Dict[str, float]] = None


class AnnotateRequest(BaseModel):
    profile: Optional[Dict[str, float]] = None
    demo_id: Optional[str] = None


class DemoProfileEntry(BaseModel):
    id: str
    label: str
    tissue: Optional[str] = None
    source: Optional[str] = None
    ground_truth_auc: Optional[Dict[str, float]] = None


class DemoProfile(BaseModel):
    id: str
    label: str
    name: Optional[str] = None
    cosmic_id: Optional[str] = None
    tissue: Optional[str] = None
    source: Optional[str] = None
    note: Optional[str] = None
    n_genes: Optional[int] = None
    ground_truth_auc: Optional[Dict[str, float]] = None
    profile: Dict[str, float]


# --- Agent (live /ask) models ---


class AskRequest(BaseModel):
    query: str


class AskStep(BaseModel):
    tool: str
    args: Dict[str, Any]
    result_summary: str


class AskResult(BaseModel):
    query: str
    steps: List[AskStep]
    answer: str
    used_llm: bool
    llm_error: Optional[str] = None


# --- Interpret (apply the interpretation key to an actual result) models ---


class InterpretRequest(BaseModel):
    kind: str  # "drug" | "cancer_type" | "annotation"
    ref: str  # drug/cancer-type slug, or demo profile id


class InterpretResult(BaseModel):
    kind: str
    ref: str
    interpretation: str
    used_llm: bool
    llm_error: Optional[str] = None


class Root(BaseModel):
    name: str
    docs: str
    drugs: int
    cancer_types: int
    annotate: bool
    demo_profiles: int


class Health(BaseModel):
    status: str


# ---------------------------------------------------------------------------
# App
# ---------------------------------------------------------------------------

app = FastAPI(
    title="Strata Report API",
    version="1.0.0",
    description=(
        "Serves precomputed, standardized cross-indication methylation->drug-response "
        "reports (drug + cancer-type) so an agent can call them. The OpenAPI schema is "
        "the standardized contract; report shapes are produced by scripts/build_reports.py."
    ),
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # demo: open; tighten to ["https://methylkb.vercel.app"] for prod
    allow_credentials=False,
    allow_methods=["GET", "POST"],
    allow_headers=["*"],
)


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text())


def _load_manifest() -> Dict[str, Any]:
    path = REPORTS_DIR / "manifest.json"
    if not path.exists():
        # Empty-but-valid manifest so the service is healthy before the batch finishes.
        return {"drugs": [], "cancer_types": [], "gates": []}
    return _read_json(path)


def _load_demo_index() -> List[Dict[str, Any]]:
    path = DEMO_DIR / "index.json"
    if not path.exists():
        return []
    return _read_json(path).get("profiles", [])


def _load_demo_profile(demo_id: str) -> Optional[Dict[str, Any]]:
    path = DEMO_DIR / f"{demo_id}.json"
    if not path.exists():
        return None
    return _read_json(path)


@app.get("/", response_model=Root)
def root() -> Root:
    m = _load_manifest()
    annotate_ok = (REPORTS_DIR / "annotation_model.json").exists()
    return Root(
        name="Strata Report API",
        docs="/docs",
        drugs=len(m.get("drugs", [])),
        cancer_types=len(m.get("cancer_types", [])),
        annotate=annotate_ok,
        demo_profiles=len(_load_demo_index()),
    )


@app.get("/health", response_model=Health)
def health() -> Health:
    return Health(status="ok")


@app.get("/manifest", response_model=Manifest)
def manifest() -> Manifest:
    return Manifest(**_load_manifest())


@app.get("/drugs", response_model=List[DrugManifestEntry])
def drugs() -> List[DrugManifestEntry]:
    return [DrugManifestEntry(**d) for d in _load_manifest().get("drugs", [])]


@app.get("/cancer_types", response_model=List[CancerTypeManifestEntry])
def cancer_types() -> List[CancerTypeManifestEntry]:
    return [CancerTypeManifestEntry(**c) for c in _load_manifest().get("cancer_types", [])]


@app.get("/report/drug/{slug}", response_model=DrugReport)
def drug_report(slug: str) -> DrugReport:
    path = REPORTS_DIR / "drug" / f"{slug}.json"
    if not path.exists():
        raise HTTPException(status_code=404, detail=f"drug report not found: {slug}")
    return DrugReport(**_read_json(path))


@app.get("/report/cancer_type/{slug}", response_model=CancerTypeReport)
def cancer_type_report(slug: str) -> CancerTypeReport:
    path = REPORTS_DIR / "cancer_type" / f"{slug}.json"
    if not path.exists():
        raise HTTPException(status_code=404, detail=f"cancer_type report not found: {slug}")
    return CancerTypeReport(**_read_json(path))


# ---------------------------------------------------------------------------
# Patient annotation (the scorer)
# ---------------------------------------------------------------------------


@app.get("/demo_profiles", response_model=List[DemoProfileEntry])
def demo_profiles() -> List[DemoProfileEntry]:
    return [DemoProfileEntry(**p) for p in _load_demo_index()]


@app.get("/demo_profile/{demo_id}", response_model=DemoProfile)
def demo_profile(demo_id: str) -> DemoProfile:
    prof = _load_demo_profile(demo_id)
    if prof is None:
        raise HTTPException(status_code=404, detail=f"demo profile not found: {demo_id}")
    return DemoProfile(**prof)


@app.post("/annotate", response_model=AnnotationResult)
def annotate(req: AnnotateRequest) -> AnnotationResult:
    if req.demo_id is not None:
        prof = _load_demo_profile(req.demo_id)
        if prof is None:
            raise HTTPException(status_code=404, detail=f"demo profile not found: {req.demo_id}")
        out = annotate_profile(prof["profile"], load_annotation_model())
        out["source_label"] = prof.get("label")
        out["tissue"] = prof.get("tissue")
        out["ground_truth_auc"] = prof.get("ground_truth_auc")
        return AnnotationResult(**out)
    if req.profile is not None:
        out = annotate_profile(req.profile, load_annotation_model())
        return AnnotationResult(**out)
    raise HTTPException(status_code=422, detail="provide either 'profile' or 'demo_id'")


# ---------------------------------------------------------------------------
# Live agent (the full loop)
# ---------------------------------------------------------------------------


@app.post("/ask", response_model=AskResult)
def ask(req: AskRequest) -> AskResult:
    return AskResult(**run_ask(req.query))


# ---------------------------------------------------------------------------
# Interpret a result (applies the fixed interpretation key to the real numbers)
# ---------------------------------------------------------------------------


@app.post("/interpret", response_model=InterpretResult)
def interpret(req: InterpretRequest) -> InterpretResult:
    if req.kind not in {"drug", "cancer_type", "annotation"}:
        raise HTTPException(
            status_code=422,
            detail="kind must be one of: drug, cancer_type, annotation",
        )
    return InterpretResult(**run_interpret(req.kind, req.ref))
