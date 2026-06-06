# MethylGKB · Strata — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A demo where an agent does unsupervised methylation subgroup discovery on GDSC, recovers the MGMT→temozolomide responder subgroup as a positive control, surfaces a novel methylation-defined responder subgroup, anchors it in TCGA-LAML patients, and writes a commercial-development opportunity brief.

**Architecture:** A deterministic Python engine (`strata/`) — data loaders → feature selection → clustering → response association → driver annotation → figures — exposed both as a CLI (rung A, no LLM) and as Claude tool-use functions (rung B, the agent). Output is rendered into the existing MethylGKB `site/`. Nothing existing in the repo is removed.

**Tech Stack:** Python 3.12 via `uv`; pandas, numpy, scikit-learn, scipy/statsmodels, lifelines, umap-learn, matplotlib, plotly; Anthropic SDK (Claude tool-use) for the agent layer.

**Build order is strict: Phase 0 (MGMT gate) → Phase A (deterministic floor) → Phase B (agent + patient) → Phase C (stretch).** Phase A alone is a complete, demoable result. Do not start Phase B until A is green.

---

## File Structure

```
strata/
  __init__.py
  config.py            # paths, constants (DATA_DIR, N_TOP_FEATURES, RANDOM_STATE)
  data/
    __init__.py
    gdsc.py            # download + load GDSC methylation, drug response, line annotations
    tcga_laml.py       # download + load TCGA-LAML methylation + clinical/survival
  engine/
    __init__.py
    features.py        # select_variable_features
    cluster.py         # cluster_samples
    associate.py       # AssociationResult, test_drug_association, score_responder, test_survival_association
    annotate.py        # differential_cpgs, annotate_cpgs
  figures.py           # plot_umap, plot_driver_heatmap, plot_auc_box, plot_km
  agent/
    __init__.py
    tools.py           # TOOLS schema + dispatch wrapping engine functions
    loop.py            # run_discovery_agent (Claude tool-use orchestration)
    brief.py           # write_brief (templated for rung A; agent-written for rung B)
  cli.py               # rung-A end-to-end deterministic pipeline
tests/
  conftest.py          # synthetic fixtures
  test_gdsc.py
  test_tcga.py
  test_features.py
  test_cluster.py
  test_associate.py
  test_annotate.py
  test_agent_tools.py
data/                  # cached downloads — gitignored
notes/
  GDSC_SCHEMA.md       # actual file paths + column names discovered in Task 1
  MGMT_GATE.md         # the go/no-go verdict from Task 2
```

**Interface contracts (consistent names used throughout):**
- `load_gdsc_methylation() -> pd.DataFrame` — rows = COSMIC_ID (str), cols = CpG-island id (str), values = beta [0,1].
- `load_gdsc_drug_response() -> pd.DataFrame` — long: columns `[COSMIC_ID, drug_name, ln_ic50, auc]`.
- `load_gdsc_annotations() -> pd.DataFrame` — index COSMIC_ID, cols include `cell_line_name`, `tissue`.
- `select_variable_features(meth, n_top=2000) -> list[str]`
- `cluster_samples(meth, features, k=None, random_state=0) -> pd.Series` — index = COSMIC_ID, value = int cluster label.
- `AssociationResult` dataclass: `effect_size: float, pvalue: float, responder_label: int, grade: str` (grade ∈ {`++`,`+`,`neutral`,`-`,`--`}).
- `test_drug_association(clusters: pd.Series, drug_auc: pd.Series) -> AssociationResult`
- `score_responder(auc_by_cluster: dict[int,float], baseline: float) -> dict[int,str]`
- `test_survival_association(clusters, time: pd.Series, event: pd.Series) -> dict` — `{logrank_p, km: {label: (timeline, survival)}}`
- `differential_cpgs(meth, clusters, target_label) -> pd.DataFrame` — cols `[cpg, delta_beta, pvalue]`, sorted by |delta_beta|.
- `annotate_cpgs(cpgs: list[str]) -> pd.DataFrame` — cols `[cpg, gene, region]`.

---

## Phase 0 — The MGMT gate (go/no-go)

### Task 0: Project scaffold

**Files:**
- Create: `pyproject.toml`, `strata/__init__.py`, `strata/config.py`, `tests/conftest.py`
- Modify: `.gitignore`

- [ ] **Step 1: Init uv project and add deps**

```bash
cd ~/Desktop/Hack
uv init --python 3.12 --no-readme --name methylgkb-strata
uv add pandas numpy scikit-learn scipy statsmodels lifelines umap-learn matplotlib plotly anthropic requests
uv add --dev pytest
```

- [ ] **Step 2: Add data cache to .gitignore**

Append to `.gitignore`:
```
# Strata data cache + outputs
data/
strata/**/__pycache__/
.venv/
outputs/
```

- [ ] **Step 3: Write `strata/config.py`**

```python
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
OUTPUT_DIR = ROOT / "outputs"
DATA_DIR.mkdir(exist_ok=True)
OUTPUT_DIR.mkdir(exist_ok=True)

N_TOP_FEATURES = 2000
RANDOM_STATE = 0
```

- [ ] **Step 4: Empty package + conftest**

Create `strata/__init__.py` (empty) and `tests/conftest.py`:
```python
import numpy as np
import pandas as pd
import pytest

@pytest.fixture
def synthetic_meth():
    """50 samples x 100 CpGs; two latent groups differing at 10 CpGs."""
    rng = np.random.default_rng(0)
    n, p = 50, 100
    base = rng.uniform(0.1, 0.9, size=(n, p))
    group = np.array([0] * 25 + [1] * 25)
    base[group == 1, :10] += 0.4  # group-1 hypermethylated at first 10 CpGs
    base = np.clip(base, 0, 1)
    samples = [f"S{i}" for i in range(n)]
    cpgs = [f"cg{i:03d}" for i in range(p)]
    return pd.DataFrame(base, index=samples, columns=cpgs), pd.Series(group, index=samples)
```

- [ ] **Step 5: Verify env**

Run: `uv run python -c "import pandas, sklearn, lifelines, umap, anthropic; print('ok')"`
Expected: `ok`

- [ ] **Step 6: Commit**

```bash
git add pyproject.toml uv.lock strata/ tests/conftest.py .gitignore
git commit -m "chore: scaffold strata package + deps"
```

---

### Task 1: GDSC data acquisition + schema discovery

GDSC file layouts change; **discover the real schema, don't assume it.** This task produces working loaders and a `notes/GDSC_SCHEMA.md` recording the actual files/columns used.

**Files:**
- Create: `strata/data/__init__.py`, `strata/data/gdsc.py`, `notes/GDSC_SCHEMA.md`, `tests/test_gdsc.py`

- [ ] **Step 1: Download GDSC bulk files**

Fetch into `data/gdsc/` from `https://www.cancerrxgene.org/downloads/bulk_download`:
- Drug response: GDSC2 fitted dose-response (xlsx) — has `COSMIC_ID`, `DRUG_NAME`, `LN_IC50`, `AUC`.
- Methylation: the preprocessed methylation (CpG-island beta) matrix.
- Cell-line details (annotations) xlsx — `COSMIC_ID`, `Sample Name`, `GDSC Tissue descriptor`.

```bash
mkdir -p data/gdsc
# Use the actual current URLs from the bulk_download page; verify each downloads (>1MB) before proceeding.
```

- [ ] **Step 2: Inspect and record the real schema**

Run a throwaway inspection (`uv run python -c "import pandas as pd; print(pd.read_excel('data/gdsc/<file>').columns.tolist())"`) for each file. Write the actual filenames, sheet names, and column names into `notes/GDSC_SCHEMA.md`. **All later code references the names recorded here.**

- [ ] **Step 3: Write loaders in `strata/data/gdsc.py`**

Implement `load_gdsc_methylation()`, `load_gdsc_drug_response()`, `load_gdsc_annotations()` to the contracts in File Structure, reading the files/columns recorded in `notes/GDSC_SCHEMA.md`. Coerce `COSMIC_ID` to `str`, methylation oriented samples-as-rows, drug response in the long format `[COSMIC_ID, drug_name, ln_ic50, auc]`.

- [ ] **Step 4: Write the smoke test `tests/test_gdsc.py`**

```python
from strata.data import gdsc

def test_methylation_shape():
    m = gdsc.load_gdsc_methylation()
    assert m.shape[0] > 500 and m.shape[1] > 1000
    assert m.index.dtype == object               # COSMIC_ID as str
    assert ((m.values >= 0) & (m.values <= 1)).mean() > 0.99  # betas

def test_drug_response_columns():
    d = gdsc.load_gdsc_drug_response()
    assert {"COSMIC_ID", "drug_name", "ln_ic50", "auc"} <= set(d.columns)
    assert (d["drug_name"].str.lower() == "temozolomide").any()
```

- [ ] **Step 5: Run tests**

Run: `uv run pytest tests/test_gdsc.py -v`
Expected: PASS (if `temozolomide` assertion fails, check exact drug-name casing/spelling in `GDSC_SCHEMA.md` and adjust the loader to normalize names).

- [ ] **Step 6: Commit**

```bash
git add strata/data/ notes/GDSC_SCHEMA.md tests/test_gdsc.py
git commit -m "feat: GDSC loaders + schema notes"
```

---

### Task 2: MGMT positive-control verification — THE GATE

Confirm the demo's load-bearing claim is real: an unsupervised methylation split of cell lines separates temozolomide responders, driven by MGMT. **If this fails, stop and pick a fallback methylation→drug pair before building Phase A.**

**Files:**
- Create: `notes/MGMT_GATE.md`, `scripts/mgmt_gate.py`

- [ ] **Step 1: Write `scripts/mgmt_gate.py`**

```python
"""MGMT positive-control gate: does MGMT methylation track temozolomide response in GDSC?"""
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
from strata.data import gdsc

meth = gdsc.load_gdsc_methylation()
drug = gdsc.load_gdsc_drug_response()
tmz = (drug[drug["drug_name"].str.lower() == "temozolomide"]
       .set_index("COSMIC_ID")["auc"])

# Locate MGMT feature(s) — column names are island/gene-tagged; record the match in notes.
mgmt_cols = [c for c in meth.columns if "MGMT" in c.upper()]
print("MGMT feature columns:", mgmt_cols)
assert mgmt_cols, "No MGMT methylation feature found — see fallback in MGMT_GATE.md"

mgmt = meth[mgmt_cols].mean(axis=1)
common = mgmt.index.intersection(tmz.index)
mgmt, tmz_c = mgmt.loc[common], tmz.loc[common]

hi = mgmt > mgmt.median()           # hypermethylated MGMT (silenced) -> expect LOWER AUC (more sensitive)
u, p = mannwhitneyu(tmz_c[hi], tmz_c[~hi], alternative="less")
print(f"n={len(common)}  median AUC hi-MGMT={tmz_c[hi].median():.3f}  "
      f"lo-MGMT={tmz_c[~hi].median():.3f}  Mann-Whitney p={p:.4g}")
```

- [ ] **Step 2: Run the gate**

Run: `uv run python scripts/mgmt_gate.py`
Expected: MGMT columns found; hi-MGMT median AUC **lower** than lo-MGMT; `p < 0.05`.

- [ ] **Step 3: Record the verdict in `notes/MGMT_GATE.md`**

Record: the MGMT feature column name(s), n, both medians, the p-value, and **PASS/FALLBACK**. If FALLBACK: scan for any feature×drug pair with a strong methylation→AUC association (e.g. iterate top-variable CpGs against a few drugs) and record the chosen positive control instead. The rest of the plan refers to "the positive control" — swap MGMT/temozolomide for the fallback pair throughout if needed.

- [ ] **Step 4: Commit**

```bash
git add scripts/mgmt_gate.py notes/MGMT_GATE.md
git commit -m "feat: MGMT positive-control gate + verdict"
```

---

## Phase A — Deterministic engine (the demoable floor)

### Task 3: Feature selection

**Files:** Create `strata/engine/__init__.py`, `strata/engine/features.py`, `tests/test_features.py`

- [ ] **Step 1: Failing test**

```python
from strata.engine.features import select_variable_features

def test_selects_most_variable(synthetic_meth):
    meth, _ = synthetic_meth
    feats = select_variable_features(meth, n_top=10)
    assert len(feats) == 10
    assert set(feats) == {f"cg{i:03d}" for i in range(10)}  # the 10 differential CpGs are most variable
```

- [ ] **Step 2: Run, expect fail** — `uv run pytest tests/test_features.py -v` → FAIL (module missing).

- [ ] **Step 3: Implement**

```python
import pandas as pd

def select_variable_features(meth: pd.DataFrame, n_top: int = 2000) -> list[str]:
    variances = meth.var(axis=0).sort_values(ascending=False)
    return variances.head(n_top).index.tolist()
```

- [ ] **Step 4: Run, expect pass.**

- [ ] **Step 5: Commit** — `git add strata/engine tests/test_features.py && git commit -m "feat: variable-feature selection"`

---

### Task 4: Clustering

**Files:** Create `strata/engine/cluster.py`, `tests/test_cluster.py`

- [ ] **Step 1: Failing test**

```python
from strata.engine.features import select_variable_features
from strata.engine.cluster import cluster_samples
from sklearn.metrics import adjusted_rand_score

def test_recovers_latent_groups(synthetic_meth):
    meth, truth = synthetic_meth
    feats = select_variable_features(meth, n_top=10)
    labels = cluster_samples(meth, feats, k=2, random_state=0)
    assert adjusted_rand_score(truth, labels.loc[truth.index]) > 0.8
```

- [ ] **Step 2: Run, expect fail.**

- [ ] **Step 3: Implement** (auto-pick k by silhouette when `k is None`):

```python
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

def cluster_samples(meth: pd.DataFrame, features: list[str], k: int | None = None,
                    random_state: int = 0) -> pd.Series:
    X = meth[features].values
    if k is None:
        best_k, best_s = 2, -1.0
        for cand in range(2, 7):
            lab = KMeans(n_clusters=cand, random_state=random_state, n_init=10).fit_predict(X)
            s = silhouette_score(X, lab)
            if s > best_s:
                best_k, best_s = cand, s
        k = best_k
    labels = KMeans(n_clusters=k, random_state=random_state, n_init=10).fit_predict(X)
    return pd.Series(labels, index=meth.index, name="cluster")
```

- [ ] **Step 4: Run, expect pass.**

- [ ] **Step 5: Commit** — `git commit -am "feat: methylation clustering with auto-k"`

---

### Task 5: Response association + responder grading

**Files:** Create `strata/engine/associate.py`, `tests/test_associate.py`

- [ ] **Step 1: Failing test**

```python
import pandas as pd
from strata.engine.associate import test_drug_association, AssociationResult

def test_finds_responder_cluster():
    # cluster 0 sensitive (low AUC), cluster 1 resistant (high AUC)
    clusters = pd.Series([0,0,0,1,1,1], index=list("abcdef"))
    auc = pd.Series([0.1,0.15,0.12, 0.8,0.85,0.82], index=list("abcdef"))
    res = test_drug_association(clusters, auc)
    assert isinstance(res, AssociationResult)
    assert res.responder_label == 0
    assert res.pvalue < 0.05
    assert res.grade in {"++", "+"}
```

- [ ] **Step 2: Run, expect fail.**

- [ ] **Step 3: Implement**

```python
from dataclasses import dataclass
import numpy as np, pandas as pd
from scipy.stats import kruskal

@dataclass
class AssociationResult:
    effect_size: float
    pvalue: float
    responder_label: int
    grade: str

def score_responder(auc_by_cluster: dict[int, float], baseline: float) -> dict[int, str]:
    grades = {}
    for label, auc in auc_by_cluster.items():
        delta = baseline - auc            # positive delta => more sensitive than baseline
        if delta > 0.20:   grades[label] = "++"
        elif delta > 0.07: grades[label] = "+"
        elif delta < -0.20: grades[label] = "--"
        elif delta < -0.07: grades[label] = "-"
        else: grades[label] = "neutral"
    return grades

def test_drug_association(clusters: pd.Series, drug_auc: pd.Series) -> AssociationResult:
    common = clusters.index.intersection(drug_auc.index)
    clusters, drug_auc = clusters.loc[common], drug_auc.loc[common]
    groups = [drug_auc[clusters == c].values for c in sorted(clusters.unique())]
    stat, p = kruskal(*groups)
    means = {int(c): float(drug_auc[clusters == c].mean()) for c in clusters.unique()}
    baseline = float(drug_auc.median())
    grades = score_responder(means, baseline)
    responder = min(means, key=means.get)          # lowest AUC = most sensitive
    effect = baseline - means[responder]
    return AssociationResult(effect_size=effect, pvalue=float(p),
                             responder_label=responder, grade=grades[responder])
```

- [ ] **Step 4: Run, expect pass.**

- [ ] **Step 5: Add survival association (for Phase B / TCGA) in the same file**

```python
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

def test_survival_association(clusters: pd.Series, time: pd.Series, event: pd.Series) -> dict:
    common = clusters.index.intersection(time.index).intersection(event.index)
    clusters, time, event = clusters.loc[common], time.loc[common], event.loc[common]
    lr = multivariate_logrank_test(time, clusters, event)
    km = {}
    for c in sorted(clusters.unique()):
        m = clusters == c
        kmf = KaplanMeierFitter().fit(time[m], event[m])
        km[int(c)] = (kmf.timeline.tolist(), kmf.survival_function_.iloc[:, 0].tolist())
    return {"logrank_p": float(lr.p_value), "km": km}
```

- [ ] **Step 6: Commit** — `git commit -am "feat: drug-response association, responder grading, survival test"`

---

### Task 6: Driver-CpG identification + gene annotation

**Files:** Create `strata/engine/annotate.py`, `tests/test_annotate.py`

- [ ] **Step 1: Failing test**

```python
import pandas as pd
from strata.engine.annotate import differential_cpgs

def test_finds_driver_cpgs(synthetic_meth):
    meth, truth = synthetic_meth
    diff = differential_cpgs(meth, truth, target_label=1)
    top10 = set(diff.head(10)["cpg"])
    assert top10 == {f"cg{i:03d}" for i in range(10)}   # the planted differential CpGs
```

- [ ] **Step 2: Run, expect fail.**

- [ ] **Step 3: Implement**

```python
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu

def differential_cpgs(meth: pd.DataFrame, clusters: pd.Series, target_label: int) -> pd.DataFrame:
    common = meth.index.intersection(clusters.index)
    meth, clusters = meth.loc[common], clusters.loc[common]
    in_grp = clusters == target_label
    rows = []
    for cpg in meth.columns:
        a, b = meth.loc[in_grp, cpg], meth.loc[~in_grp, cpg]
        try:
            _, p = mannwhitneyu(a, b)
        except ValueError:
            p = 1.0
        rows.append((cpg, float(a.mean() - b.mean()), float(p)))
    return (pd.DataFrame(rows, columns=["cpg", "delta_beta", "pvalue"])
            .reindex(columns=["cpg", "delta_beta", "pvalue"])
            .sort_values("delta_beta", key=lambda s: s.abs(), ascending=False)
            .reset_index(drop=True))

def annotate_cpgs(cpgs: list[str]) -> pd.DataFrame:
    """Map CpG-island feature ids to gene symbols.
    GDSC island feature ids embed the gene symbol — parse it; record the exact format in GDSC_SCHEMA.md.
    """
    rows = [(c, _gene_from_island_id(c), "promoter/island") for c in cpgs]
    return pd.DataFrame(rows, columns=["cpg", "gene", "region"])

def _gene_from_island_id(island_id: str) -> str:
    # GDSC island ids are gene-tagged (e.g. "MGMT_chr10_..."). Adjust to the real format from Task 1.
    return island_id.split("_")[0]
```

- [ ] **Step 4: Run, expect pass** (the synthetic `cg###` ids have no gene; `differential_cpgs` is what the test checks — `annotate_cpgs` is exercised on real GDSC ids in Task 8).

- [ ] **Step 5: Commit** — `git commit -am "feat: differential CpG detection + gene annotation"`

---

### Task 7: Figures

**Files:** Create `strata/figures.py`

- [ ] **Step 1: Implement the four plotters**

```python
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd, umap
from strata.config import OUTPUT_DIR

def plot_umap(meth, features, clusters, response=None, fname="umap.png"):
    emb = umap.UMAP(random_state=0).fit_transform(meth[features].values)
    fig, ax = plt.subplots(figsize=(5, 4))
    c = response.loc[meth.index] if response is not None else clusters.loc[meth.index]
    sc = ax.scatter(emb[:, 0], emb[:, 1], c=c.values, cmap="viridis", s=12)
    fig.colorbar(sc); ax.set_title("Methylation UMAP"); fig.tight_layout()
    fig.savefig(OUTPUT_DIR / fname, dpi=150); plt.close(fig); return OUTPUT_DIR / fname

def plot_auc_box(clusters, auc, fname="auc_box.png"):
    df = pd.DataFrame({"cluster": clusters, "auc": auc}).dropna()
    fig, ax = plt.subplots(figsize=(5, 4))
    df.boxplot(column="auc", by="cluster", ax=ax)
    ax.set_title("Drug AUC by subgroup"); fig.suptitle(""); fig.tight_layout()
    fig.savefig(OUTPUT_DIR / fname, dpi=150); plt.close(fig); return OUTPUT_DIR / fname

def plot_driver_heatmap(meth, clusters, driver_cpgs, fname="heatmap.png"):
    order = clusters.sort_values().index
    sub = meth.loc[order, driver_cpgs]
    fig, ax = plt.subplots(figsize=(6, 4))
    im = ax.imshow(sub.values, aspect="auto", cmap="RdBu_r", vmin=0, vmax=1)
    fig.colorbar(im); ax.set_title("Driver CpG methylation"); fig.tight_layout()
    fig.savefig(OUTPUT_DIR / fname, dpi=150); plt.close(fig); return OUTPUT_DIR / fname

def plot_km(survival_result, fname="km.png"):
    fig, ax = plt.subplots(figsize=(5, 4))
    for label, (t, s) in survival_result["km"].items():
        ax.step(t, s, where="post", label=f"subgroup {label}")
    ax.set_xlabel("time"); ax.set_ylabel("survival"); ax.legend()
    ax.set_title(f"KM (log-rank p={survival_result['logrank_p']:.3g})"); fig.tight_layout()
    fig.savefig(OUTPUT_DIR / fname, dpi=150); plt.close(fig); return OUTPUT_DIR / fname
```

- [ ] **Step 2: Smoke-check** — `uv run python -c "from strata import figures; print('import ok')"`

- [ ] **Step 3: Commit** — `git add strata/figures.py && git commit -m "feat: demo figures (umap, auc box, heatmap, km)"`

---

### Task 8: Deterministic CLI — rung A end-to-end

This produces the full result with no LLM: MGMT recovery + one novel responder subgroup + figures + a templated brief. **After this task a complete demo exists.**

**Files:** Create `strata/agent/__init__.py`, `strata/agent/brief.py`, `strata/cli.py`

- [ ] **Step 1: Templated brief in `strata/agent/brief.py`**

```python
def write_brief(drug, assoc, driver_genes, n_samples) -> str:
    return f"""# Commercial-Development Opportunity — {drug}

**Discovered subgroup:** cluster {assoc.responder_label} ({assoc.grade} responder),
n={n_samples}, AUC effect={assoc.effect_size:.3f}, p={assoc.pvalue:.3g}.

**Methylation signature:** driven by hypermethylation at {", ".join(driver_genes[:5])}.

**Why it's invisible to current taxonomy:** the stratum is defined by a methylation
state, not by histology or a single mutation — it cuts across the indication.

**Commercial implication:** a methylation companion-marker could define a responder
indication for {drug}, raising trial enrichment and peak-year revenue while reducing
opportunity-identification time and cost.
"""
```

- [ ] **Step 2: CLI in `strata/cli.py`**

```python
"""Rung-A deterministic discovery pipeline on GDSC."""
import argparse, pandas as pd
from strata.data import gdsc
from strata.engine.features import select_variable_features
from strata.engine.cluster import cluster_samples
from strata.engine.associate import test_drug_association
from strata.engine.annotate import differential_cpgs, annotate_cpgs
from strata import figures
from strata.agent.brief import write_brief
from strata.config import N_TOP_FEATURES, OUTPUT_DIR

def run(drug_name: str):
    meth = gdsc.load_gdsc_methylation()
    drug = gdsc.load_gdsc_drug_response()
    auc = drug[drug["drug_name"].str.lower() == drug_name.lower()].set_index("COSMIC_ID")["auc"]
    common = meth.index.intersection(auc.index)
    meth, auc = meth.loc[common], auc.loc[common]

    feats = select_variable_features(meth, N_TOP_FEATURES)
    clusters = cluster_samples(meth, feats, k=None)
    assoc = test_drug_association(clusters, auc)
    drivers = differential_cpgs(meth, clusters, assoc.responder_label).head(20)
    genes = annotate_cpgs(drivers["cpg"].tolist())["gene"].unique().tolist()

    figures.plot_umap(meth, feats, clusters, response=auc, fname=f"{drug_name}_umap.png")
    figures.plot_auc_box(clusters, auc, fname=f"{drug_name}_auc.png")
    figures.plot_driver_heatmap(meth, clusters, drivers["cpg"].tolist(), fname=f"{drug_name}_heatmap.png")

    brief = write_brief(drug_name, assoc, genes, len(common))
    (OUTPUT_DIR / f"{drug_name}_brief.md").write_text(brief)
    print(brief)
    return assoc

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--drug", default="Temozolomide")
    run(ap.parse_args().drug)
```

- [ ] **Step 3: Run the positive control** — `uv run python -m strata.cli --drug Temozolomide`
Expected: responder cluster identified with `p < 0.05`; MGMT among driver genes; figures + brief written to `outputs/`.

- [ ] **Step 4: Scan for a novel hit** — run across several other GDSC drugs; pick one with a significant, biologically interpretable responder subgroup (record the choice). This is the "novel" subgroup for the demo.

```bash
for d in Trametinib Olaparib Venetoclax Dasatinib Nutlin-3a; do uv run python -m strata.cli --drug "$d"; done
```

- [ ] **Step 5: Commit** — `git add strata/cli.py strata/agent/brief.py && git commit -m "feat: rung-A deterministic discovery pipeline (MGMT recovery + novel hit)"`

---

## Phase B — Agent + patient vignette

### Task 9: TCGA-LAML patient vignette

**Files:** Create `strata/data/tcga_laml.py`, `tests/test_tcga.py`

- [ ] **Step 1: Download TCGA-LAML methylation + clinical**

From the GDC portal (project `TCGA-LAML`): harmonized 450K beta matrix + clinical (overall survival time, vital status). Save under `data/tcga_laml/`. Record exact files in `notes/GDSC_SCHEMA.md` (append a TCGA section).

- [ ] **Step 2: Loader + smoke test**

```python
# strata/data/tcga_laml.py
import pandas as pd
from strata.config import DATA_DIR

def load_tcga_methylation() -> pd.DataFrame:
    # rows = patient barcode, cols = probe id, values = beta
    ...

def load_tcga_clinical() -> pd.DataFrame:
    # index = patient barcode; cols = ['os_time', 'os_event'] (event=1 dead, 0 censored)
    ...
```

```python
# tests/test_tcga.py
from strata.data import tcga_laml
def test_clinical_has_survival():
    c = tcga_laml.load_tcga_clinical()
    assert {"os_time", "os_event"} <= set(c.columns)
    assert c["os_event"].isin([0, 1]).all()
```

- [ ] **Step 3: Run test, expect pass.**

- [ ] **Step 4: Survival run** — small script: cluster TCGA-LAML methylation → `test_survival_association` → `plot_km`. Confirm a subgroup split with `logrank_p` reported.

- [ ] **Step 5: Commit** — `git commit -am "feat: TCGA-LAML loaders + survival vignette"`

---

### Task 10: Engine functions as Claude tools

**Files:** Create `strata/agent/tools.py`, `tests/test_agent_tools.py`

- [ ] **Step 1: Failing test for dispatch**

```python
from strata.agent.tools import dispatch_tool
import pandas as pd

def test_dispatch_cluster(monkeypatch, synthetic_meth):
    meth, _ = synthetic_meth
    monkeypatch.setattr("strata.agent.tools._load_meth", lambda dataset: meth)
    out = dispatch_tool("cluster", {"dataset": "synthetic", "k": 2})
    assert "n_clusters" in out and out["n_clusters"] == 2
```

- [ ] **Step 2: Run, expect fail.**

- [ ] **Step 3: Implement `tools.py`** — a `TOOLS` list (Anthropic tool schemas for `select_features`, `cluster`, `associate_drug`, `find_drivers`) and a `dispatch_tool(name, args) -> dict` that calls the engine and returns JSON-serializable summaries (cluster sizes, association p/effect/grade, top driver genes). Keep a module-level cache of the loaded matrices keyed by dataset name.

- [ ] **Step 4: Run, expect pass.**

- [ ] **Step 5: Commit** — `git commit -am "feat: engine functions exposed as Claude tools"`

---

### Task 11: The discovery agent loop

**Files:** Create `strata/agent/loop.py`; extend `strata/agent/brief.py` with `write_brief_llm`

- [ ] **Step 1: Implement `run_discovery_agent(dataset, goal)`** — an Anthropic tool-use loop (model `claude-opus-4-8`, prompt caching on the system prompt + tool schemas). System prompt encodes the protocol: select features → cluster → test association → if no subgroup grades `+`/`++` with `p<0.05`, re-cluster with different k/feature count → once found, get drivers → call `write_brief_llm`. The agent decides params and judges signal-vs-noise against the grading criteria.

- [ ] **Step 2: Manual run** — `uv run python -c "from strata.agent.loop import run_discovery_agent; print(run_discovery_agent('gdsc','find a methylation-defined responder subgroup'))"`
Expected: the agent iterates tool calls and returns a brief naming a responder subgroup + drivers. Verify it recovers MGMT/temozolomide when pointed there.

- [ ] **Step 3: Commit** — `git commit -am "feat: agent-driven discovery loop + LLM brief"`

---

### Task 12: Wire the demo output into the site

**Files:** Modify `site/report.html` (or add `site/strata.html`); add a small JSON dump of results

- [ ] **Step 1:** Have the CLI/agent also write `outputs/result.json` (responder subgroup, grade, p, drivers, figure paths, brief markdown).
- [ ] **Step 2:** Extend the existing MethylGKB site with a Strata view that renders `result.json` — the figures + the brief — without removing existing pages.
- [ ] **Step 3:** Manual check in browser; **record the backup demo video at this rung.**
- [ ] **Step 4: Commit** — `git commit -am "feat: surface Strata results in MethylGKB site"`

---

## Phase C — Stretch (only if Friday has room)

### Task 13: Grounded annotation + scale framing

- [ ] **Step 1:** In `annotate_cpgs` (or a new `annotate_grounded`), optionally enrich driver genes by querying the existing PolymerGenomicsAPI / methylation KB for pathway/biophysics context; fall back silently to the parsed gene symbol if unavailable.
- [ ] **Step 2:** Add one slide/section framing scale: clustering + association over millions of methylomes is GPU-shaped (cuML), and the locked controlled-access data is where this runs next — the bridge to the platform horizon (existing `VISION.md`).
- [ ] **Step 3: Commit** — `git commit -am "feat: grounded annotation + scale framing (stretch)"`

---

## Self-Review notes (author check)

- **Spec coverage:** §3 datasets → Tasks 1, 2, 9; §4 run-of-show → Tasks 8 (MGMT + novel), 9 (patient), 11–12 (agent + brief + close); §5 components 1–5 → Tasks 1/9, 3, 4–5, 6, 10–12; §6 ladder (Day-1/A/B/C) → Phase 0 / A / B / C; §9 keep-everything → no deletions anywhere, site extended not replaced.
- **Placeholders:** data-download steps (Tasks 1, 9) intentionally defer exact URLs/columns to discovery-and-record, because GDSC/GDC layouts must be read live — this is a deliberate discovery step with a recorded artifact (`notes/`), not a hand-wave. All engine/agent code is concrete.
- **Type consistency:** `AssociationResult(effect_size, pvalue, responder_label, grade)`, `cluster_samples → pd.Series`, `differential_cpgs → [cpg, delta_beta, pvalue]`, `annotate_cpgs → [cpg, gene, region]` used consistently across Tasks 5/6/8/10/11.
