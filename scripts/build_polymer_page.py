"""Build site/strata.html as an interactive tutorial from polymer_demo.ipynb.

Reads the notebook, pulls each code cell's source + captured stdout verbatim,
references the figures already extracted into ``site/polymer/`` (PNG per cell),
and renders a static, design-system-consistent tutorial page. Interaction
(scrollspy stage rail, "Run tutorial" walkthrough, code toggles, syntax
highlighting) lives in ``site/polymer.js`` + ``site/polymer.css``.

Run:  python scripts/build_polymer_page.py
"""

from __future__ import annotations

import html
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
NB = ROOT / "polymer_demo.ipynb"
OUT = ROOT / "site" / "strata.html"

# cell index -> (figure filename in site/polymer/, caption)
FIGS = {
    5: ("training_curve.png", "Stage 2 · training — one encoder across all 286 drugs (train vs val loss)"),
    6: ("heldout_fit.png", "Stage 2 · held-out fit: predicted vs true AUC, and the residual distribution"),
    8: ("drug_ranking.png", "Stage 3 · top-15 ranked drugs for the held-out patient (green = better than average)"),
    10: ("genes_of_interest.png", "Stage 4 · the patient's most hyper- and hypo-methylated genes vs cohort"),
    12: ("driver_genes.png", "Stage 5 · driver genes for the chosen drug (green drives sensitivity, red drives resistance)"),
}

# Tutorial structure: which notebook cells belong to each stage + the prose.
STAGES = [
    {
        "n": 1,
        "marker": "§02 — Stage 1 · Load &amp; split",
        "title": "Load the methylome, drug response, and build a patient-blind split.",
        "lede": [
            "Three inputs: <strong>gene-level methylation</strong> (the predictor), "
            "<strong>promoter methylation</strong> (for the hyper/hypo annotation), and "
            "<strong>GDSC2 drug AUCs</strong> (the outcome). We split on <strong>cell lines</strong> — "
            "so every test \u201cpatient\u201d is a sample the model has never seen.",
        ],
        "cells": [1, 3],
        "metrics": None,
    },
    {
        "n": 2,
        "marker": "§03 — Stage 2 · Train the encoder",
        "title": "Pretrain one methylation encoder across all 286 drugs.",
        "lede": [
            "Top-variance genes are selected on <strong>train lines only</strong> and z-scored with "
            "train statistics — no leakage. Each drug is a learned embedding, so "
            "<code>predicted AUC = drug_bias + \u27e8encoder(methylome), drug_embedding\u27e9</code>. "
            "We watch the train/val curves, then check the fit on held-out lines against a "
            "per-drug-mean baseline.",
        ],
        "cells": [5, 6],
        "metrics": [
            ("blue", "r = 0.857", "Held-out Pearson correlation, predicted vs true AUC"),
            ("teal", "\u221220.9%", "RMSE vs the per-drug-mean baseline (0.0746 vs 0.0943)"),
            ("blue", "\u03c1 = +0.489", "Within-patient personalization (de-meaned Spearman)"),
        ],
    },
    {
        "n": 3,
        "marker": "§04 — Stage 3 · Rank every drug",
        "title": "Score all 286 drugs for one held-out patient.",
        "lede": [
            "Pick a held-out cell line as a stand-in <strong>patient</strong>, push its methylome "
            "through the encoder once, and dot it against every drug embedding. Rank by predicted "
            "AUC (lower = more sensitive). Green bars are drugs predicted <em>better for this "
            "patient than average</em> — that is personalization, not just raw potency.",
        ],
        "cells": [8],
        "metrics": None,
    },
    {
        "n": 4,
        "marker": "§05 — Stage 4 · Genes of interest",
        "title": "Which genes are hyper- or hypo-methylated in this patient?",
        "lede": [
            "The patient's z-score for each panel gene <em>is</em> its standardized input, so it "
            "directly says how <strong>hyper-</strong> (z&nbsp;&gt;&nbsp;0) or <strong>hypo-</strong> "
            "(z&nbsp;&lt;&nbsp;0) methylated that gene is versus the population. We line up the "
            "promoter beta for the same genes as the focal view.",
        ],
        "cells": [10],
        "metrics": None,
    },
    {
        "n": 5,
        "marker": "§06 — Stage 5 · Why this drug?",
        "title": "Driver genes = gradient attribution × methylation deviation.",
        "lede": [
            "For the top pick we take the <strong>gradient</strong> of its predicted AUC with "
            "respect to the input methylation and multiply by the patient's deviation. That product "
            "is each gene's first-order contribution: <strong>negative</strong> means its "
            "(hyper/hypo)methylation <strong>drives sensitivity</strong>, positive drives "
            "resistance — an auditable, gene-level rationale tied directly to the prediction.",
        ],
        "cells": [12],
        "metrics": None,
    },
]


def cell_source(cell) -> str:
    return "".join(cell["source"]).rstrip("\n")


def cell_stdout(cell) -> str:
    chunks = []
    for o in cell.get("outputs", []):
        ot = o.get("output_type")
        if ot == "stream":
            chunks.append("".join(o.get("text", [])))
        elif ot in ("execute_result", "display_data"):
            data = o.get("data", {})
            txt = data.get("text/plain")
            if txt and "Figure" not in "".join(txt):
                chunks.append("".join(txt))
    return "".join(chunks).rstrip("\n")


def render_cell(idx: int, cell, exec_count) -> str:
    code = html.escape(cell_source(cell))
    stdout = cell_stdout(cell)
    fig = FIGS.get(idx)
    parts = [
        '<div class="nb-cell">',
        '  <div class="nb-toolbar">',
        '    <span class="nb-lang">python</span>',
        f'    <span class="nb-in">In [{exec_count}]</span>',
        '    <span class="nb-actions">',
        '      <button class="nb-btn" data-act="copy" type="button">copy</button>',
        '      <button class="nb-btn" data-act="toggle" type="button">hide code</button>',
        '    </span>',
        '  </div>',
        f'  <pre class="nb-code"><code>{code}</code></pre>',
    ]
    if stdout or fig:
        parts.append('  <div class="nb-out">')
        if stdout:
            parts.append('    <div class="nb-out-label">stdout</div>')
            parts.append(f'    <pre class="nb-stdout">{html.escape(stdout)}</pre>')
        if fig:
            fn, cap = fig
            parts.append('    <figure class="nb-fig">')
            parts.append(f'      <img src="polymer/{fn}" alt="{html.escape(cap)}" loading="lazy" />')
            parts.append(f'      <figcaption>{cap}</figcaption>')
            parts.append('    </figure>')
        parts.append('  </div>')
    parts.append('</div>')
    return "\n".join("      " + p for p in parts)


def render_metrics(metrics) -> str:
    cells = []
    for cls, val, label in metrics:
        cells.append(
            '        <div class="poly-metric">\n'
            f'          <div class="pm-val {cls}">{val}</div>\n'
            f'          <div class="pm-l">{label}</div>\n'
            '        </div>'
        )
    return '      <div class="poly-metrics">\n' + "\n".join(cells) + "\n      </div>"


def main() -> None:
    nb = json.loads(NB.read_text())
    cells = nb["cells"]

    # stage rail chips
    chip_html = []
    chip_labels = ["Load", "Train", "Rank", "Genes", "Why"]
    for i, lbl in enumerate(chip_labels, start=1):
        chip_html.append(
            f'    <button class="sr-chip" data-target="stage{i}" type="button">'
            f'<span class="sr-n"><span class="sr-n-txt">{i}</span></span>'
            f'<span class="sr-label">{lbl}</span></button>'
        )
    chips = "\n".join(chip_html)

    # stages
    stage_html = []
    for st in STAGES:
        body = []
        body.append('<section class="wrap section stage" id="stage{}">'.format(st["n"]))
        body.append(f'  <p class="section-marker">{st["marker"]}</p>')
        body.append(f'  <h2 class="h2" style="max-width:30ch;">{st["title"]}</h2>')
        for p in st["lede"]:
            body.append(f'  <p class="stage-lede">{p}</p>')
        if st["metrics"]:
            body.append(render_metrics(st["metrics"]))
        for ci in st["cells"]:
            exec_count = cells[ci].get("execution_count") or ""
            body.append(render_cell(ci, cells[ci], exec_count))
        body.append('</section>')
        stage_html.append("\n".join(body))
    stages_joined = "\n\n".join(stage_html)

    page = TEMPLATE.replace("{{CHIPS}}", chips).replace("{{STAGES}}", stages_joined)
    OUT.write_text(page)
    print(f"wrote {OUT}  ({len(page):,} bytes)")


TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width, initial-scale=1.0" />
<title>Strata — methylation → drug response, with the genes that drive it</title>
<meta name="description" content="An interactive walkthrough of Polymer: given a methylome, predict the response to every drug and explain why — the hyper/hypo-methylated genes that drive the prediction. Rendered from the polymer_demo notebook." />
<link rel="preconnect" href="https://fonts.googleapis.com" />
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin />
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500;700&display=swap" rel="stylesheet" />
<link rel="stylesheet" href="styles.css" />
<link rel="stylesheet" href="polymer.css" />
</head>
<body>

<header class="site-head">
  <a class="brand" href="index.html">METHYLGKB</a>
  <span class="brand-divider"></span>
  <span class="brand-sub">by Polymer Genomics</span>
  <nav class="site-nav">
    <a class="nav-link" href="index.html">Overview</a>
    <a class="nav-link" href="science.html">The Science</a>
    <a class="nav-link" href="explore.html">Explore</a>
    <a class="nav-link is-active" href="strata.html">Strata</a>
    <a class="nav-link" href="strata-v2.html">Strata v2</a>
  </nav>
</header>

<main>

  <!-- HERO -->
  <section class="wrap hero" style="padding-bottom:36px;">
    <p class="section-marker">§01 — Polymer · the per-patient engine</p>
    <h1 style="font-size:clamp(34px,5vw,52px);max-width:22ch;">
      Predict every drug for a methylome — <span class="em">and the genes that drive it.</span>
    </h1>
    <p class="lede">
      <strong>The product, in one notebook.</strong> Given a tumor or cell-line
      <strong>methylome</strong>, Polymer predicts the <strong>response (AUC) to every drug</strong>
      and ranks them — then explains <em>why</em>: the genes that are hyper- or hypo-methylated in
      this sample and push the predicted response up or down. <strong>Research-use today</strong>;
      a methodology demo on open GDSC data, not a clinical device.
    </p>

    <!-- THE PIPELINE (the loop, named explicitly) -->
    <p class="poly-pipe">
      <span class="pp-label">the pipeline</span>
      methylome
      <span class="pp-arrow">→</span> <code>encoder</code> (pretrained across all drugs)
      <span class="pp-arrow">→</span> rank 286 drugs
      <span class="pp-arrow">→</span> hyper/hypo genes
      <span class="pp-arrow">→</span> driver genes (why)
    </p>

    <div class="poly-runctl">
      <button class="btn btn-primary" id="runTutorial" type="button">▶ Run tutorial</button>
      <span class="rc-hint">Five stages, played top to bottom. Or scroll — every cell shows real code + output.</span>
    </div>

    <div class="poly-note" style="border-left-color:var(--primary);">
      <strong>Honest framing.</strong> Trained on GDSC <strong>cell lines</strong>, not patients;
      gene-level betas average over promoter CpGs; attribution is first-order (local). The
      <strong>method</strong>, not the clinical claim, is the product.
    </div>
  </section>

  <!-- STICKY STAGE RAIL -->
  <div class="stage-rail">
    <span class="sr-progress"></span>
    <div class="stage-rail-inner" role="tablist" aria-label="Tutorial stages">
{{CHIPS}}
    </div>
  </div>

{{STAGES}}

  <!-- CLOSING: demo script -->
  <section class="wrap section" id="closing">
    <p class="section-marker">§07 — Demo script · what to say</p>
    <h2 class="h2">The four-beat walkthrough.</h2>
    <div class="poly-script">
      <div class="ps-item"><span class="ps-n">1</span><p class="ps-t"><strong>Stage 2</strong> — &ldquo;We pretrain one methylation encoder across all 286 drugs. Here&rsquo;s the training curve; it generalizes to held-out patients, beating the per-drug-mean baseline.&rdquo;</p></div>
      <div class="ps-item"><span class="ps-n">2</span><p class="ps-t"><strong>Stage 3</strong> — &ldquo;Feed a new patient&rsquo;s methylome → Polymer ranks every drug. Green bars are drugs it thinks are <em>better for this patient than average</em> — personalization, not just potency.&rdquo;</p></div>
      <div class="ps-item"><span class="ps-n">3</span><p class="ps-t"><strong>Stage 4</strong> — &ldquo;These are the genes hyper/hypo-methylated in this patient versus the population.&rdquo;</p></div>
      <div class="ps-item"><span class="ps-n">4</span><p class="ps-t"><strong>Stage 5</strong> — &ldquo;And here&rsquo;s <em>why</em> the top drug was chosen: the specific (hyper/hypo)methylated genes driving the predicted sensitivity — an auditable, gene-level rationale.&rdquo;</p></div>
    </div>

    <div class="poly-note">
      <strong>Caveats to own on stage:</strong> GDSC cell lines (not patients); gene-level betas;
      attribution is first-order (local). <strong>Re-run knobs:</strong> set
      <code>PATIENT = cand[k]</code> (Stage&nbsp;3) for another patient; <code>TARGET_DRUG = "..."</code>
      (Stage&nbsp;5) to explain any drug; raise <code>EPOCHS</code> for a tighter fit.
    </div>

    <div class="hero-actions" style="margin-top:28px;">
      <a class="btn btn-secondary" href="science.html">The evidence</a>
      <a class="btn btn-secondary" href="explore.html">Explore the data</a>
    </div>
  </section>

</main>

<footer class="site-foot">
  <div class="foot-links">
    <a href="index.html">Overview</a>
    <span class="foot-dot">·</span>
    <a href="science.html">The Science</a>
    <span class="foot-dot">·</span>
    <a href="explore.html">Explore</a>
    <span class="foot-dot">·</span>
    <a href="strata.html">Strata</a>
  </div>
  <div class="foot-meta">
    <a class="ruo" href="#" title="Not for clinical or diagnostic use. Prototype. Data provided as-is.">Research use only</a>
    <span class="foot-dot">·</span>
    <span class="copy">© 2026 Polymer Genomics</span>
  </div>
</footer>

<script src="polymer.js"></script>
</body>
</html>
"""


if __name__ == "__main__":
    main()
