"""
Polymer — methylation-guided therapy selection (demo UI).

Paste a methylome (gene -> beta) and get a ranked list of drugs with predicted
response. Runs the trained recommender in-process.

    streamlit run polymer_app.py
"""

import numpy as np
import pandas as pd
import streamlit as st

from polymer_infer import Polymer

st.set_page_config(page_title="Polymer", layout="wide")

# ---------------------------------------------------------------- styling ----
st.markdown(
    """
    <style>
      .stApp { background: #0b1020; color: #e7ecf5; }
      .poly-hero { background: linear-gradient(120deg,#5b8cff 0%,#9b5bff 60%,#ff5bd0 100%);
                   padding: 26px 30px; border-radius: 16px; margin-bottom: 8px; }
      .poly-hero h1 { color:#fff; font-size: 40px; margin:0; letter-spacing:.5px; }
      .poly-hero p  { color:#eef2ff; font-size:17px; margin:6px 0 0; opacity:.95; }
      .poly-pill { display:inline-block; background:#141b33; border:1px solid #2a3870;
                   border-radius:10px; padding:10px 14px; margin:6px 8px 0 0; }
      .poly-pill b { color:#9bb4ff; font-size:20px; }
      .poly-pill span { color:#aab3c9; font-size:12px; display:block; }
      .headline { background:#101733; border:1px solid #2a3870; border-radius:14px;
                  padding:18px 22px; margin-top:6px; }
      .headline h2 { margin:0; color:#7CFFB2; font-size:26px; }
      .headline p { margin:4px 0 0; color:#aab3c9; }
    </style>
    """,
    unsafe_allow_html=True,
)


@st.cache_resource
def load_engine():
    return Polymer()


poly = load_engine()
m = poly.metrics

st.markdown(
    """
    <div class="poly-hero">
      <h1>Polymer</h1>
      <p>Methylation-guided therapy selection — paste a methylome, rank every drug.</p>
    </div>
    """,
    unsafe_allow_html=True,
)

c1, c2, c3, c4 = st.columns(4)
c1.markdown(f'<div class="poly-pill"><b>{m.get("r",0):.2f}</b>'
            f'<span>Pearson r (held-out lines)</span></div>', unsafe_allow_html=True)
c2.markdown(f'<div class="poly-pill"><b>{m.get("demeaned_rho",0):.2f}</b>'
            f'<span>personalization &rho; (de-meaned)</span></div>', unsafe_allow_html=True)
c3.markdown(f'<div class="poly-pill"><b>{len(poly.drug_names)}</b>'
            f'<span>drugs scored</span></div>', unsafe_allow_html=True)
c4.markdown(f'<div class="poly-pill"><b>{len(poly.genes):,}</b>'
            f'<span>methylation features</span></div>', unsafe_allow_html=True)

st.write("")
left, right = st.columns([0.42, 0.58], gap="large")

# ------------------------------------------------------------------ input ----
with left:
    st.subheader("1 · Methylome input")
    labels = [e["label"] for e in poly.examples]
    pick = st.selectbox("Load an example patient (held-out, never seen in training)",
                        ["— paste your own below —"] + labels)
    if st.button("Load example into box", use_container_width=True) and pick in labels:
        ex = next(e for e in poly.examples if e["label"] == pick)
        st.session_state["methyl_text"] = "\n".join(f"{g}\t{b}" for g, b in ex["betas"].items())
        st.session_state["true_top"] = ex.get("true_top_sensitive", [])
        st.rerun()

    st.session_state.setdefault("methyl_text", "")
    text = st.text_area(
        "Gene → beta table  (one `GENE<tab or comma>beta` per line; beta in 0–1)",
        key="methyl_text", height=320,
        placeholder="CDKN2A\t0.83\nMGMT\t0.10\nMTAP\t0.47\n...",
    )
    go = st.button("Predict drug response", type="primary", use_container_width=True)

# ------------------------------------------------------------------ output ---
with right:
    st.subheader("2 · Predicted response (ranked)")
    if not go:
        st.info("Load an example patient or paste a methylome, then click "
                "**Predict drug response**.")
    else:
        methylome = poly.parse_table(text)
        if not methylome:
            st.error("Couldn't parse any `gene  beta` rows. Check the format.")
            st.stop()
        rows, matched = poly.recommend(methylome, top_n=20)
        df = pd.DataFrame(rows)
        df["sensitivity"] = (1.0 - df["pred_auc"]).clip(0, 1)

        top = rows[0]
        st.markdown(
            f"""<div class="headline">
              <h2>Top recommendation: {top['drug']}</h2>
              <p>predicted AUC <b>{top['pred_auc']:.3f}</b>
                 &nbsp;·&nbsp; {top['delta_vs_avg']:+.3f} vs the average cell line
                 &nbsp;·&nbsp; target: {top['target'] or 'n/a'}
                 &nbsp;·&nbsp; {top['pathway'] or ''}</p>
            </div>""",
            unsafe_allow_html=True,
        )
        st.caption(f"Matched {matched:,} / {len(poly.genes):,} panel genes. "
                   f"Lower AUC = more sensitive. **vs avg < 0** = more sensitive than "
                   f"the typical line (the personalized signal).")

        show = df[["drug", "pred_auc", "sensitivity", "delta_vs_avg", "target", "pathway"]]
        st.dataframe(
            show, hide_index=True, use_container_width=True, height=430,
            column_config={
                "drug": st.column_config.TextColumn("Drug"),
                "pred_auc": st.column_config.NumberColumn("Pred AUC", format="%.3f"),
                "sensitivity": st.column_config.ProgressColumn(
                    "Predicted sensitivity", min_value=0.0, max_value=1.0, format=" "),
                "delta_vs_avg": st.column_config.NumberColumn("vs avg", format="%+.3f"),
                "target": st.column_config.TextColumn("Target"),
                "pathway": st.column_config.TextColumn("Pathway"),
            },
        )

        if st.session_state.get("true_top"):
            hits = [d for d in [r["drug"] for r in rows[:5]]
                    if d in set(st.session_state["true_top"])]
            st.success(f"Ground-truth check (this held-out line's true most-sensitive "
                       f"drugs): {', '.join(st.session_state['true_top'])}. "
                       f"Polymer's top-5 caught {len(hits)}: {', '.join(hits) or 'none'}.")

with st.expander("How Polymer works  ·  method & honest caveats"):
    st.markdown(
        f"""
**Model.** A single methylation **encoder** was pretrained across **all
{len(poly.drug_names)} drugs** on {m.get('n_train_lines','?')} training cell lines. For a methylome it
computes an embedding, then scores every drug as
`drug_bias + ⟨encoder(methylome), drug_embedding⟩`. The `drug_bias` captures a drug's
average potency; the dot product is the **personalized** interaction.

**Validation (cell-line-blind).** Pearson r **{m.get('r',0):.2f}** vs a per-drug-mean
baseline of {m.get('base_r',0):.2f}. After removing each drug's mean — which strips out
everything a one-size-fits-all model knows — the within-patient ranking still holds at
**ρ = {m.get('demeaned_rho',0):.2f}**. That residual is the real personalization.

**Caveats (we say these on stage).** Trained on cancer **cell lines** (GDSC), not patients —
clinical use needs patient methylomes. Features are **gene-level** betas (promoter-CpG
resolution is averaged out). Drugs are all seen in training (not cold-start).
        """
    )
