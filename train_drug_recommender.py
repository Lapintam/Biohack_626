"""
"Given my methylome, which drug would work for me?"

A drug-response recommender that pretrains a single METHYLATION ENCODER across
ALL drugs at once, then ranks every drug for a cell line by predicted response.

Architecture (content-based matrix factorization):

    cell_emb = Encoder(methylation)            # shared across every drug
    pred_AUC(line, drug) = drug_bias[drug] + <cell_emb(line), drug_emb[drug]>

The `drug_bias` term absorbs each drug's average potency (the part that is the
same for everyone); the dot product is the PERSONALIZED interaction -- the only
part that can say a drug is differentially good *for this methylome*.

Evaluation is cell-line-blind (held-out lines) and reports three things:
  1. Global RMSE / Pearson r vs the per-drug-mean baseline.
  2. PERSONALIZATION: within-line Spearman between predicted and true AUC across
     a line's drugs, both raw and AFTER removing each drug's mean (the honest
     test -- a per-drug-mean model scores exactly 0 on the de-meaned version).
  3. Top-k recovery: of the model's top-5 recommended (most-sensitive) drugs for
     a held-out line, how many fall in that line's true top-10 -- vs the
     everyone-gets-the-same-list baseline.

Saves the pretrained encoder to `outputs/methylation_encoder.pt`.

Usage:  python train_drug_recommender.py
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

import torch
import torch.nn as nn

from strata.data import gdsc
from strata.config import OUTPUT_DIR

SEED = 1337
N_GENES = 3000
EMB = 32
HID = 256
EPOCHS = 40
BATCH = 4096
LR = 2e-3
VAL_FRAC, TEST_FRAC = 0.15, 0.15
TOPK_REC, TOPK_TRUE = 5, 10


def main():
    torch.manual_seed(SEED)
    rng = np.random.default_rng(SEED)

    # ---------- data ----------
    meth = gdsc.load_gdsc_methylation()
    drug = gdsc.load_gdsc_drug_response()
    anno = gdsc.load_gdsc_annotations()

    lines = meth.index.to_numpy()
    line_pos = {c: i for i, c in enumerate(lines)}

    d = drug[drug["COSMIC_ID"].isin(set(lines))].dropna(subset=["auc"]).copy()
    drug_names = sorted(d["drug_name"].unique())
    drug_pos = {dn: i for i, dn in enumerate(drug_names)}
    n_lines, n_drugs = len(lines), len(drug_names)

    li = d["COSMIC_ID"].map(line_pos).to_numpy()
    di = d["drug_name"].map(drug_pos).to_numpy()
    y = d["auc"].to_numpy(dtype=np.float64)
    print(f"{len(y):,} (line,drug) pairs | {n_lines} lines | {n_drugs} drugs")

    # ---------- cell-line-blind split (on lines) ----------
    order = rng.permutation(n_lines)
    n_te = int(n_lines * TEST_FRAC); n_va = int(n_lines * VAL_FRAC)
    test_lines = set(order[:n_te]); val_lines = set(order[n_te:n_te + n_va])
    split = np.where(np.isin(li, list(test_lines)), 2,
                     np.where(np.isin(li, list(val_lines)), 1, 0))
    tr, va, te = (split == 0), (split == 1), (split == 2)
    train_line_mask = ~(np.isin(np.arange(n_lines), list(test_lines | val_lines)))
    print(f"lines  train={train_line_mask.sum()} val={len(val_lines)} test={len(test_lines)}")

    # ---------- features: top-variance genes on TRAIN lines, z-scored on train ----------
    Mfull = meth.to_numpy(dtype=np.float32)
    var = Mfull[train_line_mask].var(0)
    panel = np.argsort(var)[::-1][:N_GENES]
    M = Mfull[:, panel]
    mu = M[train_line_mask].mean(0); sd = M[train_line_mask].std(0); sd[sd == 0] = 1.0
    Mz = torch.tensor((M - mu) / sd)                       # (n_lines, N_GENES)

    # standardize target on train pairs
    ymu, ysd = y[tr].mean(), y[tr].std() + 1e-8
    yz = torch.tensor((y - ymu) / ysd, dtype=torch.float32)
    li_t = torch.tensor(li); di_t = torch.tensor(di)

    # ---------- baseline: per-drug mean (fit on train) ----------
    drug_mean = pd.Series(y[tr]).groupby(pd.Series(di[tr])).mean()
    glob = y[tr].mean()
    base_pred = np.array([drug_mean.get(dd, glob) for dd in di])

    def rmse_r(mask, pred):
        yt = y[mask]
        rmse = float(np.sqrt(np.mean((yt - pred[mask]) ** 2)))
        r = float(pearsonr(yt, pred[mask])[0])
        return rmse, r

    b_rmse, b_r = rmse_r(te, base_pred)

    # ---------- model ----------
    class Recommender(nn.Module):
        def __init__(self):
            super().__init__()
            self.encoder = nn.Sequential(
                nn.Linear(N_GENES, HID), nn.ReLU(), nn.BatchNorm1d(HID), nn.Dropout(0.3),
                nn.Linear(HID, EMB),
            )
            self.drug_emb = nn.Embedding(n_drugs, EMB)
            self.drug_bias = nn.Embedding(n_drugs, 1)
            nn.init.normal_(self.drug_emb.weight, std=0.1)
            nn.init.zeros_(self.drug_bias.weight)

        def encode(self, line_idx):
            return self.encoder(Mz[line_idx])

        def forward(self, line_idx, drug_idx):
            ce = self.encode(line_idx)
            de = self.drug_emb(drug_idx)
            return (self.drug_bias(drug_idx).squeeze(1) + (ce * de).sum(1))

    net = Recommender()
    opt = torch.optim.Adam(net.parameters(), lr=LR, weight_decay=1e-5)
    loss_fn = nn.MSELoss()

    tr_idx = np.where(tr)[0]

    @torch.no_grad()
    def predict(mask_idx):
        net.eval()
        out = np.zeros(len(mask_idx), dtype=np.float64)
        for s in range(0, len(mask_idx), 16384):
            b = mask_idx[s:s + 16384]
            p = net(li_t[b], di_t[b]).numpy()
            out[s:s + len(b)] = p
        return out * ysd + ymu

    va_idx = np.where(va)[0]; te_idx = np.where(te)[0]
    best_val, best_state = float("inf"), None
    for ep in range(1, EPOCHS + 1):
        net.train()
        perm = tr_idx[rng.permutation(len(tr_idx))]
        for s in range(0, len(perm), BATCH):
            b = perm[s:s + BATCH]
            opt.zero_grad()
            pred = net(li_t[b], di_t[b])
            loss = loss_fn(pred, yz[b])
            loss.backward(); opt.step()
        vp = predict(va_idx); v_rmse = float(np.sqrt(np.mean((y[va] - vp) ** 2)))
        if v_rmse < best_val - 1e-5:
            best_val = v_rmse; best_state = {k: v.clone() for k, v in net.state_dict().items()}
        if ep % 5 == 0 or ep == 1:
            print(f"epoch {ep:2d}  val RMSE={v_rmse:.4f}  (baseline {b_rmse:.4f})")
    if best_state: net.load_state_dict(best_state)

    # ---------- 1. global metrics ----------
    te_pred_full = np.full(len(y), np.nan); te_pred_full[te_idx] = predict(te_idx)
    m_rmse, m_r = rmse_r(te, te_pred_full)
    print("\n================ TEST (held-out cell lines) ================")
    print(f"per-drug-mean baseline   RMSE={b_rmse:.4f}  r={b_r:.3f}")
    print(f"recommender (this model) RMSE={m_rmse:.4f}  r={m_r:.3f}  "
          f"({(b_rmse-m_rmse)/b_rmse*100:+.1f}% RMSE)")

    # ---------- 2. personalization: within-line ranking ----------
    dmean_arr = np.array([drug_mean.get(dd, glob) for dd in di])
    raw_rhos, dem_rhos = [], []
    for L in test_lines:
        sel = te_idx[li[te_idx] == L]
        if len(sel) < 10:
            continue
        yt, yp = y[sel], te_pred_full[sel]
        raw_rhos.append(spearmanr(yp, yt).correlation)
        # de-mean each drug -> personalized residual ranking
        dm = dmean_arr[sel]
        dem_rhos.append(spearmanr(yp - dm, yt - dm).correlation)
    print(f"\npersonalization (within-line Spearman across a line's drugs):")
    print(f"  raw         mean rho = {np.nanmean(raw_rhos):+.3f}  (per-drug-mean model ~ tissue-level)")
    print(f"  de-meaned   mean rho = {np.nanmean(dem_rhos):+.3f}  "
          f"(HONEST personalization; per-drug-mean baseline = 0 by construction)")

    # ---------- 3. top-k recovery ----------
    hits_model, hits_base = [], []
    for L in test_lines:
        sel = te_idx[li[te_idx] == L]
        if len(sel) < TOPK_TRUE + 1:
            continue
        drugs_here = di[sel]
        true_rank = drugs_here[np.argsort(y[sel])]               # most-sensitive first
        true_top = set(true_rank[:TOPK_TRUE].tolist())
        model_top = set(drugs_here[np.argsort(te_pred_full[sel])][:TOPK_REC].tolist())
        base_top = set(drugs_here[np.argsort(dmean_arr[sel])][:TOPK_REC].tolist())
        hits_model.append(len(model_top & true_top))
        hits_base.append(len(base_top & true_top))
    print(f"\ntop-{TOPK_REC} recommendation recovery (hits in true top-{TOPK_TRUE} most-sensitive):")
    print(f"  recommender        = {np.mean(hits_model):.2f} / {TOPK_REC}")
    print(f"  same-list baseline = {np.mean(hits_base):.2f} / {TOPK_REC}   (per-drug-mean)")

    # ---------- demo: what drug for a held-out line? ----------
    # pick the held-out line with the most measured drugs (skip lines with no pairs)
    cand = [(L, int((li[te_idx] == L).sum())) for L in test_lines]
    demoL = max(cand, key=lambda t: t[1])[0]
    sel = te_idx[li[te_idx] == demoL]
    cid = lines[demoL]
    ct = anno.reindex([cid])["cancer_type"].iloc[0] if cid in anno.index else "?"
    order_p = np.argsort(te_pred_full[sel])[:TOPK_REC]
    print(f"\n--- DEMO: recommendations for held-out line {cid} ({ct}) ---")
    print(f"{'drug':<22}{'pred AUC':>10}{'true AUC':>10}")
    for k in order_p:
        print(f"{drug_names[di[sel][k]]:<22}{te_pred_full[sel][k]:>10.3f}{y[sel][k]:>10.3f}")

    # ---------- save pretrained encoder (encoder only; legacy) ----------
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    path = os.path.join(OUTPUT_DIR, "methylation_encoder.pt")
    torch.save({"encoder": net.encoder.state_dict(), "panel": panel,
                "mu": mu, "sd": sd, "ymu": ymu, "ysd": ysd,
                "drug_names": drug_names, "n_genes": N_GENES, "emb": EMB}, path)
    print(f"\nsaved pretrained methylation encoder -> {path}")

    # ---------- save FULL Polymer artifact (everything inference needs) ----------
    panel_genes = [str(meth.columns[j]) for j in panel]

    # per-drug average AUC (raw units), keyed by drug name -> the "everyone" prior
    drug_mean_by_name = {drug_names[int(idx)]: float(v) for idx, v in drug_mean.items()}

    # drug target / pathway for display (from the merged table if present)
    drug_meta = {dn: {"target": "", "pathway": ""} for dn in drug_names}
    merged_csv = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gdsc2_merged.csv")
    try:
        mm = pd.read_csv(merged_csv, usecols=["DRUG_NAME", "PUTATIVE_TARGET", "PATHWAY_NAME"])
        first = mm.groupby("DRUG_NAME").first()
        for dn in drug_names:
            if dn in first.index:
                drug_meta[dn] = {"target": str(first.loc[dn, "PUTATIVE_TARGET"]),
                                 "pathway": str(first.loc[dn, "PATHWAY_NAME"])}
    except Exception as e:
        print(f"(drug metadata unavailable: {e})")

    # a few held-out example "patients" (raw panel betas) for the demo's Load Example
    examples = []
    seen_ct = set()
    for L, npairs in sorted(cand, key=lambda t: -t[1]):
        cidL = lines[L]
        ctL = anno.reindex([cidL])["cancer_type"].iloc[0] if cidL in anno.index else "Unknown"
        if ctL in seen_ct or npairs < 20:
            continue
        seen_ct.add(ctL)
        selL = te_idx[li[te_idx] == L]
        true_top = [drug_names[int(di[selL][k])] for k in np.argsort(y[selL])[:5]]
        betas = {g: round(float(M[L, gi]), 4) for gi, g in enumerate(panel_genes)}
        examples.append({"label": f"{ctL} — {cidL}", "cancer_type": str(ctL),
                         "betas": betas, "true_top_sensitive": true_top})
        if len(examples) >= 5:
            break

    poly = os.path.join(OUTPUT_DIR, "polymer_model.pt")
    torch.save({
        "model_state": net.state_dict(),
        "panel_genes": panel_genes,
        "mu": np.asarray(mu), "sd": np.asarray(sd), "ymu": float(ymu), "ysd": float(ysd),
        "drug_names": drug_names, "drug_meta": drug_meta, "drug_mean": drug_mean_by_name,
        "n_genes": N_GENES, "emb": EMB, "hid": HID, "n_drugs": n_drugs,
        "examples": examples,
        "metrics": {"rmse": float(m_rmse), "r": float(m_r),
                    "base_rmse": float(b_rmse), "base_r": float(b_r),
                    "demeaned_rho": float(np.nanmean(dem_rhos)),
                    "n_train_lines": int(train_line_mask.sum())},
    }, poly)
    print(f"saved FULL Polymer artifact -> {poly}  ({len(examples)} example patients)")


if __name__ == "__main__":
    main()
