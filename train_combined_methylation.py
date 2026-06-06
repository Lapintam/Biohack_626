"""
Does FUSING gene-level + promoter methylation help predict drug response?

Apples-to-apples test: restrict to the cell lines present in BOTH matrices, use
ONE cell-line-blind split (seed 1337), and train the SAME recommender
architecture on three feature sets:

    gene_only  : top-K variance gene-level betas
    prom_only  : top-K variance promoter betas
    combined   : [gene_only | prom_only]  (concatenated)

Because the cohort and split are identical across all three, any difference is
purely the feature set. Reported metrics match the Polymer benchmark, headlined
by personalized_R^2 (R^2 on per-drug-mean residuals; per-drug-mean baseline = 0).

Usage:  python train_combined_methylation.py [--k 3000] [--epochs 40]
"""

import argparse
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

import torch
import torch.nn as nn

from strata.data import gdsc

SEED = 1337
EMB, HID, BATCH, LR = 32, 256, 4096, 2e-3
VAL_FRAC = TEST_FRAC = 0.15


def metrics(y, p, drug_idx, drug_mean_tr, lines_te, li):
    rmse = float(np.sqrt(np.mean((y - p) ** 2)))
    r = float(pearsonr(y, p)[0])
    ss = lambda a, b: 1 - np.sum((a - b) ** 2) / np.sum((a - a.mean()) ** 2)
    r2 = float(ss(y, p))
    dm = np.array([drug_mean_tr.get(d, y.mean()) for d in drug_idx])
    pers_r2 = float(ss(y - dm, p - dm))
    # within-line de-meaned Spearman
    rhos = []
    for L in set(lines_te):
        m = li == L
        if m.sum() < 10:
            continue
        yy, pp, dd = y[m], p[m], dm[m]
        if np.ptp(pp - dd) and np.ptp(yy - dd):
            rhos.append(spearmanr(pp - dd, yy - dd).correlation)
    return {"RMSE": rmse, "R2": r2, "Pearson_r": r,
            "personalized_R2": pers_r2,
            "within_line_rho_demeaned": float(np.nanmean(rhos))}


def build_feats(M, train_mask, k):
    """Top-k variance columns (train-selected), NaN->train col-mean, z-scored on train."""
    M = M.astype(np.float64).copy()
    if np.isnan(M).any():
        cm = np.nanmean(M[train_mask], axis=0)
        cm = np.where(np.isfinite(cm), cm, 0.0)
        r, c = np.where(np.isnan(M))
        M[r, c] = np.take(cm, c)
    var = M[train_mask].var(0)
    panel = np.argsort(var)[::-1][:k]
    X = M[:, panel]
    mu = X[train_mask].mean(0); sd = X[train_mask].std(0); sd[sd == 0] = 1.0
    return ((X - mu) / sd).astype(np.float32)


class Rec(nn.Module):
    def __init__(self, n_feat, n_drugs):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(n_feat, HID), nn.ReLU(), nn.BatchNorm1d(HID), nn.Dropout(0.3),
            nn.Linear(HID, EMB))
        self.drug_emb = nn.Embedding(n_drugs, EMB)
        self.drug_bias = nn.Embedding(n_drugs, 1)
        nn.init.normal_(self.drug_emb.weight, std=0.1); nn.init.zeros_(self.drug_bias.weight)

    def forward(self, Xz, li, di):
        ce = self.encoder(Xz[li])
        return self.drug_bias(di).squeeze(1) + (ce * self.drug_emb(di)).sum(1)


class TwoTower(nn.Module):
    """Per-modality encoders fused by a learned per-dimension gate.
    fused = g * gene_emb + (1 - g) * prom_emb, with g in (0,1)^EMB learned from both."""
    def __init__(self, n_gene, n_prom, n_drugs):
        super().__init__()
        def enc(n):
            return nn.Sequential(nn.Linear(n, HID), nn.ReLU(), nn.BatchNorm1d(HID),
                                 nn.Dropout(0.3), nn.Linear(HID, EMB))
        self.gene_enc = enc(n_gene)
        self.prom_enc = enc(n_prom)
        self.gate = nn.Linear(2 * EMB, EMB)
        self.drug_emb = nn.Embedding(n_drugs, EMB)
        self.drug_bias = nn.Embedding(n_drugs, 1)
        nn.init.normal_(self.drug_emb.weight, std=0.1); nn.init.zeros_(self.drug_bias.weight)

    def fuse(self, Fg, Fp):
        ge, pe = self.gene_enc(Fg), self.prom_enc(Fp)
        g = torch.sigmoid(self.gate(torch.cat([ge, pe], dim=1)))
        return g * ge + (1 - g) * pe, g

    def forward(self, Fg, Fp, li, di):
        fused, _ = self.fuse(Fg[li], Fp[li])
        return self.drug_bias(di).squeeze(1) + (fused * self.drug_emb(di)).sum(1)


def train_eval_twotower(Fg, Fp, li, di, y, tr, va, te, n_drugs, epochs):
    torch.manual_seed(SEED); rng = np.random.default_rng(SEED)
    Fg_t, Fp_t = torch.tensor(Fg), torch.tensor(Fp)
    li_t, di_t = torch.tensor(li), torch.tensor(di)
    ymu, ysd = y[tr].mean(), y[tr].std() + 1e-8
    yz = torch.tensor((y - ymu) / ysd, dtype=torch.float32)
    net = TwoTower(Fg.shape[1], Fp.shape[1], n_drugs)
    opt = torch.optim.Adam(net.parameters(), lr=LR, weight_decay=1e-5)
    loss_fn = nn.MSELoss()
    tr_idx = np.where(tr)[0]

    @torch.no_grad()
    def predict(mask):
        net.eval()
        return net(Fg_t, Fp_t, li_t[mask], di_t[mask]).numpy() * ysd + ymu

    best, best_state = 1e9, None
    for ep in range(epochs):
        net.train(); perm = tr_idx[rng.permutation(len(tr_idx))]
        for s in range(0, len(perm), BATCH):
            b = torch.tensor(perm[s:s + BATCH])
            opt.zero_grad()
            loss = loss_fn(net(Fg_t, Fp_t, li_t[b], di_t[b]), yz[b]); loss.backward(); opt.step()
        vp = predict(torch.tensor(np.where(va)[0]))
        vr = float(np.sqrt(np.mean((y[va] - vp) ** 2)))
        if vr < best - 1e-5:
            best, best_state = vr, {k: v.clone() for k, v in net.state_dict().items()}
    if best_state:
        net.load_state_dict(best_state)
    # mean gate on test lines: how much weight goes to the GENE tower (0..1)
    net.eval()
    with torch.no_grad():
        _, g = net.fuse(Fg_t, Fp_t)
        gate_gene = float(g.mean().item())
    return predict(torch.tensor(np.where(te)[0])), gate_gene


def train_eval(Xz, li, di, y, tr, va, te, n_drugs, epochs):
    torch.manual_seed(SEED); rng = np.random.default_rng(SEED)
    Xz_t = torch.tensor(Xz); li_t = torch.tensor(li); di_t = torch.tensor(di)
    ymu, ysd = y[tr].mean(), y[tr].std() + 1e-8
    yz = torch.tensor((y - ymu) / ysd, dtype=torch.float32)
    net = Rec(Xz.shape[1], n_drugs)
    opt = torch.optim.Adam(net.parameters(), lr=LR, weight_decay=1e-5)
    loss_fn = nn.MSELoss()
    tr_idx = np.where(tr)[0]

    @torch.no_grad()
    def predict(mask):
        net.eval()
        return net(Xz_t, li_t[mask], di_t[mask]).numpy() * ysd + ymu

    best, best_state = 1e9, None
    for ep in range(epochs):
        net.train(); perm = tr_idx[rng.permutation(len(tr_idx))]
        for s in range(0, len(perm), BATCH):
            b = torch.tensor(perm[s:s + BATCH])
            opt.zero_grad()
            loss = loss_fn(net(Xz_t, li_t[b], di_t[b]), yz[b]); loss.backward(); opt.step()
        vp = predict(torch.tensor(np.where(va)[0]))
        vr = float(np.sqrt(np.mean((y[va] - vp) ** 2)))
        if vr < best - 1e-5:
            best, best_state = vr, {k: v.clone() for k, v in net.state_dict().items()}
    if best_state:
        net.load_state_dict(best_state)
    return predict(torch.tensor(np.where(te)[0]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=3000, help="top-variance features per modality")
    ap.add_argument("--epochs", type=int, default=40)
    args = ap.parse_args()

    gene = gdsc.load_gdsc_methylation_genelevel(); gene.index = gene.index.astype(str)
    prom = gdsc.load_gdsc_methylation_promoter();  prom.index = prom.index.astype(str)
    drug = gdsc.load_gdsc_drug_response()

    common = sorted(set(gene.index) & set(prom.index))
    drug = drug[drug["COSMIC_ID"].astype(str).isin(common)].dropna(subset=["auc"])
    drug = drug.groupby(["COSMIC_ID", "drug_name"], as_index=False)["auc"].mean()
    lines = sorted(set(drug["COSMIC_ID"].astype(str)))   # lines with BOTH meth + drug
    line_pos = {c: i for i, c in enumerate(lines)}
    drug_names = sorted(drug["drug_name"].unique())
    drug_pos = {d: i for i, d in enumerate(drug_names)}
    n_lines, n_drugs = len(lines), len(drug_names)

    li = drug["COSMIC_ID"].astype(str).map(line_pos).to_numpy()
    di = drug["drug_name"].map(drug_pos).to_numpy()
    y = drug["auc"].to_numpy(float)
    print(f"common cohort: {n_lines} lines (gene+promoter+drug) | {n_drugs} drugs | {len(y):,} pairs")

    rng = np.random.default_rng(SEED)
    order = rng.permutation(n_lines)
    n_te, n_va = int(n_lines * TEST_FRAC), int(n_lines * VAL_FRAC)
    test_l = set(order[:n_te]); val_l = set(order[n_te:n_te + n_va])
    tr = ~np.isin(li, list(test_l | val_l)); va = np.isin(li, list(val_l)); te = np.isin(li, list(test_l))
    train_line_mask = ~np.isin(np.arange(n_lines), list(test_l | val_l))
    print(f"split: train={train_line_mask.sum()} lines | val={len(val_l)} | test={len(test_l)}")

    drug_mean_tr = pd.Series(y[tr]).groupby(pd.Series(di[tr])).mean().to_dict()

    Mgene = gene.loc[lines].to_numpy()
    Mprom = prom.loc[lines].to_numpy()
    Fg = build_feats(Mgene, train_line_mask, args.k)
    Fp = build_feats(Mprom, train_line_mask, args.k)
    feats = {
        "gene_only": Fg,
        "prom_only": Fp,
        "combined":  np.hstack([Fg, Fp]),
    }

    rows = []
    for name, Xz in feats.items():
        pred = train_eval(Xz, li, di, y, tr, va, te, n_drugs, args.epochs)
        m = metrics(y[te], pred, di[te], drug_mean_tr, li[te], li[te])
        m["model"] = f"{name} ({Xz.shape[1]}f)"
        rows.append(m); print(f"  trained {name:12s} -> personalized_R2={m['personalized_R2']:.4f}")

    # gated two-tower late fusion (separate encoders per modality)
    pred, gate_gene = train_eval_twotower(Fg, Fp, li, di, y, tr, va, te, n_drugs, args.epochs)
    m = metrics(y[te], pred, di[te], drug_mean_tr, li[te], li[te])
    m["model"] = "two_tower (gated)"
    rows.append(m)
    print(f"  trained two_tower    -> personalized_R2={m['personalized_R2']:.4f} "
          f"(mean gate->gene = {gate_gene:.2f})")

    lb = pd.DataFrame(rows).set_index("model")[
        ["RMSE", "R2", "Pearson_r", "personalized_R2", "within_line_rho_demeaned"]]
    print("\n===== gene vs promoter vs concat vs TWO-TOWER (identical cohort + split) =====\n")
    print(lb.round(4).to_string())
    best = lb["personalized_R2"].idxmax()
    print(f"\nbest by personalized_R2: {best}")
    print(f"two-tower gate leaned {gate_gene*100:.0f}% toward the gene-level modality.")


if __name__ == "__main__":
    main()
