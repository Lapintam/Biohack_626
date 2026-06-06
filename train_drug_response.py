"""
Train a neural net to predict drug response (LN_IC50) across cell lines, using
the merged GDSC2 table (`gdsc2_merged.csv`).

Goal: generalize to *held-out* cell lines -- i.e. predict response for lines the
model never saw in training. This is the "cell-line-blind" regime and is the
honest test of whether the cell-line metadata + drug identity carry signal
beyond the trivial per-drug average.

Features in this table are coarse (tissue, cancer_type, MSI, mutational_burden,
ploidy, age) + drug identity. There is NO per-gene expression here, so the bar
to beat is the per-drug-mean baseline. We report RMSE and Pearson r on a
held-out set of cell lines, against two baselines.

Model: learned embeddings for the drug + categorical cell-line metadata,
concatenated with standardized numeric features, fed to a small MLP.

Run:  python train_drug_response.py
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

import torch
import torch.nn as nn

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "gdsc2_merged.csv")

SEED = 1337
LABEL = "LN_IC50"
GROUP = "SANGER_MODEL_ID"          # split on this -> cell-line-blind
DRUG = "DRUG_ID"
CAT_COLS = ["DRUG_ID", "tissue", "cancer_type", "tissue_status",
            "growth_properties", "gender", "ethnicity", "msi_status",
            "smoking_status"]
NUM_COLS = ["mutational_burden", "ploidy_wes", "age_at_sampling"]

VAL_FRAC, TEST_FRAC = 0.15, 0.15
EPOCHS, BATCH, LR, PATIENCE = 120, 1024, 2e-3, 12


def metrics(y_true, y_pred):
    rmse = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))
    r = float(pearsonr(y_true, y_pred)[0]) if len(y_true) > 2 else float("nan")
    return rmse, r


def group_split(groups, seed=SEED):
    """Partition unique groups (cell lines) into train/val/test index arrays."""
    rng = np.random.default_rng(seed)
    uniq = np.array(sorted(pd.unique(groups)))
    rng.shuffle(uniq)
    n = len(uniq)
    n_test = int(n * TEST_FRAC)
    n_val = int(n * VAL_FRAC)
    test_g = set(uniq[:n_test])
    val_g = set(uniq[n_test:n_test + n_val])
    g = np.asarray(groups)
    te = np.where(np.isin(g, list(test_g)))[0]
    va = np.where(np.isin(g, list(val_g)))[0]
    tr = np.where(~np.isin(g, list(test_g | val_g)))[0]
    return tr, va, te


# ---------------------------------------------------------------------------
# Baselines (fit on train only)
# ---------------------------------------------------------------------------
def baseline_per_drug(df, tr, te):
    dm = df.iloc[tr].groupby(DRUG)[LABEL].mean()
    glob = df.iloc[tr][LABEL].mean()
    pred = df.iloc[te][DRUG].map(dm).fillna(glob).to_numpy()
    return metrics(df.iloc[te][LABEL].to_numpy(), pred)


def baseline_per_drug_tissue(df, tr, te):
    dt = df.iloc[tr].groupby([DRUG, "tissue"])[LABEL].mean()
    dm = df.iloc[tr].groupby(DRUG)[LABEL].mean()
    glob = df.iloc[tr][LABEL].mean()
    sub = df.iloc[te]
    pred = np.array([
        dt.get((d, t), dm.get(d, glob))
        for d, t in zip(sub[DRUG].to_numpy(), sub["tissue"].to_numpy())
    ])
    return metrics(sub[LABEL].to_numpy(), pred)


# ---------------------------------------------------------------------------
# Encoders fit on TRAIN only (no leakage)
# ---------------------------------------------------------------------------
def encode_categoricals(df, tr):
    """Return integer code matrix and per-column cardinalities.
    Index 0 is reserved for 'unknown / unseen at train time'."""
    codes = np.zeros((len(df), len(CAT_COLS)), dtype=np.int64)
    cards = []
    for j, c in enumerate(CAT_COLS):
        train_vals = pd.unique(df.iloc[tr][c])
        vocab = {v: i + 1 for i, v in enumerate(sorted(train_vals, key=str))}
        codes[:, j] = df[c].map(vocab).fillna(0).astype(np.int64).to_numpy()
        cards.append(len(vocab) + 1)
    return codes, cards


def encode_numerics(df, tr):
    X = df[NUM_COLS].to_numpy(dtype=np.float32)
    mu = X[tr].mean(axis=0)
    sd = X[tr].std(axis=0)
    sd[sd == 0] = 1.0
    return (X - mu) / sd


class DrugResponseNet(nn.Module):
    def __init__(self, cards, n_num):
        super().__init__()
        # embedding width grows with cardinality, capped
        dims = [min(50, max(4, int(round(1.6 * c ** 0.56)))) for c in cards]
        self.embs = nn.ModuleList(
            [nn.Embedding(c, d) for c, d in zip(cards, dims)]
        )
        emb_total = sum(dims)
        self.body = nn.Sequential(
            nn.Linear(emb_total + n_num, 256), nn.ReLU(), nn.BatchNorm1d(256), nn.Dropout(0.3),
            nn.Linear(256, 128), nn.ReLU(), nn.BatchNorm1d(128), nn.Dropout(0.2),
            nn.Linear(128, 1),
        )

    def forward(self, xc, xn):
        e = [emb(xc[:, j]) for j, emb in enumerate(self.embs)]
        h = torch.cat(e + [xn], dim=1)
        return self.body(h)


def main():
    torch.manual_seed(SEED)
    np.random.seed(SEED)

    df = pd.read_csv(CSV)
    df = df.dropna(subset=[LABEL]).reset_index(drop=True)

    tr, va, te = group_split(df[GROUP].to_numpy())
    print(f"rows: {len(df):,} | cell lines: {df[GROUP].nunique()} | drugs: {df[DRUG].nunique()}")
    print(f"split (cell-line-blind)  train={len(tr):,}  val={len(va):,}  test={len(te):,}")
    print(f"  train lines={df.iloc[tr][GROUP].nunique()}  "
          f"val lines={df.iloc[va][GROUP].nunique()}  test lines={df.iloc[te][GROUP].nunique()}")

    # ----- baselines -----
    b1 = baseline_per_drug(df, tr, te)
    b2 = baseline_per_drug_tissue(df, tr, te)
    print(f"\nbaseline per-drug         RMSE={b1[0]:.3f}  r={b1[1]:.3f}")
    print(f"baseline per-drug+tissue  RMSE={b2[0]:.3f}  r={b2[1]:.3f}")

    # ----- features -----
    Xc, cards = encode_categoricals(df, tr)
    Xn = encode_numerics(df, tr)
    y = df[LABEL].to_numpy(dtype=np.float32)
    ymu, ysd = y[tr].mean(), y[tr].std() + 1e-8

    dev = "cuda" if torch.cuda.is_available() else "cpu"
    Xc_t = torch.tensor(Xc, device=dev)
    Xn_t = torch.tensor(Xn, device=dev)
    y_t = torch.tensor((y - ymu) / ysd, device=dev).unsqueeze(1)
    tr_t = torch.tensor(tr, device=dev)

    net = DrugResponseNet(cards, len(NUM_COLS)).to(dev)
    opt = torch.optim.Adam(net.parameters(), lr=LR, weight_decay=1e-5)
    sched = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, factor=0.5, patience=4)
    loss_fn = nn.MSELoss()

    def predict(idx):
        net.eval()
        with torch.no_grad():
            i = torch.tensor(idx, device=dev)
            p = net(Xc_t[i], Xn_t[i]).cpu().numpy().ravel()
        return p * ysd + ymu

    best_val, best_state, bad = float("inf"), None, 0
    for ep in range(1, EPOCHS + 1):
        net.train()
        perm = tr_t[torch.randperm(len(tr_t), device=dev)]
        for i in range(0, len(perm), BATCH):
            b = perm[i:i + BATCH]
            opt.zero_grad()
            loss = loss_fn(net(Xc_t[b], Xn_t[b]), y_t[b])
            loss.backward()
            opt.step()
        val_rmse, val_r = metrics(y[va], predict(va))
        sched.step(val_rmse)
        if val_rmse < best_val - 1e-4:
            best_val, best_state, bad = val_rmse, {k: v.detach().clone() for k, v in net.state_dict().items()}, 0
        else:
            bad += 1
        if ep % 10 == 0 or ep == 1:
            print(f"epoch {ep:3d}  val RMSE={val_rmse:.3f}  r={val_r:.3f}  (best={best_val:.3f})")
        if bad >= PATIENCE:
            print(f"early stop at epoch {ep} (no val improvement for {PATIENCE})")
            break

    if best_state is not None:
        net.load_state_dict(best_state)

    nn_rmse, nn_r = metrics(y[te], predict(te))
    print("\n================ TEST (held-out cell lines) ================")
    print(f"baseline per-drug         RMSE={b1[0]:.3f}  r={b1[1]:.3f}")
    print(f"baseline per-drug+tissue  RMSE={b2[0]:.3f}  r={b2[1]:.3f}")
    print(f"neural net (embed MLP)    RMSE={nn_rmse:.3f}  r={nn_r:.3f}")
    delta = (b1[0] - nn_rmse) / b1[0] * 100
    print(f"\nNN vs per-drug baseline: {delta:+.1f}% RMSE  "
          f"({'beats' if nn_rmse < b1[0] else 'does NOT beat'} baseline)")

    torch.save(net.state_dict(), os.path.join(HERE, "drug_response_net.pt"))
    print("saved drug_response_net.pt")


if __name__ == "__main__":
    main()
