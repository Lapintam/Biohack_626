"""
Polymer — inference engine.

Loads the trained recommender (`outputs/polymer_model.pt`) and turns a pasted
methylome (gene -> beta) into a ranked list of drugs with predicted AUC.

    from polymer_infer import Polymer
    poly = Polymer()
    methylome = poly.parse_table("CDKN2A\t0.83\nMGMT\t0.10\n...")
    rows, n_matched = poly.recommend(methylome, top_n=15)

Lower predicted AUC = more sensitive. `delta_vs_avg` < 0 means the model thinks
this methylome is MORE sensitive to that drug than the average cell line — the
personalized signal.
"""

import os
import numpy as np
import torch
import torch.nn as nn

ARTIFACT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "outputs", "polymer_model.pt")


class _Recommender(nn.Module):
    """Must mirror the architecture trained in train_drug_recommender.py."""
    def __init__(self, n_genes, hid, emb, n_drugs):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(n_genes, hid), nn.ReLU(), nn.BatchNorm1d(hid), nn.Dropout(0.3),
            nn.Linear(hid, emb),
        )
        self.drug_emb = nn.Embedding(n_drugs, emb)
        self.drug_bias = nn.Embedding(n_drugs, 1)

    def encode(self, x):
        return self.encoder(x)


class Polymer:
    def __init__(self, path=ARTIFACT):
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"{path} not found. Train it first: python train_drug_recommender.py")
        ck = torch.load(path, weights_only=False)
        self.genes = list(ck["panel_genes"])
        self.gene_idx = {g: i for i, g in enumerate(self.genes)}
        self.mu = np.asarray(ck["mu"], dtype=np.float64)
        self.sd = np.asarray(ck["sd"], dtype=np.float64)
        self.ymu, self.ysd = float(ck["ymu"]), float(ck["ysd"])
        self.drug_names = list(ck["drug_names"])
        self.drug_meta = ck["drug_meta"]
        self.drug_mean = ck["drug_mean"]
        self.examples = ck.get("examples", [])
        self.metrics = ck.get("metrics", {})

        net = _Recommender(ck["n_genes"], ck["hid"], ck["emb"], ck["n_drugs"])
        net.load_state_dict(ck["model_state"])
        net.eval()
        self.net = net
        with torch.no_grad():
            self.D = net.drug_emb.weight.detach().numpy()             # (n_drugs, emb)
            self.bias = net.drug_bias.weight.detach().numpy().ravel()  # (n_drugs,)

    # ---- input parsing ----------------------------------------------------
    @staticmethod
    def parse_table(text: str) -> dict:
        """Parse pasted 'gene<sep>beta' lines. Accepts tab, comma, or whitespace.
        Skips a header row and any unparseable lines."""
        out = {}
        for raw in text.strip().splitlines():
            line = raw.strip()
            if not line:
                continue
            parts = None
            for sep in ("\t", ",", None):
                p = line.split(sep) if sep else line.split()
                if len(p) >= 2:
                    parts = p
                    break
            if not parts:
                continue
            gene = parts[0].strip().strip('"')
            try:
                out[gene] = float(parts[1])
            except ValueError:
                continue  # header like "gene,beta"
        return out

    # ---- featurize + score ------------------------------------------------
    def vectorize(self, methylome: dict):
        """Align a {gene: beta} dict to the panel; missing genes default to the
        training mean (z=0). Returns (z-scored vector, n_genes_matched)."""
        x = self.mu.copy()
        matched = 0
        for g, b in methylome.items():
            j = self.gene_idx.get(g)
            if j is not None and np.isfinite(b):
                x[j] = b
                matched += 1
        z = (x - self.mu) / self.sd
        return z.astype(np.float32), matched

    @torch.no_grad()
    def recommend(self, methylome: dict, top_n: int | None = None):
        z, matched = self.vectorize(methylome)
        emb = self.net.encode(torch.tensor(z).unsqueeze(0)).numpy().ravel()
        scores_std = self.bias + self.D @ emb              # standardized AUC
        auc = scores_std * self.ysd + self.ymu             # back to AUC units
        rows = []
        for i, dn in enumerate(self.drug_names):
            mean = self.drug_mean.get(dn, float("nan"))
            meta = self.drug_meta.get(dn, {})
            rows.append({
                "drug": dn,
                "pred_auc": float(auc[i]),
                "avg_auc": float(mean) if mean == mean else float("nan"),
                "delta_vs_avg": float(auc[i] - mean) if mean == mean else float("nan"),
                "target": meta.get("target", ""),
                "pathway": meta.get("pathway", ""),
            })
        rows.sort(key=lambda r: r["pred_auc"])             # most sensitive first
        if top_n:
            rows = rows[:top_n]
        return rows, matched


if __name__ == "__main__":
    poly = Polymer()
    print(f"Polymer loaded: {len(poly.genes)} panel genes, {len(poly.drug_names)} drugs")
    print(f"model metrics: {poly.metrics}")
    ex = poly.examples[0]
    print(f"\nExample patient: {ex['label']}  (true most-sensitive: {ex['true_top_sensitive']})")
    rows, matched = poly.recommend(ex["betas"], top_n=10)
    print(f"matched {matched}/{len(poly.genes)} panel genes\n")
    print(f"{'drug':<24}{'pred AUC':>9}{'vs avg':>9}  target")
    for r in rows:
        print(f"{r['drug']:<24}{r['pred_auc']:>9.3f}{r['delta_vs_avg']:>+9.3f}  {r['target']}")
