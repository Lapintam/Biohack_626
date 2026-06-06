import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import umap
from strata.config import OUTPUT_DIR


def plot_umap(meth, features, clusters, response=None, fname="umap.png"):
    emb = umap.UMAP(random_state=0).fit_transform(meth[features].values)
    fig, ax = plt.subplots(figsize=(5, 4))
    c = response.loc[meth.index] if response is not None else clusters.loc[meth.index]
    sc = ax.scatter(emb[:, 0], emb[:, 1], c=c.values, cmap="viridis", s=12)
    fig.colorbar(sc)
    ax.set_title("Methylation UMAP")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / fname, dpi=150)
    plt.close(fig)
    return OUTPUT_DIR / fname


def plot_auc_box(clusters, auc, fname="auc_box.png"):
    df = pd.DataFrame({"cluster": clusters, "auc": auc}).dropna()
    fig, ax = plt.subplots(figsize=(5, 4))
    df.boxplot(column="auc", by="cluster", ax=ax)
    ax.set_title("Drug AUC by subgroup")
    fig.suptitle("")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / fname, dpi=150)
    plt.close(fig)
    return OUTPUT_DIR / fname


def plot_driver_heatmap(meth, clusters, driver_cpgs, fname="heatmap.png"):
    order = clusters.sort_values().index
    sub = meth.loc[order, driver_cpgs]
    fig, ax = plt.subplots(figsize=(6, 4))
    im = ax.imshow(sub.values, aspect="auto", cmap="RdBu_r", vmin=0, vmax=1)
    fig.colorbar(im)
    ax.set_title("Driver CpG methylation")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / fname, dpi=150)
    plt.close(fig)
    return OUTPUT_DIR / fname


def plot_km(survival_result, fname="km.png"):
    fig, ax = plt.subplots(figsize=(5, 4))
    for label, (t, s) in survival_result["km"].items():
        ax.step(t, s, where="post", label=f"subgroup {label}")
    ax.set_xlabel("time")
    ax.set_ylabel("survival")
    ax.legend()
    ax.set_title(f"KM (log-rank p={survival_result['logrank_p']:.3g})")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / fname, dpi=150)
    plt.close(fig)
    return OUTPUT_DIR / fname
