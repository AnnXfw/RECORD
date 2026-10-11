#!/usr/bin/env python3
"""绘制随机置换零分布与真实BG4重叠数量的比较图。"""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_null(path):
    values = []
    with open(path, encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            values.append(int(row["BG4_overlap_peaks"]))
    return np.asarray(values, dtype=float)


def read_observed(path):
    with open(path, encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["g4_set"] == "BG4":
                return {
                    "target": row["target"],
                    "observed": int(row["observed_overlap"]),
                    "null_mean": float(row["null_mean"]),
                    "fold": float(row["fold_enrichment"]),
                    "p": float(row["empirical_p"]),
                    "nperm": int(row["permutations"]),
                }
    raise ValueError("汇总表中没有BG4记录")


def draw(output, null, result):
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 9,
        "axes.linewidth": 0.9,
        "svg.fonttype": "none",
    })
    fig, ax = plt.subplots(figsize=(4.2, 4.5))

    # 使用确定性抖动，SVG每次重绘保持一致
    rng = np.random.default_rng(20260723)
    jitter = rng.uniform(-0.12, 0.12, size=null.size)
    null_plot = null + 1.0
    observed_plot = result["observed"] + 1.0

    ax.scatter(
        jitter, null_plot,
        s=7, color="#A6A6A6", alpha=0.35,
        edgecolors="none", rasterized=False,
    )
    parts = ax.violinplot(
        [null_plot], positions=[0], widths=0.55,
        showmeans=False, showmedians=True, showextrema=False,
    )
    for body in parts["bodies"]:
        body.set_facecolor("white")
        body.set_edgecolor("black")
        body.set_alpha(0.85)
        body.set_linewidth(0.8)
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(1.0)

    ax.scatter(
        [1], [observed_plot],
        marker="D", s=50,
        facecolor="#C51B7D", edgecolor="black",
        linewidth=0.7, zorder=5,
    )
    ax.plot(
        [0.18, 0.82],
        [max(null_plot) * 1.25, observed_plot / 1.25],
        color="#737373", linestyle=(0, (3, 3)), linewidth=0.8,
    )

    ax.set_yscale("log")
    ax.set_xlim(-0.55, 1.55)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Randomized", "Observed"])
    ax.set_ylabel("BG4-overlapping peaks (log scale)")
    ax.text(
        0.97, 0.95,
        f"Fold enrichment = {result['fold']:.1f}×\n"
        f"Empirical P = {result['p']:.4g}\n"
        f"{result['nperm']:,} permutations",
        transform=ax.transAxes,
        ha="right", va="top", fontsize=8,
    )
    ax.tick_params(direction="out", length=3, width=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)
    fig.subplots_adjust(left=0.22, right=0.97, bottom=0.14, top=0.97)
    fig.savefig(output, format="svg", bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--permutations", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    null = read_null(args.permutations)
    result = read_observed(args.summary)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    draw(args.output, null, result)
    print(
        f"{result['target']}: observed={result['observed']}, "
        f"null_mean={result['null_mean']}, fold={result['fold']}, "
        f"p={result['p']}"
    )


if __name__ == "__main__":
    main()
