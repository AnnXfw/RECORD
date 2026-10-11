#!/usr/bin/env python3
"""绘制BG4阳性蛋白峰的三重复log2(Rep/Ctrl)逐峰热图。"""

import argparse
import csv
import gzip
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm


def load_matrix(path):
    metadata = []
    values = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        first = handle.readline().rstrip()
        if not first.startswith("@"):
            raise ValueError("不是有效的deepTools矩阵")
        header = json.loads(first[1:])
        boundaries = header["sample_boundaries"]
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) <= 6:
                continue
            metadata.append(fields[:6])
            row = []
            for item in fields[6:]:
                try:
                    row.append(float(item))
                except ValueError:
                    row.append(np.nan)
            values.append(row)

    full = np.asarray(values, dtype=float)
    if full.size == 0:
        raise ValueError("矩阵为空")
    samples = np.stack([
        full[:, start:end]
        for start, end in zip(boundaries[:-1], boundaries[1:])
    ])
    if samples.shape[0] != 3:
        raise ValueError(f"预期3个重复，实际为{samples.shape[0]}")
    return metadata, samples, header


def write_sorted_tsv(path, metadata, samples, order):
    fields = [
        "rank", "chrom", "start", "end", "name", "score", "strand",
        "rep1_mean_log2ratio", "rep2_mean_log2ratio",
        "rep3_mean_log2ratio", "three_rep_mean_log2ratio",
    ]
    per_rep = np.nanmean(samples, axis=2)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(fields)
        for rank, idx in enumerate(order, start=1):
            rep_values = per_rep[:, idx]
            writer.writerow([
                rank, *metadata[idx],
                *[f"{value:.8g}" for value in rep_values],
                f"{np.nanmean(rep_values):.8g}",
            ])


def draw(path, target, samples, order, upstream, downstream):
    sorted_samples = samples[:, order, :]
    finite = sorted_samples[np.isfinite(sorted_samples)]
    if finite.size == 0:
        raise ValueError("没有有限信号值")
    # 对称色标，以98百分位抑制少数极端值对主体颜色的影响
    limit = float(np.nanpercentile(np.abs(finite), 98))
    if not np.isfinite(limit) or limit <= 0:
        limit = 1.0
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 9,
        "axes.linewidth": 0.8,
        "svg.fonttype": "none",
    })
    fig = plt.figure(figsize=(7.4, 7.0))
    grid = fig.add_gridspec(
        2, 3,
        height_ratios=[1.15, 4.3],
        hspace=0.07, wspace=0.12,
    )
    profile_axes = [fig.add_subplot(grid[0, idx]) for idx in range(3)]
    heatmap_axes = [
        fig.add_subplot(grid[1, idx], sharey=None if idx == 0 else None)
        for idx in range(3)
    ]
    x = np.linspace(-upstream, downstream, sorted_samples.shape[2], endpoint=False)
    profile_colors = ["#222222", "#3056D3", "#E53935"]
    mean_profiles = np.nanmean(sorted_samples, axis=1)
    profile_min = float(np.nanmin(mean_profiles))
    profile_max = float(np.nanmax(mean_profiles))
    profile_span = max(profile_max - profile_min, 0.5)

    image = None
    for idx, (profile_ax, heatmap_ax) in enumerate(
        zip(profile_axes, heatmap_axes)
    ):
        # 上方平均曲线与下方同一重复热图严格对应
        profile_ax.plot(
            x, mean_profiles[idx],
            color=profile_colors[idx], linewidth=1.35,
        )
        profile_ax.axhline(0, color="#BDBDBD", linewidth=0.6)
        profile_ax.axvline(
            0, color="black", linestyle=(0, (3, 3)), linewidth=0.55
        )
        profile_ax.set_xlim(-upstream, downstream)
        profile_ax.set_ylim(
            profile_min - 0.08 * profile_span,
            profile_max + 0.10 * profile_span,
        )
        profile_ax.set_title(f"Rep {idx + 1}/Ctrl {idx + 1}", fontsize=9)
        profile_ax.tick_params(
            axis="x", labelbottom=False, length=0,
            direction="out", width=0.7,
        )
        profile_ax.tick_params(axis="y", direction="out", length=2.5, width=0.7)
        profile_ax.spines["top"].set_visible(False)
        profile_ax.spines["right"].set_visible(False)
        if idx == 0:
            profile_ax.set_ylabel("Mean\nlog2(Rep/Ctrl)", fontsize=8)
        else:
            profile_ax.tick_params(labelleft=False)
            profile_ax.spines["left"].set_visible(False)

        image = heatmap_ax.imshow(
            sorted_samples[idx],
            aspect="auto",
            interpolation="nearest",
            cmap="RdBu_r",
            norm=norm,
            extent=[-upstream, downstream, sorted_samples.shape[1], 0],
            rasterized=True,
        )
        heatmap_ax.axvline(
            0, color="black", linestyle=(0, (3, 3)), linewidth=0.55
        )
        heatmap_ax.set_xlim(-upstream, downstream)
        heatmap_ax.set_xticks([-upstream, 0, downstream])
        heatmap_ax.set_xticklabels([
            f"-{upstream // 1000} kb",
            "BG4 center",
            f"+{downstream // 1000} kb",
        ], rotation=45, ha="right")
        heatmap_ax.tick_params(direction="out", length=2.5, width=0.7)
        heatmap_ax.spines["top"].set_visible(False)
        heatmap_ax.spines["right"].set_visible(False)
        if idx > 0:
            heatmap_ax.tick_params(labelleft=False)
            heatmap_ax.spines["left"].set_visible(False)
    heatmap_axes[0].set_ylabel(
        f"{target} BG4-positive peaks\n(n={sorted_samples.shape[1]:,})"
    )

    colorbar = fig.colorbar(
        image, ax=heatmap_axes, fraction=0.025, pad=0.025,
        label="log2(Rep/Ctrl)",
    )
    colorbar.outline.set_linewidth(0.7)
    fig.subplots_adjust(left=0.14, right=0.88, bottom=0.15, top=0.95)
    fig.savefig(path, format="svg", bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    metadata, samples, header = load_matrix(args.matrix)

    # 每个峰先在每个重复内求窗口均值，再跨三个重复求均值排序
    peak_score = np.nanmean(np.nanmean(samples, axis=2), axis=0)
    peak_score = np.nan_to_num(peak_score, nan=-np.inf)
    order = np.argsort(peak_score)[::-1]

    tsv = out / f"{args.target}.BG4_positive_heatmap.sorted_peaks.tsv"
    svg = out / f"{args.target}.BG4_positive_heatmap.svg"
    write_sorted_tsv(tsv, metadata, samples, order)
    draw(
        svg, args.target, samples, order,
        header["upstream"][0], header["downstream"][0],
    )
    print(
        f"{args.target}: peaks={samples.shape[1]}, bins={samples.shape[2]}, "
        f"SVG={svg}, TSV={tsv}"
    )


if __name__ == "__main__":
    main()
