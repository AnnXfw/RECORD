#!/usr/bin/env python3
"""比较BG4阳性与BG4阴性蛋白峰的CUT&Tag富集强度。"""

import argparse
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyBigWig
from scipy.stats import mannwhitneyu


def read_bed(path):
    records = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                continue
            if end > start:
                records.append((fields[0], start, end))
    return records


def bw_mean(handle, chrom, start, end):
    chroms = handle.chroms()
    if chrom not in chroms:
        return np.nan
    start = max(0, start)
    end = min(end, chroms[chrom])
    if end <= start:
        return np.nan
    value = handle.stats(chrom, start, end, type="mean", exact=False)[0]
    if value is None or not math.isfinite(value):
        return np.nan
    return float(value)


def score_group(records, handles, target, group):
    rows = []
    for chrom, start, end in records:
        rep_values = [
            bw_mean(handle, chrom, start, end)
            for handle in handles
        ]
        finite = [value for value in rep_values if math.isfinite(value)]
        mean_value = float(np.mean(finite)) if finite else np.nan
        rows.append({
            "target": target,
            "group": group,
            "chrom": chrom,
            "start": start,
            "end": end,
            "rep1_log2ratio": rep_values[0],
            "rep2_log2ratio": rep_values[1],
            "rep3_log2ratio": rep_values[2],
            "mean_log2ratio": mean_value,
        })
    return rows


def fmt(value):
    return "NA" if not math.isfinite(value) else f"{value:.8g}"


def write_tsv(path, rows):
    columns = [
        "target", "group", "chrom", "start", "end",
        "rep1_log2ratio", "rep2_log2ratio", "rep3_log2ratio",
        "mean_log2ratio",
    ]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=columns)
        writer.writeheader()
        for row in rows:
            output = dict(row)
            for key in columns[5:]:
                output[key] = fmt(output[key])
            writer.writerow(output)


def p_label(p):
    if p < 1e-4:
        return "P < 0.0001"
    return f"P = {p:.4f}"


def draw(path, target, positive, negative, pvalue):
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 9,
        "axes.linewidth": 0.9,
        "svg.fonttype": "none",
    })
    fig, ax = plt.subplots(figsize=(4.5, 4.7))
    data = [negative, positive]
    positions = [0, 1]

    violins = ax.violinplot(
        data, positions=positions, widths=0.72,
        showmeans=False, showmedians=False, showextrema=False,
    )
    violin_styles = [
        ("white", "black"),
        ("#51458A", "black"),
    ]
    for body, (face, edge) in zip(violins["bodies"], violin_styles):
        body.set_facecolor(face)
        body.set_edgecolor(edge)
        body.set_linewidth(0.9)
        body.set_alpha(0.9)

    box = ax.boxplot(
        data, positions=positions, widths=0.22,
        patch_artist=True, showfliers=False,
        medianprops={"color": "white", "linewidth": 1.2},
        whiskerprops={"color": "black", "linewidth": 0.8},
        capprops={"color": "black", "linewidth": 0.8},
        boxprops={"edgecolor": "black", "linewidth": 0.8},
    )
    box["boxes"][0].set_facecolor("#7A7A7A")
    box["boxes"][1].set_facecolor("#2B225F")

    # 显示逐峰数据；大组仅确定性抽取最多600个点，避免SVG过重
    rng = np.random.default_rng(20260723)
    for position, values, color in zip(
        positions, data, ["#737373", "#6A51A3"]
    ):
        shown = values
        if values.size > 600:
            shown = rng.choice(values, size=600, replace=False)
        jitter = rng.uniform(-0.18, 0.18, size=shown.size)
        ax.scatter(
            position + jitter, shown,
            s=7, color=color, alpha=0.28,
            edgecolors="none", rasterized=False,
        )

    ymin = min(float(np.min(negative)), float(np.min(positive)))
    ymax = max(float(np.max(negative)), float(np.max(positive)))
    span = max(ymax - ymin, 1.0)
    bracket_y = ymax + 0.08 * span
    ax.plot(
        [0, 0, 1, 1],
        [bracket_y - 0.02 * span, bracket_y, bracket_y, bracket_y - 0.02 * span],
        color="black", linewidth=0.8,
    )
    ax.text(0.5, bracket_y + 0.025 * span, p_label(pvalue),
            ha="center", va="bottom", fontsize=9)
    ax.set_ylim(ymin - 0.08 * span, bracket_y + 0.16 * span)

    ax.set_xticks(positions)
    ax.set_xticklabels([
        f"BG4−\n(n={negative.size:,})",
        f"BG4+\n(n={positive.size:,})",
    ])
    ax.set_ylabel(f"{target} mean log2(Rep/Ctrl) signal")
    ax.tick_params(direction="out", length=3, width=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)
    fig.subplots_adjust(left=0.22, right=0.97, bottom=0.16, top=0.95)
    fig.savefig(path, format="svg", bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", required=True)
    parser.add_argument("--with-bg4", required=True)
    parser.add_argument("--without-bg4", required=True)
    parser.add_argument("--bigwig", nargs=3, required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    handles = [pyBigWig.open(path) for path in args.bigwig]
    try:
        positive_rows = score_group(
            read_bed(args.with_bg4), handles, args.target, "BG4_positive"
        )
        negative_rows = score_group(
            read_bed(args.without_bg4), handles, args.target, "BG4_negative"
        )
    finally:
        for handle in handles:
            handle.close()

    all_rows = negative_rows + positive_rows
    tsv = out / f"{args.target}.BG4_positive_negative.per_peak.tsv"
    write_tsv(tsv, all_rows)

    positive = np.asarray([
        row["mean_log2ratio"] for row in positive_rows
        if math.isfinite(row["mean_log2ratio"])
    ])
    negative = np.asarray([
        row["mean_log2ratio"] for row in negative_rows
        if math.isfinite(row["mean_log2ratio"])
    ])
    if positive.size == 0 or negative.size == 0:
        raise SystemExit("至少一个比较组没有有效信号")

    test = mannwhitneyu(
        positive, negative, alternative="two-sided", method="auto"
    )
    svg = out / f"{args.target}.BG4_positive_negative.signal.svg"
    draw(svg, args.target, positive, negative, test.pvalue)

    summary = out / f"{args.target}.BG4_positive_negative.summary.tsv"
    with open(summary, "w", encoding="utf-8") as handle:
        handle.write(
            "target\tBG4_positive_n\tBG4_negative_n\t"
            "BG4_positive_median\tBG4_negative_median\t"
            "median_difference\tmannwhitney_U\tp_value\n"
        )
        handle.write(
            f"{args.target}\t{positive.size}\t{negative.size}\t"
            f"{np.median(positive):.8g}\t{np.median(negative):.8g}\t"
            f"{np.median(positive)-np.median(negative):.8g}\t"
            f"{test.statistic:.8g}\t{test.pvalue:.8g}\n"
        )
    print(summary.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
