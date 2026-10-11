#!/usr/bin/env python3
"""绘制目标蛋白CUT&Tag在BG4阳性峰中心的Nature Chemistry风格曲线。"""

import argparse
import gzip
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load_profiles(path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        first = handle.readline().rstrip()
        if not first.startswith("@"):
            raise ValueError("不是deepTools矩阵")
        header = json.loads(first[1:])
        boundaries = header["sample_boundaries"]
        rows = []
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) > 6:
                rows.append([
                    float(value) if value.lower() != "nan" else np.nan
                    for value in fields[6:]
                ])
    matrix = np.asarray(rows, dtype=float)
    if matrix.size == 0:
        raise ValueError("矩阵为空")
    profiles = np.asarray([
        np.nanmean(matrix[:, start:end], axis=0)
        for start, end in zip(boundaries[:-1], boundaries[1:])
    ])
    if profiles.shape[0] != 6:
        raise ValueError(f"预期6个样本，实际为{profiles.shape[0]}个")
    return profiles, header, matrix.shape[0]


def draw(matrix_path, output, target):
    profiles, header, region_count = load_profiles(matrix_path)
    upstream = header["upstream"][0]
    downstream = header["downstream"][0]
    x = np.linspace(-upstream, downstream, profiles.shape[1], endpoint=False)
    keep = (x >= -400) & (x <= 400)

    labels = ["Rep 1", "Rep 2", "Rep 3", "Ctrl 1", "Ctrl 2", "Ctrl 3"]
    colors = ["#222222", "#3056D3", "#E53935", "#8C8C8C", "#9ECAE1", "#E6B800"]
    widths = [1.5, 1.5, 1.5, 1.0, 1.0, 1.0]

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 9,
        "axes.linewidth": 0.9,
        "svg.fonttype": "none",
    })
    fig, ax = plt.subplots(figsize=(5.0, 4.2))
    for profile, label, color, width in zip(profiles, labels, colors, widths):
        ax.plot(x[keep], profile[keep], color=color, linewidth=width, label=label)

    ax.axvline(0, color="black", linestyle=(0, (5, 4)), linewidth=0.8)
    ax.set_xlim(-400, 400)
    ax.set_xticks([-400, 0, 400])
    ax.set_xticklabels(["-0.4 kb", "BG4 peak center", "+0.4 kb"])
    ax.set_ylabel(f"{target} CUT&Tag signal (RPKM)")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.tick_params(direction="out", length=3, width=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)
    fig.subplots_adjust(left=0.18, right=0.98, bottom=0.18, top=0.97)
    fig.savefig(output, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"{target}: {region_count} BG4-positive peaks; output={output}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    draw(args.matrix, args.output, args.target)


if __name__ == "__main__":
    main()
