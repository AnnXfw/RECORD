#!/usr/bin/env python3
"""从deepTools矩阵绘制Nature Chemistry风格的BG4/Input平均信号曲线。"""

import argparse
import gzip
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load_matrix(path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header_line = handle.readline().rstrip()
        if not header_line.startswith("@"):
            raise ValueError("不是有效的deepTools矩阵")
        header = json.loads(header_line[1:])
        boundaries = header["sample_boundaries"]
        labels = header["sample_labels"]
        rows = []
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) <= 6:
                continue
            values = []
            for item in fields[6:]:
                try:
                    value = float(item)
                except ValueError:
                    value = np.nan
                values.append(value)
            rows.append(values)

    matrix = np.asarray(rows, dtype=float)
    if matrix.size == 0:
        raise ValueError("矩阵中没有有效区域")
    profiles = []
    for start, end in zip(boundaries[:-1], boundaries[1:]):
        profiles.append(np.nanmean(matrix[:, start:end], axis=0))
    return labels, np.asarray(profiles), header


def draw(path, matrix_path, target):
    labels, profiles, header = load_matrix(matrix_path)
    expected = [
        "BG4_R1a", "BG4_R1c", "BG4_R2c", "BG4_R2b",
        "BG4_R2a", "Input_R1", "Input_R2",
    ]
    positions = {label: labels.index(label) for label in expected}

    bg4_rep1 = np.nanmean(
        profiles[[positions["BG4_R1a"], positions["BG4_R1c"]]], axis=0
    )
    bg4_rep2 = np.nanmean(
        profiles[[
            positions["BG4_R2c"],
            positions["BG4_R2b"],
            positions["BG4_R2a"],
        ]],
        axis=0,
    )
    input1 = profiles[positions["Input_R1"]]
    input2 = profiles[positions["Input_R2"]]

    upstream = header["upstream"][0]
    downstream = header["downstream"][0]
    x = np.linspace(-upstream, downstream, bg4_rep1.size, endpoint=False)

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 9,
        "axes.linewidth": 0.9,
        "svg.fonttype": "none",
    })
    fig, ax = plt.subplots(figsize=(5.2, 4.3))
    ax.plot(x, bg4_rep1, color="#202A87", linewidth=1.6, label="BG4 Rep 1")
    ax.plot(x, bg4_rep2, color="#D73027", linewidth=1.6, label="BG4 Rep 2")
    ax.plot(x, input1, color="#777777", linewidth=1.15, label="Input 1")
    ax.plot(x, input2, color="#BDBDBD", linewidth=1.15, label="Input 2")
    ax.axvline(0, color="black", linestyle=(0, (5, 4)), linewidth=0.8)

    ax.set_xlim(-upstream, downstream)
    ax.set_xticks([-2000, 0, 2000])
    ax.set_xticklabels(["-2 kb", "Peak center", "+2 kb"])
    ax.set_xlabel("")
    ax.set_ylabel("Published BG4 ChIP-seq signal")
    ax.legend(frameon=False, loc="upper right", fontsize=8)
    ax.tick_params(direction="out", length=3, width=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)
    fig.subplots_adjust(left=0.17, right=0.98, bottom=0.17, top=0.97)
    fig.savefig(path, format="svg", bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    draw(args.output, args.matrix, args.target)
    print(f"完成：{args.output}")


if __name__ == "__main__":
    main()
