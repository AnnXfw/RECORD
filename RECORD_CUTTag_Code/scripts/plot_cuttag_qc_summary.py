#!/usr/bin/env python3
"""根据QC汇总表绘制每个靶点的深度、比对率和FRiP图。"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_table(path):
    grouped = defaultdict(list)
    with open(path, encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            row["replicate"] = int(row["replicate"])
            row["mapped_reads"] = int(row["mapped_reads"])
            row["alignment_rate"] = float(row["alignment_rate"])
            row["frip_percent"] = float(row["frip_percent"])
            grouped[row["target"]].append(row)
    return grouped


def draw_target(path, target, rows):
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 9,
        "axes.linewidth": 0.9,
        "hatch.linewidth": 0.8,
        "svg.fonttype": "none",
    })
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.8))
    metrics = [
        ("mapped_reads", "Mapped reads (million)", 1e6),
        ("alignment_rate", "Overall alignment rate (%)", 1.0),
        ("frip_percent", "FRiP in final consensus peaks (%)", 1.0),
    ]
    reps = np.arange(1, 4)
    width = 0.34

    by_condition = {
        condition: sorted(
            [row for row in rows if row["condition"] == condition],
            key=lambda row: row["replicate"],
        )
        for condition in ["experiment", "control"]
    }
    if any(len(values) != 3 for values in by_condition.values()):
        raise ValueError(f"{target}实验或对照不是3个重复")

    for ax, (metric, ylabel, divisor) in zip(axes, metrics):
        experiment = [
            row[metric] / divisor for row in by_condition["experiment"]
        ]
        control = [
            row[metric] / divisor for row in by_condition["control"]
        ]
        ax.bar(
            reps - width / 2, experiment, width,
            label="Experiment", facecolor="black",
            edgecolor="black", linewidth=0.8,
        )
        ax.bar(
            reps + width / 2, control, width,
            label="Control", facecolor="white",
            edgecolor="black", hatch="////", linewidth=0.8,
        )
        ax.set_xticks(reps)
        ax.set_xticklabels(["Rep 1", "Rep 2", "Rep 3"])
        ax.set_ylabel(ylabel)
        ax.tick_params(direction="out", length=3, width=0.8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(False)
        if metric == "alignment_rate":
            lower = max(0, min(experiment + control) - 10)
            ax.set_ylim(lower, 100)
        else:
            ax.set_ylim(bottom=0)

    axes[2].legend(frameon=False, fontsize=8, loc="upper right")
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.18, top=0.96, wspace=0.40)
    fig.savefig(path, format="svg", bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--target", choices=["KLF16", "HEY1"])
    args = parser.parse_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    grouped = read_table(args.input)
    targets = [args.target] if args.target else ["KLF16", "HEY1"]
    for target in targets:
        draw_target(out / f"{target}.QC_summary.svg", target, grouped[target])
        print(f"{target}: {len(grouped[target])} samples")


if __name__ == "__main__":
    main()
