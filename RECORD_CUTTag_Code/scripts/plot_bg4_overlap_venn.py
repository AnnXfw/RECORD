#!/usr/bin/env python3
"""按Nature Chemistry参考风格绘制蛋白峰与BG4峰的重叠示意图。"""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


def read_bg4_result(path):
    with open(path, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["g4_set"] == "BG4":
                return {
                    "target": row["target"],
                    "peak_count": int(row["peak_count"]),
                    "overlap": int(row["observed_overlap"]),
                    "percent": float(row["overlap_percent"]),
                }
    raise ValueError(f"结果表中找不到BG4记录：{path}")


def count_bed(path):
    count = 0
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.strip() and not line.startswith(("#", "track", "browser")):
                count += 1
    return count


def draw(path, result, bg4_count):
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 10,
        "svg.fonttype": "none",
    })
    fig, ax = plt.subplots(figsize=(5.4, 4.1))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7)
    ax.set_aspect("equal")
    ax.axis("off")

    # 参考论文：蛋白峰为深灰色，BG4为蓝色，重叠数字放白色小框
    left = Circle(
        (4.15, 3.65), 2.25,
        facecolor="#4A4748", edgecolor="none", alpha=0.96,
    )
    right = Circle(
        (6.05, 3.65), 1.75,
        facecolor="#3F5798", edgecolor="none", alpha=0.93,
    )
    ax.add_patch(left)
    ax.add_patch(right)

    ax.text(
        1.65, 0.55,
        f"{result['target']} CUT&Tag\n({result['peak_count']:,})",
        ha="center", va="center", fontsize=10,
    )
    ax.text(
        8.25, 0.55,
        f"BG4 ChIP\n({bg4_count:,})",
        ha="center", va="center", fontsize=10,
    )
    ax.text(
        5.25, 3.65,
        f"{result['overlap']:,}\n({result['percent']:.1f}%)",
        ha="center", va="center", fontsize=9,
        bbox={
            "boxstyle": "square,pad=0.28",
            "facecolor": "white",
            "edgecolor": "none",
        },
    )

    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.03, top=0.98)
    fig.savefig(path, format="svg", transparent=True, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", required=True)
    parser.add_argument("--bg4-bed", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    result = read_bg4_result(args.summary)
    bg4_count = count_bed(args.bg4_bed)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    draw(args.output, result, bg4_count)
    print(
        f"{result['target']}: {result['overlap']}/{result['peak_count']} "
        f"({result['percent']:.1f}%) target peaks overlap "
        f"{bg4_count} published BG4 peaks"
    )


if __name__ == "__main__":
    main()
