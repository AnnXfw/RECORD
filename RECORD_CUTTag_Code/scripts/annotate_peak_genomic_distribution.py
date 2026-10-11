#!/usr/bin/env python3
"""注释全部蛋白峰和全部BG4峰，并绘制Nature Chemistry风格分组柱状图。"""

import argparse
import bisect
import gzip
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


CATEGORY_ORDER = [
    "-1 kb to TSS",
    "5' UTR",
    "Exon",
    "Intron",
    "3' UTR",
    "+1 kb to TES",
    "Intergenic",
]


def attributes(text):
    return dict(re.findall(r'(\S+)\s+"([^"]+)"', text))


def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [tuple(item) for item in merged]


def intersect_interval(interval, bounds):
    start = max(interval[0], bounds[0])
    end = min(interval[1], bounds[1])
    return (start, end) if end > start else None


def subtract_interval(interval, masks):
    pieces = [interval]
    for mstart, mend in masks:
        updated = []
        for start, end in pieces:
            if mend <= start or mstart >= end:
                updated.append((start, end))
            else:
                if mstart > start:
                    updated.append((start, mstart))
                if mend < end:
                    updated.append((mend, end))
        pieces = updated
    return pieces


def parse_gtf(path):
    transcripts = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attr_text = fields
            if not re.fullmatch(r"chr(?:[0-9]+|X|Y)", chrom):
                continue
            attr = attributes(attr_text)
            tid = attr.get("transcript_id")
            if not tid:
                continue
            # GTF为1-based闭区间，转换为0-based半开区间
            start0, end0 = int(start) - 1, int(end)
            tx = transcripts.setdefault(tid, {
                "chrom": chrom, "strand": strand,
                "start": start0, "end": end0,
                "exons": [], "cds": [],
            })
            tx["start"] = min(tx["start"], start0)
            tx["end"] = max(tx["end"], end0)
            if feature == "exon":
                tx["exons"].append((start0, end0))
            elif feature == "CDS":
                tx["cds"].append((start0, end0))
    return transcripts


def build_features(transcripts):
    features = defaultdict(lambda: defaultdict(list))
    for tx in transcripts.values():
        chrom, strand = tx["chrom"], tx["strand"]
        tx_start, tx_end = tx["start"], tx["end"]
        exons = merge_intervals(tx["exons"])
        cds = merge_intervals(tx["cds"])
        if not exons:
            continue

        if strand == "+":
            tss, tes = tx_start, tx_end
            promoter = (max(0, tss - 1000), tss)
            downstream = (tes, tes + 1000)
        else:
            tss, tes = tx_end, tx_start
            promoter = (tss, tss + 1000)
            downstream = (max(0, tes - 1000), tes)

        if promoter[1] > promoter[0]:
            features[chrom]["-1 kb to TSS"].append(promoter)
        if downstream[1] > downstream[0]:
            features[chrom]["+1 kb to TES"].append(downstream)

        utr5, utr3 = [], []
        if cds:
            cds_start = min(item[0] for item in cds)
            cds_end = max(item[1] for item in cds)
            for exon in exons:
                left = intersect_interval(exon, (0, cds_start))
                right = intersect_interval(exon, (cds_end, 10**12))
                if strand == "+":
                    if left:
                        utr5.append(left)
                    if right:
                        utr3.append(right)
                else:
                    if left:
                        utr3.append(left)
                    if right:
                        utr5.append(right)

        features[chrom]["5' UTR"].extend(utr5)
        features[chrom]["3' UTR"].extend(utr3)

        utr_masks = merge_intervals(utr5 + utr3)
        for exon in exons:
            features[chrom]["Exon"].extend(subtract_interval(exon, utr_masks))

        # 内含子为转录本范围扣除外显子后的部分
        features[chrom]["Intron"].extend(
            subtract_interval((tx_start, tx_end), exons)
        )

    for chrom in features:
        for category in list(features[chrom]):
            features[chrom][category] = merge_intervals(
                features[chrom][category]
            )
    return features


class IntervalLookup:
    def __init__(self, intervals):
        self.intervals = intervals
        self.starts = [item[0] for item in intervals]

    def contains(self, position):
        idx = bisect.bisect_right(self.starts, position) - 1
        return idx >= 0 and self.intervals[idx][1] > position


def build_lookup(features):
    return {
        chrom: {
            category: IntervalLookup(intervals)
            for category, intervals in categories.items()
        }
        for chrom, categories in features.items()
    }


def read_peak_intervals(path):
    peaks = []
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
            peaks.append((fields[0], start, end))
    return peaks


def peak_centers(peaks):
    return [(chrom, (start + end) // 2) for chrom, start, end in peaks]


def annotate(peaks, lookup):
    # UTR属于外显子子区，因此优先于普通Exon；启动子优先级最高
    priority = [
        "-1 kb to TSS", "5' UTR", "3' UTR",
        "Exon", "Intron", "+1 kb to TES",
    ]
    counts = {category: 0 for category in CATEGORY_ORDER}
    for chrom, center in peaks:
        assigned = "Intergenic"
        chromosome = lookup.get(chrom, {})
        for category in priority:
            checker = chromosome.get(category)
            if checker and checker.contains(center):
                assigned = category
                break
        counts[assigned] += 1
    return counts


def write_table(path, target, all_counts):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("target\tpeak_set\tcategory\tcount\tpercent\n")
        for peak_set, counts in all_counts.items():
            total = sum(counts.values())
            for category in CATEGORY_ORDER:
                count = counts[category]
                percent = 100.0 * count / total if total else 0
                handle.write(
                    f"{target}\t{peak_set}\t{category}\t"
                    f"{count}\t{percent:.4f}\n"
                )


def plot_single_svg(path, target, all_counts):
    x = np.arange(len(CATEGORY_ORDER))

    # 采用Nature Chemistry参考图的紧凑黑白风格
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 9,
        "axes.linewidth": 0.9,
        "hatch.linewidth": 0.8,
        "svg.fonttype": "none",
    })
    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    width = 0.34
    plot_sets = [
        (target, -width / 2, "black", "black", None),
        ("BG4", width / 2, "white", "black", None),
    ]
    for peak_set, offset, facecolor, edgecolor, hatch in plot_sets:
        counts = all_counts[peak_set]
        total = sum(counts.values())
        values = [
            100.0 * counts[category] / total if total else 0
            for category in CATEGORY_ORDER
        ]
        ax.bar(
            x + offset, values, width,
            label=peak_set,
            facecolor=facecolor,
            edgecolor=edgecolor,
            hatch=hatch,
            linewidth=0.8,
        )

    ax.set_ylabel("Proportion (%)", fontsize=10)
    ax.set_xlabel("")
    ax.set_xticks(x)
    ax.set_xticklabels(CATEGORY_ORDER, rotation=70, ha="right", va="top")
    ax.legend(
        frameon=False, loc="upper right", handlelength=1.5,
        borderaxespad=0.2, labelspacing=0.35,
    )
    # 两个独立面板固定相同纵轴，避免视觉比例误导
    ax.set_ylim(0, 100)
    ax.set_yticks(np.arange(0, 101, 20))
    ax.tick_params(direction="out", length=3, width=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)
    fig.subplots_adjust(left=0.14, right=0.98, bottom=0.34, top=0.97)
    fig.savefig(path, format="svg", bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--klf16-peaks", required=True)
    parser.add_argument("--hey1-peaks", required=True)
    parser.add_argument("--bg4-peaks", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    transcripts = parse_gtf(args.gtf)
    features = build_features(transcripts)
    lookup = build_lookup(features)

    protein_files = {
        "KLF16": args.klf16_peaks,
        "HEY1": args.hey1_peaks,
    }
    bg4_peaks = read_peak_intervals(args.bg4_peaks)
    bg4_counts = annotate(peak_centers(bg4_peaks), lookup)

    for target, path in protein_files.items():
        protein_peaks = read_peak_intervals(path)
        all_counts = {
            target: annotate(peak_centers(protein_peaks), lookup),
            "BG4": bg4_counts,
        }
        write_table(
            out / f"{target}.genomic_distribution.tsv",
            target, all_counts,
        )
        plot_single_svg(
            out / f"{target}.genomic_distribution.svg",
            target, all_counts,
        )
        print(
            target,
            "protein_peaks=", len(protein_peaks),
            "bg4_peaks=", len(bg4_peaks),
        )


if __name__ == "__main__":
    main()
