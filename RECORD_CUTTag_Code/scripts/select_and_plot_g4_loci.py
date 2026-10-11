#!/usr/bin/env python3
"""筛选KLF16/HEY1与BG4共定位的代表峰，并生成pyGenomeTracks配置。"""

import argparse
import configparser
import math
import os
import statistics
import subprocess
from pathlib import Path

import pyBigWig


def read_bed(path):
    records = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 3 or not fields[0].startswith("chr"):
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                continue
            if end > start:
                records.append((fields[0], start, end, fields))
    return records


def region_mean(bw, chrom, start, end):
    chroms = bw.chroms()
    if chrom not in chroms:
        return 0.0
    start = max(0, start)
    end = min(end, chroms[chrom])
    if end <= start:
        return 0.0
    value = bw.stats(chrom, start, end, type="mean", exact=False)[0]
    return 0.0 if value is None or not math.isfinite(value) else float(value)


def robust_z(values):
    median = statistics.median(values)
    deviations = [abs(v - median) for v in values]
    mad = statistics.median(deviations)
    if mad == 0:
        sd = statistics.pstdev(values)
        return [(v - median) / sd if sd else 0.0 for v in values]
    return [0.67448975 * (v - median) / mad for v in values]


def write_track_config(path, target, ratio_bws, bg4_bw, input_bw,
                       target_peaks, bg4_bed, potential_bed, gtf):
    colors = ["#253494", "#2C7FB8", "#41B6C4"]
    cfg = configparser.ConfigParser()
    cfg.optionxform = str

    for idx, (bw, color) in enumerate(zip(ratio_bws, colors), start=1):
        section = f"{target}_ratio_rep{idx}"
        cfg[section] = {
            "file": str(bw),
            "title": f"{target} Rep{idx}/Ctrl{idx} log2",
            "file_type": "bigwig",
            "color": color,
            "height": "2.2",
            "number_of_bins": "900",
            "nans_to_zeros": "true",
            "min_value": "-1",
            "summary_method": "mean",
        }

    # pyGenomeTracks要求特殊轨道名以“spacer”开头且后接空格编号
    cfg["spacer"] = {"height": "0.25"}
    cfg["published_bg4"] = {
        "file": str(bg4_bw),
        "title": "Published BG4 ChIP (K562)",
        "file_type": "bigwig",
        "color": "#6A3D9A",
        "height": "2.4",
        "number_of_bins": "900",
        "nans_to_zeros": "true",
        "min_value": "0",
        "summary_method": "mean",
    }
    cfg["published_input"] = {
        "file": str(input_bw),
        "title": "Published Input (K562)",
        "file_type": "bigwig",
        "color": "#9E9E9E",
        "height": "1.8",
        "number_of_bins": "900",
        "nans_to_zeros": "true",
        "min_value": "0",
        "summary_method": "mean",
    }

    cfg["target_peaks"] = {
        "file": str(target_peaks),
        "title": f"{target} enriched peaks",
        "file_type": "bed",
        "color": "#111111",
        "height": "0.8",
        "display": "collapsed",
        "labels": "false",
    }
    cfg["bg4_peaks"] = {
        "file": str(bg4_bed),
        "title": "High-confidence BG4 peaks",
        "file_type": "bed",
        "color": "#7B3294",
        "height": "0.8",
        "display": "collapsed",
        "labels": "false",
    }
    cfg["potential_g4"] = {
        "file": str(potential_bed),
        "title": "Potential G4",
        "file_type": "bed",
        "color": "#BDBDBD",
        "height": "0.8",
        "display": "collapsed",
        "labels": "false",
    }
    cfg["genes"] = {
        "file": str(gtf),
        "title": "RefSeq genes (hg19)",
        "file_type": "gtf",
        "prefered_name": "gene_name",
        "merge_transcripts": "true",
        "color": "#252525",
        "height": "3.0",
        "fontsize": "9",
        "labels": "true",
        "style": "UCSC",
    }
    cfg["x-axis"] = {"where": "bottom"}

    with open(path, "w", encoding="utf-8") as handle:
        cfg.write(handle, space_around_delimiters=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", required=True, choices=["KLF16", "HEY1"])
    parser.add_argument("--root", required=True)
    parser.add_argument("--top-n", type=int, default=3)
    parser.add_argument("--window", type=int, default=20000)
    args = parser.parse_args()

    root = Path(args.root)
    results = root / "results"
    reference = root / "reference"
    bg4_dir = reference / "bg4_k562_gse107690"

    if args.target == "KLF16":
        label = "KLF16_strict_2of3"
        target_peaks = results / "consensus/KLF16/KLF16.consensus.2of3.bed"
    else:
        label = "HEY1_pooled_default_2of3"
        target_peaks = (
            results
            / "diagnostic_peak_calling/HEY1/default/HEY1.default.consensus.2of3.bed"
        )

    candidates = (
        results / f"g4_enrichment/{label}/{label}_with_BG4.bed"
    )
    ratio_bws = [
        results / f"g4_signal_redraw/{label}/{args.target}.Rep{i}_vs_Ctrl{i}.log2ratio.bw"
        for i in range(1, 4)
    ]
    bg4_bw = (
        bg4_dir
        / "GSM2876090_B_REP1.SLX-12319.K562_asynch_a_701_504.clean.bw"
    )
    input_bw = (
        bg4_dir
        / "GSM2876092_B_REP1.SLX-12319.K562_asynch_IP_704_504.clean.bw"
    )
    bg4_bed = reference / "G4.bed"
    potential_bed = reference / "Potential_G4.bed"
    gtf = reference / "hg19.refGene.gtf.gz"
    out = results / f"locus_tracks/{args.target}"
    out.mkdir(parents=True, exist_ok=True)

    required = [
        candidates, target_peaks, *ratio_bws, bg4_bw, input_bw,
        bg4_bed, potential_bed, gtf,
    ]
    for path in required:
        if not path.exists() or path.stat().st_size == 0:
            raise SystemExit(f"缺少或为空：{path}")

    records = read_bed(candidates)
    if not records:
        raise SystemExit(f"没有可用候选峰：{candidates}")

    ratio_handles = [pyBigWig.open(str(path)) for path in ratio_bws]
    bg4_handle = pyBigWig.open(str(bg4_bw))

    scored = []
    for chrom, start, end, fields in records:
        center = (start + end) // 2
        score_start, score_end = center - 250, center + 250
        ratios = [
            region_mean(handle, chrom, score_start, score_end)
            for handle in ratio_handles
        ]
        bg4_signal = region_mean(
            bg4_handle, chrom, score_start, score_end
        )
        scored.append({
            "chrom": chrom,
            "start": start,
            "end": end,
            "center": center,
            "ratio_mean": statistics.mean(ratios),
            "ratio_reps": ratios,
            "bg4_signal": bg4_signal,
            "fields": fields,
        })

    for handle in ratio_handles:
        handle.close()
    bg4_handle.close()

    ratio_z = robust_z([row["ratio_mean"] for row in scored])
    bg4_z = robust_z([math.log1p(row["bg4_signal"]) for row in scored])
    for row, rz, bz in zip(scored, ratio_z, bg4_z):
        # 两类证据等权，且略微惩罚重复间方向不一致
        consistency = sum(value > 0 for value in row["ratio_reps"]) / 3.0
        row["combined_score"] = rz + bz + consistency

    scored.sort(key=lambda row: row["combined_score"], reverse=True)

    # 避免三张图选到同一局部区域；中心至少相隔100 kb
    selected = []
    for row in scored:
        if all(
            row["chrom"] != old["chrom"]
            or abs(row["center"] - old["center"]) >= 100000
            for old in selected
        ):
            selected.append(row)
        if len(selected) == args.top_n:
            break

    config_path = out / f"{args.target}.tracks.ini"
    write_track_config(
        config_path, args.target, ratio_bws, bg4_bw, input_bw,
        target_peaks, bg4_bed, potential_bed, gtf,
    )

    summary_path = out / f"{args.target}.selected_loci.tsv"
    with open(summary_path, "w", encoding="utf-8") as handle:
        handle.write(
            "rank\tchrom\tpeak_start\tpeak_end\tplot_start\tplot_end\t"
            "rep1_log2\trep2_log2\trep3_log2\tmean_log2\t"
            "BG4_signal\tcombined_score\tSVG\n"
        )
        for rank, row in enumerate(selected, start=1):
            plot_start = max(0, row["center"] - args.window)
            plot_end = row["center"] + args.window
            svg = out / f"{args.target}.locus{rank}.{row['chrom']}_{plot_start}_{plot_end}.svg"
            values = [
                rank, row["chrom"], row["start"], row["end"],
                plot_start, plot_end,
                *[f"{v:.6g}" for v in row["ratio_reps"]],
                f"{row['ratio_mean']:.6g}",
                f"{row['bg4_signal']:.6g}",
                f"{row['combined_score']:.6g}",
                svg.name,
            ]
            handle.write("\t".join(map(str, values)) + "\n")

            subprocess.run([
                "pyGenomeTracks",
                "--tracks", str(config_path),
                "--region", f"{row['chrom']}:{plot_start}-{plot_end}",
                "--outFileName", str(svg),
                "--title", f"{args.target} / published BG4 co-localized locus {rank}",
                "--width", "38",
                "--trackLabelFraction", "0.18",
                "--dpi", "300",
            ], check=True)

    print(f"已生成{len(selected)}个代表位点：{out}")


if __name__ == "__main__":
    main()
