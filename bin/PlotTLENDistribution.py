#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Plot insert-size distribution and write TLEN summary statistics."""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot sampled positive TLEN values and summarize insert sizes."
    )
    parser.add_argument("--Input", required=True, help="Input TLEN TSV (absolute path).")
    parser.add_argument("--Output", required=True, help="Output PDF path (absolute path).")
    parser.add_argument("--Summary", required=True, help="Output summary TSV path (absolute path).")
    parser.add_argument("--Sample", required=True, help="Sample identifier.")
    return parser.parse_args()


def calculate_mad(values, median_value):
    return float(np.median(np.abs(values - median_value)))


def require_absolute_path(path, label):
    if not os.path.isabs(path):
        raise ValueError(f"[ERROR] {label} path must be absolute: {path}")


def format_value(value):
    if value is None or not np.isfinite(value):
        return "NA"
    return f"{value:.3f}"


def write_summary(path, sample, stats):
    columns = [
        "SampleID",
        "Total_Read_Pairs",
        "Mean_Insert_Size",
        "SD_Insert_Size",
        "Median_Insert_Size",
        "MAD",
        "MAD_Scale",
        "Robust_SD",
        "N_MAD",
        "Abnormal_TLEN_Threshold",
        "Q1",
        "Q3",
        "P95",
        "P99",
        "Plot_Upper_P99_5",
    ]

    row = [
        sample,
        str(stats["count"]),
        format_value(stats.get("mean")),
        format_value(stats.get("sd")),
        format_value(stats.get("median")),
        format_value(stats.get("mad")),
        format_value(stats.get("mad_scale")),
        format_value(stats.get("robust_sd")),
        str(stats.get("n_mad", 3)),
        format_value(stats.get("threshold")),
        format_value(stats.get("q1")),
        format_value(stats.get("q3")),
        format_value(stats.get("p95")),
        format_value(stats.get("p99")),
        format_value(stats.get("plot_upper")),
    ]

    with open(path, "w", encoding="utf-8") as out:
        out.write("\t".join(columns) + "\n")
        out.write("\t".join(row) + "\n")


def main():
    args = parse_args()
    for path, label in (
        (args.Input, "Input"),
        (args.Output, "Output"),
        (args.Summary, "Summary"),
    ):
        require_absolute_path(path, label)

    if not os.path.isfile(args.Input):
        raise FileNotFoundError(f"[ERROR] Input file does not exist: {args.Input}")

    df = pd.read_csv(args.Input, sep="\t")
    if "TLEN" not in df.columns:
        raise ValueError("[ERROR] Input file must contain a column named 'TLEN'.")

    values = pd.to_numeric(df["TLEN"], errors="coerce").dropna().to_numpy(dtype=float)
    values = values[values > 0]

    os.makedirs(os.path.dirname(args.Output), exist_ok=True)
    os.makedirs(os.path.dirname(args.Summary), exist_ok=True)

    if values.size == 0:
        stats = {"count": 0, "n_mad": 3}
        fig, ax = plt.subplots(figsize=(9, 5.5))
        ax.text(0.5, 0.5, "No valid positive TLEN values", ha="center", va="center",
                transform=ax.transAxes)
        ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(args.Output, dpi=300)
        plt.close(fig)
        write_summary(args.Summary, args.Sample, stats)
        print("[WARNING] No valid positive TLEN values were found.")
        print(f"[INFO] Output figure: {args.Output}")
        print(f"[INFO] Summary file: {args.Summary}")
        return

    count = int(values.size)
    mean_value = float(np.mean(values))
    sd_value = float(np.std(values, ddof=1)) if count > 1 else 0.0
    median_value = float(np.median(values))
    mad = calculate_mad(values, median_value)
    mad_scale = 1.4826
    robust_sd = mad_scale * mad
    n_mad = 3
    threshold = median_value + n_mad * robust_sd
    q1, q3, p95, p99, plot_upper = [
        float(x) for x in np.percentile(values, [25, 75, 95, 99, 99.5])
    ]
    values_for_plot = values[values <= plot_upper]

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.hist(values_for_plot, bins=100, alpha=0.85, edgecolor="black", linewidth=0.3)
    ax.axvline(median_value, linestyle="--", linewidth=1.5,
               label=f"Median = {median_value:.1f} bp")
    ax.axvline(threshold, linestyle="--", linewidth=1.5,
               label=f"Diagnostic threshold = {threshold:.1f} bp")
    ax.set_xlabel("Insert Size / TLEN (bp)")
    ax.set_ylabel("Read pair count")
    ax.set_title(
        f"{args.Sample} | N = {count:,} | Median = {median_value:.1f} bp | "
        f"MAD = {mad:.1f} | Robust SD = {robust_sd:.1f}",
        fontsize=9,
        pad=6,
    )
    ax.legend(loc="upper right", frameon=True)
    fig.tight_layout()
    fig.savefig(args.Output, dpi=300)
    plt.close(fig)

    stats = {
        "count": count,
        "mean": mean_value,
        "sd": sd_value,
        "median": median_value,
        "mad": mad,
        "mad_scale": mad_scale,
        "robust_sd": robust_sd,
        "n_mad": n_mad,
        "threshold": threshold,
        "q1": q1,
        "q3": q3,
        "p95": p95,
        "p99": p99,
        "plot_upper": plot_upper,
    }
    write_summary(args.Summary, args.Sample, stats)

    print("[INFO] TLEN distribution analysis completed.")
    print(f"[INFO] Sample: {args.Sample}")
    print(f"[INFO] Total read pairs: {count:,}")
    print(f"[INFO] Median insert size: {median_value:.3f}")
    print(f"[INFO] Diagnostic TLEN threshold: {threshold:.3f}")
    print(f"[INFO] Output figure: {args.Output}")
    print(f"[INFO] Summary file: {args.Summary}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)
