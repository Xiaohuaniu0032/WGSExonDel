#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
HCMExonDel multi-track deletion evidence plot
Chinese color-rendered version

Panels:
A. Gene and exon structure
B. Split-read evidence
C. Base-level depth ratios and CUSUM-supported intervals
D. Discordant read-pair evidence
E. Integrated HCMExonDel results

Example:

python PlotHCMExonDelMultiTrack.py \
    /path/to/HCMExonDel/results/Sample_01 \
    --font /path/to/fonts/NotoSansCJKsc-Regular.otf
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import font_manager

from matplotlib.patches import Rectangle, PathPatch, Patch
from matplotlib.path import Path as MplPath
from matplotlib.ticker import FuncFormatter


# ============================================================
# Default parameters
# ============================================================

DEFAULT_GTF = None


# ============================================================
# Color palette
# ============================================================

COLORS = {

    # Text and axes
    "text": "#222222",
    "axis": "#444444",
    "gene_line": "#666666",

    # Candidate deletion interval
    "candidate_fill": "#EAF2F8",
    "candidate_edge": "#8FA3B5",

    # Panel A
    "affected_exon": "#2F5D8A",

    # Panel B
    "split": "#2A9D8F",
    "split_dark": "#1F7F73",

    # Panel C
    "depth": "#264653",
    "ref_line": "#8A8A8A",
    "threshold_line": "#A5A5A5",

    # Panel D
    "discordant": "#7A5195",
    "discordant_dark": "#5B3F78",

    # Panel E
    "high": "#2A9D8F",
    "moderate": "#D99036",
    "low": "#B56576",
}


# Colors for multiple CUSUM intervals
# Use darker colors to keep white interval labels readable.
INTERVAL_COLORS = [
    "#C94F3D",
    "#D77A2B",
    "#A77C20",
    "#4F8060",
    "#7562A7",
    "#457B8C",
]


# ============================================================
# Matplotlib default style
# ============================================================

plt.rcParams["font.size"] = 11
plt.rcParams["axes.titlesize"] = 11
plt.rcParams["axes.labelsize"] = 11

plt.rcParams["xtick.labelsize"] = 10
plt.rcParams["ytick.labelsize"] = 10

plt.rcParams["legend.fontsize"] = 9

plt.rcParams["axes.linewidth"] = 0.9

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

plt.rcParams["axes.unicode_minus"] = False


# ============================================================
# General utility functions
# ============================================================

def die(msg):

    print(
        f"[ERROR] {msg}",
        file=sys.stderr
    )

    sys.exit(1)


def info(msg):

    print(
        f"[INFO] {msg}",
        file=sys.stderr
    )


# ============================================================
# Chinese font setup
# ============================================================

def setup_chinese_font(font_arg=None):
    """
    Configure a Chinese font.

    Providing a font file is recommended:

    --font /path/to/NotoSansCJKsc-Regular.otf
    """

    if font_arg:

        font_path = Path(
            font_arg
        )

        if font_path.exists():

            font_manager.fontManager.addfont(
                str(font_path)
            )

            prop = font_manager.FontProperties(
                fname=str(font_path)
            )

            font_name = prop.get_name()

            plt.rcParams["font.family"] = (
                "sans-serif"
            )

            plt.rcParams[
                "font.sans-serif"
            ] = [
                font_name
            ]

            plt.rcParams[
                "axes.unicode_minus"
            ] = False

            plt.rcParams[
                "pdf.fonttype"
            ] = 42

            plt.rcParams[
                "ps.fonttype"
            ] = 42

            # Convert SVG text to paths
            # to avoid missing-font issues on other systems.
            plt.rcParams[
                "svg.fonttype"
            ] = "path"

            info(
                f"使用中文字体文件: "
                f"{font_path}"
            )

            info(
                f"Matplotlib 字体名称: "
                f"{font_name}"
            )

            return prop

        # If the value is not a file path, treat it as a font family name.
        installed_names = {
            f.name
            for f in
            font_manager.fontManager.ttflist
        }

        if font_arg in installed_names:

            plt.rcParams[
                "font.family"
            ] = "sans-serif"

            plt.rcParams[
                "font.sans-serif"
            ] = [
                font_arg
            ]

            plt.rcParams[
                "axes.unicode_minus"
            ] = False

            info(
                f"使用中文字体: "
                f"{font_arg}"
            )

            return None

        die(
            "指定的中文字体不存在: "
            f"{font_arg}"
        )

    # --------------------------------------------------------
    # Search common Chinese fonts automatically.
    # --------------------------------------------------------

    preferred = [
        "Noto Sans CJK SC",
        "Source Han Sans SC",
        "Source Han Sans CN",
        "Microsoft YaHei",
        "SimHei",
        "WenQuanYi Micro Hei",
        "Noto Sans CJK JP",
    ]

    installed_names = {
        f.name
        for f in
        font_manager.fontManager.ttflist
    }

    for name in preferred:

        if name in installed_names:

            plt.rcParams[
                "font.family"
            ] = "sans-serif"

            plt.rcParams[
                "font.sans-serif"
            ] = [
                name
            ]

            info(
                f"自动选择中文字体: "
                f"{name}"
            )

            return None

    die(
        "没有找到可用中文字体。\n"
        "请通过 --font 指定中文字体文件，例如：\n"
        "--font "
        "/path/to/fonts/NotoSansCJKsc-Regular.otf"
    )


# ============================================================
# String utilities
# ============================================================

def strip_version(tx):

    if pd.isna(tx):
        return ""

    return re.sub(
        r"\.\d+$",
        "",
        str(tx)
    )


def format_coord(x):

    return f"{int(x):,}"


# ============================================================
# File discovery
# ============================================================

def find_exact_or_recursive(
    root,
    exact_path,
    patterns
):

    p = (
        root
        / exact_path
    )

    if p.exists():
        return p

    hits = []

    for pattern in patterns:

        hits.extend(
            root.rglob(
                pattern
            )
        )

    hits = sorted(
        set(hits)
    )

    if len(hits) == 0:

        return None

    if len(hits) > 1:

        info(
            "发现多个匹配文件:"
        )

        for x in hits:

            info(
                f"    {x}"
            )

        info(
            f"使用: {hits[0]}"
        )

    return hits[0]


def read_tsv(path):

    if path is None:

        return None

    return pd.read_csv(
        path,
        sep="\t",
        low_memory=False
    )


# ============================================================
# GTF
# ============================================================

def parse_gtf_attributes(attr):

    result = {}

    for item in attr.strip().split(";"):

        item = item.strip()

        if not item:
            continue

        m = re.match(
            r'(\S+)\s+"([^"]+)"',
            item
        )

        if m:

            result[
                m.group(1)
            ] = m.group(2)

    return result


def parse_gtf_exons(
    gtf_file,
    gene,
    transcript
):

    tx_base = strip_version(
        transcript
    )

    exons = []
    transcript_rows = []

    with open(
        gtf_file,
        "r",
        encoding="utf-8"
    ) as fh:

        for line in fh:

            if line.startswith("#"):
                continue

            fields = (
                line.rstrip("\n")
                .split("\t")
            )

            if len(fields) < 9:
                continue

            (
                chrom,
                source,
                feature,
                start,
                end,
                score,
                strand,
                frame,
                attr
            ) = fields

            if feature not in {
                "transcript",
                "exon"
            }:
                continue

            a = parse_gtf_attributes(
                attr
            )

            gene_name = a.get(
                "gene_name",
                a.get(
                    "gene_id",
                    ""
                )
            )

            if gene_name != gene:
                continue

            tid = a.get(
                "transcript_id",
                ""
            )

            if (
                tx_base
                and
                strip_version(tid)
                != tx_base
            ):
                continue

            rec = {

                "chrom": chrom,

                "start": int(
                    start
                ),

                "end": int(
                    end
                ),

                "strand": strand,

                "transcript_id": tid,

                "exon_number": a.get(
                    "exon_number",
                    ""
                ),

                "exon_id": a.get(
                    "exon_id",
                    ""
                ),
            }

            if feature == "transcript":

                transcript_rows.append(
                    rec
                )

            elif feature == "exon":

                exons.append(
                    rec
                )

    if not exons:

        die(
            "GTF 中未找到目标转录本"
            "外显子注释: "
            f"{gene}/{transcript}"
        )

    exon_df = (
        pd.DataFrame(
            exons
        )
        .sort_values(
            [
                "start",
                "end"
            ]
        )
        .reset_index(
            drop=True
        )
    )

    actual_tx = (
        exon_df.iloc[0][
            "transcript_id"
        ]
    )

    if transcript_rows:

        tr = (
            transcript_rows[0]
        )

        tx_start = tr[
            "start"
        ]

        tx_end = tr[
            "end"
        ]

        strand = tr[
            "strand"
        ]

        chrom = tr[
            "chrom"
        ]

    else:

        tx_start = int(
            exon_df[
                "start"
            ].min()
        )

        tx_end = int(
            exon_df[
                "end"
            ].max()
        )

        strand = (
            exon_df.iloc[0][
                "strand"
            ]
        )

        chrom = (
            exon_df.iloc[0][
                "chrom"
            ]
        )

    transcript_info = {

        "chrom": chrom,

        "start": tx_start,

        "end": tx_end,

        "strand": strand,

        "transcript_id": (
            actual_tx
        ),
    }

    return (
        exon_df,
        transcript_info
    )


def exon_number_from_annotation(
    text
):

    if pd.isna(
        text
    ):

        return set()

    return set(
        re.findall(
            r"EX(\d+)",
            str(text)
        )
    )


# ============================================================
# Evidence
# ============================================================

def evidence_rank(x):

    rank = {
        "High": 3,
        "Moderate": 2,
        "Low": 1,
    }

    return rank.get(
        str(x),
        0
    )


def evidence_level_cn(
    level
):

    mapping = {
        "High": "高",
        "Moderate": "中",
        "Low": "低",
    }

    return mapping.get(
        str(level),
        str(level)
    )


# ============================================================
# HCMExonDel result-file layout
# ============================================================

def discover_files(
    input_dir,
    sample,
    gene
):

    files = {}

    files["candidate"] = (
        find_exact_or_recursive(
            input_dir,

            Path(
                "S6_Report"
            )
            / "annotated_candidates.tsv",

            [
                "annotated_candidates.tsv"
            ]
        )
    )

    files["merged_candidate"] = (
        find_exact_or_recursive(
            input_dir,
            Path("S5_Candidates") / "merged_candidates.tsv",
            ["merged_candidates.tsv"]
        )
    )

    files["split_cluster"] = (
        find_exact_or_recursive(
            input_dir,

            Path(
                "S3_SplitReads"
            )
            / "split_clusters.tsv",

            [
                "split_clusters.tsv"
            ]
        )
    )

    files["split_events"] = (
        find_exact_or_recursive(
            input_dir,
            Path("S3_SplitReads") / "split_deletion_events.tsv",
            ["split_deletion_events.tsv"]
        )
    )

    files["split_assignments"] = (
        find_exact_or_recursive(
            input_dir,
            Path("S3_SplitReads") / "split_cluster_assignments.tsv",
            ["split_cluster_assignments.tsv"]
        )
    )

    files["dr_cluster"] = (
        find_exact_or_recursive(
            input_dir,

            Path(
                "S4_DiscordantReads"
            )
            / "discordant_clusters.tsv",

            ["discordant_clusters.tsv"]
        )
    )

    files["dr_pairs"] = (
        find_exact_or_recursive(
            input_dir,

            Path(
                "S4_DiscordantReads"
            )
            / "discordant_pairs.tsv",

            ["discordant_pairs.tsv"]
        )
    )

    files["dr_assignments"] = (
        find_exact_or_recursive(
            input_dir,
            Path("S4_DiscordantReads") / "discordant_cluster_assignments.tsv",
            ["discordant_cluster_assignments.tsv"]
        )
    )

    depth_root = (
        input_dir
        / "S2_Depth"
    )

    depth_hits = sorted(
        depth_root.rglob(
            f"{sample}_{gene}_"
            "*.base_depth_ratio.tsv"
        )
    )

    cusum_hits = sorted(
        depth_root.rglob(
            f"{sample}_{gene}_"
            "*.cusum.tsv"
        )
    )

    files[
        "depth_ratio"
    ] = (
        depth_hits[0]
        if depth_hits
        else None
    )

    files[
        "cusum"
    ] = (
        cusum_hits[0]
        if cusum_hits
        else None
    )

    return files


# ============================================================
# Candidate
# ============================================================

def select_candidate(
    candidate_df,
    gene,
    candidate_id=None
):

    df = candidate_df.copy()

    if candidate_id:
        tmp = df[
            df["Candidate_ID"].astype(str).eq(candidate_id)
        ].copy()
        if tmp.empty:
            die(f"未找到候选事件 Candidate_ID: {candidate_id}")
        return tmp.iloc[0]

    if "Gene" not in df.columns:
        die("最终报告缺少 Gene 字段")

    gene_mask = df["Gene"].fillna("").astype(str).apply(
        lambda value: gene in [item.strip() for item in value.split(",")]
    )
    tmp = df[gene_mask].copy()

    if tmp.empty:
        die(f"未找到注释到基因 {gene} 的候选事件")

    tmp["_rank"] = tmp["Evidence_Level"].map(evidence_rank)
    sort_cols = ["_rank"]
    if "Affected_Exon_Count" in tmp.columns:
        sort_cols.append("Affected_Exon_Count")
    if "Split_Reads" in tmp.columns:
        sort_cols.append("Split_Reads")
    if "Size_bp" in tmp.columns:
        sort_cols.append("Size_bp")

    return tmp.sort_values(
        sort_cols,
        ascending=[False] * len(sort_cols)
    ).iloc[0]


# ============================================================
# Shared plotting elements
# ============================================================

def add_candidate_span(
    ax,
    start,
    end
):

    ax.axvspan(
        start,
        end,

        facecolor=COLORS[
            "candidate_fill"
        ],

        alpha=0.62,

        zorder=0
    )


def add_candidate_boundaries(
    ax,
    start,
    end
):

    ax.axvline(
        start,

        color=COLORS[
            "candidate_edge"
        ],

        linestyle="--",

        linewidth=0.9,

        zorder=1
    )

    ax.axvline(
        end,

        color=COLORS[
            "candidate_edge"
        ],

        linestyle="--",

        linewidth=0.9,

        zorder=1
    )


def format_genomic_axis(
    ax
):

    ax.xaxis.set_major_formatter(
        FuncFormatter(
            lambda x, pos:
            f"{int(x):,}"
        )
    )


# ============================================================
# Panel A
# ============================================================

def draw_gene_track(
    ax,
    exon_df,
    transcript_info,
    affected_exons,
    x_min,
    x_max,
    gene
):

    ax.set_xlim(
        x_min,
        x_max
    )

    ax.set_ylim(
        0,
        1
    )

    y = 0.43

    strand_cn = (
        "正链"
        if (
            transcript_info[
                "strand"
            ] == "+"
        )
        else "负链"
    )

    # --------------------------------------------------------
    # Transcript information
    # Place the label at the upper left to avoid overlap with exon labels.
    # --------------------------------------------------------

    ax.text(
        0.01,
        0.96,

        (
            f"{gene} / "
            f"{transcript_info['transcript_id']} "
            f"（{strand_cn}）"
        ),

        transform=ax.transAxes,

        ha="left",
        va="top",

        fontsize=11,

        color=COLORS[
            "text"
        ]
    )

    # --------------------------------------------------------
    # Gene backbone
    # --------------------------------------------------------

    ax.plot(
        [
            max(
                transcript_info[
                    "start"
                ],
                x_min
            ),
            min(
                transcript_info[
                    "end"
                ],
                x_max
            )
        ],

        [
            y,
            y
        ],

        color=COLORS[
            "gene_line"
        ],

        linewidth=1.2
    )

    # --------------------------------------------------------
    # Exons within the displayed interval
    # --------------------------------------------------------

    local = exon_df[
        (
            exon_df[
                "end"
            ] >= x_min
        )
        &
        (
            exon_df[
                "start"
            ] <= x_max
        )
    ].copy()

    for _, row in (
        local.iterrows()
    ):

        start = max(
            int(
                row[
                    "start"
                ]
            ),
            x_min
        )

        end = min(
            int(
                row[
                    "end"
                ]
            ),
            x_max
        )

        ex_num = str(
            row[
                "exon_number"
            ]
        )

        is_affected = (
            ex_num
            in affected_exons
        )

        rect = Rectangle(
            (
                start,
                y - 0.14
            ),

            max(
                1,
                end - start + 1
            ),

            0.28,

            facecolor=(
                COLORS[
                    "affected_exon"
                ]
                if is_affected
                else "white"
            ),

            edgecolor=(
                COLORS[
                    "affected_exon"
                ]
                if is_affected
                else COLORS[
                    "text"
                ]
            ),

            linewidth=1.2
        )

        ax.add_patch(
            rect
        )

        ax.text(
            (
                start
                + end
            )
            / 2,

            y + 0.19,

            f"EX{ex_num}",

            ha="center",
            va="bottom",

            fontsize=9,

            fontweight="bold",

            color=COLORS[
                "text"
            ]
        )

    # --------------------------------------------------------
    # Transcript direction
    # --------------------------------------------------------

    if (
        transcript_info[
            "strand"
        ] == "+"
    ):

        ax.annotate(
            "",

            xy=(
                x_max,
                0.13
            ),

            xytext=(
                x_min,
                0.13
            ),

            arrowprops=dict(
                arrowstyle="->",
                linewidth=1.0,
                color=COLORS[
                    "axis"
                ]
            )
        )

    else:

        ax.annotate(
            "",

            xy=(
                x_min,
                0.13
            ),

            xytext=(
                x_max,
                0.13
            ),

            arrowprops=dict(
                arrowstyle="->",
                linewidth=1.0,
                color=COLORS[
                    "axis"
                ]
            )
        )

    ax.set_ylabel(
        "基因"
    )

    ax.set_yticks([])

    ax.set_title(
        "A  基因及外显子结构",

        loc="left",

        fontsize=12,

        fontweight="bold",

        color=COLORS[
            "text"
        ]
    )

    ax.spines[
        "left"
    ].set_visible(
        False
    )

    ax.spines[
        "right"
    ].set_visible(
        False
    )

    ax.spines[
        "top"
    ].set_visible(
        False
    )


# ============================================================
# Panel B
# ============================================================

def draw_split_reads(
    ax,
    split_cluster_row,
    split_raw_df,
    x_min,
    x_max
):

    cluster_id = str(
        split_cluster_row[
            "Cluster_ID"
        ]
    )

    df = split_raw_df[
        split_raw_df[
            "Cluster_ID"
        ]
        .astype(str)
        .eq(cluster_id)
    ].copy()

    df = (
        df.reset_index(
            drop=True
        )
    )

    n = len(
        df
    )

    if n == 0:

        die(
            "未找到所选 Cluster "
            "对应的分裂读段"
        )

    # --------------------------------------------------------
    # Retain all supporting reads.
    # --------------------------------------------------------

    for idx, (_, row) in enumerate(
        df.iterrows(),
        start=1
    ):

        y = (
            n
            - idx
            + 1
        )

        p1 = (
            int(
                row[
                    "Primary_Start"
                ]
            ),
            int(
                row[
                    "Primary_End"
                ]
            )
        )

        p2 = (
            int(
                row[
                    "SA_Start"
                ]
            ),
            int(
                row[
                    "SA_End"
                ]
            )
        )

        blocks = sorted(
            [
                p1,
                p2
            ],
            key=lambda x:
            x[0]
        )

        left = (
            blocks[0]
        )

        right = (
            blocks[1]
        )

        # Connect the deletion gap.
        ax.plot(
            [
                left[1],
                right[0]
            ],

            [
                y,
                y
            ],

            color=COLORS[
                "split"
            ],

            linewidth=1.0,

            alpha=0.95
        )

        # Alignment blocks on both sides.
        for s, e in blocks:

            ax.add_patch(
                Rectangle(
                    (
                        s,
                        y - 0.22
                    ),

                    max(
                        1,
                        e - s + 1
                    ),

                    0.44,

                    facecolor=COLORS[
                        "split_dark"
                    ],

                    edgecolor=COLORS[
                        "split_dark"
                    ],

                    linewidth=0.5
                )
            )

    ax.set_xlim(
        x_min,
        x_max
    )

    ax.set_ylim(
        0,
        n + 2
    )

    ax.set_ylabel(
        "分裂读段"
    )

    ax.set_yticks([])

    ax.set_title(
        (
            "B  分裂读段证据："
            f"{cluster_id}，"
            "支持数 = "
            f"{int(split_cluster_row['Support_Reads'])}"
        ),

        loc="left",

        fontsize=12,

        fontweight="bold",

        color=COLORS[
            "text"
        ]
    )

    ax.spines[
        "right"
    ].set_visible(
        False
    )

    ax.spines[
        "top"
    ].set_visible(
        False
    )


# ============================================================
# Panel C
# ============================================================

def draw_depth_track(
    ax,
    depth_df,
    cusum_df,
    candidate_start,
    candidate_end,
    x_min,
    x_max
):

    # --------------------------------------------------------
    # Raw base-level depth ratios
    # --------------------------------------------------------

    df = depth_df[
        (
            depth_df[
                "POS"
            ] >= x_min
        )
        &
        (
            depth_df[
                "POS"
            ] <= x_max
        )
    ].copy()

    if df.empty:

        die(
            "绘图区间内没有"
            "深度比值数据"
        )

    ratio_col = (
        "Ratio_By_Median"
    )

    if (
        ratio_col
        not in df.columns
    ):

        die(
            f"{ratio_col} 不存在于 "
            "base_depth_ratio.tsv"
        )

    # --------------------------------------------------------
    # Use the original values directly.
    # Do not apply a rolling median.
    # --------------------------------------------------------

    ax.plot(
        df[
            "POS"
        ],

        df[
            ratio_col
        ],

        color=COLORS[
            "depth"
        ],

        linewidth=0.85,

        zorder=3
    )

    # --------------------------------------------------------
    # Ratio = 1.0
    # --------------------------------------------------------

    ax.axhline(
        1.0,

        color=COLORS[
            "ref_line"
        ],

        linewidth=0.9,

        linestyle="-",

        zorder=1
    )

    # --------------------------------------------------------
    # Ratio = 0.65
    # --------------------------------------------------------

    ax.axhline(
        0.65,

        color=COLORS[
            "threshold_line"
        ],

        linewidth=1.0,

        linestyle="--",

        zorder=1
    )

    # --------------------------------------------------------
    # Select CUSUM intervals overlapping the candidate deletion.
    # --------------------------------------------------------

    local_cusum = (
        cusum_df[
            (
                cusum_df[
                    "END"
                ] >= candidate_start
            )
            &
            (
                cusum_df[
                    "START"
                ] <= candidate_end
            )
        ]
        .copy()
        .sort_values(
            [
                "START",
                "END"
            ]
        )
        .reset_index(
            drop=True
        )
    )

    # --------------------------------------------------------
    # Y axis
    # --------------------------------------------------------

    ratio_values = (
        pd.to_numeric(
            df[
                ratio_col
            ],
            errors="coerce"
        )
    )

    q995 = (
        np.nanquantile(
            ratio_values,
            0.995
        )
    )

    ymax = max(
        1.45,
        min(
            1.80,
            q995 * 1.06
        )
    )

    ax.set_ylim(
        0,
        ymax
    )

    ax.set_xlim(
        x_min,
        x_max
    )

    # --------------------------------------------------------
    # CUSUM interval bar
    # --------------------------------------------------------

    bar_y = (
        0.035
    )

    bar_h = (
        0.125
    )

    legend_handles = []

    for i, (_, row) in enumerate(
        local_cusum.iterrows(),
        start=1
    ):

        s = int(
            row[
                "START"
            ]
        )

        e = int(
            row[
                "END"
            ]
        )

        mid = (
            s + e
        ) / 2

        mean_ratio = float(
            row[
                "MEAN_RATIO"
            ]
        )

        interval_color = (
            INTERVAL_COLORS[
                (
                    i - 1
                )
                % len(
                    INTERVAL_COLORS
                )
            ]
        )

        # ----------------------------------------------------
        # Interval blocks
        # ----------------------------------------------------

        rect = Rectangle(
            (
                s,
                bar_y
            ),

            e - s + 1,

            bar_h,

            facecolor=interval_color,

            edgecolor=(
                interval_color
            ),

            linewidth=1.0,

            zorder=5
        )

        ax.add_patch(
            rect
        )

        # ----------------------------------------------------
        # Show only numeric labels inside intervals.
        # Use white labels consistently.
        # ----------------------------------------------------

        ax.text(
            mid,

            bar_y
            + bar_h / 2,

            str(i),

            ha="center",
            va="center",

            fontsize=10,

            fontweight="bold",

            color="white",

            zorder=7
        )

        # ----------------------------------------------------
        # Legend on the right
        # ----------------------------------------------------

        legend_label = (
            f"区间 {i}\n"
            f"{format_coord(s)}–"
            f"{format_coord(e)}\n"
            "平均深度比 = "
            f"{mean_ratio:.3f}"
        )

        legend_handles.append(
            Patch(
                facecolor=(
                    interval_color
                ),

                edgecolor=(
                    interval_color
                ),

                label=(
                    legend_label
                )
            )
        )

    # --------------------------------------------------------
    # Reference value
    # --------------------------------------------------------

    label_x = (
        x_min
        + (
            x_max
            - x_min
        )
        * 0.002
    )

    ax.text(
        label_x,

        1.0
        + ymax * 0.010,

        "1.00",

        fontsize=9,

        va="bottom",
        ha="left",

        color=COLORS[
            "ref_line"
        ]
    )

    ax.text(
        label_x,

        0.65
        + ymax * 0.010,

        "0.65",

        fontsize=9,

        va="bottom",
        ha="left",

        color=COLORS[
            "ref_line"
        ]
    )

    # --------------------------------------------------------
    # CUSUM legend
    # --------------------------------------------------------

    if legend_handles:

        legend = ax.legend(
            handles=(
                legend_handles
            ),

            title=(
                "CUSUM 区间"
            ),

            loc=(
                "upper left"
            ),

            bbox_to_anchor=(
                1.015,
                1.00
            ),

            borderaxespad=0.0,

            frameon=True,

            facecolor="white",

            edgecolor="#666666",

            framealpha=1.0,

            fontsize=9,

            title_fontsize=10,

            labelspacing=1.05,

            handlelength=2.0,

            handleheight=1.4,

            borderpad=0.8
        )

        legend.get_title().set_fontweight(
            "bold"
        )

    ax.set_ylabel(
        "深度比值"
    )

    ax.set_title(
        (
            "C  碱基层面深度比值及 "
            "CUSUM 支持区间"
        ),

        loc="left",

        fontsize=12,

        fontweight="bold",

        color=COLORS[
            "text"
        ]
    )

    ax.spines[
        "right"
    ].set_visible(
        False
    )

    ax.spines[
        "top"
    ].set_visible(
        False
    )


# ============================================================
# Panel D
# ============================================================

def make_arc(
    x1,
    x2,
    y_base,
    height
):

    mid = (
        x1
        + x2
    ) / 2

    vertices = [

        (
            x1,
            y_base
        ),

        (
            mid,
            height
        ),

        (
            x2,
            y_base
        )
    ]

    codes = [

        MplPath.MOVETO,

        MplPath.CURVE3,

        MplPath.CURVE3
    ]

    return MplPath(
        vertices,
        codes
    )


def draw_discordant(
    ax,
    dr_cluster_row,
    dr_pairs_df,
    x_min,
    x_max
):

    if dr_cluster_row is None:
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(0, 1.16)
        ax.set_ylabel("异常配对")
        ax.set_yticks([])
        ax.set_title(
            "D  异常配对读段证据：无重叠的 PASS 聚类",
            loc="left",
            fontsize=12,
            fontweight="bold",
            color=COLORS["text"]
        )
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        return

    cluster_id = str(
        dr_cluster_row[
            "Cluster_ID"
        ]
    )

    df = dr_pairs_df[
        dr_pairs_df[
            "Cluster_ID"
        ]
        .astype(str)
        == cluster_id
    ].copy()

    df = (
        df.reset_index(
            drop=True
        )
    )

    n = len(
        df
    )

    if n == 0:

        die(
            "未找到异常配对"
            "支持读段"
        )

    # --------------------------------------------------------
    # Retain all supporting read pairs.
    # --------------------------------------------------------

    for i, (_, row) in enumerate(
        df.iterrows()
    ):

        x1 = int(
            row[
                "Left_Breakpoint"
            ]
        )

        x2 = int(
            row[
                "Right_Breakpoint"
            ]
        )

        # Use different arc heights to make overlapping links visible.
        height = (
            0.34
            + (
                i % 12
            ) * 0.035
        )

        path = make_arc(
            x1,
            x2,

            y_base=0.08,

            height=height
        )

        patch = PathPatch(
            path,

            fill=False,

            linewidth=0.85,

            edgecolor=COLORS[
                "discordant"
            ],

            alpha=0.90
        )

        ax.add_patch(
            patch
        )

        # breakpoint marks
        ax.plot(
            [
                x1,
                x2
            ],

            [
                0.08,
                0.08
            ],

            linestyle="None",

            marker="|",

            markersize=6,

            color=COLORS[
                "discordant_dark"
            ]
        )

    # --------------------------------------------------------
    # DR Cluster
    # --------------------------------------------------------

    cstart = int(
        dr_cluster_row[
            "Start"
        ]
    )

    cend = int(
        dr_cluster_row[
            "End"
        ]
    )

    ax.plot(
        [
            cstart,
            cend
        ],

        [
            0.91,
            0.91
        ],

        color=COLORS[
            "discordant_dark"
        ],

        linewidth=2.5,

        solid_capstyle="butt"
    )

    ax.text(
        (
            cstart
            + cend
        ) / 2,

        0.97,

        cluster_id,

        ha="center",

        va="bottom",

        fontsize=9,

        color=COLORS[
            "discordant_dark"
        ]
    )

    ax.set_xlim(
        x_min,
        x_max
    )

    ax.set_ylim(
        0,
        1.16
    )

    ax.set_ylabel(
        "异常配对"
    )

    ax.set_yticks([])

    median_insert = (
        dr_cluster_row.get(
            "Median_Insert_Size",
            np.nan
        )
    )

    title = (
        "D  异常配对读段证据："
        f"{cluster_id}，"
        "支持数 = "
        f"{int(dr_cluster_row['Support_Pairs'])}"
    )

    if not pd.isna(
        median_insert
    ):

        title += (
            "，中位插入片段长度 = "
            f"{float(median_insert):.0f} bp"
        )

    ax.set_title(
        title,

        loc="left",

        fontsize=12,

        fontweight="bold",

        color=COLORS[
            "text"
        ]
    )

    ax.spines[
        "right"
    ].set_visible(
        False
    )

    ax.spines[
        "top"
    ].set_visible(
        False
    )


# ============================================================
# Panel E
# ============================================================

def draw_summary(
    ax,
    candidate,
    dr_cluster
):

    ax.axis(
        "off"
    )

    chrom = str(
        candidate[
            "Chrom"
        ]
    )

    start = int(
        candidate[
            "Start"
        ]
    )

    end = int(
        candidate[
            "End"
        ]
    )

    size = int(
        candidate[
            "Size_bp"
        ]
    )

    affected = sorted(
        exon_number_from_annotation(
            candidate.get(
                "Affected_Exons",
                ""
            )
        ),

        key=lambda x:
        int(x)
    )

    affected_text = (
        ",".join(
            f"EX{x}"
            for x in affected
        )
        if affected
        else "NA"
    )

    split_count = int(
        candidate.get(
            "Split_Reads",
            0
        )
    )

    depth_fraction = float(
        candidate.get(
            "Depth_Coverage_Fraction",
            np.nan
        )
    )

    dr_fraction = float(
        candidate.get(
            "Discordant_Overlap_Fraction",
            np.nan
        )
    )

    dr_count_value = candidate.get("Discordant_Pairs", 0)
    dr_count = int(dr_count_value) if str(dr_count_value).isdigit() else 0

    level = str(
        candidate.get(
            "Evidence_Level",
            "NA"
        )
    )

    level_cn = (
        evidence_level_cn(
            level
        )
    )

    # --------------------------------------------------------
    # Panel E title
    # Do not repeat the transcript identifier.
    # --------------------------------------------------------

    ax.set_title(
        "E  HCMExonDel 综合检测结果",

        loc="left",

        fontsize=12,

        fontweight="bold",

        color=COLORS[
            "text"
        ]
    )

    # --------------------------------------------------------
    # Information-box layout
    # --------------------------------------------------------

    boxes = [

        # Candidate
        (
            0.01,
            0.13,
            0.33,
            0.72
        ),

        # Split
        (
            0.36,
            0.13,
            0.14,
            0.72
        ),

        # Depth
        (
            0.52,
            0.13,
            0.16,
            0.72
        ),

        # Discordant
        (
            0.70,
            0.13,
            0.16,
            0.72
        ),

        # Evidence
        (
            0.88,
            0.13,
            0.10,
            0.72
        ),
    ]

    # ========================================================
    # Candidate box
    # ========================================================

    x, y, w, h = (
        boxes[0]
    )

    ax.add_patch(
        Rectangle(
            (
                x,
                y
            ),

            w,
            h,

            transform=ax.transAxes,

            facecolor="white",

            edgecolor=COLORS[
                "axis"
            ],

            linewidth=1.2
        )
    )

    ax.text(
        x + 0.015,

        y + 0.55,

        "候选缺失",

        transform=ax.transAxes,

        ha="left",

        va="center",

        fontsize=10,

        fontweight="bold",

        color=COLORS[
            "text"
        ]
    )

    ax.text(
        x + 0.015,

        y + 0.32,

        (
            f"{chrom}:"
            f"{format_coord(start)}-"
            f"{format_coord(end)}"
        ),

        transform=ax.transAxes,

        ha="left",

        va="center",

        fontsize=9,

        color=COLORS[
            "text"
        ]
    )

    ax.text(
        x + 0.015,

        y + 0.12,

        (
            "缺失长度: "
            f"{format_coord(size)} bp   "
            "受影响外显子: "
            f"{affected_text}"
        ),

        transform=ax.transAxes,

        ha="left",

        va="center",

        fontsize=9,

        color=COLORS[
            "text"
        ]
    )

    # ========================================================
    # Split box
    # ========================================================

    x, y, w, h = (
        boxes[1]
    )

    ax.add_patch(
        Rectangle(
            (
                x,
                y
            ),

            w,
            h,

            transform=ax.transAxes,

            facecolor="#F1FAF8",

            edgecolor=COLORS[
                "split"
            ],

            linewidth=1.2
        )
    )

    ax.text(
        x + w / 2,

        y + 0.53,

        "分裂读段",

        transform=ax.transAxes,

        ha="center",

        va="center",

        fontsize=10,

        fontweight="bold",

        color=COLORS[
            "split_dark"
        ]
    )

    ax.text(
        x + w / 2,

        y + 0.20,

        f"{split_count}",

        transform=ax.transAxes,

        ha="center",

        va="center",

        fontsize=14,

        fontweight="bold",

        color=COLORS[
            "split_dark"
        ]
    )

    # ========================================================
    # Depth box
    # ========================================================

    x, y, w, h = (
        boxes[2]
    )

    ax.add_patch(
        Rectangle(
            (
                x,
                y
            ),

            w,
            h,

            transform=ax.transAxes,

            facecolor="#F2F6F7",

            edgecolor=COLORS[
                "depth"
            ],

            linewidth=1.2
        )
    )

    ax.text(
        x + w / 2,

        y + 0.53,

        "深度支持",

        transform=ax.transAxes,

        ha="center",

        va="center",

        fontsize=10,

        fontweight="bold",

        color=COLORS[
            "depth"
        ]
    )

    ax.text(
        x + w / 2,

        y + 0.20,

        f"{depth_fraction:.2%}",

        transform=ax.transAxes,

        ha="center",

        va="center",

        fontsize=13,

        fontweight="bold",

        color=COLORS[
            "depth"
        ]
    )

    # ========================================================
    # Discordant box
    # ========================================================

    x, y, w, h = (
        boxes[3]
    )

    ax.add_patch(
        Rectangle(
            (
                x,
                y
            ),

            w,
            h,

            transform=ax.transAxes,

            facecolor="#F7F3FA",

            edgecolor=COLORS[
                "discordant"
            ],

            linewidth=1.2
        )
    )

    ax.text(
        x + w / 2,

        y + 0.57,

        "异常配对",

        transform=ax.transAxes,

        ha="center",

        va="center",

        fontsize=10,

        fontweight="bold",

        color=COLORS[
            "discordant_dark"
        ]
    )

    ax.text(
        x + w / 2,

        y + 0.31,

        f"{dr_count} 对",

        transform=ax.transAxes,

        ha="center",

        va="center",

        fontsize=11,

        fontweight="bold",

        color=COLORS[
            "discordant_dark"
        ]
    )

    ax.text(
        x + w / 2,

        y + 0.11,

        f"{dr_fraction:.2%}",

        transform=ax.transAxes,

        ha="center",

        va="center",

        fontsize=10,

        color=COLORS[
            "discordant_dark"
        ]
    )

    # ========================================================
    # Evidence box
    # ========================================================

    x, y, w, h = (
        boxes[4]
    )

    level_style = {

        "High": (
            COLORS[
                "high"
            ],
            "white"
        ),

        "Moderate": (
            COLORS[
                "moderate"
            ],
            "white"
        ),

        "Low": (
            COLORS[
                "low"
            ],
            "white"
        ),
    }

    face, text_color = (
        level_style.get(
            level,
            (
                "#777777",
                "white"
            )
        )
    )

    ax.add_patch(
        Rectangle(
            (
                x,
                y
            ),

            w,
            h,

            transform=ax.transAxes,

            facecolor=face,

            edgecolor=face,

            linewidth=1.2
        )
    )

    ax.text(
        x + w / 2,

        y + 0.54,

        "证据等级",

        transform=ax.transAxes,

        ha="center",

        va="center",

        fontsize=10,

        fontweight="bold",

        color=text_color
    )

    ax.text(
        x + w / 2,

        y + 0.20,

        level_cn,

        transform=ax.transAxes,

        ha="center",

        va="center",

        fontsize=14,

        fontweight="bold",

        color=text_color
    )


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Render the Chinese color version of "
            "HCMExonDel multi-track deletion evidence plot"
        )
    )

    parser.add_argument(
        "input_dir",

        help=(
            "HCMExonDel sample result directory"
        )
    )

    parser.add_argument(
        "--gene",

        default="FHOD3",

        help=(
            "Target gene "
            "[默认: FHOD3]"
        )
    )

    parser.add_argument(
        "--gtf",

        required=True,

        help=(
            "GTF annotation file. Provide an absolute path."
        )
    )

    parser.add_argument(
        "--candidate-id",
        default=None,
        help="Optional integrated Candidate_ID from S6_Report"
    )

    parser.add_argument(
        "--plot-flank",

        type=int,

        default=5000,

        help=(
            "Flanking region length around the candidate deletion (bp) "
            "[默认: 5000]"
        )
    )

    parser.add_argument(
        "--font",

        default=None,

        help=(
            "Chinese font file path, for example "
            "/path/to/NotoSansCJKsc-Regular.otf"
        )
    )

    parser.add_argument(
        "--output-prefix",

        default=None,

        help=(
            "Output file prefix"
        )
    )

    parser.add_argument(
        "--dpi",

        type=int,

        default=300
    )

    args = (
        parser.parse_args()
    )

    # ========================================================
    # Chinese font setup
    # ========================================================

    setup_chinese_font(
        args.font
    )

    # ========================================================
    # Input directory
    # ========================================================

    input_dir = Path(
        args.input_dir
    ).resolve()

    if not input_dir.exists():

        die(
            "Input directory does not exist: "
            f"{input_dir}"
        )

    gtf_file = Path(
        args.gtf
    )

    if not gtf_file.exists():

        die(
            "GTF file does not exist: "
            f"{gtf_file}"
        )

    sample = (
        input_dir.name
    )

    gene = (
        args.gene
    )

    info(
        f"样本    : {sample}"
    )

    info(
        f"基因    : {gene}"
    )

    info(
        f"输入目录: {input_dir}"
    )

    # ========================================================
    # Input files
    # ========================================================

    files = discover_files(
        input_dir,
        sample,
        gene
    )

    required = [

        "candidate",

        "merged_candidate",

        "split_cluster",

        "split_events",

        "split_assignments",

        "dr_cluster",

        "dr_pairs",

        "dr_assignments",

        "depth_ratio",

        "cusum"
    ]

    for key in required:

        if (
            files[
                key
            ] is None
        ):

            die(
                "Required file is missing: "
                f"{key}"
            )

        info(
            f"{key:14s}: "
            f"{files[key]}"
        )

    # ========================================================
    # Read input data
    # ========================================================

    candidate_df = read_tsv(
        files[
            "candidate"
        ]
    )

    merged_candidate_df = read_tsv(
        files[
            "merged_candidate"
        ]
    )

    sr_cluster_df = read_tsv(
        files[
            "split_cluster"
        ]
    )

    sr_event_df = read_tsv(files["split_events"])
    sr_assignment_df = read_tsv(files["split_assignments"])
    sr_raw_df = sr_event_df.merge(
        sr_assignment_df,
        on="Event_ID",
        how="inner"
    )

    dr_cluster_df = read_tsv(
        files[
            "dr_cluster"
        ]
    )

    dr_pair_df = read_tsv(files["dr_pairs"])
    dr_assignment_df = read_tsv(files["dr_assignments"])
    dr_pairs_df = dr_pair_df.merge(
        dr_assignment_df,
        on="Pair_ID",
        how="inner"
    )

    depth_df = read_tsv(
        files[
            "depth_ratio"
        ]
    )

    cusum_df = read_tsv(
        files[
            "cusum"
        ]
    )

    # ========================================================
    # Candidate
    # ========================================================

    candidate = (
        select_candidate(
            candidate_df,
            gene,
            args.candidate_id
        )
    )

    chrom = str(
        candidate[
            "Chrom"
        ]
    )

    cand_start = int(
        candidate[
            "Start"
        ]
    )

    cand_end = int(
        candidate[
            "End"
        ]
    )

    candidate_id = str(candidate["Candidate_ID"])
    merged_row = merged_candidate_df[
        merged_candidate_df["Candidate_ID"].astype(str).eq(candidate_id)
    ]
    if merged_row.empty:
        die(f"S5 merged candidate not found for Candidate_ID: {candidate_id}")

    cluster_id = str(merged_row.iloc[0]["Split_Cluster_ID"])
    gene_names = [item.strip() for item in str(candidate.get("Gene", "")).split(",")]
    transcript_values = [item.strip() for item in str(candidate.get("Transcript", "")).split(",")]
    transcript = ""
    if gene in gene_names:
        gene_index = gene_names.index(gene)
        if gene_index < len(transcript_values):
            transcript = transcript_values[gene_index].split(";")[0]
    if not transcript and transcript_values:
        transcript = transcript_values[0].split(";")[0]

    if not transcript or transcript == "NA":
        die("Candidate event lacks a valid MANE Transcript")

    info(
        "Selected candidate event: "
        f"{candidate_id} / {cluster_id} "
        f"{chrom}:{cand_start}-{cand_end}"
    )

    # ========================================================
    # Split cluster
    # ========================================================

    sr_tmp = sr_cluster_df[
        sr_cluster_df[
            "Cluster_ID"
        ]
        .astype(str)
        == cluster_id
    ]

    if sr_tmp.empty:

        die(
            "Split-read cluster not found: "
            f"{cluster_id}"
        )

    sr_cluster = (
        sr_tmp.iloc[0]
    )

    # ========================================================
    # Discordant cluster
    # ========================================================

    has_discordant_support = str(candidate.get("Discordant_Support", "No")).lower() in {
        "yes", "1", "true", "pass", "supported"
    }
    dr_candidates = dr_cluster_df.copy()
    if has_discordant_support:
        dr_candidates = dr_candidates[
            dr_candidates["Cluster_Status"].astype(str).str.lower().eq("pass")
            & dr_candidates["Chrom"].astype(str).eq(chrom)
        ].copy()
    else:
        dr_candidates = dr_candidates.iloc[0:0].copy()

    if not dr_candidates.empty:
        overlap_start = dr_candidates["Start"].astype(int).clip(lower=cand_start)
        overlap_end = dr_candidates["End"].astype(int).clip(upper=cand_end)
        overlap_bases = (overlap_end - overlap_start + 1).clip(lower=0)
        dr_candidates["_candidate_overlap_fraction"] = overlap_bases / (cand_end - cand_start + 1)
        dr_candidates["_support_pairs_numeric"] = pd.to_numeric(
            dr_candidates["Support_Pairs"], errors="coerce"
        ).fillna(0)
        dr_candidates = dr_candidates[dr_candidates["_candidate_overlap_fraction"] > 0]

    if dr_candidates.empty:
        dr_cluster = None
        dr_pairs_for_plot = dr_pairs_df.iloc[0:0].copy()
    else:
        dr_cluster = dr_candidates.sort_values(
            ["_candidate_overlap_fraction", "_support_pairs_numeric"],
            ascending=[False, False]
        ).iloc[0]
        dr_id = str(dr_cluster["Cluster_ID"])
        dr_pairs_for_plot = dr_pairs_df[
            dr_pairs_df["Cluster_ID"].astype(str).eq(dr_id)
        ].copy()

    # ========================================================
    # GTF
    # ========================================================

    (
        exon_df,
        transcript_info
    ) = parse_gtf_exons(
        gtf_file,
        gene,
        transcript
    )

    affected_exons = (
        exon_number_from_annotation(
            candidate.get(
                "Affected_Exons",
                ""
            )
        )
    )

    # ========================================================
    # Plot interval
    # ========================================================

    x_min = max(
        1,

        cand_start
        - args.plot_flank
    )

    x_max = (
        cand_end
        + args.plot_flank
    )

    # ========================================================
    # Figure
    # ========================================================

    fig = plt.figure(
        figsize=(
            15.5,
            10.8
        ),

        constrained_layout=False,

        facecolor="white"
    )

    gs = fig.add_gridspec(

        nrows=5,

        ncols=1,

        height_ratios=[
            1.20,
            2.25,
            2.85,
            2.15,
            1.50
        ],

        hspace=0.53
    )

    ax_gene = fig.add_subplot(
        gs[0]
    )

    ax_split = fig.add_subplot(
        gs[1],

        sharex=ax_gene
    )

    ax_depth = fig.add_subplot(
        gs[2],

        sharex=ax_gene
    )

    ax_dr = fig.add_subplot(
        gs[3],

        sharex=ax_gene
    )

    ax_summary = fig.add_subplot(
        gs[4]
    )

    # ========================================================
    # Candidate interval background
    # ========================================================

    for ax in [

        ax_gene,

        ax_split,

        ax_depth,

        ax_dr
    ]:

        add_candidate_span(
            ax,
            cand_start,
            cand_end
        )

        add_candidate_boundaries(
            ax,
            cand_start,
            cand_end
        )

    # ========================================================
    # A
    # ========================================================

    draw_gene_track(
        ax_gene,

        exon_df,

        transcript_info,

        affected_exons,

        x_min,

        x_max,

        gene
    )

    # ========================================================
    # B
    # ========================================================

    draw_split_reads(
        ax_split,

        sr_cluster,

        sr_raw_df,

        x_min,

        x_max
    )

    # ========================================================
    # C
    # ========================================================

    draw_depth_track(
        ax_depth,

        depth_df,

        cusum_df,

        cand_start,

        cand_end,

        x_min,

        x_max
    )

    # ========================================================
    # D
    # ========================================================

    draw_discordant(
        ax_dr,

        dr_cluster,

        dr_pairs_for_plot,

        x_min,

        x_max
    )

    # ========================================================
    # E
    # ========================================================

    draw_summary(
        ax_summary,

        candidate,

        dr_cluster
    )

    # ========================================================
    # X axis
    # ========================================================

    for ax in [

        ax_gene,

        ax_split,

        ax_depth
    ]:

        ax.tick_params(
            axis="x",

            labelbottom=False
        )

    format_genomic_axis(
        ax_dr
    )

    ax_dr.set_xlabel(
        (
            f"{chrom} 基因组坐标"
            "（GRCh37/hg19）"
        ),

        fontsize=11,

        color=COLORS[
            "text"
        ]
    )

    # ========================================================
    # Overall title
    #
    # Do not use a colon.
    # ========================================================

    fig.suptitle(
        (
            f"{sample}样本 "
            f"{gene} 多信号缺失证据图"
        ),

        fontsize=16,

        fontweight="bold",

        y=0.985,

        color=COLORS[
            "text"
        ]
    )

    # Reserve space on the right for the CUSUM legend.
    fig.subplots_adjust(
        left=0.075,
        right=0.82,
        top=0.945,
        bottom=0.06
    )

    # ========================================================
    # Output
    # ========================================================

    if args.output_prefix:

        out_prefix = Path(
            args.output_prefix
        )

    else:

        outdir = (
            input_dir
            / "S6_Report"
        )

        outdir.mkdir(
            parents=True,
            exist_ok=True
        )

        out_prefix = (
            outdir
            /
            (
                f"{sample}."
                f"{gene}."
                "multisignal.zh_color"
            )
        )

    png = (
        str(
            out_prefix
        )
        + ".png"
    )

    pdf = (
        str(
            out_prefix
        )
        + ".pdf"
    )

    svg = (
        str(
            out_prefix
        )
        + ".svg"
    )

    fig.savefig(
        png,

        dpi=args.dpi,

        bbox_inches="tight",

        facecolor="white"
    )

    fig.savefig(
        pdf,

        bbox_inches="tight",

        facecolor="white"
    )

    fig.savefig(
        svg,

        bbox_inches="tight",

        facecolor="white"
    )

    plt.close(
        fig
    )

    info(
        "绘图完成"
    )

    info(
        f"PNG: {png}"
    )

    info(
        f"PDF: {pdf}"
    )

    info(
        f"SVG: {svg}"
    )


if __name__ == "__main__":

    main()

