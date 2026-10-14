import hail as hl
import numpy as np
import pandas as pd

from bokeh.plotting import figure, show, output_file
from bokeh.models import (
    ColumnDataSource, 
    Range1d, 
    CustomJSTickFormatter, 
    HoverTool, 
    Span, 
    Label,
    BoxAnnotation
)
from bokeh.layouts import column


# Define global substitution categories and harmonized color palette
CATEGORIES = [
    'A->C / T->G (Tv)',
    'A->G / T->C (Ti)',
    'A->T / T->A (Tv)',
    'G->A / C->T (Ti)',
    'G->C / C->G (Tv)',
    'G->T / C->A (Tv)'
]

# Color gradient (Bottom to Top): Red -> Orange -> Orange-Yellow -> Yellow -> Green -> Blue
COLORS = [
    "#d73027",  # 1. Red (A->C / T->G)
    "#f46d43",  # 2. Orange (A->G / T->C)
    "#fdae61",  # 3. Orange-Yellow (A->T / T->A)
    "#fee08b",  # 4. Yellow (G->A / C->T)
    "#a6d96a",  # 5. Soft Green (G->C / C->G)
    "#4575b4"   # 6. Muted Blue (G->T / C->A)
]


def extract_titv_summary_from_hail(mt: hl.MatrixTable) -> pd.DataFrame:
    """
    Extracts Ti/Tv ratios from mt.sample_qc and aggregates 6-class mutation spectrum proportions per sample.
    """
    # 1. Automatically run sample_qc if it hasn't been executed on this MatrixTable
    if 'sample_qc' not in mt.col:
        mt = hl.sample_qc(mt)

    # 2. Filter to biallelic SNPs
    mt_snps = mt.filter_rows(
        (hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1])
    )

    ref = mt_snps.alleles[0]
    alt = mt_snps.alleles[1]

    # 3. Annotate variants with 6 substitution category strings
    mt_snps = mt_snps.annotate_rows(
        sub_cat=hl.case()
        .when(((ref == 'A') & (alt == 'C')) | ((ref == 'T') & (alt == 'G')), 'A->C / T->G (Tv)')
        .when(((ref == 'A') & (alt == 'G')) | ((ref == 'T') & (alt == 'C')), 'A->G / T->C (Ti)')
        .when(((ref == 'A') & (alt == 'T')) | ((ref == 'T') & (alt == 'A')), 'A->T / T->A (Tv)')
        .when(((ref == 'G') & (alt == 'A')) | ((ref == 'C') & (alt == 'T')), 'G->A / C->T (Ti)')
        .when(((ref == 'G') & (alt == 'C')) | ((ref == 'C') & (alt == 'G')), 'G->C / C->G (Tv)')
        .when(((ref == 'G') & (alt == 'T')) | ((ref == 'C') & (alt == 'A')), 'G->T / C->A (Tv)')
        .or_missing()
    )

    # 4. Pull r_ti_tv from mt.sample_qc and count substitution categories
    sample_qc_mt = mt_snps.annotate_cols(
        ti_tv_ratio=mt_snps.sample_qc.r_ti_tv,
        cat_counts=hl.agg.filter(
            mt_snps.GT.is_non_ref(),
            hl.agg.counter(mt_snps.sub_cat)
        )
    )

    # 5. Collect column annotations into Pandas
    df_raw = sample_qc_mt.cols().to_pandas()

    df_formatted = pd.DataFrame()
    df_formatted['sample'] = df_raw['s']
    df_formatted['ti_tv_ratio'] = df_raw['ti_tv_ratio']

    # Extract substitution category counts
    for cat in CATEGORIES:
        df_formatted[cat] = df_raw['cat_counts'].apply(
            lambda d: d.get(cat, 0) if isinstance(d, dict) else 0
        )

    # Convert counts to 100% proportions
    total_snps = df_formatted[CATEGORIES].sum(axis=1)
    for cat in CATEGORIES:
        df_formatted[cat] = df_formatted[cat] / total_snps.replace(0, 1.0)

    return df_formatted


def create_titv_qc_plot(
    df: pd.DataFrame, 
    initial_visible: int = 30, 
    exome_min: float = 2.8, 
    exome_max: float = 3.3,
    exome_ideal: float = 3.0
):
    """
    Creates a 2-panel Bokeh QC plot for Exome Sequencing data:
      1. Top: Histogram of cohort Ti/Tv distribution with shaded target zone and 3.0 ideal value.
      2. Bottom: Sample-level 100% stacked bar chart with matching visual legend order.
    """
    num_samples = len(df)

    df = df.copy()
    df['x'] = np.arange(1, num_samples + 1)
    sample_names = df['sample'].tolist()

    source = ColumnDataSource(df)

    # JS Formatter for sample labels on tick marks
    js_formatter = CustomJSTickFormatter(args=dict(names=sample_names), code="""
        var index = Math.round(tick) - 1;
        if (index >= 0 && index < names.length) {
            return names[index];
        }
        return '';
    """)

    
    # PANEL 1: Cohort-Wide Ti/Tv Distribution (Histogram)
    
    counts, edges = np.histogram(df['ti_tv_ratio'].dropna(), bins=35)
    hist_source = ColumnDataSource(data=dict(
        top=counts,
        bottom=np.zeros_like(counts),
        left=edges[:-1],
        right=edges[1:]
    ))

    p_hist = figure(
        height=260,
        width=1050,
        title="Sample Ti/Tv Ratio Distribution",
        tools="pan,wheel_zoom,reset,save",
        toolbar_location="right",
        x_axis_label="Ti/Tv Ratio",
        y_axis_label="Frequency"
    )

    p_hist.quad(
        bottom='bottom', top='top', left='left', right='right',
        source=hist_source, 
        fill_color="#6b6eb4", 
        line_color="#3d3f73", 
        line_width=0.8,
        alpha=0.85
    )

    # Shaded target range (2.8 - 3.3)
    acceptable_zone = BoxAnnotation(
        left=exome_min, right=exome_max, 
        fill_color="#c2e6c2", fill_alpha=0.45, 
        line_color="#2ca02c", line_dash="dashed", line_width=1.5
    )
    p_hist.add_layout(acceptable_zone)

    # Range boundaries (2.8 and 3.3)
    line_min = Span(location=exome_min, dimension='height', line_color='#2ca02c', line_dash='dashed', line_width=1.5)
    line_max = Span(location=exome_max, dimension='height', line_color='#2ca02c', line_dash='dashed', line_width=1.5)
    p_hist.add_layout(line_min)
    p_hist.add_layout(line_max)

    # Ideal baseline (3.0)
    line_ideal = Span(location=exome_ideal, dimension='height', line_color='#d95f02', line_dash='dashed', line_width=2.0)
    p_hist.add_layout(line_ideal)

    # Annotation Callouts
    y_max = max(counts) if len(counts) > 0 and max(counts) > 0 else 1
    target_label = Label(
        x=exome_min + 0.02, y=y_max * 0.90, 
        text=f"Expected Exome Range ({exome_min} - {exome_max})", 
        text_color="#1b7837", text_font_size="9.5pt", text_font_style="bold"
    )
    ideal_label = Label(
        x=exome_ideal + 0.01, y=y_max * 0.75, 
        text=f"Ideal Baseline ({exome_ideal})", 
        text_color="#d95f02", text_font_size="9pt", text_font_style="bold"
    )
    p_hist.add_layout(target_label)
    p_hist.add_layout(ideal_label)

    # Hover tool
    hist_hover = HoverTool(tooltips=[
        ("Ti/Tv Range", "@left{0.00} to @right{0.00}"),
        ("Sample Count", "@top")
    ])
    p_hist.add_tools(hist_hover)

    p_hist.yaxis.axis_label_text_font_style = "italic"
    p_hist.title.text_font_style = "bold"


    # PANEL 2: Sample Mutational Spectrum (Pan/Drag enabled)

    p_detail = figure(
        x_range=Range1d(0.5, initial_visible + 0.5, bounds=(0.5, num_samples + 0.5)),
        y_range=(0, 1),
        height=380,
        width=1050,
        tools="xpan,xwheel_zoom,reset,save",
        active_drag="xpan",
        active_scroll="xwheel_zoom",
        toolbar_location="right",
        title="Sample Mutational Spectrum (Click & Drag horizontally to scroll through cohort)",
        y_axis_label="Proportion"
    )

    bars = p_detail.vbar_stack(
        stackers=CATEGORIES,
        x='x',
        width=0.85,
        color=COLORS,
        source=source,
        legend_label=CATEGORIES
    )

    p_detail.xaxis.formatter = js_formatter
    p_detail.xaxis.major_label_orientation = 1.2
    p_detail.xgrid.grid_line_color = None

    p_detail.yaxis.axis_label_text_font_style = "italic"
    p_detail.title.text_font_style = "bold"

    # Match legend order to visual bar stack (Blue on top, Red on bottom)
    legend = p_detail.legend[0]
    legend.items = legend.items[::-1]
    
    p_detail.add_layout(legend, 'right')
    p_detail.legend.background_fill_alpha = 0.0
    p_detail.legend.border_line_alpha = 0.0

    # Hover tool
    detail_hover = HoverTool(
        renderers=bars,
        tooltips=[
            ("Sample ID", "@sample"),
            ("Ti/Tv Ratio", "@ti_tv_ratio{0.00}"),
            ("Category", "$name"),
            ("Proportion", "@$name{0.000}")
        ]
    )
    p_detail.add_tools(detail_hover)

    return column(p_hist, p_detail)


if __name__ == "__main__":
    # Initialize Hail session
    hl.init()

    # Load MatrixTable
    mt = hl.read_matrix_table("synthetic_exome.mt")

    # Extract summary metrics using precomputed sample_qc fields
    df_summary = extract_titv_summary_from_hail(mt)

    # Build Bokeh layout
    qc_plot_layout = create_titv_qc_plot(
        df=df_summary, 
        exome_min=2.8, 
        exome_max=3.3, 
        exome_ideal=3.0
    )

    # Export interactive HTML widget
    output_file("titv_qc_report.html")
    show(qc_plot_layout)