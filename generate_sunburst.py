import hail as hl
import plotly.graph_objects as go

# 1. Initialize Hail Session

hl.init(default_reference="GRCh38", log="/tmp/hail_sunburst.log")

# 2. Define Metric Aggregation Helper

def execute_filter_step(mt, step_name, short_code, keep_expr):
    """
    Filters MatrixTable rows based on keep_expr and calculates per-step 
    SNV vs Indel metrics for both Rows (variants) and Entries (genotypes).
    """
    removed_mt = mt.filter_rows(keep_expr, keep=False)
    mt_filtered = mt.filter_rows(keep_expr, keep=True)

    # hl.is_snp takes ref allele (alleles[0]) and alt allele (alleles[1])
    is_snv = hl.is_snp(removed_mt.alleles[0], removed_mt.alleles[1])

    removed_rows = removed_mt.count_rows()
    removed_entries = removed_mt.aggregate_entries(hl.agg.count()) if removed_rows > 0 else 0

    snv_rows = removed_mt.aggregate_rows(hl.agg.count_where(is_snv)) if removed_rows > 0 else 0
    snv_entries = removed_mt.aggregate_entries(hl.agg.count_where(is_snv)) if removed_entries > 0 else 0

    indel_rows = removed_rows - snv_rows
    indel_entries = removed_entries - snv_entries

    total_rows = mt_filtered.count_rows()
    total_entries = mt_filtered.aggregate_entries(hl.agg.count())

    metrics = {
        "step_name": step_name,
        "short_code": short_code,
        "removed_rows": removed_rows,
        "total_rows": total_rows,
        "removed_entries": removed_entries,
        "total_entries": total_entries,
        "snv_rows": snv_rows,
        "snv_entries": snv_entries,
        "indel_rows": indel_rows,
        "indel_entries": indel_entries,
    }

    return mt_filtered, metrics

# 3. Plotly Sunburst Generation Function

def create_progressive_sunburst(filter_history, final_mt, current_step_index):
    initial_total_rows = filter_history[0]["total_rows"] + filter_history[0]["removed_rows"]
    initial_total_entries = filter_history[0]["total_entries"] + filter_history[0]["removed_entries"]

    if current_step_index < 0:
        cumulative_removed_rows = 0
        cumulative_removed_entries = 0
    else:
        cumulative_removed_rows = sum(filter_history[i]["removed_rows"] for i in range(current_step_index + 1))
        cumulative_removed_entries = sum(filter_history[i]["removed_entries"] for i in range(current_step_index + 1))

    retained_rows = initial_total_rows - cumulative_removed_rows
    retained_entries = initial_total_entries - cumulative_removed_entries

    # Retained SNV vs Indel metrics
    if retained_rows > 0:
        is_snv_expr = hl.is_snp(final_mt.alleles[0], final_mt.alleles[1])
        retained_snv_rows = final_mt.aggregate_rows(hl.agg.count_where(is_snv_expr))
        retained_snv_entries = final_mt.aggregate_entries(hl.agg.count_where(is_snv_expr))
        retained_indel_rows = retained_rows - retained_snv_rows
        retained_indel_entries = retained_entries - retained_snv_entries
    else:
        retained_snv_rows = 0
        retained_snv_entries = 0
        retained_indel_rows = 0
        retained_indel_entries = 0

    ids = ["total_initial", "retained", "total_removed"]
    labels = ["Total Initial", "Retained", "Total Removed"]
    parents = ["", "total_initial", "total_initial"]
    values = [initial_total_entries, retained_entries, cumulative_removed_entries]

    customdata = [
        [initial_total_rows, "Total Initial Dataset"],
        [retained_rows, "Retained Data"],
        [cumulative_removed_rows, "Total Removed Data"]
    ]

    colors = ["#ECEFF1", "#2ECC71", "#E53935"]

    # Retained Sub-layer: Blue Tones
    retained_snv_color = "#1565C0"   # Dark Blue
    retained_indel_color = "#64B5F6" # Light Blue

    if retained_entries > 0:
        ids.extend(["retained_snv", "retained_indel"])
        labels.extend(["SNV", "INDEL"])
        parents.extend(["retained", "retained"])
        values.extend([retained_snv_entries, retained_indel_entries])
        customdata.extend([
            [retained_snv_rows, "Retained Data - SNVs"],
            [retained_indel_rows, "Retained Data - Indels"]
        ])
        colors.extend([retained_snv_color, retained_indel_color])

    # Removed Steps Palette
    snv_orange = "#FB8C00"   # Medium-Dark Orange
    indel_orange = "#FFB74D" # Soft Light Orange
    past_reds = ["#FFCDD2", "#EF9A9A", "#E57373", "#EF5350", "#E53935", "#D32F2F", "#C62828"]
    current_red = "#FF0000"

    # Step 2.1 Calculated Consequences Parent Node
    step2_1_1 = filter_history[0] if current_step_index >= 0 else None
    step2_1_2 = filter_history[1] if current_step_index >= 1 else None

    if step2_1_1:
        step2_1_entries = step2_1_1["removed_entries"] + (step2_1_2["removed_entries"] if step2_1_2 else 0)
        step2_1_rows = step2_1_1["removed_rows"] + (step2_1_2["removed_rows"] if step2_1_2 else 0)

        ids.append("step_2_1")
        labels.append("2.1<br>IMPACT")
        parents.append("total_removed")
        values.append(step2_1_entries)
        customdata.append([step2_1_rows, "2.1 Calculated Consequences"])
        colors.append(current_red if current_step_index in [0, 1] else past_reds[0])

        # Substep 2.1.1: LOW
        ids.append("step_2_1_1")
        labels.append("2.1.1<br>LOW")
        parents.append("step_2_1")
        values.append(step2_1_1["removed_entries"])
        customdata.append([step2_1_1["removed_rows"], step2_1_1["step_name"]])
        colors.append("#FFCDD2")

        ids.extend(["step_2_1_1_snv", "step_2_1_1_indel"])
        labels.extend(["SNV", "INDEL"])
        parents.extend(["step_2_1_1", "step_2_1_1"])
        values.extend([step2_1_1["snv_entries"], step2_1_1["indel_entries"]])
        customdata.extend([
            [step2_1_1["snv_rows"], f"{step2_1_1['step_name']} - SNVs"],
            [step2_1_1["indel_rows"], f"{step2_1_1['step_name']} - Indels"]
        ])
        colors.extend([snv_orange, indel_orange])

        # Substep 2.1.2: MODIFIER
        if step2_1_2:
            ids.append("step_2_1_2")
            labels.append("2.1.2<br>MODIFIER")
            parents.append("step_2_1")
            values.append(step2_1_2["removed_entries"])
            customdata.append([step2_1_2["removed_rows"], step2_1_2["step_name"]])
            colors.append("#EF9A9A")

            ids.extend(["step_2_1_2_snv", "step_2_1_2_indel"])
            labels.extend(["SNV", "INDEL"])
            parents.extend(["step_2_1_2", "step_2_1_2"])
            values.extend([step2_1_2["snv_entries"], step2_1_2["indel_entries"]])
            customdata.extend([
                [step2_1_2["snv_rows"], f"{step2_1_2['step_name']} - SNVs"],
                [step2_1_2["indel_rows"], f"{step2_1_2['step_name']} - Indels"]
            ])
            colors.extend([snv_orange, indel_orange])

    # Main Steps: 2.2 through 2.8
    for idx in range(2, current_step_index + 1):
        step = filter_history[idx]
        step_name = step["step_name"]
        rem_rows = step["removed_rows"]
        rem_entries = step["removed_entries"]

        step_id = f"step_{idx}"
        slice_label = step.get("short_code", "<br>".join(step_name.split()))

        ids.append(step_id)
        labels.append(slice_label)
        parents.append("total_removed")
        values.append(rem_entries)
        customdata.append([rem_rows, step_name])

        step_color = current_red if idx == current_step_index else past_reds[min(idx - 1, len(past_reds) - 1)]
        colors.append(step_color)

        snv_rows = step["snv_rows"]
        snv_entries = step["snv_entries"]
        indel_rows = step["indel_rows"]
        indel_entries = step["indel_entries"]

        ids.extend([f"{step_id}_snv", f"{step_id}_indel"])
        labels.extend(["SNV", "INDEL"])
        parents.extend([step_id, step_id])
        values.extend([snv_entries, indel_entries])
        customdata.extend([
            [snv_rows, f"{step_name} - SNVs"],
            [indel_rows, f"{step_name} - Indels"]
        ])
        colors.extend([snv_orange, indel_orange])

    customdata_formatted = []
    for idx, cd in enumerate(customdata):
        row_cnt = cd[0]
        step_name = cd[1]
        row_pct = (row_cnt / initial_total_rows) if initial_total_rows > 0 else 0
        entry_cnt = values[idx]
        entry_pct = (entry_cnt / initial_total_entries) if initial_total_entries > 0 else 0
        customdata_formatted.append([row_cnt, row_pct, step_name, entry_pct])

    hovertemplate = (
        "<span style='color: #00A8E8; font-weight: bold;'>Step:</span> "
        "<span style='color: #333333;'>%{customdata[2]}</span><br>"
        "<span style='color: #00A8E8; font-weight: bold;'>Rows:</span> "
        "<span style='color: #333333;'>%{customdata[0]:,} (%{customdata[1]:.2%})</span><br>"
        "<span style='color: #00A8E8; font-weight: bold;'>Entries:</span> "
        "<span style='color: #333333;'>%{value:,} (%{customdata[3]:.2%})</span>"
        "<extra></extra>"
    )

    fig = go.Figure(go.Sunburst(
        ids=ids,
        labels=labels,
        parents=parents,
        values=values,
        customdata=customdata_formatted,
        branchvalues="total",
        sort=False,
        insidetextorientation="horizontal",
        marker=dict(colors=colors, line=dict(color="#FFFFFF", width=1.5)),
        textfont=dict(size=10),
        hovertemplate=hovertemplate,
        hoverlabel=dict(bgcolor="#FFFFFF", bordercolor="#D3D3D3", font_size=12, font_family="Arial, sans-serif")
    ))

    title_text = "<b>Filtering Summary: Baseline</b>" if current_step_index < 0 else f"<b>Filtering Summary: Step {filter_history[current_step_index]['step_name'].split()[0]}</b><br><sup>{filter_history[current_step_index]['step_name']}</sup>"
    fig.update_layout(title=dict(text=title_text, x=0.02, y=0.96), margin=dict(t=80, l=10, r=10, b=10), width=650, height=650)

    return fig



# 4. Pipeline Execution

if __name__ == "__main__":
    mt = hl.read_matrix_table("synthetic_exome.mt")

    filter_history = []

    # Substep 2.1.1: LOW Impact Consequences Filter
    mt, step2_1_1 = execute_filter_step(
        mt,
        step_name="2.1.1 LOW Impact Filter",
        short_code="2.1.1<br>LOW",
        keep_expr=(mt.vep.impact != "LOW")
    )
    filter_history.append(step2_1_1)

    # Substep 2.1.2: MODIFIER Impact Consequences Filter
    mt, step2_1_2 = execute_filter_step(
        mt,
        step_name="2.1.2 MODIFIER Impact Filter",
        short_code="2.1.2<br>MODIFIER",
        keep_expr=(mt.vep.impact != "MODIFIER")
    )
    filter_history.append(step2_1_2)

    # Step 2.2: African Allele Frequency Filter
    mt, step2_2 = execute_filter_step(
        mt,
        step_name="2.2 African Allele Frequency Filter",
        short_code="2.2<br>AFR AF",
        keep_expr=(mt.info.AF_afr[0] < 0.05) | hl.is_missing(mt.info.AF_afr)
    )
    filter_history.append(step2_2)

    # Step 2.3: Cohort Allele Frequency Filter
    mt, step2_3 = execute_filter_step(
        mt,
        step_name="2.3 Allele Frequency Filter",
        short_code="2.3<br>Cohort AF",
        keep_expr=(mt.info.AF[0] < 0.05)
    )
    filter_history.append(step2_3)

    # Step 2.4: Quality & Allele Balance
    ad_sum = hl.sum(mt.AD)
    ab = hl.or_missing(ad_sum > 0, mt.AD[1] / ad_sum)
    is_bad_hom_ref = mt.GT.is_hom_ref() & (ab > 0.1)
    is_low_qual_dp = (mt.GQ < 10) | (mt.DP < 20)
    mt = mt.filter_entries(hl.is_defined(mt.GT) & ~(is_bad_hom_ref | is_low_qual_dp))

    mt, step2_4 = execute_filter_step(
        mt,
        step_name="2.4 Quality & Allele Balance",
        short_code="2.4<br>Qual/AB",
        keep_expr=(hl.agg.count_where(hl.is_defined(mt.GT)) > 0)
    )
    filter_history.append(step2_4)

    # Step 2.5: Homozygous Reference Filter
    mt, step2_5 = execute_filter_step(
        mt,
        step_name="2.5 Homozygous Reference",
        short_code="2.5<br>Homo Ref",
        keep_expr=(hl.agg.count_where(mt.GT.is_non_ref()) > 0)
    )
    filter_history.append(step2_5)

    # Step 2.6: Prediction Scores
    if 'CADD_phred' in mt.info.dtype.fields and 'SIFT_score' in mt.info.dtype.fields:
        prediction_keep_expr = (mt.info.CADD_phred >= 15.0) | (mt.info.SIFT_score < 0.05)
    else:
        print("Notice: CADD_phred / SIFT_score not present in MatrixTable schema. Bypassing Step 2.6 for test run.")
        prediction_keep_expr = hl.literal(True)

    mt, step2_6 = execute_filter_step(
        mt,
        step_name="2.6 Prediction Scores",
        short_code="2.6<br>Scores",
        keep_expr=prediction_keep_expr
    )
    filter_history.append(step2_6)

    # Step 2.7: Empty Variants/Sample Filter
    mt = mt.filter_cols(hl.agg.count_where(mt.GT.is_non_ref()) > 0)
    mt, step2_7 = execute_filter_step(
        mt,
        step_name="2.7 Empty",
        short_code="2.7<br>Empty",
        keep_expr=(hl.agg.count_where(mt.GT.is_non_ref()) > 0)
    )
    filter_history.append(step2_7)

    # Step 2.8: Gene Panel Subset Filter
    mt, step2_8 = execute_filter_step(
        mt,
        step_name="2.8 Candidate Gene Panel Filter",
        short_code="2.8<br>Genes",
        keep_expr=hl.literal(["GENE_1", "GENE_2", "GENE_3", "GENE_4", "GENE_5"]).contains(mt.vep.gene) == False
    )
    filter_history.append(step2_8)

    
    # 5. Export Visualization
   
    fig = create_progressive_sunburst(filter_history, final_mt=mt, current_step_index= 8)
    # Index: -1 = baseline, 0 = step 2.1.1 and 8 = step 2.8
    output_html = "filtering_sunburst.html"
    fig.write_html(output_html)
    print(f"Success! Chart saved to {output_html}")