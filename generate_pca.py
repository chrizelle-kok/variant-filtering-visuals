# generate_pca.py
import hail as hl
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

# 1. Initialize Hail

hl.init(default_reference="GRCh38", log="/tmp/hail_pca.log")

# Change how axis is displayed

def get_simplified_ticks(series):
    min_val = series.min()
    max_val = series.max()
    
    vals = []
    texts = []
    
    # Place '-' midpoint in the negative domain if negative values exist
    if min_val < 0:
        vals.append(min_val * 0.5)
        texts.append("-")
        
    # Place '0' at origin
    vals.append(0.0)
    texts.append("0")
    
    # Place '+' midpoint in the positive domain if positive values exist
    if max_val > 0:
        vals.append(max_val * 0.5)
        texts.append("+")
        
    return vals, texts

# 2. PCA Execution & 3D Plot Generation

def run_pca(mt_path="synthetic_exome.mt", k=10):

    mt = hl.read_matrix_table(mt_path)

    eigenvalues, scores_table, _ = hl.hwe_normalized_pca(mt.GT, k=k)

    # Convert Hail Scores Table to a Pandas DataFrame
    pca_df = scores_table.to_pandas()

    # Extract PC1, PC2, PC3 into separate columns
    pca_df["PC1"] = pca_df["scores"].apply(lambda x: x[0])
    pca_df["PC2"] = pca_df["scores"].apply(lambda x: x[1])
    pca_df["PC3"] = pca_df["scores"].apply(lambda x: x[2])

    # Check if sample population/ancestry metadata exists in mt.cols()
    col_fields = mt.col.dtype.fields
    if "pop" in col_fields:
        pop_df = mt.cols().select("pop").to_pandas()
        pca_df = pca_df.merge(pop_df, on="s", how="left")
        color_col = "pop"
    elif "ancestry" in col_fields:
        pop_df = mt.cols().select("ancestry").to_pandas()
        pca_df = pca_df.merge(pop_df, on="s", how="left")
        color_col = "ancestry"
    else:
        pca_df["Group"] = "Cohort"
        color_col = "Group"

    # Calculate percentage of variance explained
    total_top_k_eval = sum(eigenvalues) if sum(eigenvalues) > 0 else 1.0
    var_explained = [(e / total_top_k_eval) * 100 for e in eigenvalues[:3]]
    cum_var_explained = sum(var_explained[:3])

    title_text = (
        f"<b>3D Principal Component Analysis</b><br>"
        f"<sup>Total Variance Explained (PC1–PC3): {cum_var_explained:.2f}%</sup>"
    )

    # Build Plotly Express 3D Scatter Plot
    fig = px.scatter_3d(
        pca_df,
        x="PC1",
        y="PC2",
        z="PC3",
        color=color_col,
        custom_data=["s", color_col] if color_col != "Group" else ["s"],
        title=title_text
    )

    # Hover template
    if color_col != "Group":
        hovertemplate = (
            "<span style='color: #00A8E8; font-weight: bold;'>Sample ID:</span> "
            "<span style='color: #333333;'>%{customdata[0]}</span><br>"
            "<span style='color: #00A8E8; font-weight: bold;'>Group:</span> "
            "<span style='color: #333333;'>%{customdata[1]}</span><br>"
            "<span style='color: #00A8E8; font-weight: bold;'>PC1:</span> "
            "<span style='color: #333333;'>%{x:.4f}</span><br>"
            "<span style='color: #00A8E8; font-weight: bold;'>PC2:</span> "
            "<span style='color: #333333;'>%{y:.4f}</span><br>"
            "<span style='color: #00A8E8; font-weight: bold;'>PC3:</span> "
            "<span style='color: #333333;'>%{z:.4f}</span>"
            "<extra></extra>"
        )
    else:
        hovertemplate = (
            "<span style='color: #00A8E8; font-weight: bold;'>Sample ID:</span> "
            "<span style='color: #333333;'>%{customdata[0]}</span><br>"
            "<span style='color: #00A8E8; font-weight: bold;'>PC1:</span> "
            "<span style='color: #333333;'>%{x:.4f}</span><br>"
            "<span style='color: #00A8E8; font-weight: bold;'>PC2:</span> "
            "<span style='color: #333333;'>%{y:.4f}</span><br>"
            "<span style='color: #00A8E8; font-weight: bold;'>PC3:</span> "
            "<span style='color: #333333;'>%{z:.4f}</span>"
            "<extra></extra>"
        )

    # Apply styling parameters
    fig.update_traces(
        hovertemplate=hovertemplate,
        hoverlabel=dict(bgcolor="#FFFFFF", bordercolor="#D3D3D3", font_size=12, font_family="Arial, sans-serif"),
        marker=dict(size=5, opacity=0.85)
    )

    # Axis tick calculations for simplified ('-', '0', '+') display
    x_vals, x_texts = get_simplified_ticks(pca_df["PC1"])
    y_vals, y_texts = get_simplified_ticks(pca_df["PC2"])
    z_vals, z_texts = get_simplified_ticks(pca_df["PC3"])

    # Configure 3D Scene Axes and Layout
    fig.update_scenes(
        xaxis=dict(
            title=f"PC1 ({var_explained[0]:.2f}%)",
            showspikes=False,
            zeroline=False,
            tickmode="array",
            tickvals=x_vals,
            ticktext=x_texts,
        ),
        yaxis=dict(
            title=f"PC2 ({var_explained[1]:.2f}%)",
            showspikes=False,
            zeroline=False,
            tickmode="array",
            tickvals=y_vals,
            ticktext=y_texts,
        ),
        zaxis=dict(
            title=f"PC3 ({var_explained[2]:.2f}%)",
            showspikes=False,
            zeroline=False,
            tickmode="array",
            tickvals=z_vals,
            ticktext=z_texts,
        )
    )

    # Add true 3D central axes passing through origin (0, 0, 0)
    x_min, x_max = pca_df["PC1"].min(), pca_df["PC1"].max()
    y_min, y_max = pca_df["PC2"].min(), pca_df["PC2"].max()
    z_min, z_max = pca_df["PC3"].min(), pca_df["PC3"].max()

    # X-axis center line (y=0, z=0)
    fig.add_trace(go.Scatter3d(
        x=[x_min, x_max], y=[0, 0], z=[0, 0],
        mode="lines",
        line=dict(color="black", width=3),
        hoverinfo="none",
        showlegend=False
    ))

    # Y-axis center line (x=0, z=0)
    fig.add_trace(go.Scatter3d(
        x=[0, 0], y=[y_min, y_max], z=[0, 0],
        mode="lines",
        line=dict(color="black", width=3),
        hoverinfo="none",
        showlegend=False
    ))

    # Z-axis center line (x=0, y=0)
    fig.add_trace(go.Scatter3d(
        x=[0, 0], y=[0, 0], z=[z_min, z_max],
        mode="lines",
        line=dict(color="black", width=3),
        hoverinfo="none",
        showlegend=False
    ))

    fig.update_layout(
        title=dict(x=0.02, y=0.95),
        margin=dict(l=0, r=0, b=0, t=60),
        height=650
    )

    return fig

if __name__ == "__main__":
    fig = run_pca("synthetic_exome.mt")
    output_html = "pca_3d_report.html"
    fig.write_html(output_html)
    print(f"Success! 3D PCA plot saved to {output_html}")