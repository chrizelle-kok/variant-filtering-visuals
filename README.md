# variant-filtering-visuals


This repository contains Python scripts to generate independant visualizations in the form of interactive `.html` plot outputs for evaluating human germline Next Generation Sequencing (NGS) data in variant filtering workflows.

The scripts extract and process variant data loaded into a [Hail](https://hail.is/) MatrixTable (converting relevant metrics into [pandas](https://pandas.pydata.org/) DataFrames for downstream processing) and output standalone, interactive figures generated with [Plotly](https://plotly.com/python/) and [Bokeh](https://docs.bokeh.org/).


## 📁 Repository Structure

```text
├── synthetic_exomes.mt/     # Pre-generated Hail MatrixTable test input (from syntetic exomes)
├── data/
│   └── my_dataset.mt/       # hail MatrixTable folder containing sample QC metrics
├── generate_pca.py          # Generates interactive PCA plot
├── generate_titv_plot.py    # Generates Ti/Tv ratio plots
├── generate_sunburst.py     # Generates variant consequence sunburst chart
├── pca_3d_report.html       # Pre-rendered interactive PCA output
├── titv_qc_report.html      # Pre-rendered interactive Ti/Tv output
├── filtering_sunburst.html  # Pre-rendered interactive sunburst output
└── README.md                # Project documentation
```

## 🛠️ Environment & Prerequisites

### Dependencies
Install the required packages using `pip`:

```bash
pip install hail pandas bokeh plotly
```
*Note: Hail requires Python 3.8+ and Java 8 or 11.*

## 🚀 Workflow & Execution

Each script runs independently and processes the included synthetic_exomes.mt MatrixTable. Where applicable, metrics are exported from Hail into pandas DataFrames to streamline data manipulation prior to rendering.

To execute or re-generate any visualization, run the respective script:

```bash
# Run PCA analysis
python generate_pca.py

# Run Transistion/Transversion analysis
python generate_titv_plot.py

# Run variant retained/removed counts
python generate_sunburst.py
```
Each script outputs or overwrites its corresponding `.html` file. Open any generated `.html` file in your browser to inspect the interactive figure.

## 👤 Author & Acknowledgments
* **Contributor:** Chrizelle Kok

* **Project Role:** Designed visualization scripts for variant filtering outputs for integration into the Quarto reporting framework.

* **Main Workflow Repository:** Developed under the supervision of [@phelelani](https://github.com/phelelani). This project complements the main variant filtering workflow hosted at [nf-variantfilt](https://github.com/phelelani/nf-variantfilt).
* **Disclaimer:**  **Google Gemini** was utilized to assist in writing and refining portions of the Python code for data processing and HTML visualization generation as well as repository setup instructions. All generated code was independently executed, tested, and validated by the author.