# Multi-Gene Correlations to Signature

A containerized R tool for computing Pearson correlations between a gene panel and a signature score across user-defined sample categories. Generates per-category bar plots, scatter plots, and an aggregate heatmap.

## Example Outputs

Representative run (TCGA PRAD; categories TP vs NT) with default plotting parameters:

- Heatmap (correlations by gene and category):
  ![Correlation heatmap](docs/images/heatmap.png)
- Per-category bar plot:
  ![Bar plot](docs/images/barplot_all.png)
- Per-category scatter (signature vs gene expression):
  ![Scatter plot](docs/images/scatter_all.png)

These screenshots use metadata values `TP` and `NT`; match your own metadata values in `--categories` (and optionally `--rename_categories` if you want display labels to differ from underlying values).

## Features

- Correlate gene expression against a signature score
- Per-category correlation analysis with statistical testing (FDR correction) and scatter facets annotated with Pearson r/p values
- Automated gene symbol updates via L2P
- Flexible input formats (CSV/TSV auto-detection)
- Comprehensive outputs: tables, bar plots, heatmap, scatter plots
- Optional clinical correlation module to relate gene/signature scores to numeric or categorical metadata
- Command-line interface for batch processing
- Available as a Docker/OCI image (GHCR). Singularity/Apptainer is supported by pulling the Docker image (singularity pull … docker://…) to produce a .sif.

## Container Images

**Docker:**
```bash
docker pull ghcr.io/nidap-community/multi-gene-correlations:latest
```

**Singularity:**
```bash
singularity pull multi-gene-correlations.sif docker://ghcr.io/nidap-community/multi-gene-correlations:latest
```

## Quick Start

### Docker

```bash
docker run --rm \
  -v /path/to/data:/data \
  -v /path/to/output:/output \
  ghcr.io/nidap-community/multi-gene-correlations:latest \
  --counts /data/counts.tsv \
  --metadata /data/metadata.csv \
  --gene_column Gene \
  --sample_column SampleID \
  --genes "GATM,GAMT,CKM,SLC6A8" \
  --category_column Tissue \
  --categories "Normal,Tumor,all" \
  --signature_name "MySignature" \
  --signature_genes "PTEN,SPRY2,CKM" \
  --output_dir /output
```

### Singularity

```bash
singularity run \
  --bind /path/to/data:/data \
  --bind /path/to/output:/output \
  multi-gene-correlations.sif \
  --counts /data/counts.tsv \
  --metadata /data/metadata.csv \
  --gene_column Gene \
  --sample_column SampleID \
  --genes "GATM,GAMT,CKM,SLC6A8" \
  --category_column Tissue \
  --categories "Normal,Tumor,all" \
  --signature_name "MySignature" \
  --signature_genes "PTEN,SPRY2,CKM" \
  --output_dir /output
```

## Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--counts` | Path to normalized counts file (genes × samples). CSV or TSV format. |
| `--metadata` | Path to sample metadata file with category assignments. |
| `--gene_column` | Column name in counts file containing gene identifiers. |
| `--sample_column` | Column name in metadata containing sample identifiers. |
| `--genes` | Comma-separated list of genes to correlate (or path to file with one gene per line). |
| `--category_column` | Column name in metadata defining sample categories. |
| `--categories` | Comma-separated list of category values to analyze; include the literal lowercase `all` when you want aggregate plots/tables spanning every requested category. |
| `--signature_genes` | Comma-separated list of genes defining the signature score (averaged). |

## Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--signature_name` | "Signature" | Display name for the signature in plots. |
| `--samples` | (all) | Comma-separated list to restrict analysis to specific samples. |
| `--output_dir` | `/output` | Base directory for all output files. |
| `--barplots_dir` | `<output_dir>/barplots` | Directory for correlation bar plots. |
| `--heatmap_dir` | `<output_dir>/heatmap` | Directory for the correlation heatmap. |
| `--scatter_dir` | `<output_dir>/scatter` | Directory for scatter plots. |
| `--tables_dir` | `<output_dir>/tables` | Directory for correlation tables. |
| `--margin_adjust` | `-11` | Heatmap row label margin adjustment (cm). |
| `--show_row_dendrogram` | `false` | Display row dendrogram in heatmap (`true`/`false`). |
| `--heatmap_width` | `6` | Width (cm) per heatmap column. |
| `--add_text_to_heatmap` | `true` | Print correlation values in heatmap cells. |
| `--rename_categories` | (none) | Comma-separated display labels for categories (in order). |
| `--sum_duplicates` | `false` | Sum duplicate gene rows instead of selecting maximum. |
| `--clinical_columns` | (none) | Comma-separated metadata columns to correlate against; accepts CSV file paths. |
| `--clinical_use_signature` | `true` | Correlate the computed signature score against the requested clinical columns. |
| `--clinical_use_genes` | `false` | Correlate each gene in `--genes` with the requested clinical columns (may generate large tables). |
| `--clinical_box_colors` | (auto) | Comma-separated color codes (hex or names) for categorical clinical box plots; defaults to 20 distinct colors. |
 
### Clinical correlations

Enable clinical associations by pointing `--clinical_columns` to one or more metadata columns:

```
--clinical_columns "Stage,PSA" --clinical_use_signature true --clinical_use_genes true
```

- Numeric columns use Pearson correlations (signature vs measurement and per-gene vs measurement).
- Categorical columns are expanded into one-hot indicators per level (point-biserial correlations). Signature-level box plots + bar plots, as well as **gene vs numeric clinical scatter plots with Pearson r and p annotations**, are saved under `output/clinical/`.
- Customize categorical box-plot colors via `--clinical_box_colors`. Provide at least as many colors as levels (hex like `#1f77b4` or names); otherwise a default palette of 20 distinct colors is used.
- Signature vs categorical clinical box plots annotate each level with the exact Pearson r and p values derived from the indicator correlations so reviewers can match the boxes to the tables.
- Numeric clinical scatter plots annotate each gene facet and the signature plot with the exact Pearson r and p values to mirror the table outputs.
- Results are written to `tables/clinical_signature_correlations.csv` (signature) and `tables/clinical_gene_correlations.csv` (genes, optional).

## Bind mount layout

The container follows the Container-as-a-Function model and expects explicit bind mounts:

- `/data`: read-only directory containing counts, metadata, and optional gene lists.
- `/output`: writable directory where the pipeline writes `tables/`, `barplots/`, `scatter/`, `heatmap/`, and `clinical/` artifacts.

Always map host paths into those locations. The helper script below enforces the layout automatically.

### Singularity helper script

Use [container/run_clinical_plots.sh](container/run_clinical_plots.sh) to standardize Singularity invocations:

```bash
DATA_DIR=$PWD/test_data \
OUTPUT_DIR=$PWD/test_output \
container/run_clinical_plots.sh \
  ~/images/multi-gene-correlations.sif \
  --counts /data/counts.tsv \
  --metadata /data/metadata.csv \
  --gene_column Gene \
  --sample_column SampleID \
  --genes "GATM,GAMT,CKM,SLC6A8" \
  --category_column Tissue \
  --categories "Normal,Tumor,all" \
  --signature_name Creatine \
  --signature_genes "PTEN,SPRY2" \
  --clinical_columns "Stage,PSA" \
  --output_dir /output
```

`DATA_DIR`, `OUTPUT_DIR`, and optional `EXTRA_BINDS` environment variables control the mounts. The script also propagates `R_LIBS_USER=/usr/local/lib/R/site-library` so Singularity runs resolve pre-installed R packages even when `HOME` is remapped.

## Input File Formats

### Counts File
Tab- or comma-separated file with genes in rows and samples in columns:

```
Gene    Sample1    Sample2    Sample3
PTEN    123.4      456.7      789.0
CKM     234.5      567.8      890.1
```

### Metadata File
Tab- or comma-separated file with sample annotations:

```
SampleID,Tissue,Patient
Sample1,Normal,P001
Sample2,Tumor,P001
Sample3,Tumor,P002
```

## Output Files

The tool generates the following outputs:

```
output/
├── tables/
│   ├── pathways.csv                      # Gene vs signature correlations
│   ├── clinical_signature_correlations.csv   # Signature vs clinical metadata (optional)
│   └── clinical_gene_correlations.csv        # Gene panel vs clinical metadata (optional)
├── metadata/
│   └── run_parameters.json   # Run configuration and gene mapping details
├── barplots/
|   ├── barplot_Category1.png
|   ├── barplot_Category2.png
|   └── barplot_all.png (if the literal lowercase `all` was included in --categories)
├── heatmap/
│   └── heatmap.png           # Combined correlation heatmap
├── scatter/
|   ├── scatter_Category1.png # Signature vs gene expression with r/p annotations per facet
|   ├── scatter_Category2.png
|   └── scatter_all.png (if the literal lowercase `all` was included in --categories)
└── clinical/ (optional)
  ├── clinical_scatter_PSA.png           # Signature vs numeric clinical variable
  ├── clinical_scatter_signature_vs_PSA.png # Signature-specific scatter (numeric clinical)
  ├── clinical_scatter_genes_vs_PSA.png      # Gene expression vs numeric clinical with r/p labels
  ├── clinical_box_Stage.png             # Signature vs categorical clinical variable
  ├── clinical_bar_signature_vs_clinical.png  # Bar plot of signature vs all clinical targets
  └── clinical_bar_genes_vs_Stage_Stage1.png   # Gene panel correlations per clinical target (one file per target)
```

## Gene List Files

Instead of comma-separated lists, you can provide file paths:

**genes.txt:**
```
GATM
GAMT
CKM
SLC6A8
```

**Usage:**
```bash
docker run --rm \
  -v /path/to/data:/data \
  -v /path/to/output:/output \
  ghcr.io/nidap-community/multi-gene-correlations:latest \
  --counts /data/counts.tsv \
  --metadata /data/metadata.csv \
  --genes /data/genes.txt \
  --signature_genes /data/signature.txt \
  ...
```

## Advanced Usage

### Custom Output Directories

Specify separate directories for each output type:

```bash
docker run --rm \
  -v /path/to/data:/data \
  -v /path/to/output:/output \
  ghcr.io/nidap-community/multi-gene-correlations:latest \
  --counts /data/counts.tsv \
  --metadata /data/metadata.csv \
  --output_dir /output \
  --barplots_dir /output/plots/bars \
  --heatmap_dir /output/plots/heatmap \
  --tables_dir /output/results \
  ...
```

### Debugging with Shell Access

**Docker:**
```bash
docker run --rm -it \
  -v /path/to/data:/data \
  --entrypoint bash \
  ghcr.io/nidap-community/multi-gene-correlations:latest
```

**Singularity:**
```bash
singularity shell \
  --bind /path/to/data:/data \
  multi-gene-correlations.sif
```

### Running on HPC with Singularity

Example SLURM job script:

```bash
#!/bin/bash
#SBATCH --job-name=correlations
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00

module load singularity

singularity run \
  --bind /data/project:/data \
  --bind /scratch/output:/output \
  /path/to/multi-gene-correlations.sif \
  --counts /data/counts.tsv \
  --metadata /data/metadata.csv \
  --gene_column Gene \
  --sample_column SampleID \
  --genes "GATM,GAMT,CKM,SLC6A8,AR,PTEN,SPRY2" \
  --category_column Condition \
  --categories "Control,Treatment,all" \
  --signature_name "Creatine" \
  --signature_genes "PTEN,SPRY2,CKM" \
  --output_dir /output
```

## Local R Workflow (renv)

If you prefer to run the analysis scripts directly on your workstation (outside Docker/Singularity), a pinned `renv.lock` is now provided.

1. Install `renv` once: `Rscript -e "install.packages('renv', repos='https://cran.r-project.org')"`
2. Restore the project library: `Rscript -e "renv::restore(prompt = FALSE)"`
3. For each interactive R session, run `source('renv/activate.R')` from the repo root to place the renv library first in `.libPaths()`.
4. Invoke `Rscript entrypoint.R ...` (same arguments shown above) or call `Multi_Gene_Correlations_to_Signature()` directly after `source('renv/activate.R'); source('Multi_Gene_Correlations_to_Signature.R')`.

Notes:
- The repo does **not** auto-activate renv to avoid interfering with container runs; always source `renv/activate.R` explicitly when working locally.
- `renv.lock` targets R 4.5.x; use the closest release you have available (matching the container whenever possible).
- After changing R package dependencies, run `renv::snapshot()` so the lockfile stays in sync.
- You can mirror the container smoke test locally via `Rscript entrypoint.R --counts test_data/counts.tsv ... --output_dir test_output/local_smoke` once renv is restored.

## Building from Source

### Clone Repository

```bash
git clone https://github.com/NIDAP-Community/multi-gene-correlations.git
cd multi-gene-correlations
```

### Build Docker Image

```bash
docker build -t multi-gene-correlations:latest -f container/Dockerfile .
# On Apple Silicon, add: --platform linux/amd64
```

### Convert to Singularity

```bash
singularity build multi-gene-correlations.sif docker-daemon://multi-gene-correlations:latest
```

## Requirements

- **R version:** 4.5.2
- **Key R packages:** dplyr, tidyr, stringr, ggplot2, ComplexHeatmap, circlize, RColorBrewer, l2p, l2psupp
- **System:** Linux container runtime (Docker or Singularity)

## Citation

This tool implements the NIDAP "Multi-Gene Correlations to Signature" template (version 14).

Template ID: `ri.vector.main.template.363cef88-16c4-433b-a709-58aa605f0958`

## License

This software is provided as-is for research use.

## Support

For issues or questions, please open an issue at:
https://github.com/NIDAP-Community/multi-gene-correlations/issues

## Version

**Current Version:** 1.0.1  
**Base Image:** rocker/r-ver:4.5.2  
**Template Version:** 14
