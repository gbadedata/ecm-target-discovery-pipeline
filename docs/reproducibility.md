# Reproducibility

## How to Reproduce This Analysis

### Prerequisites

- R >= 4.4 with Bioconductor packages (see `install_r_packages.R`)
- Python 3.12 with conda (see `environment.yml`)
- Nextflow >= 24.0 (optional, for full workflow orchestration)
- Git

### Step 1: Clone and set up environments

```bash
git clone https://github.com/gbadedata/ecm-target-discovery-pipeline.git
cd ecm-target-discovery-pipeline

# Python environment
conda env create -f environment.yml
conda activate ecm-pipeline

# R packages
Rscript install_r_packages.R
```

### Step 2: Download data

Follow instructions in `metadata/data_sources.md` to download CPTAC PDAC data. Place files in `data/raw/`.

### Step 3: Validate data

```bash
python src/python/validate_data.py
```

This checks file integrity against `metadata/checksums.txt` and validates the sample manifest.

### Step 4: Run the analysis

**Option A: Nextflow (recommended)**
```bash
nextflow run main.nf -params-file params.yaml
```

**Option B: Step by step**
Each analysis step can be run independently. See the numbered scripts in `src/r/` and `src/python/`.

### Step 5: Generate reports

```bash
quarto render reports/technical_report.qmd
```

## Reproducibility Controls

| Control | Implementation |
|---|---|
| Random seeds | Set in `params.yaml` (42), applied in every R and Python script |
| Package versions | Locked via `environment.yml` (Python) and `install_r_packages.R` (R) |
| Parameters | All thresholds and settings in `params.yaml`, not hard-coded |
| Data provenance | Sources, URLs, and checksums in `metadata/` |
| Environment | Conda environment for Python, R package versions logged |
| Workflow | Nextflow execution report captures runtime, resources, and versions |
