# AML Biomarker Discovery Framework

A reproducible RNA-seq analysis framework for identifying diagnostic and prognostic biomarkers in Acute Myeloid Leukemia using differential expression analysis, machine learning, functional enrichment, and external validation.

---

## Overview

The AML Biomarker Discovery Framework is a modular R based pipeline designed to perform end to end transcriptomic analysis from raw featureCounts output to publication ready figures and biomarker validation.

The framework emphasizes reproducibility, modularity, automation, and publication quality outputs while following best practices for computational biology research.

---

## Features

- Automated featureCounts merging
- Metadata generation and validation
- Differential expression analysis using DESeq2
- Variance Stabilizing Transformation
- Principal Component Analysis
- Volcano plot generation
- Heatmap visualization
- Batch effect assessment
- Biomarker candidate selection
- Functional enrichment analysis
- ROC curve analysis
- Random Forest validation
- External validation using TCGA
- Publication ready figures
- Automated analysis reports

---

## Pipeline Workflow

```text
featureCounts
      │
      ▼
Metadata Preparation
      │
      ▼
DESeq2 Differential Expression
      │
      ▼
Variance Stabilizing Transformation
      │
      ▼
PCA
Volcano Plot
Heatmap
      │
      ▼
Batch Assessment
      │
      ▼
Biomarker Selection
      │
      ▼
Functional Enrichment
      │
      ▼
ROC Analysis
      │
      ▼
Random Forest Validation
      │
      ▼
External Validation
      │
      ▼
Publication Figures
      │
      ▼
Publication Report
```

---

## Project Structure

```text
AML_Project/

├── config.R
├── run_pipeline.R
├── run_validation.R
├── scripts/
├── functions/
├── data/
├── results/
├── README.md
├── LICENSE
├── CHANGELOG.md
├── CITATION.cff
```

---

## Installation

Clone the repository

```bash
git clone https://github.com/isreal001-cmd/AML-Biomarker-Discovery-Pipeline.git
```

Open the project in RStudio.

Run the complete analysis

```r
source("run_pipeline.R")
```

Run only the validation pipeline

```r
source("run_validation.R")
```

---

## Outputs

The framework automatically generates

- Differential expression tables
- PCA plots
- Volcano plots
- Heatmaps
- Functional enrichment results
- ROC curves
- Biomarker validation results
- Publication ready figures
- Pipeline logs
- Session information

---

## Reproducibility

The framework records

- Pipeline version
- R version
- Package versions
- Session information
- Processing logs

to ensure complete reproducibility.

---

## Citation

If you use this framework in your research, please cite the accompanying manuscript after publication.

Citation information will be updated upon publication.

---

## License

This project is distributed under the MIT License.

---

## Author

**Isreal Oluwafemi Abiodun**

Bioinformatician | Computational Biology | Cancer Genomics | RNA-seq Analysis

---

## Version

Current Release

**Version 1.0.0**