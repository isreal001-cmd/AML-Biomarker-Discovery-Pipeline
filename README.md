# AML Biomarker Discovery Pipeline

A reproducible RNA sequencing framework for discovering and validating diagnostic biomarkers in **Acute Myeloid Leukemia (AML)** using differential expression analysis, machine learning, pathway enrichment, and external validation.

---
  
  ## Overview
  
  The AML Biomarker Discovery Pipeline is a modular bioinformatics framework developed to perform end to end RNA sequencing analysis from raw featureCounts output to externally validated biomarker candidates.

The framework combines statistical analysis, machine learning, visualization, and external validation into a reproducible workflow designed for research and publication.

---
  
  ## Key Features
  
  ### RNA Sequencing Analysis
  
  * FeatureCounts count matrix merging
* Sample metadata generation
* Differential expression analysis using DESeq2
* Variance Stabilizing Transformation (VST)

### Exploratory Data Analysis

* Principal Component Analysis (PCA)
* Volcano plot generation
* Heatmap visualization
* Batch effect assessment

### Biomarker Discovery

* Differential gene prioritization
* Functional enrichment analysis
* ROC curve analysis
* Biomarker expression visualization
* Random Forest based feature importance

### External Validation

* Validation dataset preparation
* TCGA RNA sequencing integration
* Independent biomarker validation
* Validation visualizations

### Reproducibility

* Automated package checking
* Global configuration management
* Session information logging
* Publication summary generation
* Modular pipeline execution
* Validation specific runner

---
  
  # Pipeline Workflow
  
  ```text
FeatureCounts Output

│

▼

Merge Count Matrix

│

▼

Prepare Metadata

│

▼

Differential Expression Analysis

│

▼

Variance Stabilizing Transformation

│

▼

Quality Assessment

│

▼

Biomarker Discovery

│

▼

Functional Enrichment

│

▼

Machine Learning Validation

│

▼

External Validation (TCGA)

│

▼

Publication Ready Outputs
```

---
  
  # Repository Structure
  
  ```text
AML_Project/
  
  ├── data/
  │
├── results/
  │
├── scripts/
  │
├── config.R
│
├── plotting_functions.R
│
├── validation_functions.R
│
├── run_pipeline.R
│
├── run_validation.R
│
└── README.md
```

---
  
  # Pipeline Modules
  
  | Step | Module                                |
  | ---- | ------------------------------------- |
  | 00   | Package Installation and Verification |
  | 01   | Merge featureCounts Output            |
  | 02   | Metadata Preparation                  |
  | 03   | Differential Expression Analysis      |
  | 04   | VST Normalization                     |
  | 05   | Principal Component Analysis          |
  | 06   | Volcano Plot                          |
  | 07   | Heatmap                               |
  | 08   | Batch Assessment                      |
  | 09   | Biomarker Candidate Selection         |
  | 10   | Functional Enrichment                 |
  | 11   | ROC Analysis                          |
  | 12   | Biomarker Expression                  |
  | 13   | Random Forest Validation              |
  | 14A  | Validation Dataset Preparation        |
  | 14B  | TCGA Data Processing                  |
  | 15   | External Validation                   |
  | 16   | Validation Visualization              |
  | 17   | Session Information                   |
  | 18   | Publication Summary                   |
  
  ---
  
  # Installation
  
  Clone the repository

```bash
git clone https://github.com/YOUR_USERNAME/AML_Biomarker_Discovery_Pipeline.git
```

Open the project in RStudio.

---
  
  # Running the Complete Pipeline
  
  ```r
source("run_pipeline.R")
```

This executes every analysis module from preprocessing through external validation.

---
  
  # Running Only External Validation
  
  ```r
source("run_validation.R")
```

This executes only the validation workflow.

---
  
  # Pipeline Output
  
  The framework automatically generates

* Differential expression results
* Normalized expression matrices
* PCA figures
* Volcano plots
* Heatmaps
* Functional enrichment analysis
* ROC analysis
* Random Forest results
* External validation outputs
* Publication summary
* Session information
* Pipeline logs

---
  
  # Reproducibility
  
  To ensure reproducible research, every execution records

* R session information
* Installed package versions
* Pipeline version
* Execution logs

---
  
  # Current Version
  
  **Version:** 1.0.0

Current release includes

* Complete RNA sequencing workflow
* Automated pipeline execution
* External TCGA validation
* Publication summary generation
* Reproducible project configuration

---
  
  # Future Development
  
  Planned improvements include

* GEO validation module
* GTEx validation
* User supplied dataset validation
* Interactive HTML reports
* Shiny dashboard
* Multi cohort consensus biomarker analysis

---
  
  # Citation
  
  A formal software citation will be provided through the repository using a `CITATION.cff` file.

---
  
  # License
  
  This project will be distributed under the MIT License.

---
  
  # Author
  
  **Isreal Oluwafemi Abiodun**
  
  Bioinformatician

Computational Biology

Cancer Genomics

RNA Sequencing Analysis

Machine Learning

---
  
  If you use this framework in your research, please consider citing the repository once the official release becomes available.
