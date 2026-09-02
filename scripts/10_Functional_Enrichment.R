##############################################################
# AML Biomarker Discovery Pipeline
#
# Script:
# 10_Functional_Enrichment.R
#
# Purpose:
# Perform Gene Ontology and KEGG enrichment analysis
# of significant AML biomarker candidates.
#
# Author:
# Isreal Oluwafemi Abiodun
#
# Version:
# 1.1.0
##############################################################


##############################################################
# Load Configuration
##############################################################

source("config.R")

suppressPackageStartupMessages({
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(DOSE)
  
})


##############################################################
# Stage Header
##############################################################

cat("\n")
cat("=========================================\n")
cat("Stage 10 : Functional Enrichment Analysis\n")
cat("=========================================\n\n")

start_time <- Sys.time()


##############################################################
# Directories
##############################################################

ENRICHMENT_DIR <- file.path(
  RESULTS_DIR,
  "enrichment"
)

PLOT_DIR <- file.path(
  ENRICHMENT_DIR,
  "plots"
)

TABLE_DIR <- file.path(
  ENRICHMENT_DIR,
  "tables"
)

dir.create(
  ENRICHMENT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  PLOT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


##############################################################
# Input File
##############################################################

BIOMARKER_FILE <- file.path(
  RESULTS_DIR,
  "biomarkers",
  "Ranked_Biomarkers.csv"
)


##############################################################
# Validate Input
##############################################################

if (!file.exists(BIOMARKER_FILE)) {
  
  stop(
    "Ranked biomarkers not found.\n",
    "Expected file:\n",
    BIOMARKER_FILE
  )
  
}


##############################################################
# Load Biomarkers
##############################################################

cat("Loading ranked biomarkers...\n")

biomarkers <- readr::read_csv(
  BIOMARKER_FILE,
  show_col_types = FALSE
)

cat(
  "Genes Loaded :",
  nrow(biomarkers),
  "\n\n"
)


##############################################################
# Remove Missing Gene IDs
##############################################################

biomarkers <- biomarkers |>
  dplyr::filter(
    !is.na(GeneID),
    GeneID != ""
  )


##############################################################
# Remove Ensembl Version Numbers
##############################################################

biomarkers$GeneID <- sub(
  "\\..*$",
  "",
  biomarkers$GeneID
)


##############################################################
# Remove Duplicate Ensembl IDs
##############################################################

biomarkers <- biomarkers |>
  dplyr::distinct(
    GeneID,
    .keep_all = TRUE
  )


cat(
  "Genes Remaining :",
  nrow(biomarkers),
  "\n\n"
)


##############################################################
# Convert Ensembl IDs to Entrez IDs
##############################################################

cat("Converting Ensembl IDs to Entrez IDs...\n")

gene_map <- clusterProfiler::bitr(
  biomarkers$GeneID,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)


##############################################################
# Remove Duplicate Mappings
##############################################################

gene_map <- gene_map |>
  dplyr::distinct(
    ENSEMBL,
    .keep_all = TRUE
  )


##############################################################
# Create Entrez Gene List
##############################################################

gene_ids <- unique(
  gene_map$ENTREZID
)

cat(
  "Mapped genes :",
  length(gene_ids),
  "\n\n"
)


##############################################################
# Mapping Statistics
##############################################################

input_gene_count <- nrow(biomarkers)

mapped_gene_count <- length(
  unique(gene_map$ENSEMBL)
)

mapping_rate <- (
  mapped_gene_count /
    input_gene_count
) * 100

unmapped_gene_count <- (
  input_gene_count -
    mapped_gene_count
)


cat(
  "Mapping Rate :",
  round(mapping_rate, 2),
  "%\n"
)

cat(
  "Unmapped Genes :",
  unmapped_gene_count,
  "\n\n"
)


##############################################################
# Save Gene Mapping
##############################################################

readr::write_csv(
  gene_map,
  file.path(
    TABLE_DIR,
    "Ensembl_to_Entrez_Mapping.csv"
  )
)


##############################################################
# GO Biological Process
##############################################################

cat(
  "Running GO Biological Process...\n"
)

go_bp <- clusterProfiler::enrichGO(
  
  gene = gene_ids,
  
  OrgDb = org.Hs.eg.db,
  
  ont = "BP",
  
  keyType = "ENTREZID",
  
  pAdjustMethod = "BH",
  
  pvalueCutoff = 0.05,
  
  qvalueCutoff = 0.20,
  
  readable = TRUE
  
)


##############################################################
# GO Cellular Component
##############################################################

cat(
  "Running GO Cellular Component...\n"
)

go_cc <- clusterProfiler::enrichGO(
  
  gene = gene_ids,
  
  OrgDb = org.Hs.eg.db,
  
  ont = "CC",
  
  keyType = "ENTREZID",
  
  pAdjustMethod = "BH",
  
  pvalueCutoff = 0.05,
  
  qvalueCutoff = 0.20,
  
  readable = TRUE
  
)


##############################################################
# GO Molecular Function
##############################################################

cat(
  "Running GO Molecular Function...\n"
)

go_mf <- clusterProfiler::enrichGO(
  
  gene = gene_ids,
  
  OrgDb = org.Hs.eg.db,
  
  ont = "MF",
  
  keyType = "ENTREZID",
  
  pAdjustMethod = "BH",
  
  pvalueCutoff = 0.05,
  
  qvalueCutoff = 0.20,
  
  readable = TRUE
  
)


##############################################################
# KEGG Pathway Analysis
##############################################################

cat(
  "Running KEGG pathway analysis...\n"
)

kegg <- tryCatch(
  
  clusterProfiler::enrichKEGG(
    
    gene = gene_ids,
    
    organism = "hsa",
    
    pvalueCutoff = 0.05,
    
    pAdjustMethod = "BH"
    
  ),
  
  error = function(e) {
    
    cat(
      "KEGG analysis unavailable:\n",
      conditionMessage(e),
      "\n"
    )
    
    NULL
    
  }
  
)


##############################################################
# Save Complete GO Results
##############################################################

if (!is.null(go_bp) && nrow(go_bp) > 0) {
  
  readr::write_csv(
    
    as.data.frame(go_bp),
    
    file.path(
      ENRICHMENT_DIR,
      "GO_BP.csv"
    )
    
  )
  
}


if (!is.null(go_cc) && nrow(go_cc) > 0) {
  
  readr::write_csv(
    
    as.data.frame(go_cc),
    
    file.path(
      ENRICHMENT_DIR,
      "GO_CC.csv"
    )
    
  )
  
}


if (!is.null(go_mf) && nrow(go_mf) > 0) {
  
  readr::write_csv(
    
    as.data.frame(go_mf),
    
    file.path(
      ENRICHMENT_DIR,
      "GO_MF.csv"
    )
    
  )
  
}


##############################################################
# Save KEGG Results
##############################################################

if (!is.null(kegg) && nrow(kegg) > 0) {
  
  readr::write_csv(
    
    as.data.frame(kegg),
    
    file.path(
      ENRICHMENT_DIR,
      "KEGG.csv"
    )
    
  )
  
}


##############################################################
# Save RDS Objects
##############################################################

saveRDS(
  
  go_bp,
  
  file.path(
    ENRICHMENT_DIR,
    "GO_BP.rds"
  )
  
)


saveRDS(
  
  go_cc,
  
  file.path(
    ENRICHMENT_DIR,
    "GO_CC.rds"
  )
  
)


saveRDS(
  
  go_mf,
  
  file.path(
    ENRICHMENT_DIR,
    "GO_MF.rds"
  )
  
)


if (!is.null(kegg)) {
  
  saveRDS(
    
    kegg,
    
    file.path(
      ENRICHMENT_DIR,
      "KEGG.rds"
    )
    
  )
  
}


##############################################################
# Save Top 30 Enrichment Results
##############################################################

save_top_results <- function(
    enrichment_object,
    filename,
    n = 30
) {
  
  if (
    !is.null(enrichment_object) &&
    nrow(enrichment_object) > 0
  ) {
    
    result_table <- as.data.frame(
      enrichment_object
    )
    
    result_table <- result_table |>
      dplyr::arrange(
        p.adjust
      ) |>
      dplyr::slice_head(
        n = n
      )
    
    readr::write_csv(
      result_table,
      file.path(
        TABLE_DIR,
        filename
      )
    )
    
  }
  
}


save_top_results(
  go_bp,
  "Top30_GO_BP.csv"
)


save_top_results(
  go_cc,
  "Top30_GO_CC.csv"
)


save_top_results(
  go_mf,
  "Top30_GO_MF.csv"
)


save_top_results(
  kegg,
  "Top30_KEGG.csv"
)


##############################################################
# Generate GO BP Dotplot
##############################################################

if (
  !is.null(go_bp) &&
  nrow(go_bp) > 0
) {
  
  png(
    
    filename = file.path(
      PLOT_DIR,
      "GO_BP_Dotplot.png"
    ),
    
    width = 1800,
    height = 1400,
    res = 200
    
  )
  
  print(
    
    enrichplot::dotplot(
      go_bp,
      showCategory = 20
    )
    
  )
  
  dev.off()
  
}


##############################################################
# Generate GO CC Dotplot
##############################################################

if (
  !is.null(go_cc) &&
  nrow(go_cc) > 0
) {
  
  png(
    
    filename = file.path(
      PLOT_DIR,
      "GO_CC_Dotplot.png"
    ),
    
    width = 1800,
    height = 1400,
    res = 200
    
  )
  
  print(
    
    enrichplot::dotplot(
      go_cc,
      showCategory = 20
    )
    
  )
  
  dev.off()
  
}


##############################################################
# Generate GO MF Dotplot
##############################################################

if (
  !is.null(go_mf) &&
  nrow(go_mf) > 0
) {
  
  png(
    
    filename = file.path(
      PLOT_DIR,
      "GO_MF_Dotplot.png"
    ),
    
    width = 1800,
    height = 1400,
    res = 200
    
  )
  
  print(
    
    enrichplot::dotplot(
      go_mf,
      showCategory = 20
    )
    
  )
  
  dev.off()
  
}


##############################################################
# Generate KEGG Dotplot
##############################################################

if (
  !is.null(kegg) &&
  nrow(kegg) > 0
) {
  
  png(
    
    filename = file.path(
      PLOT_DIR,
      "KEGG_Dotplot.png"
    ),
    
    width = 1800,
    height = 1400,
    res = 200
    
  )
  
  print(
    
    enrichplot::dotplot(
      kegg,
      showCategory = 20
    )
    
  )
  
  dev.off()
  
}


##############################################################
# Save Summary Report
##############################################################

SUMMARY_FILE <- file.path(
  
  ENRICHMENT_DIR,
  
  "Stage10_Enrichment_Summary.txt"
  
)


summary_lines <- c(
  
  "AML BIOMARKER DISCOVERY PIPELINE",
  
  "Stage 10 : Functional Enrichment Analysis",
  
  "",
  
  "=========================================",
  
  "SUMMARY",
  
  "=========================================",
  
  "",
  
  paste(
    "Input genes:",
    input_gene_count
  ),
  
  paste(
    "Mapped genes:",
    mapped_gene_count
  ),
  
  paste(
    "Unmapped genes:",
    unmapped_gene_count
  ),
  
  paste(
    "Mapping rate:",
    round(mapping_rate, 2),
    "%"
  ),
  
  "",
  
  paste(
    "GO Biological Process terms:",
    ifelse(
      is.null(go_bp),
      0,
      nrow(go_bp)
    )
  ),
  
  paste(
    "GO Cellular Component terms:",
    ifelse(
      is.null(go_cc),
      0,
      nrow(go_cc)
    )
  ),
  
  paste(
    "GO Molecular Function terms:",
    ifelse(
      is.null(go_mf),
      0,
      nrow(go_mf)
    )
  ),
  
  paste(
    "KEGG pathway terms:",
    ifelse(
      is.null(kegg),
      0,
      nrow(kegg)
    )
  ),
  
  "",
  
  "=========================================",
  
  "ANALYSIS PARAMETERS",
  
  "=========================================",
  
  "",
  
  "GO p value cutoff: 0.05",
  
  "GO q value cutoff: 0.20",
  
  "Multiple testing correction: Benjamini Hochberg",
  
  "KEGG p value cutoff: 0.05",
  
  "KEGG multiple testing correction: Benjamini Hochberg",
  
  "",
  
  "=========================================",
  
  "INTERPRETATION",
  
  "=========================================",
  
  "",
  
  "The analysis identified significantly enriched",
  
  "Gene Ontology categories and KEGG pathways",
  
  "among the AML biomarker gene set.",
  
  "",
  
  "Enrichment results provide functional context",
  
  "for the candidate biomarkers and support biological",
  
  "interpretation of the molecular processes associated",
  
  "with the identified AML biomarker candidates.",
  
  "",
  
  "Ensembl identifiers that could not be mapped",
  
  "to Entrez identifiers were excluded from the",
  
  "downstream enrichment analysis.",
  
  "",
  
  "=========================================",
  
  "END OF STAGE 10",
  
  "========================================="
  
)


writeLines(
  
  summary_lines,
  
  SUMMARY_FILE
  
)


##############################################################
# Final Summary
##############################################################

end_time <- Sys.time()

cat("\n")
cat("=========================================\n")
cat("Stage 10 Completed Successfully\n")
cat("=========================================\n\n")

cat(
  "Input Genes       :",
  input_gene_count,
  "\n"
)

cat(
  "Mapped Genes      :",
  mapped_gene_count,
  "\n"
)

cat(
  "Mapping Rate      :",
  round(mapping_rate, 2),
  "%\n"
)

cat(
  "Unmapped Genes    :",
  unmapped_gene_count,
  "\n\n"
)

cat(
  "GO BP Terms       :",
  ifelse(
    is.null(go_bp),
    0,
    nrow(go_bp)
  ),
  "\n"
)

cat(
  "GO CC Terms       :",
  ifelse(
    is.null(go_cc),
    0,
    nrow(go_cc)
  ),
  "\n"
)

cat(
  "GO MF Terms       :",
  ifelse(
    is.null(go_mf),
    0,
    nrow(go_mf)
  ),
  "\n"
)

cat(
  "KEGG Terms        :",
  ifelse(
    is.null(kegg),
    0,
    nrow(kegg)
  ),
  "\n\n"
)

cat(
  "Results Directory :",
  ENRICHMENT_DIR,
  "\n"
)

cat(
  "Plots Directory   :",
  PLOT_DIR,
  "\n"
)

cat(
  "Tables Directory  :",
  TABLE_DIR,
  "\n"
)

cat(
  "Summary Report    :",
  SUMMARY_FILE,
  "\n\n"
)

cat(
  "Time Elapsed      :",
  round(
    end_time - start_time,
    2
  ),
  "\n\n"
)

cat(
  "Stage 10 is now complete.\n\n"
)