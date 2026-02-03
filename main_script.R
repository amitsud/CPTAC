# ==========================================================================
# Name: CPTAC Multi-Omic Integration & Correlation Pipeline
# Description: A comprehensive R pipeline for integrating and analyzing 
#              Pancancer Proteogenomic data (CNV, RNA-seq, and Proteomics).
#              Includes data cleaning, covariate regression (Age, Sex, 
#              Purity, Cohort), KDE normalization, and gene-level 
#              multi-omic correlation analysis.
# Author: Amit Sud
# Date: February 2026
# ==========================================================================


# 1. LIBRARIES ------------------------------------------------------------

library(tidyverse)
library(vroom)
library(edgeR)
library(limma)
library(data.table)
library(ggExtra)
library(patchwork)
library(pheatmap)
library(broom)
library(stats)
library(matrixStats)
library(scales)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dndscv)
library(R.utils)
library(ggpubr)
library(boot)
library(mclust)
library(ppcor)
library(lme4)
library(sva)
library(cowplot)
library(RColorBrewer)
library(ggrepel)

# 2. REFERENCE DATA & SAMPLE DESCRIPTIONS ---------------------------------

# Path setup
base_path <- ""

# Gene references
ensembl_reference <- vroom(file.path(base_path, "reference/symbol_gene_mart_export.txt"))
cptac_reference <- vroom(file.path(base_path, "proteome/Proteome_BCM_GENCODE_v34_harmonized_v1/README/Gene_annotation_and_representable_isoform_mapping_table.txt"))

setDT(cptac_reference)

# Remove duplicates in cptac_reference by taking the first occurrence for each gene
cptac_reference_unique <- cptac_reference[!duplicated(gene), .(gene, gene_name_BCM_refined)]

# Load sample details and filter for Tumor only
sample_descriptions <- vroom(file.path(base_path, "rnaseq/RNA_Broad_v1/sample_descriptions.tsv"))
sample_descriptions <- sample_descriptions[sample_descriptions$tissue_type == "Tumor", ]
sample_descriptions <- sample_descriptions[sample_descriptions$project_id != "CPTAC3_polyA", ]

setDT(sample_descriptions)

# Standardize cohort names
sample_descriptions[, cohort := case_when(
  cohort == "Ovary" ~ "OV",
  cohort == "Colon" ~ "COAD",
  TRUE ~ cohort
)]

# 3. FILE SYSTEM INVENTORY ------------------------------------------------

# List all relevant files recursively
file_list <- list.files(
  path = base_path, 
  pattern = "\\.(txt|tsv|csv|gz)$", 
  full.names = TRUE, 
  recursive = TRUE
)

files_df <- data.frame(
  file_path = file_list,
  file_name = basename(file_list),
  stringsAsFactors = FALSE
)

# Extract metadata from paths
files_df <- files_df %>%
  mutate(
    modality = str_match(file_path, paste0(base_path, "/([^/]+)/"))[,2],
    cancer_type = str_extract(file_path, "(?i)BRCA|ccRCC|COAD|GBM|HGSC|HNSCC|LSCC|Lung_Squamous_Cell_Carcinoma|LUAD|PDAC|OV|UCEC|PanCan|ALL"),
    center = str_extract(file_path, "(?i)BCM|Broad|WashU|UMich")
  )

# 4. READ CLINICAL DATA ---------------------------------------------------

cancer_types <- c("BRCA", "ccRCC", "COAD", "GBM", "HGSC", "HNSCC", "LSCC", "LUAD", "PDAC", "UCEC")
input_dir_meta <- file.path(base_path, "clinical/Clinical_meta_data_v1/")

meta_data_list <- list()

for (cancer in cancer_types) {
  file_name_meta <- paste0(input_dir_meta, cancer, "_meta.txt")
  
  if (!file.exists(file_name_meta)) {
    cat(paste("File not found for cancer type:", cancer, "\n"))
    next
  }
  
  meta_data <- fread(file_name_meta)
  meta_data <- meta_data[-1, ] # Remove second row (often headers/units)
  meta_data[, cohort := cancer]
  meta_data_list[[cancer]] <- meta_data
  
  cat(paste("Processed meta data for:", cancer, "\n"))
}

final_merged_data_meta <- rbindlist(meta_data_list, use.names = TRUE, fill = TRUE)
final_merged_data_meta <- final_merged_data_meta %>% filter(Tumor == "Yes")

# Organize Columns
setnames(final_merged_data_meta, old = "idx", new = "sample_id")
setcolorder(final_merged_data_meta, c("cohort", setdiff(names(final_merged_data_meta), "cohort")))

# Recode Stage to Numeric
final_merged_data_meta[, Stage := fcase(
  Stage == "Stage I", 1,
  Stage == "Stage II", 2,
  Stage == "Stage III", 3,
  Stage == "Stage IV", 4,
  default = NA_real_
)]

# 5. READ CNV DATA --------------------------------------------------------

input_dir_cnv <- file.path(base_path, "wes/CNV_BCM_v1/")
final_merged_data_cnv <- NULL

for (cancer in cancer_types) {
  file_name_cnv <- paste0(input_dir_cnv, cancer, "_WES_CNV_gene_gistic_level.txt")
  
  if (!file.exists(file_name_cnv)) {
    cat(paste("File not found for cancer type:", cancer, "\n"))
    next
  }
  
  wes_data_cnv <- fread(file_name_cnv)
  wes_data_cnv <- merge(wes_data_cnv, cptac_reference_unique, by.x = "idx", by.y = "gene", all.x = TRUE)
  
  cancer_columns_cnv <- setdiff(names(wes_data_cnv), c("idx", "gene_name_BCM_refined"))
  wes_data_cnv <- wes_data_cnv[, c("idx", "gene_name_BCM_refined", cancer_columns_cnv), with = FALSE]
  
  if (is.null(final_merged_data_cnv)) {
    final_merged_data_cnv <- wes_data_cnv
  } else {
    final_merged_data_cnv <- merge(final_merged_data_cnv, wes_data_cnv, by = c("idx", "gene_name_BCM_refined"), all = TRUE)
  }
  cat(paste("Merged cnv wes data for:", cancer, "\n"))
}

# Check matches
cnv_columns <- setdiff(names(final_merged_data_cnv), c("idx", "gene_name_BCM_refined"))
final_merged_data_meta[, matched_cnv := sample_id %in% cnv_columns]

# 6. READ TRANSCRIPTOME DATA ----------------------------------------------

input_dir_rnaseq <- file.path(base_path, "rnaseq/RNA_BCM_v1/")
final_merged_data_rnaseq <- NULL

for (cancer in cancer_types) {
  file_name_rnaseq <- paste0(input_dir_rnaseq, cancer, "_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt")
  
  if (!file.exists(file_name_rnaseq)) {
    cat(paste("File not found for cancer type:", cancer, "\n"))
    next
  }
  
  rnaseq_data <- fread(file_name_rnaseq)
  rnaseq_data <- merge(rnaseq_data, cptac_reference_unique, by.x = "idx", by.y = "gene", all.x = TRUE)
  
  cancer_columns_rnaseq <- setdiff(names(rnaseq_data), c("idx", "gene_name_BCM_refined"))
  rnaseq_data <- rnaseq_data[, c("idx", "gene_name_BCM_refined", cancer_columns_rnaseq), with = FALSE]
  
  if (is.null(final_merged_data_rnaseq)) {
    final_merged_data_rnaseq <- rnaseq_data
  } else {
    final_merged_data_rnaseq <- merge(final_merged_data_rnaseq, rnaseq_data, by = c("idx", "gene_name_BCM_refined"), all = TRUE)
  }
  cat(paste("Merged rnaseq data for:", cancer, "\n"))
}

rna_columns <- setdiff(names(final_merged_data_rnaseq), c("idx", "gene_name_BCM_refined"))
final_merged_data_meta[, matched_rnaseq := sample_id %in% rna_columns]

# 7. READ PROTEOME DATA ---------------------------------------------------

cancer_types_prot <- c("BRCA", "ccRCC", "COAD", "GBM", "OV", "HNSCC", "LSCC", "LUAD", "PDAC", "UCEC")
input_dir_proteome <- file.path(base_path, "proteome/Proteome_BCM_GENCODE_v34_harmonized_v1/")
file_suffix_proteome <- "_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"

final_merged_data_proteome <- NULL

for (cancer in cancer_types_prot) {
  file_name_proteome <- paste0(input_dir_proteome, cancer, file_suffix_proteome)
  
  if (!file.exists(file_name_proteome)) {
    cat(paste("File not found for cancer type:", cancer, "\n"))
    next
  }
  
  proteome_data <- fread(file_name_proteome)
  proteome_data <- merge(proteome_data, cptac_reference_unique[, .(gene, gene_name_BCM_refined)], by.x = "idx", by.y = "gene", all.x = TRUE)
  
  cancer_columns_proteome <- setdiff(names(proteome_data), c("idx", "gene_name_BCM_refined"))
  proteome_data <- proteome_data[, c("idx", "gene_name_BCM_refined", cancer_columns_proteome), with = FALSE]
  
  if (is.null(final_merged_data_proteome)) {
    final_merged_data_proteome <- proteome_data
  } else {
    final_merged_data_proteome <- merge(final_merged_data_proteome, proteome_data, by = c("idx", "gene_name_BCM_refined"), all = TRUE)
  }
  cat(paste("Merged data for:", cancer, "\n"))
}

proteome_columns <- setdiff(names(final_merged_data_proteome), c("idx", "gene_name_BCM_refined"))
final_merged_data_meta[, matched_proteome := sample_id %in% proteome_columns]

# 8. GENERATE MATRICES ----------------------------------------------------

# Identify intersection of genes and samples
common_genes <- Reduce(intersect, list(
  final_merged_data_cnv$gene_name_BCM_refined,
  final_merged_data_rnaseq$gene_name_BCM_refined,
  final_merged_data_proteome$gene_name_BCM_refined
))

cnv_samples <- setdiff(names(final_merged_data_cnv), c("idx", "gene_name_BCM_refined"))
rna_samples <- setdiff(names(final_merged_data_rnaseq), c("idx", "gene_name_BCM_refined"))
proteome_samples <- setdiff(names(final_merged_data_proteome), c("idx", "gene_name_BCM_refined"))
meta_samples <- final_merged_data_meta$sample_id

common_samples <- Reduce(intersect, list(cnv_samples, rna_samples, proteome_samples, meta_samples))

# Prune and align datasets
prune_data <- function(dt, genes, samples) {
  dt_pruned <- dt[gene_name_BCM_refined %in% genes, ]
  setkey(dt_pruned, gene_name_BCM_refined)
  dt_pruned <- dt_pruned[, c("idx", "gene_name_BCM_refined", samples), with = FALSE]
  return(dt_pruned)
}

cnv_pruned <- prune_data(final_merged_data_cnv, common_genes, common_samples)
rna_pruned <- prune_data(final_merged_data_rnaseq, common_genes, common_samples)
proteome_pruned <- prune_data(final_merged_data_proteome, common_genes, common_samples)

# Create Matrices
as_aligned_matrix <- function(dt) {
  mat <- as.matrix(dt[, -c("idx", "gene_name_BCM_refined"), with = FALSE])
  rownames(mat) <- dt$gene_name_BCM_refined
  return(mat)
}

cnv_matrix <- as_aligned_matrix(cnv_pruned)
rna_matrix <- as_aligned_matrix(rna_pruned)
proteome_matrix <- as_aligned_matrix(proteome_pruned)

# Final Diagnostics
cat("\nFinal Matrix Summary:\n")
cat("Genes:", nrow(cnv_matrix), "| Samples:", ncol(cnv_matrix), "\n")








# QUALITY CONTROL & FILTERING PIPELINE ------------------------------------

# 1. INITIAL PRUNING: REMOVE ALL-ZERO ROWS/COLS ---------------------------

# Identify genes (rows) with all zeros across any matrix
genes_all_zero <- Reduce(union, list(
  rownames(rna_matrix)[apply(rna_matrix, 1, function(x) all(x == 0))],
  rownames(proteome_matrix)[apply(proteome_matrix, 1, function(x) all(x == 0))],
  rownames(cnv_matrix)[apply(cnv_matrix, 1, function(x) all(x == 0))]
))

cat("Number of genes with all zeros in any matrix:", length(genes_all_zero), "\n")

# Identify samples (columns) with all zeros across any matrix
samples_all_zero <- Reduce(union, list(
  colnames(rna_matrix)[apply(rna_matrix, 2, function(x) all(x == 0))],
  colnames(proteome_matrix)[apply(proteome_matrix, 2, function(x) all(x == 0))],
  colnames(cnv_matrix)[apply(cnv_matrix, 2, function(x) all(x == 0))]
))

cat("Number of samples with all zeros in any matrix:", length(samples_all_zero), "\n")

# Apply initial pruning
rna_matrix_pruned      <- rna_matrix[!rownames(rna_matrix) %in% genes_all_zero, !colnames(rna_matrix) %in% samples_all_zero]
proteome_matrix_pruned <- proteome_matrix[!rownames(proteome_matrix) %in% genes_all_zero, !colnames(proteome_matrix) %in% samples_all_zero]
cnv_matrix_pruned      <- cnv_matrix[!rownames(cnv_matrix) %in% genes_all_zero, !colnames(cnv_matrix) %in% samples_all_zero]

cat("Dimensions after initial pruning:\n")
cat("RNA:", dim(rna_matrix_pruned), "| Proteome:", dim(proteome_matrix_pruned), "| CNV:", dim(cnv_matrix_pruned), "\n")

# 2. DATA CONSISTENCY: JACCARD INDEX --------------------------------------

# Convert to binary presence/absence to check data overlap per sample
proteome_binary <- ifelse(!is.na(proteome_matrix_pruned), 1, 0)
rna_binary      <- ifelse(!is.na(rna_matrix_pruned), 1, 0)

if (!all(rownames(proteome_binary) == rownames(rna_binary)) || !all(colnames(proteome_binary) == colnames(rna_binary))) {
  stop("Matrices are not aligned. Please align before calculating Jaccard index.")
}

# Calculate Jaccard index for each sample
jaccard_indices <- sapply(1:ncol(proteome_binary), function(i) {
  intersect_count <- sum(proteome_binary[, i] & rna_binary[, i])
  union_count     <- sum(proteome_binary[, i] | rna_binary[, i])
  if (union_count == 0) return(NA)
  return(intersect_count / union_count)
})

jaccard_df <- data.frame(
  sample_id = colnames(proteome_matrix_pruned),
  jaccard_index = jaccard_indices
)

# Merge with metadata for cohort-level visualization
jaccard_df <- merge(jaccard_df, final_merged_data_meta[, .(sample_id, cohort)], by = "sample_id", all.x = TRUE)

# Plot Jaccard Index by Cohort
cohort_order <- jaccard_df %>%
  group_by(cohort) %>%
  summarize(median_jaccard = median(jaccard_index, na.rm = TRUE)) %>%
  arrange(median_jaccard) %>%
  pull(cohort)

jaccard_df$cohort <- factor(jaccard_df$cohort, levels = cohort_order)

ggplot(jaccard_df, aes(x = cohort, y = jaccard_index)) +
  geom_boxplot(fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(color = "blue", alpha = 0.5, size = 1.5, width = 0.2) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    panel.grid = element_blank()
  ) +
  labs(x = "Cohort", y = "Jaccard Index", title = "Jaccard Index by Cohort (Proteome vs Transcriptome)") +
  coord_cartesian(ylim = c(0, 1))

# 3. HELPER FUNCTIONS FOR ADVANCED QC -------------------------------------

# Identify pure samples based on WES threshold
get_pure_samples <- function(meta_features, purity_threshold = 0.2, matrix_colnames) {
  pure_samples <- meta_features %>%
    filter(!is.na(WES_purity) & WES_purity >= purity_threshold) %>%
    pull(sample_id)
  
  aligned_samples <- intersect(pure_samples, matrix_colnames)
  cat("Matched pure samples (purity >= ", purity_threshold, "):", length(aligned_samples), "\n")
  return(aligned_samples)
}

# Filter genes by detection rate
filter_genes_detected <- function(data_matrix, detection_threshold = 0.5) {
  detection_rate <- rowMeans(!is.na(data_matrix) & data_matrix != 0, na.rm = TRUE)
  genes_to_keep  <- detection_rate >= detection_threshold
  cat("Genes removed (detected in <", detection_threshold * 100, "% of samples):", sum(!genes_to_keep), "\n")
  return(data_matrix[genes_to_keep, , drop = FALSE])
}

# Filter genes by percentile of expression
filter_low_percentile <- function(data_matrix, percentile = 0.005) {
  medians <- apply(data_matrix, 1, median, na.rm = TRUE)
  threshold <- quantile(medians, percentile, na.rm = TRUE)
  genes_to_keep <- medians > threshold
  cat("Threshold (", percentile * 100, " percentile):", threshold, "\n")
  cat("Items below threshold removed:", sum(!genes_to_keep), "\n")
  return(data_matrix[genes_to_keep, , drop = FALSE])
}

# Remove rows/columns with zero variance
remove_zero_variance <- function(data_matrix) {
  zero_variance_columns <- apply(data_matrix, 2, function(x) var(x, na.rm = TRUE)) == 0
  zero_variance_rows    <- apply(data_matrix, 1, function(x) var(x, na.rm = TRUE)) == 0
  
  if (any(zero_variance_columns)) {
    cat("Zero-variance columns removed:", sum(zero_variance_columns), "\n")
    data_matrix <- data_matrix[, !zero_variance_columns, drop = FALSE]
  }
  if (any(zero_variance_rows)) {
    cat("Zero-variance rows removed:", sum(zero_variance_rows), "\n")
    data_matrix <- data_matrix[!zero_variance_rows, , drop = FALSE]
  }
  return(data_matrix)
}

# Matrix alignment function
align_matrices <- function(rna_matrix, proteome_matrix, cnv_matrix, pure_samples) {
  common_genes <- intersect(intersect(rownames(rna_matrix), rownames(proteome_matrix)), rownames(cnv_matrix))
  
  rna_matrix      <- rna_matrix[common_genes, pure_samples, drop = FALSE]
  proteome_matrix <- proteome_matrix[common_genes, pure_samples, drop = FALSE]
  cnv_matrix      <- cnv_matrix[common_genes, pure_samples, drop = FALSE]
  
  cat("Final alignment: ", length(common_genes), "common genes across", length(pure_samples), "pure samples.\n")
  return(list(rna_matrix = rna_matrix, proteome_matrix = proteome_matrix, cnv_matrix = cnv_matrix))
}

# 4. MAIN QC EXECUTION ----------------------------------------------------

process_data <- function(rna_matrix, proteome_matrix, cnv_matrix, meta_features) {
  # 1. Sample Purity Filter
  pure_samples <- get_pure_samples(meta_features, purity_threshold = 0.2, matrix_colnames = colnames(rna_matrix))
  
  rna_matrix      <- rna_matrix[, pure_samples, drop = FALSE]
  proteome_matrix <- proteome_matrix[, pure_samples, drop = FALSE]
  cnv_matrix      <- cnv_matrix[, pure_samples, drop = FALSE]
  
  # 2. RNA QC
  cat("\nProcessing RNA matrix...\n")
  rna_matrix <- rna_matrix %>%
    filter_low_percentile(percentile = 0.01) %>%
    filter_genes_detected(detection_threshold = 0.5) %>%
    remove_zero_variance()
  
  # 3. Proteome QC
  cat("\nProcessing Proteome matrix...\n")
  proteome_matrix <- proteome_matrix %>%
    filter_low_percentile(percentile = 0) %>%
    filter_genes_detected(detection_threshold = 0.5) %>%
    remove_zero_variance()
  
  # 4. Align and Impute
  aligned_matrices <- align_matrices(rna_matrix, proteome_matrix, cnv_matrix, pure_samples)
  
  # Handling missing proteomic values (Impute with 0)
  aligned_matrices$proteome_matrix[is.na(aligned_matrices$proteome_matrix)] <- 0
  
  return(aligned_matrices)
}

# Execute final pipeline
final_matrices <- process_data(rna_matrix_pruned, proteome_matrix_pruned, cnv_matrix_pruned, final_merged_data_meta)

# Extract and Verify
rna_matrix_final      <- final_matrices$rna_matrix
proteome_matrix_final <- final_matrices$proteome_matrix
cnv_matrix_final      <- final_matrices$cnv_matrix

cat("\nFinal Output Dimensions:\n")
cat("RNA:", dim(rna_matrix_final), "| Proteome:", dim(proteome_matrix_final), "| CNV:", dim(cnv_matrix_final), "\n")





# PCA & COVARIATE REGRESSION ----------------------------------------------

# 1. PCA HELPER FUNCTIONS -------------------------------------------------

perform_pca <- function(data_matrix) {
  prcomp(t(data_matrix), scale. = TRUE)
}

extract_explained_variance <- function(pca) {
  pca$sdev^2 / sum(pca$sdev^2)
}

plot_explained_variance <- function(explained_variance, title) {
  explained_df <- data.frame(
    PC = factor(paste0("PC", seq_along(explained_variance)), 
                levels = paste0("PC", seq_along(explained_variance))),
    Variance = explained_variance
  )
  ggplot(explained_df[1:10, ], aes(x = PC, y = Variance)) +
    geom_bar(stat = "identity", fill = "steelblue", color = "black") +
    theme_minimal() +
    theme(axis.line = element_line(color = "black"),
          axis.text = element_text(color = "black")) +
    labs(x = "Principal Components", y = "Explained Variance Ratio", title = title)
}

plot_pca <- function(pca, meta_data, title) {
  pca_df <- data.frame(pca$x[, 1:2], sample_id = rownames(pca$x)) %>%
    left_join(meta_data, by = "sample_id")
  
  ggplot(pca_df, aes(x = PC1, y = PC2, color = cohort)) +
    geom_point(size = 2, alpha = 0.7) +
    theme_minimal() +
    theme(axis.line = element_line(color = "black"),
          axis.text = element_text(color = "black")) +
    labs(x = "PC1", y = "PC2", title = title) +
    scale_color_discrete(name = "Cohort")
}

# 2. EXPLORATORY PCA (PRE-ADJUSTMENT) -------------------------------------

pca_transcriptome <- perform_pca(rna_matrix_final)
pca_proteome      <- perform_pca(proteome_matrix_final)

plot_explained_variance(extract_explained_variance(pca_transcriptome), "Transcriptomic PCA: Explained Variance")
plot_explained_variance(extract_explained_variance(pca_proteome), "Proteomic PCA: Explained Variance")

plot_pca(pca_transcriptome, final_merged_data_meta, "Transcriptomic PCA: First 2 PCs")
plot_pca(pca_proteome, final_merged_data_meta, "Proteomic PCA: First 2 PCs")

# 3. COVARIATE REGRESSION -------------------------------------------------

regress_out_covariates <- function(expression_matrix, meta_features) {
  adjusted_matrix <- matrix(nrow = nrow(expression_matrix), ncol = ncol(expression_matrix))
  rownames(adjusted_matrix) <- rownames(expression_matrix)
  colnames(adjusted_matrix) <- colnames(expression_matrix)
  
  for (gene in rownames(expression_matrix)) {
    model_data <- meta_features %>%
      mutate(expression = expression_matrix[gene, match(sample_id, colnames(expression_matrix))])
    
    # Fit linear model for Age, Sex, Cohort, and Purity
    model <- lm(expression ~ Age + Sex + cohort + WES_purity, data = model_data)
    
    # Extract residuals
    adj_vals <- resid(model)
    adjusted_matrix[gene, ] <- adj_vals[match(colnames(expression_matrix), model_data$sample_id)]
  }
  return(adjusted_matrix)
}

# Prepare Metadata for Regression
aligned_meta_adj <- final_merged_data_meta %>%
  filter(sample_id %in% colnames(rna_matrix_final)) %>%
  drop_na(Age, Sex, cohort, WES_purity)

# Execute Adjustment
cat("Regressing covariates from Transcriptome and Proteome...\n")
rna_matrix_adj      <- regress_out_covariates(rna_matrix_final[, aligned_meta_adj$sample_id], aligned_meta_adj)
proteome_matrix_adj <- regress_out_covariates(proteome_matrix_final[, aligned_meta_adj$sample_id], aligned_meta_adj)

# Align common genes/samples across all 3 modalities
common_samples_adj <- Reduce(intersect, list(colnames(rna_matrix_adj), colnames(proteome_matrix_adj), colnames(cnv_matrix_final)))
common_genes_adj   <- Reduce(intersect, list(rownames(rna_matrix_adj), rownames(proteome_matrix_adj), rownames(cnv_matrix_final)))

rna_matrix_final_adj      <- rna_matrix_adj[common_genes_adj, common_samples_adj]
proteome_matrix_final_adj <- proteome_matrix_adj[common_genes_adj, common_samples_adj]
cnv_matrix_final_adj      <- cnv_matrix_final[common_genes_adj, common_samples_adj]

# 4. KERNEL DENSITY ESTIMATION (KDE) NORMALIZATION -----------------------

normalize_center_kde <- function(data_matrix) {
  t(apply(data_matrix, 1, function(values) {
    density_est <- density(values, na.rm = TRUE)
    kde_mean    <- weighted.mean(density_est$x, density_est$y)
    kde_sd      <- sqrt(weighted.mean((density_est$x - kde_mean)^2, density_est$y))
    (values - kde_mean) / kde_sd
  }))
}

cat("Normalizing adjusted matrices...\n")
rna_matrix_final_norm      <- normalize_center_kde(rna_matrix_final_adj)
proteome_matrix_final_norm <- normalize_center_kde(proteome_matrix_final_adj)

# 5. VALIDATION: POST-ADJUSTMENT CORRELATION HEATMAP ----------------------

# Re-run PCA on adjusted data
pca_rna_norm  <- perform_pca(rna_matrix_final_norm)
pca_prot_norm <- perform_pca(proteome_matrix_final_norm)

# Prepare metadata for correlation check
meta_features_val <- final_merged_data_meta %>%
  filter(sample_id %in% colnames(rna_matrix_final_norm)) %>%
  transmute(
    sample_id,
    tumor  = as.numeric(factor(cohort)),
    age    = as.numeric(Age),
    sex    = as.numeric(factor(Sex)),
    purity = as.numeric(WES_purity),
    TMB    = as.numeric(TMB)
  ) %>%
  filter(!if_any(everything(), is.na))

# 1. Align BOTH objects to the common IDs
common_ids <- intersect(rownames(pca_rna_norm$x), meta_features_val$sample_id)

# Subset Meta
aligned_meta_val <- meta_features_val[match(common_ids, meta_features_val$sample_id), ]

# Subset PCA matrix
aligned_pca_x <- pca_rna_norm$x[common_ids, 1:10]

# 2. Update the helper function
calculate_correlation <- function(pcs_subset, meta_df) {
  # We use match to ensure the metadata rows are in the EXACT same order as PCA rows
  meta_ordered <- meta_df[match(rownames(pcs_subset), meta_df$sample_id), ]
  
  sapply(names(meta_ordered)[-1], function(feat) {
    # Now dimensions are guaranteed to match
    apply(pcs_subset, 2, function(pc) {
      cor(pc, as.numeric(meta_ordered[[feat]]), use = "complete.obs")
    })
  })
}

# 3. Run the correlation with the aligned PCA subset
cor_rna  <- calculate_correlation(aligned_pca_x, aligned_meta_val)
cor_prot <- calculate_correlation(pca_prot_norm$x[common_ids, 1:10], aligned_meta_val)



# Combine and Plot Heatmap
cor_heatmap_df <- bind_rows(
  as.data.frame(cor_rna)  %>% mutate(Dataset = "Transcriptome", PC = paste0("PC", 1:10)),
  as.data.frame(cor_prot) %>% mutate(Dataset = "Proteome", PC = paste0("PC", 1:10))
) %>%
  pivot_longer(cols = -c(Dataset, PC), names_to = "Feature", values_to = "Correlation") %>%
  mutate(PC = factor(PC, levels = rev(paste0("PC", 1:10))))



ggplot(cor_heatmap_df, aes(x = Feature, y = PC, fill = Correlation)) +
  geom_tile(color = "white") +
  facet_wrap(~Dataset) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Post-Adjustment Correlation: PCs vs Meta Features")





# GENE-LEVEL CORRELATIONS & MULTI-OMIC ANALYSIS ----------------------------

# 1. CALCULATE GLOBAL CORRELATIONS ----------------------------------------

# Function to calculate Pearson correlation for each gene across modalities
calculate_gene_correlations <- function(matrix1, matrix2) {
  # Ensure genes are aligned
  common_genes <- intersect(rownames(matrix1), rownames(matrix2))
  
  sapply(common_genes, function(gene) {
    cor(matrix1[gene, ], matrix2[gene, ], use = "complete.obs", method = "pearson")
  })
}

cat("Calculating multi-omic correlations...\n")

# Ensure the CNV matrix name matches correlation calls
# We use the 'adj' version created during the covariate regression step
cnv_matrix_final_adjusted_normalized <- cnv_matrix_final_adj

# Now run the correlations again
cat("Calculating multi-omic correlations...\n")

# Correlation: CNV vs Transcriptome
cnv_rna_corr  <- calculate_gene_correlations(cnv_matrix_final_adjusted_normalized, 
                                             rna_matrix_final_norm)

# Correlation: CNV vs Proteome
cnv_proteome_corr <- calculate_gene_correlations(cnv_matrix_final_adjusted_normalized, 
                                                 proteome_matrix_final_norm)

# Correlation: Transcriptome vs Proteome
transcriptome_proteome_corr <- calculate_gene_correlations(rna_matrix_final_norm, 
                                                           proteome_matrix_final_norm)


# Combine results into a data frame for plotting
correlation_summary <- data.frame(
  gene     = names(transcriptome_proteome_corr),
  CNV_RNA  = cnv_rna_corr,
  CNV_Prot = cnv_proteome_corr,
  RNA_Prot = transcriptome_proteome_corr
)



# 2. TARGETED GENE PAIR ANALYSIS BY COHORT --------------------------------

analyze_gene_pair_by_cancer <- function(expression_matrix, dataset_name, gene1 = "SEC62", gene2 = "HLA-E") {
  
  # Check genes exist
  if (!(gene1 %in% rownames(expression_matrix)) | !(gene2 %in% rownames(expression_matrix))) {
    message(paste("❌ Genes", gene1, "or", gene2, "not found in", dataset_name))
    return(NULL)
  }
  
  # Prepare plotting data
  df <- data.frame(
    sample_id = colnames(expression_matrix),
    Gene1 = expression_matrix[gene1, ],
    Gene2 = expression_matrix[gene2, ]
  ) %>%
    left_join(final_merged_data_meta[, .(sample_id, cohort)], by = "sample_id") %>%
    filter(!is.na(cohort) & is.finite(Gene1) & is.finite(Gene2))
  
  if (nrow(df) < 3) {
    message("❌ Too few samples for correlation.")
    return(NULL)
  }
  
  # Statistics
  pearson_test  <- cor.test(df$Gene1, df$Gene2, method = "pearson")
  spearman_test <- cor.test(df$Gene1, df$Gene2, method = "spearman")
  
  # Refine cohort labels for legend
  cohort_counts <- df %>%
    group_by(cohort) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(label = paste0(cohort, " (n=", n, ")"))
  
  df <- df %>% left_join(cohort_counts, by = "cohort")
  
  # Color Palette setup
  n_groups      <- nrow(cohort_counts)
  base_colors   <- brewer.pal(min(max(n_groups, 3), 8), "Set1")
  color_palette <- setNames(colorRampPalette(base_colors)(n_groups), cohort_counts$label)
  
  # Visualization
  
  p <- ggplot(df, aes(x = Gene1, y = Gene2)) +
    geom_point(aes(fill = label), color = "black", shape = 21, size = 4, stroke = 0.3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.8) +
    scale_fill_manual(values = color_palette) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.ticks   = element_line(color = "black"),
      panel.grid   = element_blank(),
      legend.position = "bottom"
    ) +
    labs(
      x = paste(gene1, "(Normalized Expression)"),
      y = paste(gene2, "(Normalized Expression)"),
      title = paste(gene1, "vs", gene2, "Correlation in", dataset_name),
      subtitle = paste0(
        "Pearson r = ", round(pearson_test$estimate, 3),
        " (p = ", format.pval(pearson_test$p.value, eps = 1e-3), ") | ",
        "Spearman ρ = ", round(spearman_test$estimate, 3), 
        " | n = ", nrow(df)
      ),
      fill = "Cancer Type"
    )
  
  print(p)
  return(list(pearson = pearson_test, spearman = spearman_test))
}

# 🔬 RUN FINAL ANALYSES
rna_analysis  <- analyze_gene_pair_by_cancer(rna_matrix_final_norm, "Transcriptome")
prot_analysis <- analyze_gene_pair_by_cancer(proteome_matrix_final_norm, "Proteome")
