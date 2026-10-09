# ############################################################
# 1. LOAD R PACKAGES
# Load the packages used to read and analyse the count data.
# ############################################################
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scDblFinder)
  library(DropletUtils)
  library(Matrix)
})

# ############################################################
# 2. PATHS AND SAMPLE CONFIGURATION
# Set the sample-specific input and output paths.
# ############################################################
sample_dir <- getwd()
sample <- basename(sample_dir)
input_h5 <- file.path(
  "/home/cv1u24/GIST/1_Raw_Processing/CellBender",
  sample,
  "CellBender_Outputs",
  paste0(sample, "_cellbender_filtered.h5")
)
output_dir <- file.path(sample_dir, "Output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(input_h5)) {
  stop(paste("Input file not found:", input_h5))
}
cat("Processing sample:", sample, "\n")
cat("Input:", input_h5, "\n")
cat("Output:", output_dir, "\n")

# ############################################################
# 3. READ THE CELLBENDER-FILTERED MATRIX
# Import the filtered counts as a SingleCellExperiment object.
# ############################################################
sce <- read10xCounts(input_h5, col.names = TRUE)

# ############################################################
# 4. RUN SCDblFinder
# Set the random seed and detect doublets with scDblFinder.
# ############################################################
set.seed(1)
sce <- scDblFinder(sce)

# ############################################################
# 5. SAVE PER-BARCODE CALLS
# Write the scDblFinder score and class for each barcode to CSV.
# ############################################################
barcode_df <- data.frame(
  sample = sample,
  barcode = colnames(sce),
  scdblfinder_score = sce$scDblFinder.score,
  scdblfinder_class = sce$scDblFinder.class,
  scdblfinder_doublet = sce$scDblFinder.class == "doublet"
)
barcode_output <- file.path(
  output_dir,
  paste0(sample, "_scdblfinder_barcode_calls.csv")
)
write.csv(barcode_df, barcode_output, row.names = FALSE)

# ############################################################
# 6. SAVE SUMMARY STATISTICS
# Summarise the number and scores of detected doublets.
# ############################################################
doublet_calls <- barcode_df$scdblfinder_doublet
scores <- barcode_df$scdblfinder_score
summary_df <- data.frame(
  sample = sample,
  n_cells = ncol(sce),
  n_genes = nrow(sce),
  scdblfinder_n_doublets = sum(doublet_calls),
  scdblfinder_pct_doublets = mean(doublet_calls) * 100,
  scdblfinder_mean_score = mean(scores, na.rm = TRUE),
  scdblfinder_median_score = median(scores, na.rm = TRUE),
  scdblfinder_max_score = max(scores, na.rm = TRUE)
)
summary_output <- file.path(
  output_dir,
  paste0(sample, "_scdblfinder_summary.csv")
)
write.csv(summary_df, summary_output, row.names = FALSE)

# ############################################################
# 7. SAVE THE R OBJECT AND REPORT OUTPUTS
# Save the annotated object and print the generated file paths.
# ############################################################
rds_output <- file.path(
  output_dir,
  paste0(sample, "_scdblfinder_calls.rds")
)
saveRDS(sce, rds_output)
cat("Saved files:\n")
cat(barcode_output, "\n")
cat(summary_output, "\n")
cat(rds_output, "\n")
cat("Done.\n")