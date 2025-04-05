#!/usr/bin/env Rscript

#install.packages("doRNG")
#install.packages("doParallel")

source("dynGENIE3.R")

# Load raw count matrix and sample metadata
raw_counts_aka <- read.delim(
  gzfile("../RawCountsStar-Aka.tsv.gz"),
  row.names = 1
)

contrastAka <- read.csv(
  "../contrastAka.csv",
  row.names = 1
)

# load regulators
regulators <- read.csv(
  "../DEseq_AkaPval1.in",
  header = F
)


# Step 1: Keep only valid samples (present in both metadata and count matrix)
valid_samples <- row.names(contrastAka)[row.names(contrastAka) %in% colnames(raw_counts_aka)]
contrastAka <- contrastAka[valid_samples, ]
raw_counts_aka_mat <- raw_counts_aka[, valid_samples]

# Step 2: Filter genes that have zero counts in >90% of replicates per timepoint
contrastAka_by_timepoint <- split(row.names(contrastAka), contrastAka$time)
filtered_genes_by_tp <- list()

for (tp in names(contrastAka_by_timepoint)) {
  samples <- contrastAka_by_timepoint[[tp]]
  submat <- raw_counts_aka_mat[, samples, drop = FALSE]
  zero_counts <- rowSums(submat == 0)
  threshold <- 0.9 * ncol(submat)
  keep_genes <- rownames(submat)[zero_counts <= threshold]
  message(tp, ": ", length(keep_genes), " genes kept")
  filtered_genes_by_tp[[tp]] <- keep_genes
}

# Step 3: Keep genes that passed in all timepoints
genes_to_keep <- Reduce(intersect, filtered_genes_by_tp)
filtered_counts <- raw_counts_aka_mat[genes_to_keep, ]

# Step 4: Compute average gene expression per timepoint
samples_by_tp <- split(row.names(contrastAka), contrastAka$time)

avg_expr_per_tp <- lapply(samples_by_tp, function(samples) {
  submat <- filtered_counts[, samples, drop = FALSE]
  rowMeans(submat)
})

# Step 5: Combine into matrix (timepoints x genes), sorted by timepoint name
TP_order <- sort(names(avg_expr_per_tp))
TS.matrix <- do.call(rbind, avg_expr_per_tp[TP_order])

# Step 6: Prepare final TS.data list for dynGENIE3
TS.data <- list(TS.matrix)

# Optional checks
cat("TS.data structure:\n")
str(TS.data)

#dyn.load("dynGENIE3_R_C_wrapper/dynGENIE3.so")
dynGENIE3_results <- dynGENIE3(TS.data = TS.data, regulators = regulators$V1)

