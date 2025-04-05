#!/usr/bin/env Rscript

# Load data from df and prepare for output
# dynGENIE3
# Polished by Kristina Gagalova

library("doRNG")
library("doParallel")

source("dynGENIE3.R")

# Read the raw counts table
data <- read.table(gzfile("../RawCountsStar-Aka.tsv.gz"), sep="\t", header=TRUE, row.names=1)

# Read the contrast file that contains the time and treatment columns per sample
contrastAka <- read.csv("../contrastAka.csv", row.names=1)

# Ensure matching samples between raw data and contrast
valid_samples <- intersect(colnames(data), rownames(contrastAka))
contrastAka <- contrastAka[valid_samples, ]
data <- data[, valid_samples]

# Map original time numbers to new labels
original_times <- sort(unique(contrastAka$time))
time_mapping <- setNames(c(0, 1, 3, 7, 9), original_times)
contrastAka$mapped_time <- time_mapping[as.character(contrastAka$time)]

# Filter genes (>90% zero counts per timepoint)
contrast_by_tp_treatment <- split(rownames(contrastAka), list(contrastAka$mapped_time, contrastAka$treatment))
filtered_genes_by_group <- lapply(contrast_by_tp_treatment, function(samples) {
  submat <- data[, samples, drop=FALSE]
  keep_genes <- rownames(submat)[rowSums(submat == 0) <= 0.9 * ncol(submat)]
  keep_genes
})

# Intersection of genes kept in all groups
genes_to_keep <- Reduce(intersect, filtered_genes_by_group)
filtered_counts <- data[genes_to_keep, ]

# Separate data by treatment, average expression per timepoint, and write to files
treatments <- unique(contrastAka$treatment)

for (treatment in treatments) {
  contrast_subset <- contrastAka[contrastAka$treatment == treatment, ]
  timepoints <- sort(unique(contrast_subset$mapped_time))
  avg_expr_matrix <- sapply(timepoints, function(tp) {
    samples <- rownames(contrast_subset[contrast_subset$mapped_time == tp, ])
    rowMeans(filtered_counts[, samples, drop=FALSE])
  })
  avg_expr_df <- data.frame(time_points = timepoints, t(avg_expr_matrix))
  colnames(avg_expr_df) <- c("time_points", genes_to_keep)
  write.table(avg_expr_df, file=paste0("averaged_expression_", treatment, ".txt"), sep="\t", quote=FALSE, row.names=FALSE)
}

# load data with function from package
TS1 <- read.expr.matrix("./averaged_expression_Pet.txt", form="rows.are.samples")
TS2 <- read.expr.matrix("./averaged_expression_Tox8F.txt", form="rows.are.samples")

regulators <- read.csv(
     "../DEseq_AkaPval1.in",
     header = F
   )

# Find intersection with regulator genes
common_genes <- Reduce(intersect, list(row.names(TS1), row.names(TS2), regulators$V1))
common_genes <- unlist(common_genes)

time.points <- list(unname(TS1[1, ]), unname(TS2[1, ]))

# Find intersection with regulator genes
common_genes <- Reduce(intersect, list(row.names(TS1), row.names(TS2), regulators$V1))
common_genes <- unlist(common_genes)

time.points <- list(unname(TS1[1, ]), unname(TS2[1, ]))

# Subset dataframes to include only intersected genes
TS1 <- TS1[common_genes, ,drop=FALSE]
TS2 <- TS2[common_genes, ,drop=FALSE]

#TS.data <- list(TS1[2:nrow(TS1),], TS2[2:nrow(TS2),])
TS.data <- list(TS1, TS2)

# dynGENIE3 run
dynGENIE3_results <- dynGENIE3(TS.data = TS.data,
                              time.points = time.points,
                              regulators = common_genes,
                               ncores=6)
