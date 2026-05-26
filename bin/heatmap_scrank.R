library(vroom)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(scales)

# Load data
perb_score <- read.table("../scripts/nf_scrank/results_gazestani/rank_scores/perbscore_all_targets.txt", header = TRUE)

# Filter out 0s (to avoid log10 issues), mark as NA
perb_score <- perb_score %>%
  mutate(perb_score = ifelse(perb_score <= 0, NA, perb_score)) %>%
  mutate(perb_score = log10(perb_score))

# Identify outliers
iqr <- IQR(perb_score$perb_score, na.rm = TRUE)
q1 <- quantile(perb_score$perb_score, 0.25, na.rm = TRUE)
q3 <- quantile(perb_score$perb_score, 0.75, na.rm = TRUE)

lower_bound <- q1 - 1.5 * iqr

# Mark outliers as NA so they show up gray
perb_score <- perb_score %>%
  mutate(perb_score = ifelse(perb_score < lower_bound, NA, perb_score))

# Rescale across non-NA values
range_vals <- range(perb_score$perb_score, na.rm = TRUE)

perb_score <- perb_score %>%
  mutate(perb_score = rescale(perb_score, to = c(0, 1), from = range_vals))

# Convert to heatmap matrix
heatmap_matrix <- acast(perb_score, cell_type ~ target, value.var = "perb_score")

# Plot
Heatmap(
  heatmap_matrix,
  name = "log10\nscaled\nperturbation\nscore",
  row_title = "Cell Types",
  column_title = "Alzheimer first trial",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_heatmap_legend = TRUE,
  na_col = "gray"  # gray for outliers
)
