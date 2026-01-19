####this script is check the baseline of differentiation of each datasets and find the xi numbers for testing local score

# Load libraries
library(tidyverse)

# Define your input folder
input_folder <- "/Users/s1950737/Desktop/ER_CMH/2_CMHoutput" #change accordingly

# Function to suggest xi based on average -log10(p)
suggest_xi_auto <- function(pvals, max_xi = 10, step = 0.5) {
  logp <- -log10(pvals)
  for (xi in seq(1, max_xi, by = step)) {
    score <- logp - xi
    if (mean(score, na.rm = TRUE) < 0) {
      return(xi)
    }
  }
  warning("⚠ No ξ found with negative mean score — using max_xi")
  return(max_xi)
}

# Create a folder for output plots
plot_dir <- "../3_diagnostics"
dir.create(plot_dir, showWarnings = FALSE)

# List all relevant files
file_list <- list.files(input_folder, pattern = "^cmh_results_.*\\.txt$", full.names = TRUE)

# Store summary results
summary_results <- data.frame()

# Process each file
for (file_path in file_list) {
  treatment_name <- str_extract(basename(file_path), "(?<=cmh_results_).*(?=\\.txt)")
  
  # Read the file
  df <- read_table(file_path, show_col_types = FALSE) %>%
    select(CHROM, POS, cmh_pval) %>%
    filter(!is.na(cmh_pval) & cmh_pval > 0 & cmh_pval <= 1)
  
  # Compute basic stats
  avg_neglog10p <- mean(-log10(df$cmh_pval))
  min_p <- min(df$cmh_pval)
  max_p <- max(df$cmh_pval)
  suggested <- suggest_xi_auto(df$cmh_pval)
  
  # Compute autocorrelation
  acf_values <- acf(-log10(df$cmh_pval), lag.max = 1, plot = FALSE)$acf
  autocorr <- acf_values[2]  # lag-1 autocorrelation
  
  # Save histogram
  hist_path <- file.path(plot_dir, paste0("hist_", treatment_name, ".png"))
  png(hist_path, width = 800, height = 600)
  hist(df$cmh_pval, breaks = 100, main = paste("Histogram of p-values:", treatment_name),
       xlab = "p-value", col = "steelblue", border = "white")
  dev.off()
  
  # Add to summary table
  summary_results <- rbind(summary_results, data.frame(
    Treatment = treatment_name,
    Mean_neglog10p = round(avg_neglog10p, 3),
    Min_p = signif(min_p, 3),
    Max_p = signif(max_p, 3),
    Autocorr_lag1 = round(autocorr, 3),
    Suggested_xi = suggested
  ))
}

# Write the summary table to a CSV
write_csv(summary_results, "lrt1_diagnostics_summary.csv")

# Show results
print(summary_results)
