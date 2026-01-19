##############################
## Summarize enhanced GO results (finalized)
## - Only classic algorithm
## - FDR (BH) and Raw p at cutoffs 0.05 and 0.01
## - Raw counts only (no extra filters)
##############################

library(dplyr)
library(readr)
library(tidyr)
library(data.table)

output_dir <- "../OUTPUT"
if (!dir.exists(output_dir)) stop("Output directory not found: ../OUTPUT")

# Detect treatments from OUTPUT (fallback to default list)
all_files <- list.files(output_dir, pattern = "_GO_enrichment_results_enhanced\\.txt$", full.names = FALSE)
if (length(all_files) > 0) {
  treatments <- unique(sub("_.*$", "", all_files))
} else {
  treatments <- c("CB", "CP", "HB", "HG", "HH")
}

ontologies <- c("BP", "MF", "CC")
cutoffs <- c(0.01, 0.05)  # finalized cutoffs
algorithms <- c("classic")

safe_nrow <- function(path) {
  if (!file.exists(path)) return(0L)
  df <- tryCatch(suppressWarnings(read.table(path, sep = "\t", header = TRUE, quote = "")),
                 error = function(e) NULL)
  if (is.null(df)) return(0L)
  nrow(df)
}

safe_unique_genes <- function(path) {
  if (!file.exists(path)) return(0L)
  df <- tryCatch(suppressWarnings(read.table(path, sep = "\t", header = TRUE, quote = "", comment.char = "")),
                 error = function(e) NULL)
  if (is.null(df)) return(0L)
  if (!"locusName" %in% names(df)) return(0L)
  length(unique(df$locusName))
}

read_table_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(suppressWarnings(read.table(path, sep = "\t", header = TRUE, quote = "", comment.char = "")),
           error = function(e) NULL)
}

summary_rows <- list()

for (tr in treatments) {
  for (ont in ontologies) {
    alg <- "classic"
    base_alg <- paste0(tr, "_", ont, "_", alg)
    for (cut in cutoffs) {
      lab <- sprintf("p%.2f", cut)

      # Raw p-value filtered enrichment and mapping (kept for reference)
      enr_raw <- file.path(output_dir, paste0(base_alg, "_GO_enrichment_results_", lab, "_enhanced.txt"))
      if (!file.exists(enr_raw)) {
        enr_raw <- file.path(output_dir, paste0(tr, "_", ont, "_GO_enrichment_results_", lab, "_enhanced.txt"))
      }
      map_raw <- file.path(output_dir, paste0(base_alg, "_significant_genes_per_GO_term_", lab, "_enhanced.txt"))
      if (!file.exists(map_raw)) {
        map_raw <- file.path(output_dir, paste0(tr, "_", ont, "_significant_genes_per_GO_term_", lab, "_enhanced.txt"))
      }

      n_terms_raw <- safe_nrow(enr_raw)
      n_genes_raw <- safe_unique_genes(map_raw)

      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        treatment = tr,
        ontology = ont,
        algorithm = alg,
        threshold_type = "raw_p",
        cutoff = cut,
        n_terms = n_terms_raw,
        n_genes = n_genes_raw,
        stringsAsFactors = FALSE
      )

      # FDR (BH) filtered enrichment and mapping (0.01 and 0.05)
      enr_fdr <- file.path(output_dir, paste0(base_alg, "_GO_enrichment_results_FDR", lab, "_enhanced.txt"))
      if (!file.exists(enr_fdr)) {
        enr_fdr <- file.path(output_dir, paste0(tr, "_", ont, "_GO_enrichment_results_FDR", lab, "_enhanced.txt"))
      }
      map_fdr <- file.path(output_dir, paste0(base_alg, "_significant_genes_per_GO_term_FDR", lab, "_enhanced.txt"))
      if (!file.exists(map_fdr)) {
        map_fdr <- file.path(output_dir, paste0(tr, "_", ont, "_significant_genes_per_GO_term_FDR", lab, "_enhanced.txt"))
      }

      n_terms_fdr <- safe_nrow(enr_fdr)
      n_genes_fdr <- safe_unique_genes(map_fdr)

      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        treatment = tr,
        ontology = ont,
        algorithm = alg,
        threshold_type = "FDR_BH",
        cutoff = cut,
        n_terms = n_terms_fdr,
        n_genes = n_genes_fdr,
        stringsAsFactors = FALSE
      )
    }
  }
}

summary_df <- dplyr::bind_rows(summary_rows) %>%
  dplyr::arrange(treatment, ontology, algorithm, threshold_type, cutoff)

outfile <- file.path(output_dir, "GO_enrichment_summary_counts_enhanced.tsv")
fwrite(summary_df, file = outfile, sep = "\t")
cat("Wrote summary counts to:", outfile, "\n")
