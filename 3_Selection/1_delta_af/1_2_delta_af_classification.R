#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

source("scripts/delta_af_common.R")

args <- list(
  freq = "inputs/raw/allele_frequencies_and_coverage.txt",
  out_dir = "outputs/△AF_classification"
)

cmd <- commandArgs(trailingOnly = TRUE)
if (length(cmd) %% 2 != 0) stop("Arguments must be provided as --key value pairs")
if (length(cmd) > 0) {
  for (i in seq(1, length(cmd), by = 2)) {
    key <- sub("^--", "", cmd[[i]])
    val <- cmd[[i + 1]]
    if (key %in% names(args)) {
      args[[key]] <- val
    } else {
      warning("Unknown option --", key)
    }
  }
}

if (!file.exists(args$freq)) stop("Frequency file not found: ", args$freq)
if (!dir.exists(args$out_dir)) dir.create(args$out_dir, recursive = TRUE)

message("Loading allele-frequency table...")
delta_dt <- load_delta_table(args$freq)

all_results <- list()

for (label in names(treatment_defs)) {
  message("Classifying ", label)
  def <- treatment_defs[[label]]
  res <- classify_treatment(delta_dt, def)
  res[, `:=`(
    treatment = label,
    cold_label = paste(def$cold, collapse = ","),
    hot_label = paste(def$hot, collapse = ",")
  )]
  setnames(res, c("cold_delta", "hot_delta"), c(paste0("AF_change_", label, "_cold"), paste0("AF_change_", label, "_hot")))
  fwrite(res, file.path(args$out_dir, paste0("classified_SNPs_", label, ".tsv")), sep = "\t")
  all_results[[label]] <- res
}

categories <- levels(all_results[[1]]$category)
summary_counts <- rbindlist(lapply(names(all_results), function(label) {
  dt <- all_results[[label]]
  counts <- dt[, .N, by = category]
  counts <- merge(data.table(category = categories), counts, by = "category", all.x = TRUE)
  counts[is.na(N), N := 0]
  counts[, `:=`(treatment = label, fraction = N / nrow(dt))]
}))
setcolorder(summary_counts, c("treatment", "category", "N", "fraction"))
fwrite(summary_counts, file.path(args$out_dir, "classification_counts.csv"))
message("Saved classification summaries to ", file.path(args$out_dir, "classification_counts.csv"))
