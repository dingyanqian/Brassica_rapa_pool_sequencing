#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

source("scripts/delta_af_common.R")

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list(
    freq = "inputs/raw/allele_frequencies_and_coverage.txt",
    iterations = 1000,
    seed = 20240906,
    out_dir = "outputs/△AF_classification"
  )
  if (length(args) %% 2 != 0) stop("Arguments must be provided as --key value pairs")
  if (length(args) > 0) {
    for (i in seq(1, length(args), by = 2)) {
      key <- sub("^--", "", args[[i]])
      val <- args[[i + 1]]
      if (!(key %in% names(opts))) {
        warning("Unknown option --", key)
        next
      }
      if (key %in% c("iterations", "seed")) {
        opts[[key]] <- as.integer(val)
      } else {
        opts[[key]] <- val
      }
    }
  }
  opts
}

args <- parse_args()
set.seed(args$seed)
if (!file.exists(args$freq)) stop("Frequency file not found: ", args$freq)
if (!dir.exists(args$out_dir)) dir.create(args$out_dir, recursive = TRUE)

delta_dt <- load_delta_table(args$freq)

assign_pairs <- function(def) {
  list(
    A = list(cold = def$cold[1], hot = def$hot[1]),
    B = list(cold = def$cold[2], hot = def$hot[2])
  )
}

permute_def <- function(def) {
  pairs <- assign_pairs(def)
  new_def <- def
  for (nm in names(pairs)) {
    if (runif(1) < 0.5) {
      tmp <- pairs[[nm]]$cold
      pairs[[nm]]$cold <- pairs[[nm]]$hot
      pairs[[nm]]$hot <- tmp
    }
  }
  new_def$cold <- c(pairs$A$cold, pairs$B$cold)
  new_def$hot <- c(pairs$A$hot, pairs$B$hot)
  new_def
}

all_categories <- c("background", "global_adaptation", "conditional_cold", "conditional_hot", "antagonistic_pleiotropy")

observed <- lapply(names(treatment_defs), function(label) {
  def <- treatment_defs[[label]]
  res <- classify_treatment(delta_dt, def)
  counts <- res[, .N, by = category]
  counts <- merge(data.table(category = all_categories), counts, by = "category", all.x = TRUE)
  counts[is.na(N), N := 0]
  counts[, `:=`(treatment = label, fraction = N / nrow(res))]
})
observed_dt <- rbindlist(observed)
observed_counts_path <- file.path(args$out_dir, "observed_counts_permutation.csv")
fwrite(observed_dt, observed_counts_path)

perm_records <- vector("list", args$iterations * length(treatment_defs))
idx <- 1

for (iter in seq_len(args$iterations)) {
  for (label in names(treatment_defs)) {
    def <- treatment_defs[[label]]
    perm_def <- permute_def(def)
    res <- classify_treatment(delta_dt, perm_def)
    counts <- res[, .N, by = category]
    counts <- merge(data.table(category = all_categories), counts, by = "category", all.x = TRUE)
    counts[is.na(N), N := 0]
    counts[, `:=`(iteration = iter, treatment = label)]
    perm_records[[idx]] <- counts
    idx <- idx + 1
  }
}

perm_dt <- rbindlist(perm_records, fill = TRUE)
perm_dt[is.na(N), N := 0]
perm_counts_path <- file.path(args$out_dir, "permutation_counts.csv")
fwrite(perm_dt, perm_counts_path)

obs_lookup <- observed_dt[, setNames(N, paste(treatment, category, sep = "|"))]
perm_dt[, observed := obs_lookup[paste(treatment, category, sep = "|")]]
perm_dt[is.na(observed), observed := 0]

summary_dt <- perm_dt[, .(mean_count = mean(N), sd_count = sd(N)), by = .(treatment, category)]
summary_dt <- merge(summary_dt, unique(observed_dt[, .(treatment, category, observed = N)]), by = c("treatment", "category"), all.x = TRUE)
pvals <- perm_dt[, .(p_value = (sum(N >= unique(observed)) + 1) / (args$iterations + 1)), by = .(treatment, category)]
summary_dt <- merge(summary_dt, pvals, by = c("treatment", "category"), all.x = TRUE)
summary_path <- file.path(args$out_dir, "permutation_summary.csv")
fwrite(summary_dt, summary_path)

message("Permutation results written to ", args$out_dir)
