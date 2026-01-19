#!/usr/bin/env Rscript

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    freq = "inputs/raw/allele_frequencies_and_coverage.txt",
    out_dir = "outputs/empirical_null",
    baseline = "CHA,CHB,HHA,HHB",
    probs = "0.95,0.99"
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    key <- substring(key, 3)
    if (i == length(args)) stop("Missing value for --", key)
    value <- args[[i + 1]]
    if (key %in% names(defaults)) {
      defaults[[key]] <- value
    } else {
      warning("Ignoring unknown option --", key)
    }
    i <- i + 2
  }
  defaults
}

opt <- parse_args()

if (!file.exists(opt$freq)) {
  stop("Allele frequency table not found at ", opt$freq)
}
if (!dir.exists(opt$out_dir)) {
  dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
}

probs <- as.numeric(strsplit(opt$probs, ",")[[1]])
probs <- probs[!is.na(probs) & probs > 0 & probs < 1]
if (length(probs) == 0) {
  stop("No valid quantile probabilities provided (use values between 0 and 1).")
}
probs <- sort(unique(probs))

message("Reading allele frequency table from ", opt$freq)
freq <- read.table(opt$freq, header = TRUE, sep = "\t", quote = "",
                   check.names = FALSE, stringsAsFactors = FALSE)

if (!"AF_pool_G1" %in% names(freq)) {
  stop("Allele frequency table must include 'AF_pool_G1'.")
}

pool_cols <- grep("^AF_pool_", names(freq), value = TRUE)
pool_ids <- sub("^AF_pool_", "", pool_cols)
pool_ids <- setdiff(pool_ids, "G1")
baseline_ids <- strsplit(opt$baseline, ",")[[1]]
baseline_ids <- trimws(baseline_ids)
baseline_ids <- baseline_ids[baseline_ids != ""]
if (length(baseline_ids) == 0) {
  stop("No baseline pools specified.")
}
missing_baseline <- setdiff(baseline_ids, pool_ids)
if (length(missing_baseline) > 0) {
  stop("Baseline pools not found in data: ", paste(missing_baseline, collapse = ", "))
}

p0 <- as.numeric(freq[["AF_pool_G1"]])
chrom <- freq[["CHROM"]]
pos <- freq[["POS"]]

get_delta <- function(pool_id) {
  col_name <- paste0("AF_pool_", pool_id)
  values <- as.numeric(freq[[col_name]])
  abs(values - p0)
}

baseline_quantiles <- lapply(baseline_ids, function(pool_id) {
  delta <- get_delta(pool_id)
  valid <- delta[!is.na(delta)]
  if (length(valid) == 0) {
    stop("No valid |ΔAF| values found for baseline pool ", pool_id)
  }
  data.frame(
    pool = pool_id,
    prob = probs,
    threshold = as.numeric(quantile(valid, probs = probs, na.rm = TRUE, names = FALSE)),
    stringsAsFactors = FALSE
  )
})

baseline_tbl <- do.call(rbind, baseline_quantiles)
write.csv(baseline_tbl, file.path(opt$out_dir, "baseline_thresholds.csv"), row.names = FALSE)

# Global baseline used for non-baseline exceedance summaries
global_thresholds <- sapply(probs, function(p) {
  max(baseline_tbl$threshold[baseline_tbl$prob == p])
})
names(global_thresholds) <- probs

non_baseline_ids <- setdiff(pool_ids, baseline_ids)

summary_list <- list()
flag_q_list <- vector("list", length(probs))
names(flag_q_list) <- paste0("q", sprintf("%02d", as.integer(probs * 100)))

for (pool_id in non_baseline_ids) {
  delta <- get_delta(pool_id)
  valid <- !is.na(delta)
  total <- sum(valid)
  if (total == 0) {
    next
  }
  pool_summary <- list(pool = pool_id, total_snps = total)
  for (i in seq_along(probs)) {
    prob <- probs[[i]]
    prob_key <- as.character(prob)
    thr <- global_thresholds[[prob_key]]
    exceed <- valid & (delta > thr)
    count_exceed <- sum(exceed)
    pool_summary[[paste0("count_above_q", sprintf("%02d", as.integer(prob * 100)))]] <- count_exceed
    pool_summary[[paste0("frac_above_q", sprintf("%02d", as.integer(prob * 100)))]] <- count_exceed / total

    if (count_exceed > 0) {
      flag_df <- data.frame(
        CHROM = chrom[exceed],
        POS = pos[exceed],
        pool = pool_id,
        abs_delta = delta[exceed],
        threshold_prob = prob,
        stringsAsFactors = FALSE
      )
      flag_q_list[[i]] <- rbind(flag_q_list[[i]], flag_df)
    }
  }
  summary_list[[length(summary_list) + 1]] <- pool_summary
}

if (length(summary_list) == 0) {
  stop("No non-baseline pools processed.")
}

summary_tbl <- do.call(rbind, lapply(summary_list, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
summary_tbl <- summary_tbl[order(summary_tbl$pool), ]
write.csv(summary_tbl, file.path(opt$out_dir, "empirical_null_summary.csv"), row.names = FALSE)

for (i in seq_along(flag_q_list)) {
  df <- flag_q_list[[i]]
  if (!is.null(df) && nrow(df) > 0) {
    prob_label <- names(flag_q_list)[[i]]
    out_path <- file.path(opt$out_dir, paste0("empirical_null_flags_", prob_label, ".csv"))
    write.csv(df, out_path, row.names = FALSE)
  }
}

message("Empirical null analysis complete. Outputs written to ", opt$out_dir)
