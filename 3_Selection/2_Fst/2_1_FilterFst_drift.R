#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list(
    freq = "inputs/raw/allele_frequencies_and_coverage.txt",
    ne = "inputs/raw/Ne.csv",
    concordant = "outputs/final_analysis/concordant_flags_final95.csv",
    out_dir = "outputs/final_analysis/fst_drift_candidates_q95",
    threshold_out = "outputs/final_analysis/fst_drift_thresholds.csv",
    reps = 5000,
    seed = 20240905
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    key <- substring(key, 3)
    if (i == length(args)) stop("Missing value for --", key)
    value <- args[[i + 1]]
    if (key %in% names(opts)) {
      if (key %in% c("reps", "seed")) {
        opts[[key]] <- as.integer(value)
      } else {
        opts[[key]] <- value
      }
    } else {
      warning("Ignoring unknown option --", key)
    }
    i <- i + 2
  }
  opts
}

args <- parse_args()
set.seed(args$seed)

if (!file.exists(args$freq)) stop("Allele frequency file not found: ", args$freq)
if (!file.exists(args$ne)) stop("Ne file not found: ", args$ne)
if (!file.exists(args$concordant)) stop("Concordant SNP file not found: ", args$concordant)

message("Reading allele frequency table from ", args$freq)
freq_dt <- fread(args$freq)

pool_columns <- grep("^AF_pool_", names(freq_dt), value = TRUE)
if (length(pool_columns) == 0) stop("No AF_pool_* columns found")

g1_column <- if ("AF_pool_G1" %in% names(freq_dt)) "AF_pool_G1" else if ("AF_G1" %in% names(freq_dt)) "AF_G1" else stop("Missing G1 column")

calc_fst <- function(p, p0) {
  denom <- (p * (1 - p0)) + (p0 * (1 - p))
  out <- (p - p0)^2 / denom
  out[denom <= 0] <- 0
  out[is.na(out)] <- 0
  out
}

message("Calculating Hudson FST vs G1")
for (pool in pool_columns) {
  fst_name <- paste0("FST_", pool, "_G1")
  freq_dt[[fst_name]] <- calc_fst(freq_dt[[pool]], freq_dt[[g1_column]])
}

ne_tab <- fread(args$ne)
if (!"TREAT" %in% names(ne_tab)) stop("Ne table must include TREAT column")

p0_values <- freq_dt[[g1_column]]
p0_values <- p0_values[!is.na(p0_values) & p0_values > 0 & p0_values < 1]
if (length(p0_values) == 0) stop("No valid baseline allele frequencies found")

drift_step <- function(p, ne) {
  if (is.na(p) || is.na(ne)) return(NA_real_)
  copies <- max(1L, round(2 * ne))
  prob <- min(max(p, 0), 1)
  rbinom(1L, copies, prob) / copies
}

temporal_fst <- function(p0, pt, ne_vals) {
  if (is.na(p0) || is.na(pt)) return(0)
  p_bar <- 0.5 * (p0 + pt)
  denom <- p_bar * (1 - p_bar)
  if (denom <= 0) return(0)
  sampling <- p0 * (1 - p0) / (2 * mean(ne_vals))
  fst <- ((pt - p0)^2 - sampling) / denom
  max(fst, 0)
}

simulate_threshold <- function(ne_vals, reps) {
  res <- numeric(reps)
  for (i in seq_len(reps)) {
    p0 <- sample(p0_values, 1)
    pt <- p0
    for (ne in ne_vals) {
      pt <- drift_step(pt, ne)
    }
    res[[i]] <- temporal_fst(p0, pt, ne_vals)
  }
  quantile(res, probs = 0.95, na.rm = TRUE)
}

message("Simulating drift-based FST thresholds (q95)")
thr_list <- list()
for (pool in sub("^AF_pool_", "", pool_columns)) {
  ne_vals <- as.numeric(ne_tab[TREAT == pool, -1, with = FALSE])
  if (length(ne_vals) == 0) stop("No Ne values for pool ", pool)
  thr_q95 <- simulate_threshold(ne_vals, args$reps)
  thr_list[[pool]] <- data.frame(pool = pool, drift_fst_q95 = thr_q95, stringsAsFactors = FALSE)
}
thr_dt <- do.call(rbind, thr_list)
fwrite(thr_dt, args$threshold_out)
message("Saved drift thresholds to ", args$threshold_out)

concordant <- fread(args$concordant)

pair_map <- list(
  CB = c("CBA", "CBB"),
  CG = c("CGA", "CGB"),
  CH = c("CHA", "CHB"),
  CP = c("CPA", "CPB"),
  HB = c("HBA", "HBB"),
  HG = c("HGA", "HGB"),
  HH = c("HHA", "HHB"),
  HP = c("HPA", "HPB")
)

pair_lookup_vec <- character()
for (label in names(pair_map)) {
  pools <- pair_map[[label]]
  pair_lookup_vec[paste(pools[1], pools[2], sep = ":")] <- label
  pair_lookup_vec[paste(pools[2], pools[1], sep = ":")] <- label
}

concordant[, treatment := pair_lookup_vec[paste(pool_a, pool_b, sep = ":")]]
concordant <- concordant[!is.na(treatment)]

freq_dt_keyed <- as.data.table(freq_dt)
setkey(freq_dt_keyed, CHROM, POS)

if (!dir.exists(args$out_dir)) dir.create(args$out_dir, recursive = TRUE)

for (label in names(pair_map)) {
  pools <- pair_map[[label]]
  fst_col_a <- paste0("FST_AF_pool_", pools[1], "_G1")
  fst_col_b <- paste0("FST_AF_pool_", pools[2], "_G1")
  if (!all(c(fst_col_a, fst_col_b) %in% names(freq_dt_keyed))) {
    warning("Missing FST columns for ", label)
    next
  }
  thr_a <- thr_dt$drift_fst_q95[thr_dt$pool == pools[1]]
  thr_b <- thr_dt$drift_fst_q95[thr_dt$pool == pools[2]]
  subset_keys <- unique(concordant[treatment == label, .(CHROM, POS)])
  if (nrow(subset_keys) == 0) {
    fwrite(data.table(CHROM = character(), POS = integer()),
           file.path(args$out_dir, paste0("SNPs_", label, ".txt")), sep = "\t")
    next
  }
  sub_dt <- freq_dt_keyed[subset_keys, .(CHROM, POS, fst_a = get(fst_col_a), fst_b = get(fst_col_b))]
  filtered <- sub_dt[fst_a > thr_a & fst_b > thr_b]
  out_path <- file.path(args$out_dir, paste0("SNPs_", label, ".txt"))
  fwrite(filtered[, .(CHROM, POS)], out_path, sep = "\t")
  message(label, ": kept ", nrow(filtered), " SNPs")
}

message("Done. Outputs in ", args$out_dir)
