#!/usr/bin/env Rscript

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    freq = "inputs/raw/allele_frequencies_and_coverage.txt",
    out_dir = "outputs/replicate_concordance",
    pairs = "CBA:CBB,CGA:CGB,CHA:CHB,CPA:CPB,HBA:HBB,HGA:HGB,HHA:HHB,HPA:HPB",
    threshold = "0.342948717948718",
    sample_n = "5000"
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

freq_path <- opt$freq
if (!file.exists(freq_path)) stop("Allele frequency table not found at ", freq_path)

out_dir <- opt$out_dir
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(out_dir, "plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

delta_threshold <- as.numeric(opt$threshold)
if (is.na(delta_threshold))
  delta_threshold <- 0
sample_n <- as.integer(opt$sample_n)
if (is.na(sample_n) || sample_n <= 0)
  sample_n <- 5000

pair_tokens <- strsplit(opt$pairs, ",")[[1]]
pair_tokens <- trimws(pair_tokens)
pairs <- lapply(pair_tokens, function(tok) {
  parts <- strsplit(tok, ":")[[1]]
  parts <- trimws(parts)
  if (length(parts) != 2) stop("Invalid pair specification: ", tok)
  parts
})

message("Reading allele frequency table from ", freq_path)
freq <- read.table(freq_path, header = TRUE, sep = "\t", quote = "",
                   check.names = FALSE, stringsAsFactors = FALSE)

if (!"AF_pool_G1" %in% names(freq)) stop("AF_pool_G1 column missing.")

p0 <- as.numeric(freq[["AF_pool_G1"]])
chrom <- freq[["CHROM"]]
pos <- freq[["POS"]]

summary_rows <- list()
flagged_list <- list()
plots_created <- 0

for (pair in pairs) {
  pool_a <- pair[[1]]
  pool_b <- pair[[2]]
  col_a <- paste0("AF_pool_", pool_a)
  col_b <- paste0("AF_pool_", pool_b)
  if (!col_a %in% names(freq) || !col_b %in% names(freq)) {
    warning("Skipping pair ", pool_a, ":", pool_b, " (missing columns)")
    next
  }
  pa <- as.numeric(freq[[col_a]])
  pb <- as.numeric(freq[[col_b]])
  valid <- !is.na(p0) & !is.na(pa) & !is.na(pb)
  if (!any(valid)) {
    warning("No overlapping SNPs for pair ", pool_a, ":", pool_b)
    next
  }
  delta_a <- pa[valid] - p0[valid]
  delta_b <- pb[valid] - p0[valid]
  sign_consistency <- mean(sign(delta_a) == sign(delta_b))
  corr <- suppressWarnings(cor(delta_a, delta_b))
  joint_exceed <- (abs(delta_a) > delta_threshold) & (abs(delta_b) > delta_threshold)
  same_sign <- sign(delta_a) == sign(delta_b)
  joint_consistent <- joint_exceed & same_sign
  frac_joint_exceed <- mean(joint_exceed)
  frac_joint_consistent <- mean(joint_consistent)

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    pair = paste0(pool_a, ":", pool_b),
    snps = length(delta_a),
    pearson_corr = ifelse(is.finite(corr), corr, NA_real_),
    sign_consistency = sign_consistency,
    frac_both_above_threshold = frac_joint_exceed,
    frac_both_above_threshold_same_sign = frac_joint_consistent,
    stringsAsFactors = FALSE
  )

  if (any(joint_consistent)) {
    flags <- data.frame(
      CHROM = chrom[valid][joint_consistent],
      POS = pos[valid][joint_consistent],
      pool_a = pool_a,
      pool_b = pool_b,
      delta_a = delta_a[joint_consistent],
      delta_b = delta_b[joint_consistent],
      abs_delta_a = abs(delta_a[joint_consistent]),
      abs_delta_b = abs(delta_b[joint_consistent]),
      stringsAsFactors = FALSE
    )
    flagged_list[[length(flagged_list) + 1]] <- flags
  }

  n_obs <- length(delta_a)
  idx <- seq_len(n_obs)
  if (n_obs > sample_n) {
    set.seed(123)
    idx <- sample(idx, sample_n)
  }
  png(file.path(plot_dir, paste0(pool_a, "_", pool_b, "_scatter.png")),
      width = 6 * 96, height = 6 * 96, res = 96)
  par(mar = c(4.5, 4.5, 2, 1))
  plot(delta_a[idx], delta_b[idx],
       pch = 16,
       col = grDevices::rgb(44, 127, 184, 120, maxColorValue = 255),
       xlab = paste0("ΔAF (", pool_a, ")"),
       ylab = paste0("ΔAF (", pool_b, ")"),
       main = paste0(pool_a, " vs ", pool_b))
  abline(0, 1, col = "#d95f0e", lwd = 2)
  abline(h = c(-delta_threshold, delta_threshold), v = c(-delta_threshold, delta_threshold),
         col = "#999999", lty = 2)
  dev.off()
  plots_created <- plots_created + 1
}

if (length(summary_rows) == 0) stop("No replicate pairs processed.")

summary_tbl <- do.call(rbind, summary_rows)
summary_tbl <- summary_tbl[order(summary_tbl$pair), ]
write.csv(summary_tbl, file.path(out_dir, "replicate_concordance_summary.csv"), row.names = FALSE)

if (length(flagged_list) > 0) {
  flagged_tbl <- do.call(rbind, flagged_list)
  write.csv(flagged_tbl, file.path(out_dir, "replicate_concordant_flags.csv"), row.names = FALSE)
}

message("Replicate concordance analysis complete. Summary written to ",
        file.path(out_dir, "replicate_concordance_summary.csv"),
        "; plots generated: ", plots_created)
