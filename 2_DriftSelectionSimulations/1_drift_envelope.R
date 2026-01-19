#!/usr/bin/env Rscript

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    ne = "inputs/raw/Ne.csv",
    freq = "inputs/raw/allele_frequencies_and_coverage.txt",
    out_dir = "outputs/drift_envelope",
    sample_n = 5000
  )
  if (length(args) == 0) return(defaults)
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
    key <- substring(key, 3)
    if (i == length(args)) {
      stop("Missing value for --", key)
    }
    value <- args[[i + 1]]
    if (key == "sample_n") {
      defaults[[key]] <- as.integer(value)
    } else if (key %in% names(defaults)) {
      defaults[[key]] <- value
    } else {
      warning("Ignoring unknown option --", key)
    }
    i <- i + 2
  }
  defaults
}

opt <- parse_args()

if (!file.exists(opt$ne)) {
  stop("Ne table not found at ", opt$ne)
}
if (!file.exists(opt$freq)) {
  stop("Allele frequency table not found at ", opt$freq)
}
if (!dir.exists(opt$out_dir)) {
  dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
}
plot_dir <- file.path(opt$out_dir, "plots")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
}

message("Reading Ne table from ", opt$ne)
ne_raw <- read.csv(opt$ne, check.names = FALSE, stringsAsFactors = FALSE)
if (!"TREAT" %in% names(ne_raw)) {
  stop("Ne table must include a 'TREAT' column.")
}
gen_cols <- grep("^G[0-9]+$", names(ne_raw), value = TRUE)
if (length(gen_cols) == 0) {
  stop("No generation columns (e.g., G1...G6) detected in Ne table.")
}

phi <- apply(ne_raw[, gen_cols, drop = FALSE], 1, function(values) {
  if (any(is.na(values) | values <= 0)) {
    return(NA_real_)
  }
  prod(1 - (1 / (2 * values)))
})
sd_multiplier <- sqrt(pmax(0, 1 - phi))
ne_params <- data.frame(pool = ne_raw$TREAT, sd_multiplier = sd_multiplier,
                        stringsAsFactors = FALSE)

message("Reading allele frequency table from ", opt$freq)
afrecords <- read.table(opt$freq, header = TRUE, sep = "\t", quote = "",
                        check.names = FALSE, stringsAsFactors = FALSE)
if (!"AF_pool_G1" %in% names(afrecords)) {
  stop("Allele frequency table must include 'AF_pool_G1'.")
}

pool_cols <- grep('^AF_pool_', names(afrecords), value = TRUE)
pool_cols <- setdiff(pool_cols, 'AF_pool_G1')
if (length(pool_cols) == 0) {
  stop('No treatment-specific AF columns detected.')
}

p0 <- afrecords[["AF_pool_G1"]]
summaries <- list()
theoretical_bounds <- list()
set.seed(42)

for (col_name in pool_cols) {
  pool_id <- sub("^AF_pool_", "", col_name)
  if (pool_id == "G1") {
    next
  }
  match_idx <- match(pool_id, ne_params$pool)
  if (is.na(match_idx)) {
    warning("Skipping pool ", pool_id, " (no matching Ne trajectory).")
    next
  }
  sd_mult <- ne_params$sd_multiplier[[match_idx]]
  if (is.na(sd_mult)) {
    warning("Skipping pool ", pool_id, " (invalid sd multiplier).")
    next
  }

  p_t <- afrecords[[col_name]]
  valid <- !is.na(p0) & !is.na(p_t)
  if (!any(valid)) {
    warning("No overlapping SNPs for pool ", pool_id)
    next
  }

  p0_valid <- p0[valid]
  abs_delta <- abs(p_t[valid] - p0_valid)
  bound <- 1.96 * sqrt(pmax(p0_valid * (1 - p0_valid), 0)) * sd_mult
  beyond <- abs_delta > bound

  summaries[[pool_id]] <- data.frame(
    pool = pool_id,
    snps = length(abs_delta),
    sd_multiplier = sd_mult,
    median_abs_delta = stats::median(abs_delta, na.rm = TRUE),
    p95_abs_delta = stats::quantile(abs_delta, 0.95, na.rm = TRUE, names = FALSE),
    median_bound = stats::median(bound, na.rm = TRUE),
    p95_bound = stats::quantile(bound, 0.95, na.rm = TRUE, names = FALSE),
    fraction_beyond = mean(beyond, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  line_p0 <- seq(0, 1, length.out = 200)
  line_bound <- 1.96 * sqrt(pmax(line_p0 * (1 - line_p0), 0)) * sd_mult
  theoretical_bounds[[pool_id]] <- data.frame(pool = pool_id, p0 = line_p0,
                                              bound = line_bound)

  sample_n <- min(opt$sample_n, length(abs_delta))
  sample_idx <- sample.int(length(abs_delta), sample_n)
  png_filename <- file.path(plot_dir, paste0("drift_envelope_", pool_id, ".png"))
  grDevices::png(png_filename, width = 6 * 96, height = 4 * 96, res = 96)
  par(mar = c(4.5, 4.5, 2, 1))
  plot(p0_valid[sample_idx], abs_delta[sample_idx], pch = 16,
       col = grDevices::rgb(44, 127, 184, 100, maxColorValue = 255),
       cex = 0.6,
       xlab = "Initial allele frequency (G1)",
       ylab = "|ΔAF|",
       main = paste0("|ΔAF| vs Neutral Envelope (", pool_id, ")"))
  lines(line_p0, line_bound, col = "#d95f0e", lwd = 2)
  suppressWarnings({
    dev.off()
  })
}

if (length(summaries) == 0) {
  stop("No pools processed. Check input files.")
}

summary_tbl <- do.call(rbind, summaries)
summary_tbl <- summary_tbl[order(summary_tbl$pool), ]
write.csv(summary_tbl, file.path(opt$out_dir, "drift_envelope_summary.csv"), row.names = FALSE)
write.csv(ne_params, file.path(opt$out_dir, "ne_sd_multipliers.csv"), row.names = FALSE)
if (length(theoretical_bounds) > 0) {
  bounds_tbl <- do.call(rbind, theoretical_bounds)
  write.csv(bounds_tbl, file.path(opt$out_dir, "theoretical_bounds.csv"), row.names = FALSE)
}

message("Drift envelope analysis complete. Outputs written to ", opt$out_dir)
