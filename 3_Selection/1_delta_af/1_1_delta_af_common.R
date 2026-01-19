library(data.table)

pool_thresholds <- c(
  CBA = 0.35, CPA = 0.35, CGA = 0.35, CHA = 0.35,
  CBB = 0.34, CPB = 0.34, CGB = 0.34, CHB = 0.34,
  HBA = 0.32, HPA = 0.32, HGA = 0.32, HHA = 0.32,
  HBB = 0.36, HPB = 0.36, HGB = 0.36, HHB = 0.36
)

treatment_defs <- list(
  CB = list(label = "CB", cold = c("CBA", "CBB"), hot = c("HBA", "HBB")),
  CP = list(label = "CP", cold = c("CPA", "CPB"), hot = c("HPA", "HPB")),
  CG = list(label = "CG", cold = c("CGA", "CGB"), hot = c("HGA", "HGB")),
  CH = list(label = "CH", cold = c("CHA", "CHB"), hot = c("HHA", "HHB"))
)

load_delta_table <- function(freq_path) {
  needed <- c("CHROM", "POS", "AF_pool_G1", paste0("AF_pool_", names(pool_thresholds)))
  dt <- fread(freq_path, select = needed, showProgress = FALSE)
  for (pool in names(pool_thresholds)) {
    col <- paste0("AF_pool_", pool)
    dt[[pool]] <- dt[[col]] - dt[["AF_pool_G1"]]
  }
  dt[, AF_pool_G1 := NULL]
  for (pool in names(pool_thresholds)) {
    dt[, paste0("AF_pool_", pool) := NULL]
  }
  dt
}

classify_treatment <- function(delta_dt, def) {
  cold_cols <- def$cold
  hot_cols <- def$hot
  cold_mat <- as.matrix(delta_dt[, ..cold_cols])
  hot_mat <- as.matrix(delta_dt[, ..hot_cols])
  cold_valid <- rowSums(is.na(cold_mat)) == 0
  hot_valid <- rowSums(is.na(hot_mat)) == 0
  cold_sig <- cold_valid
  hot_sig <- hot_valid
  for (i in seq_along(cold_cols)) {
    thr <- pool_thresholds[[cold_cols[i]]]
    cold_sig <- cold_sig & (abs(cold_mat[, i]) >= thr)
  }
  for (i in seq_along(hot_cols)) {
    thr <- pool_thresholds[[hot_cols[i]]]
    hot_sig <- hot_sig & (abs(hot_mat[, i]) >= thr)
  }
  cold_mean <- rowMeans(cold_mat, na.rm = TRUE)
  cold_mean[!cold_valid] <- NA
  hot_mean <- rowMeans(hot_mat, na.rm = TRUE)
  hot_mean[!hot_valid] <- NA
  levels <- c("background", "global_adaptation", "conditional_cold", "conditional_hot", "antagonistic_pleiotropy")
  category <- rep("background", nrow(delta_dt))
  same_dir <- (cold_mean > 0 & hot_mean > 0) | (cold_mean < 0 & hot_mean < 0)
  opp_dir <- (cold_mean > 0 & hot_mean < 0) | (cold_mean < 0 & hot_mean > 0)
  idx_global <- which(cold_sig & hot_sig & same_dir)
  idx_antag <- which(cold_sig & hot_sig & opp_dir)
  idx_cond_cold <- which(cold_sig & !hot_sig)
  idx_cond_hot <- which(!cold_sig & hot_sig)
  category[idx_global] <- "global_adaptation"
  category[idx_antag] <- "antagonistic_pleiotropy"
  category[idx_cond_cold] <- "conditional_cold"
  category[idx_cond_hot] <- "conditional_hot"
  category <- factor(category, levels = levels)
  data.table(
    CHROM = delta_dt$CHROM,
    POS = delta_dt$POS,
    cold_delta = cold_mean,
    hot_delta = hot_mean,
    cold_sig = cold_sig,
    hot_sig = hot_sig,
    category = category
  )
}

summarise_categories <- function(class_dt) {
  total <- nrow(class_dt)
  class_dt[, .N, by = category][, .(category, count = N, fraction = N / total)]
}
