# ---- Libraries and functions ----
library(data.table)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

source("scorelocalfunctions.R")  # This must be in the same directory

# ---- Parameters ----
input_folder <- "/home/yading/scratch/ER_CMH/input"
output_folder <- "/home/yading/scratch/ER_CMH/local_score_dual_xi_outputs/"
dir.create(output_folder, showWarnings = FALSE)

xi_values <- c(1, 2)
alpha_levels <- c(0.05, 0.01)

# ---- Main Loop ----
file_list <- list.files(input_folder, pattern = "^cmh_results_.*\\.txt$", full.names = TRUE)

for (file_path in file_list) {
  treatment <- gsub("^cmh_results_|\\.txt$", "", basename(file_path))
  message("\nProcessing treatment: ", treatment)
  
  # Load data
  mydata_original <- fread(file_path, select = c("CHROM", "POS", "cmh_pval"))
  colnames(mydata_original) <- c("chr", "pos", "pval")
  mydata <- mydata_original[complete.cases(mydata_original), ]
  mydata$pval[mydata$pval == 0] <- 1e-16
  setkey(mydata, chr)
  
  # Prepare genome-wide position
  Nchr <- length(unique(mydata$chr))
  chrInf <- mydata[, .(L = .N, cor = autocor(pval)), by = chr]
  setkey(chrInf, chr)
  tmp <- data.table(chr = unique(mydata$chr), S = cumsum(c(0, chrInf$L[-Nchr])))
  setkey(tmp, chr)
  mydata[tmp, posT := pos + S]
  
  for (xi in xi_values) {
    cat("  ➤ Running xi =", xi, "\n")
    
    # Compute score and Lindley process
    mydata[, score := -log10(pval) - xi]
    mydata[, lindley := lindley(score), by = chr]
    
    # Estimate Gumbel coefficients for non-uniform p-values
    message("    ... estimating Gumbel coefficients (may take time)")
    coefsG <- coefsGumb(mydata, Ls = seq(10000, 70000, 10000), nSeq = 5000)
    
    # Compute thresholds per chromosome
    chrInf[, thG05 := threshold(L, cor, coefsG$aCoef, coefsG$bCoef, 0.05)]
    chrInf[, thG01 := threshold(L, cor, coefsG$aCoef, coefsG$bCoef, 0.01)]
    
    mydata <- mydata[chrInf]
    
    # Save full data with local scores
    out_data_path <- file.path(output_folder, paste0("cmh_LocalScore_xi", xi, "_", treatment, ".txt"))
    fwrite(mydata, file = out_data_path, sep = "\t", quote = FALSE)
    
    # Detect significant regions
    sigZones05 <- mydata[, sig_sl(lindley, pos, unique(thG05)), by = chr]
    sigZones01 <- mydata[, sig_sl(lindley, pos, unique(thG01)), by = chr]
    
    # Save significant zones
    fwrite(sigZones05[peak > 0], file = file.path(output_folder, paste0("SL_xi", xi, "_signif5_", treatment, ".txt")), sep = "\t")
    fwrite(sigZones01[peak > 0], file = file.path(output_folder, paste0("SL_xi", xi, "_signif1_", treatment, ".txt")), sep = "\t")
    
    # Plot Lindley process and thresholds
    pdf(file.path(output_folder, paste0("ScoreLocalAll_xi", xi, "_", treatment, ".pdf")))
    par(mfrow = c(4,1), mar = c(5,2,1,1))
    for (g in unique(mydata$chr)) {
      plot(mydata[chr == g, pos], mydata[chr == g, lindley], type = "l",
           xlab = paste("Position chr", g), ylab = "Local Score",
           ylim = c(0, max(mydata[chr == g, max(lindley)], mydata[chr == g, unique(thG01)])))
      abline(h = mydata[chr == g, unique(thG05)], col = 'grey')
      abline(h = mydata[chr == g, unique(thG01)], col = 'grey', lty = 2)
      abline(v = sigZones01[chr == g, beg], col = 'grey', lty = 3)
      abline(v = sigZones01[chr == g, end], col = 'grey', lty = 3)
    }
    dev.off()
    
    # Manhattan plot
    manhattan_plot <- create_manhattan_plot(
      data = mydata,
      title = paste("CMH -", treatment, "- xi =", xi),
      threshold = -log10(0.05 / nrow(mydata)),
      chr_column = "chr",
      bp_column = "pos",
      p_column = "pval"
    )
    ggsave(file.path(output_folder, paste0("manhattan_plot_cmh_xi", xi, "_", treatment, ".png")),
           manhattan_plot, width = 10, height = 6, dpi = 300)
  }
}
