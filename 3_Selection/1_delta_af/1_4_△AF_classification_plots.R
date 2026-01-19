library(data.table)
library(ggplot2)

# === Step 1: Load main ΔAF data ===
df_CH <- fread("./Selection_new/1_input/Fst_all_SNPs_unfiltered/all_SNPs_CH_unfiltered.txt")
df_HH <- fread("./Selection_new/1_input/Fst_all_SNPs_unfiltered/all_SNPs_HH_unfiltered.txt")

df_CH$AF_change_CH <- (df_CH$AF_change_CHA + df_CH$AF_change_CHB)/2
df_HH$AF_change_HH <- (df_HH$AF_change_HHA + df_HH$AF_change_HHB)/2

# === Step 2: Mark ER-CMH & FST overlap SNPs (keep original triangles) ===
overlap_snps_CH <- fread("./Selection_new/2_overlap_intersect/CH_cmh_Fst_overlap_snps.txt")[, .(CHROM, POS)]
df_CH$Overlap_sig_CH <- paste(df_CH$CHROM, df_CH$POS) %in% paste(overlap_snps_CH$CHROM, overlap_snps_CH$POS)

overlap_snps_HH <- fread("./Selection_new/2_overlap_intersect/HH_cmh_Fst_overlap_snps.txt")[, .(CHROM, POS)]
df_HH$Overlap_sig_HH <- paste(df_HH$CHROM, df_HH$POS) %in% paste(overlap_snps_HH$CHROM, overlap_snps_HH$POS)

# === Step 3: Merge datasets ===
df <- merge(df_CH[, .(CHROM, POS, AF_change_CH, Overlap_sig_CH)],
            df_HH[, .(CHROM, POS, AF_change_HH, Overlap_sig_HH)],
            by = c("CHROM", "POS"))

# === Step 4: Define new significance criteria based on |ΔAF| ≥ threshold ===
# Check if |ΔAF| ≥ threshold for each temperature regime
df$sig_CH <- abs(df$AF_change_CH) >= 0.343
df$sig_HH <- abs(df$AF_change_HH) >= 0.343

# === Step 5: Categorize SNPs based on new criteria ===
df$category <- "background"

# Antagonistic pleiotropy: opposite directions AND |ΔAF| ≥ threshold in both
df$category[df$sig_CH & df$sig_HH & 
              ((df$AF_change_CH > 0 & df$AF_change_HH < 0) | 
                 (df$AF_change_CH < 0 & df$AF_change_HH > 0))] <- "antagonistic_pleiotropy"

# Global adaptation: same direction AND |ΔAF| ≥ threshold in both
df$category[df$sig_CH & df$sig_HH & 
              ((df$AF_change_CH > 0 & df$AF_change_HH > 0) | 
                 (df$AF_change_CH < 0 & df$AF_change_HH < 0))] <- "global_adaptation"

# Conditional neutrality: |ΔAF| ≥ threshold in one regime but < threshold in the other
df$category[df$sig_CH & !df$sig_HH] <- "conditional_neutrality_CH"
df$category[!df$sig_CH & df$sig_HH] <- "conditional_neutrality_HH"

# === Step 6: Add overlap information for triangle markers ===
df$has_overlap_CH <- df$Overlap_sig_CH
df$has_overlap_HH <- df$Overlap_sig_HH
df$has_any_overlap <- df$has_overlap_CH | df$has_overlap_HH

# Factor levels for plot order
df$category <- factor(df$category,
                      levels = c("background", "conditional_neutrality_CH", "conditional_neutrality_HH", 
                                 "global_adaptation", "antagonistic_pleiotropy"))

# === Step 7: Create color and shape mappings ===
category_colors <- c(
  "background" = "lightgray",
  "conditional_neutrality_CH" = "#5894C8",      # Blue for CH-specific
  "conditional_neutrality_HH" = "#EB7C6A",      # Red for HH-specific  
  "global_adaptation" = "#2E8B57",              # Dark green for same direction
  "antagonistic_pleiotropy" = "#8B008B"         # Dark magenta for opposite directions
)

# === Step 8: Plot ===
p <- ggplot(df, aes(x = AF_change_CH, y = AF_change_HH)) +
  # Plot background points first
  geom_point(data = df[category == "background"], 
             color = category_colors["background"], size = 0.6) +
  
  # Plot significant categories
  geom_point(data = df[category == "conditional_neutrality_CH"], 
             color = category_colors["conditional_neutrality_CH"], size = 1.0) +
  geom_point(data = df[category == "conditional_neutrality_HH"], 
             color = category_colors["conditional_neutrality_HH"], size = 1.0) +
  geom_point(data = df[category == "global_adaptation"], 
             color = category_colors["global_adaptation"], size = 1.2) +
  geom_point(data = df[category == "antagonistic_pleiotropy"], 
             color = category_colors["antagonistic_pleiotropy"], size = 1.2) +
  
  # Add triangle overlays for overlap SNPs
  geom_point(data = df[has_overlap_CH == TRUE], aes(x = AF_change_CH, y = AF_change_HH),
             fill = "blue", size = 4.0, shape = 24, alpha = 0.7) +
  geom_point(data = df[has_overlap_HH == TRUE], aes(x = AF_change_CH, y = AF_change_HH),
             fill = "red", size = 4.0, shape = 24, alpha = 0.7) +
  
  # Add reference lines
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.0) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.0) +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = "dotted", alpha = 0.6) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dotted", alpha = 0.6) +
  
  # Theme and labels
  theme_minimal(base_size = 12) +
  labs(
    x = "Allele Frequency Change (Cold Bumblebee)",
    y = "Allele Frequency Change (Hot Bumblebee)",
    title = "SNP Classification Based on Temperature-Dependent Selection",
    subtitle = "Circles: |ΔAF| ≥ 0.2 significance; Triangles: Overlap with other selection signals"
  ) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1))

# === Step 9: Add legend ===
# Create summary statistics
cat("\n=== SNP Classification Summary ===\n")
cat("Total SNPs:", nrow(df), "\n")
cat("Background (|ΔAF| < 0.2 in both):", sum(df$category == "background"), "\n")
cat("Conditional neutrality (CH only):", sum(df$category == "conditional_neutrality_CH"), "\n")
cat("Conditional neutrality (HH only):", sum(df$category == "conditional_neutrality_HH"), "\n")
cat("Global adaptation (same direction):", sum(df$category == "global_adaptation"), "\n")
cat("Antagonistic pleiotropy (opposite direction):", sum(df$category == "antagonistic_pleiotropy"), "\n")
cat("Overlap with CH selection signals:", sum(df$has_overlap_CH), "\n")
cat("Overlap with HH selection signals:", sum(df$has_overlap_HH), "\n")

# Save plot
ggsave("AlleleFrequencyChange_new_classification.png", p, width = 12, height = 10, dpi = 300)

# Optional: Save classified data
fwrite(df[, .(CHROM, POS, AF_change_CH, AF_change_HH, category, has_overlap_CH, has_overlap_HH)], 
       "classified_SNPs_summary.txt", sep = "\t")

