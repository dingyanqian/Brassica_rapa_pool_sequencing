##############################
## Visualize GO summary counts (enhanced outputs)
##############################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 is required. Please install.packages('ggplot2') and rerun.")
}
library(ggplot2)

outdir <- "../OUTPUT"
infile <- file.path(outdir, "GO_enrichment_summary_counts_enhanced.tsv")
if (!file.exists(infile)) stop(paste("Summary file not found:", infile))

df <- fread(infile, sep = "\t", header = TRUE)

# Normalize types and ordering
df <- df %>%
  filter(algorithm == "classic", cutoff %in% c(0.01, 0.05)) %>%
  mutate(
    treatment = factor(treatment, levels = sort(unique(treatment))),
    ontology = factor(ontology, levels = c("BP", "MF", "CC")),
    threshold_type = factor(threshold_type, levels = c("FDR_BH", "raw_p"), labels = c("FDR (BH)", "Raw p")),
    cutoff = factor(sprintf("%.2f", as.numeric(cutoff)))
  )

out_terms <- file.path(outdir, "GO_summary_terms_heatmap.pdf")
out_genes <- file.path(outdir, "GO_summary_genes_dots.pdf")
out_genes_rawp <- file.path(outdir, "GO_summary_genes_dots_rawp.pdf")
out_fdr_effect <- file.path(outdir, "GO_summary_rawp_001_vs_005_heatmap.pdf")

# 1) Heatmap: number of enriched GO terms
df_terms <- df
p_terms <- ggplot(df_terms, aes(x = cutoff, y = treatment, fill = n_terms)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = n_terms), size = 3) +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b") +
  facet_grid(threshold_type ~ ontology, scales = "free_x", space = "free_x") +
  labs(
    title = "Enriched GO terms (raw counts): FDR and Raw p",
    x = "Cutoff",
    y = "Treatment",
    fill = "# Terms"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), legend.position = "right")

ggsave(out_terms, p_terms, width = 8.5, height = 3.8, device = "pdf")
cat("Saved:", out_terms, "\n")

# 2) Dot plot: unique contributing genes (FDR overlay; RAW counts)
df_genes_fdr_005 <- df %>% filter(threshold_type == "FDR (BH)", cutoff == "0.05")
df_genes_fdr_001 <- df %>% filter(threshold_type == "FDR (BH)", cutoff == "0.01")
p_genes <- ggplot() +
  # Base: FDR 0.05 (underneath)
  geom_point(data = df_genes_fdr_005,
             aes(x = cutoff, y = treatment, size = n_genes),
             color = "#9ecae1", alpha = 0.8) +
  # Top: FDR 0.01 (on top with labels)
  geom_point(data = df_genes_fdr_001,
             aes(x = cutoff, y = treatment, size = n_genes),
             color = "#d62728", alpha = 0.95) +
  geom_text(data = df_genes_fdr_001,
            aes(x = cutoff, y = treatment, label = n_genes),
            vjust = -0.7, size = 3, color = "#d62728") +
  scale_size_continuous(range = c(2, 10)) +
  facet_wrap(~ ontology, nrow = 1) +
  labs(
    title = "Unique contributing genes (raw counts): FDR 0.05 (blue) under 0.01 (red)",
    x = "FDR cutoff",
    y = "Treatment",
    size = "# Genes"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "right")

ggsave(out_genes, p_genes, width = 8.5, height = 3.2, device = "pdf")
cat("Saved:", out_genes, "\n")

# 2b) Dot plot: unique contributing genes (Raw p overlay; RAW counts)
df_genes_raw_005 <- df %>% filter(threshold_type == "Raw p", cutoff == "0.05")
df_genes_raw_001 <- df %>% filter(threshold_type == "Raw p", cutoff == "0.01")
p_genes_raw <- ggplot() +
  geom_point(data = df_genes_raw_005,
             aes(x = cutoff, y = treatment, size = n_genes),
             color = "#9ecae1", alpha = 0.8) +
  geom_point(data = df_genes_raw_001,
             aes(x = cutoff, y = treatment, size = n_genes),
             color = "#d62728", alpha = 0.95) +
  geom_text(data = df_genes_raw_001,
            aes(x = cutoff, y = treatment, label = n_genes),
            vjust = -0.7, size = 3, color = "#d62728") +
  scale_size_continuous(range = c(2, 10)) +
  facet_wrap(~ ontology, nrow = 1) +
  labs(
    title = "Unique contributing genes (raw counts): Raw p 0.05 (blue) under 0.01 (red)",
    x = "Raw p cutoff",
    y = "Treatment",
    size = "# Genes"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "right")

ggsave(out_genes_rawp, p_genes_raw, width = 8.5, height = 3.2, device = "pdf")
cat("Saved:", out_genes_rawp, "\n")

# 3) Heatmap of difference between FDR 0.01 and 0.05 term counts
df_delta <- df %>% filter(threshold_type == "FDR (BH)") %>%
  select(treatment, ontology, cutoff, n_terms) %>%
  tidyr::pivot_wider(names_from = cutoff, values_from = n_terms, values_fill = 0) %>%
  mutate(delta_001_005 = as.numeric(`0.01`) - as.numeric(`0.05`))

p_effect <- ggplot(df_delta, aes(x = factor(1), y = treatment, fill = delta_001_005)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = delta_001_005), size = 3) +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0) +
  facet_wrap(~ ontology, nrow = 1) +
  labs(
    title = "Δ enriched terms (FDR 0.01 - 0.05)",
    x = "",
    y = "Treatment",
    fill = "Δ Terms"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), legend.position = "right",
        axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(out_fdr_effect, p_effect, width = 8.5, height = 3.2, device = "pdf")
cat("Saved:", out_fdr_effect, "\n")
