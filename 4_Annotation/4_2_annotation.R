##############################
## Enhanced GO enrichment + back-mapping
##############################

library(dplyr)
library(readr)
library(topGO)
library(clusterProfiler)
library(org.At.tair.db)
library(tidyr)
library(data.table)

# Configuration
# Evaluate multiple p-value cutoffs for GO term significance (raw p)
go_p_cutoffs <- c(0.01, 0.05, 0.10)

# Background annotation (used for GO enrichment)
background_file <- "../1_input/1_BrapaFPsc_277_v1.3.annotation_info.txt"
all_genes <- read.table(
  background_file,
  sep = "\t", header = TRUE, fill = TRUE, quote = "",
  comment.char = "", stringsAsFactors = FALSE
)

allGeneIDs <- unique(all_genes$locusName)

# Build robust gene->GO mapping: split multiple GO IDs and trim
gene2go_long <- all_genes %>%
  dplyr::select(locusName, GO) %>%
  dplyr::filter(!is.na(locusName), !is.na(GO), GO != "") %>%
  dplyr::mutate(
    GO = as.character(GO),
    GO = strsplit(GO, "[,;]")
  ) %>%
  tidyr::unnest(GO) %>%
  dplyr::mutate(GO = trimws(GO)) %>%
  dplyr::filter(GO != "")

geneID2GO <- split(gene2go_long$GO, gene2go_long$locusName)
geneID2GO <- geneID2GO[!is.na(names(geneID2GO))]

# Define treatments. Assumes for each treatment you have two files:
# 1. "../2_annotated_genes/<TREATMENT>_annotated_significant_genes.txt" (annotated table)
# 2. "../2_annotated_genes/<TREATMENT>_annotated_genes.txt" (one gene ID per line)
treatments <- c("CB", "CP", "HB", "HG", "HH")

for (treatment in treatments) {
  cat("Processing treatment:", treatment, "\n")

  # Input paths
  sig_file <- paste0("../2_annotated_genes/", treatment, "_annotated_significant_genes.txt")
  gene_file <- paste0("../2_annotated_genes/", treatment, "_annotated_genes.txt")

  ###########################
  ## Part I. Annotation processing
  ###########################

  annotated_genes <- read.table(
    sig_file, sep = "\t", header = FALSE, fill = TRUE,
    quote = "", comment.char = "", stringsAsFactors = FALSE
  )

  colnames(annotated_genes) <- c(
    "pacId", "locusName", "transcriptName", "peptideName",
    "Pfam", "Panther", "KOG", "KEGG.ec", "KO", "GO",
    "Best_hit_arabi_name", "arabi_symbol", "arabi_defline"
  )

  unique_annotated_genes <- annotated_genes %>% dplyr::distinct(locusName, .keep_all = TRUE)
  cleaned_file <- paste0(treatment, "_cleaned_annotated_genes.txt")
  write.table(unique_annotated_genes, cleaned_file, sep = "\t", row.names = FALSE, quote = FALSE)

  functional_summary <- unique_annotated_genes %>%
    dplyr::select(locusName, GO, KO, Best_hit_arabi_name, arabi_symbol, arabi_defline)
  func_summary_file <- paste0(treatment, "_functional_summary.txt")
  write.table(functional_summary, func_summary_file, sep = "\t", row.names = FALSE, quote = FALSE)

  ###########################
  ## Part II. GO enrichment for BP, MF, CC
  ###########################

  # Selected genes for this treatment (first column)
  sigGenes_data <- read.table(gene_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  sigGenes <- sigGenes_data$V1

  # Build 0/1 vector over background
  geneList <- factor(ifelse(allGeneIDs %in% sigGenes, 1, 0))
  names(geneList) <- allGeneIDs

  ontologies <- c("BP", "MF", "CC")
  for (ont in ontologies) {
    cat("  Ontology:", ont, "\n")

    # topGO object per ontology
    GOdata <- new(
      "topGOdata",
      ontology = ont,
      allGenes = geneList,
      annot = annFUN.gene2GO,
      gene2GO = geneID2GO
    )

    # Fisher's exact test; numeric p-values
    resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    pvals <- score(resultFisher)

    # Enrichment table (unfiltered) with numeric and FDR-adjusted p-values
    top_n <- length(pvals)
    go_enrichment_results <- GenTable(
      GOdata,
      classicFisher = resultFisher,
      orderBy = "classicFisher",
      topNodes = top_n
    )
    # Add numeric p-values and BH-FDR
    go_enrichment_results$classicFisher_numeric <- pvals[go_enrichment_results$GO.ID]
    pvals_fdr <- stats::p.adjust(pvals, method = "BH")
    go_enrichment_results$classicFisher_FDR_BH <- pvals_fdr[go_enrichment_results$GO.ID]
    base <- paste0(treatment, "_", ont)
    go_enrichment_file <- paste0(base, "_GO_enrichment_results_enhanced.txt")
    write.table(go_enrichment_results, go_enrichment_file, sep = "\t", row.names = FALSE, quote = FALSE)

    # For each requested cutoff, write filtered enrichment and GO->gene mapping
    for (cutoff in go_p_cutoffs) {
      cutoff_label <- sprintf("p%.2f", cutoff)  # e.g., p0.05, p0.10

      # Significant terms by numeric threshold
      sig_go_terms <- names(pvals)[pvals < cutoff]

      # Save filtered enrichment table for this cutoff
      go_enrichment_sig <- go_enrichment_results %>%
        dplyr::filter(!is.na(classicFisher_numeric), classicFisher_numeric < cutoff)
      go_enrichment_file_sig <- paste0(base, "_GO_enrichment_results_", cutoff_label, "_enhanced.txt")
      write.table(go_enrichment_sig, go_enrichment_file_sig, sep = "\t", row.names = FALSE, quote = FALSE)

      # Also provide FDR-filtered enrichment at the same cutoff
      go_enrichment_sig_fdr <- go_enrichment_results %>%
        dplyr::filter(!is.na(classicFisher_FDR_BH), classicFisher_FDR_BH < cutoff)
      go_enrichment_file_sig_fdr <- paste0(base, "_GO_enrichment_results_FDR", cutoff_label, "_enhanced.txt")
      write.table(go_enrichment_sig_fdr, go_enrichment_file_sig_fdr, sep = "\t", row.names = FALSE, quote = FALSE)

      # Back-map: genes contributing to each significant term (intersection with selected set)
      term_gene_rows <- lapply(sig_go_terms, function(go) {
        members <- genesInTerm(GOdata, go)[[1]]
        if (is.null(members)) members <- character(0)
        sig_members <- intersect(members, names(geneList)[geneList == 1])
        if (!length(sig_members)) return(NULL)
        data.frame(GO.ID = go, locusName = sig_members, stringsAsFactors = FALSE)
      })

      sig_term_genes <- dplyr::bind_rows(term_gene_rows)

      # Attach GO term label and gene annotations
      if (!is.null(sig_term_genes) && nrow(sig_term_genes) > 0) {
        sig_term_genes <- sig_term_genes %>%
          dplyr::left_join(go_enrichment_results %>% dplyr::select(GO.ID, Term), by = "GO.ID") %>%
          dplyr::left_join(functional_summary, by = "locusName") %>%
          dplyr::arrange(GO.ID, locusName)
      } else {
        sig_term_genes <- data.frame(
          GO.ID = character(0), Term = character(0), locusName = character(0),
          GO = character(0), KO = character(0), Best_hit_arabi_name = character(0),
          arabi_symbol = character(0), arabi_defline = character(0),
          stringsAsFactors = FALSE
        )
      }

      sig_genes_file <- paste0(base, "_significant_genes_per_GO_term_", cutoff_label, "_enhanced.txt")
      write.table(sig_term_genes, sig_genes_file, sep = "\t", row.names = FALSE, quote = FALSE)
      cat("    Mapped GO terms back to genes (", ont, ", ", cutoff_label, "): ", sig_genes_file, "\n", sep = "")

      # FDR-based GO->gene mapping at the same cutoff
      sig_go_terms_fdr <- names(pvals_fdr)[pvals_fdr < cutoff]
      term_gene_rows_fdr <- lapply(sig_go_terms_fdr, function(go) {
        members <- genesInTerm(GOdata, go)[[1]]
        if (is.null(members)) members <- character(0)
        sig_members <- intersect(members, names(geneList)[geneList == 1])
        if (!length(sig_members)) return(NULL)
        data.frame(GO.ID = go, locusName = sig_members, stringsAsFactors = FALSE)
      })
      sig_term_genes_fdr <- dplyr::bind_rows(term_gene_rows_fdr)
      if (!is.null(sig_term_genes_fdr) && nrow(sig_term_genes_fdr) > 0) {
        sig_term_genes_fdr <- sig_term_genes_fdr %>%
          dplyr::left_join(go_enrichment_results %>% dplyr::select(GO.ID, Term), by = "GO.ID") %>%
          dplyr::left_join(functional_summary, by = "locusName") %>%
          dplyr::arrange(GO.ID, locusName)
      } else {
        sig_term_genes_fdr <- data.frame(
          GO.ID = character(0), Term = character(0), locusName = character(0),
          GO = character(0), KO = character(0), Best_hit_arabi_name = character(0),
          arabi_symbol = character(0), arabi_defline = character(0),
          stringsAsFactors = FALSE
        )
      }
      sig_genes_file_fdr <- paste0(base, "_significant_genes_per_GO_term_FDR", cutoff_label, "_enhanced.txt")
      write.table(sig_term_genes_fdr, sig_genes_file_fdr, sep = "\t", row.names = FALSE, quote = FALSE)
      cat("    Mapped GO terms back to genes (", ont, ", FDR", cutoff_label, "): ", sig_genes_file_fdr, "\n", sep = "")
    }
  }

 
##############################
## Combined functional summary (optional)
##############################

treatments <- c("CB", "CP", "HB", "HG", "HH")

combined_list <- lapply(treatments, function(tr) {
  df <- fread(paste0(tr, "_functional_summary.txt"))
  df$treatment <- tr
  df
})

combined_df <- rbindlist(combined_list)

dup_counts <- combined_df %>%
  dplyr::group_by(locusName) %>%
  dplyr::summarise(
    treatment_count = dplyr::n_distinct(treatment),
    treatments = paste(unique(treatment), collapse = ", ")
  )

summary_table <- combined_df %>%
  dplyr::left_join(dup_counts, by = "locusName") %>%
  dplyr::arrange(dplyr::desc(treatment_count), locusName)

fwrite(summary_table, file = "Combined_functional_summary_across_treatments_enhanced.txt", sep = "\t")

cat("Summary saved as 'Combined_functional_summary_across_treatments_enhanced.txt'\n")
