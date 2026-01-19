# Brassica_rapa_pool_sequencing
# Brassica rapa Experimental Evolution – Downstream Genomic Data and Analysis Scripts

This repository contains the processed genomic-relevant input files, conda environment specifications, and complete analysis scripts used in the downstream population genomic analyses for the study:

**"Genomic responses to divergent temperature and pollinator selection in *Brassica rapa*"**  
(Ding & Schiestl, *New Phytologist*, accepeted)

Raw FASTQ sequence files are available through the European Nucleotide Archive (Project Number: PRJEB106602).

This github package includes:

- Effective population size (Ne) estimates and pool-size metadata  
- Conda environment file (`gatk_env.yml`) documenting full software versions  
- Scripts for population structure analyses, drift simulations, selection scans, and annotation  
- A reproducibility-oriented folder structure and workflow instructions  

---

## 1. Repository Structure

```
0_InputData/
    Ne.csv
    pool_sizes.txt
    gatk_env.yml

1_GeneticPopulationStructure/
    1_pi_allsites.slurm
    2_ΔAF.slurm
    3_withinpopDiversity_grenedalf.slurm.txt
    4_Heterozygosity.R
    5_PCA.R

2_DriftSelectionSimulations/
    1_drift_envelope.R
    2_simulate_neutral.py
    3_empirical_null.R

3_Selection/
    1_delta_AF/
        1_1_delta_af_common.R
        1_2_delta_af_classification.R
        1_3_delta_af_permutation.R
        1_4_△AF_classification_plots.R

    2_Fst/
        2_1_FilterFst_drift.R
        2_2_replicate_concordance.R
        2_3_flag_fst_exceeders.py

    3_ER_CMH_localscore/
        3_1_ER_cmh_modified.R
        3_2_ManhanttonPlot.R
        3_3_LocalScore_diagnostic.R
        3_4_localScore_auto_dual_xi.R
        3_5_localscore.slurm.txt
        scorelocalfunctions.R

4_Annotation/
    4_1_annotate_snps.sh
    4_2_annotation.R
    4_3_summarize_GO_results.R
    4_4_plot_GO_summary.R
```

---

## 2. Input Data Description

### **Ne.csv**
Effective population size estimates used for drift simulations and null expectation.

### **pool_sizes.txt**
Census sizes of each pool for expected drift variance computations.

### **gatk_env.yml**
Complete conda environment file listing software versions required to reproduce
SNP calling, preprocessing, and downstream analyses.

Recreate with:

```
conda env create -f gatk_env.yml
conda activate gatk_env
```
what's not in the folder but important for downstream analysis:
### **filtered_snps.vcf.gz**
Filtered SNP dataset generated from pool-seq reads using GATK HaplotypeCaller,
followed by depth/quality filtering in the publication. Used as input for all downstream analyses.

### **allele_frequencies_and_coverage.txt**
Per-SNP allele-frequency estimates and sequencing depth across all pools.
---

## 3. Analysis Modules

### **1_GeneticPopulationStructure/**
Scripts to generate allele-frequency matrices, perform PCA, and visualize population structure.

### **2_DriftSelectionSimulations/**
Scripts to simulate drift expectations (R/Python) and compare to observed data.

### **3_Selection/**
- **1_delta_AF** – allele-frequency shift analyses between ancestor and evolved pools  
- **2_Fst** – SNP-level and windowed FST differentiation  
- **3_ER_CMH_localscore** – CMH tests, ER statistics, and local score region detection  

### **4_Annotation/**
Scripts to map SNPs/regions to *Brassica rapa* genes and summarize functional categories.

---

## 4. Software Requirements

Recreate full environment using:

```
conda env create -f gatk_env.yml
conda activate gatk_env
```

Core tools included:

- GATK  
- bcftools  
- samtools  
- **R ≥ 4.0** (`tidyverse`, `data.table`, `ggplot2`, `poolfstat`, `hierfstat`, `pheatmap`)  
- **Python 3** (`pandas`, `numpy`, `scipy`)  
- bash / SLURM batch support  

---

## 5. Reproducibility

1. Start with input files in `0_InputData/`.  
2. Activate the conda environment using `gatk_env.yml`.  
3. Run modules in order: **1 → 2 → 3 → 4**.  
4. Each script indicates input/output paths and required dependencies.  
5. Outputs include PCA summaries, drift expectations, selection statistics, and annotated gene sets.  

---

## 6. Ethical and Legal Considerations

All data were generated from greenhouse-grown *Brassica rapa*.  
No human or sensitive data are included.  
Reuse permitted with citation of the associated publication.

---

## 7. Contact

**Yanqian Ding**  
yanqian.ding@systbot.uzh.ch
