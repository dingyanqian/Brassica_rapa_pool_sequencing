#####E&R-CMH script for pool-sequencing data
#############################
# Apply CMH test using ACER package
#############################

# Install and load ACER package if not already installed
#if (!require(remotes)) install.packages('remotes')
#remotes::install_github("MartaPelizzola/ACER") # if need from github
library(ACER)
library(data.table)

# Convert coverage_data to data.table for easier manipulation
coverage_data <- read.csv("allele_frequencies_and_coverage.txt")
# Read pool size information
pool_sizes <- read.table("../pool_sizes.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(pool_sizes) <- c("pool", "size")

# Extract pool information from column names
pool_info <- data.frame(
  pool = sub("^Coverage_", "", grep("^Coverage_", colnames(coverage_data), value = TRUE))
)

# Handle the special case for G1 and the other pools
# Create a custom function to parse pool names
parse_pool_info <- function(pool_name) {
  if (pool_name == "G1") {
    # Special case for G1
    return(list(
      treatment = "G1",
      replicate = 1,  # Numeric replicate
      generation = 1  # Generation 1
    ))
  } else {
    # For other pools (e.g., CBA, CBB, HBA, HBB, etc.)
    treatment <- substr(pool_name, 1, 2)  # First two characters (e.g., CB, HB)
    
    # Convert A/B to numeric replicates
    rep_char <- substr(pool_name, 3, 3)  # Third character (A or B)
    replicate <- ifelse(rep_char == "A", 1, 2)  # A = 1, B = 2
    
    # All other treatments are generation 6
    generation <- 6  # Both A and B are generation 6
    
    return(list(
      treatment = treatment,
      replicate = replicate,
      generation = generation
    ))
  }
}

# Apply the parsing function to each pool
for (i in 1:nrow(pool_info)) {
  pool_name <- pool_info$pool[i]
  info <- parse_pool_info(pool_name)
  
  pool_info$treatment[i] <- info$treatment
  pool_info$replicate[i] <- info$replicate
  pool_info$generation[i] <- info$generation
}

# Print pool information for verification
print(pool_info)

# Create a results dataframe
results <- data.frame(
  CHROM = coverage_data$CHROM,
  POS = coverage_data$POS,
  REF = coverage_data$REF,
  ALT = coverage_data$ALT
)

# Get unique treatments excluding G1 (which is a special case)
treatments <- unique(pool_info$treatment[pool_info$treatment != "G1"])

# Process each treatment separately
for (treat in treatments) {
  # Get pools for this treatment
  treat_pools <- pool_info$pool[pool_info$treatment == treat]
  
  # Skip if we don't have both replicates for this treatment
  if (length(treat_pools) < 2) {
    message("Skipping treatment ", treat, " - need at least 2 pools")
    next
  }
  
  # Convert columns to numeric to ensure proper calculations
  for (pool in c("G1", treat_pools)) {
    coverage_data[[paste0("AF_pool_", pool)]] <- as.numeric(coverage_data[[paste0("AF_pool_", pool)]])
    coverage_data[[paste0("Coverage_", pool)]] <- as.numeric(coverage_data[[paste0("Coverage_", pool)]])
  }
  
  # Remove rows with NA values in relevant columns
  cols_to_check <- c(
    paste0("AF_pool_G1"), 
    paste0("AF_pool_", treat_pools),
    paste0("Coverage_G1"), 
    paste0("Coverage_", treat_pools)
  )
  
  # Create a subset of data for this treatment with no NAs
  sub_data <- na.omit(coverage_data, cols = cols_to_check)
  
  # Filter out rows with zero coverage
  sub_data <- sub_data[
    get(paste0("Coverage_G1")) > 0 & 
      get(paste0("Coverage_", treat_pools[1])) > 0 & 
      get(paste0("Coverage_", treat_pools[2])) > 0
  ]
  
  # Skip if no data left after filtering
  if (nrow(sub_data) == 0) {
    message("No valid data for treatment ", treat, " after filtering")
    next
  }
  
  # Create allele frequency matrix
  # G1 data is duplicated for each replicate
  af_lst <- list(
    G1_rep1 = sub_data[[paste0("AF_pool_G1")]],
    G1_rep2 = sub_data[[paste0("AF_pool_G1")]],
    Gen6_rep1 = sub_data[[paste0("AF_pool_", treat_pools[1])]],
    Gen6_rep2 = sub_data[[paste0("AF_pool_", treat_pools[2])]]
  )
  af_mat <- as.matrix(do.call(cbind, af_lst))
  
  # Create coverage matrix
  cov_lst <- list(
    G1_rep1 = sub_data[[paste0("Coverage_G1")]],
    G1_rep2 = sub_data[[paste0("Coverage_G1")]],
    Gen6_rep1 = sub_data[[paste0("Coverage_", treat_pools[1])]],
    Gen6_rep2 = sub_data[[paste0("Coverage_", treat_pools[2])]]
  )
  cov_mat <- as.matrix(do.call(cbind, cov_lst))
  
  # Define generations and replicates vectors
  generations <- c(1, 1, 6, 6)  # G1, G1, Gen6, Gen6
  replicates <- c(1, 2, 1, 2)   # Rep1, Rep2, Rep1, Rep2
  
  # Assume Ne = 50 for all replicates (adjust as needed)
  Ne <- rep(50, length(unique(replicates)))
  
  # Read pool size information
  pool_sizes <- read.table("../pool_sizes.txt", header = FALSE, stringsAsFactors = FALSE)
  colnames(pool_sizes) <- c("pool", "size")
  
  # Create pool size vector for this treatment
  # Order must match the order of freq and coverage matrices: G1_rep1, G1_rep2, Gen6_rep1, Gen6_rep2
  pool_size_vec <- c(
    pool_sizes$size[pool_sizes$pool == "G1"],  # G1 for replicate 1
    pool_sizes$size[pool_sizes$pool == "G1"],  # G1 for replicate 2
    pool_sizes$size[pool_sizes$pool == treat_pools[1]],  # First replicate of treatment
    pool_sizes$size[pool_sizes$pool == treat_pools[2]]   # Second replicate of treatment
  )
  
  # Apply CMH test
  tryCatch({
    cmh_result <- ACER::adapted.cmh.test(
      freq = af_mat,
      coverage = cov_mat,
      Ne = Ne,
      gen = generations,
      repl = replicates,
      poolSize = pool_size_vec,
      RetVal = 2  # Request both test statistic and p-value in matrix form
    )
    
    # Print the structure to understand what we're getting
    print(str(cmh_result))
    
    # Create a results dataframe for this treatment
    # Handle different return formats based on RetVal parameter
    if (is.matrix(cmh_result) && ncol(cmh_result) == 2) {
      # RetVal=2 returns a matrix with test statistic in column 1 and p-value in column 2
      treatment_results <- data.frame(
        CHROM = sub_data$CHROM,
        POS = sub_data$POS,
        REF = sub_data$REF,
        ALT = sub_data$ALT,
        cmh_stat = cmh_result[,1],
        cmh_pval = cmh_result[,2]
      )
    } else if (is.vector(cmh_result)) {
      # RetVal=0 returns just p-values, RetVal=1 returns just test statistics
      # Check which one we have based on the parameter we used
      treatment_results <- data.frame(
        CHROM = sub_data$CHROM,
        POS = sub_data$POS,
        REF = sub_data$REF,
        ALT = sub_data$ALT
      )
      
      # Add the appropriate column based on what was returned
      if (RetVal == 0) {
        treatment_results$cmh_pval <- cmh_result
        # Calculate q-values
        treatment_results$cmh_qval <- p.adjust(treatment_results$cmh_pval, method = "BH")
      } else if (RetVal == 1) {
        treatment_results$cmh_stat <- cmh_result
      }
    } else {
      # Handle any other return format
      cat("Unexpected return format from adapted.cmh.test\n")
      print(cmh_result)
      
      # Create a basic results frame with just the SNP info
      treatment_results <- data.frame(
        CHROM = sub_data$CHROM,
        POS = sub_data$POS,
        REF = sub_data$REF,
        ALT = sub_data$ALT
      )
    }
    
    # Add FDR correction if we have p-values
    if ("cmh_pval" %in% colnames(treatment_results) && !exists("cmh_qval", treatment_results)) {
      treatment_results$cmh_qval <- p.adjust(treatment_results$cmh_pval, method = "BH")
    }
    
    # Save results for this treatment
    output_file <- paste0("cmh_results_", tolower(treat), ".txt")
    write.table(treatment_results, output_file, row.names = FALSE, quote = FALSE, sep = "\t")
    
    # Print summary
    cat("Treatment", treat, ":", nrow(treatment_results), "SNPs processed\n")
    if ("cmh_qval" %in% colnames(treatment_results)) {
      cat("Significant SNPs (q < 0.05):", sum(treatment_results$cmh_qval < 0.05, na.rm = TRUE), "\n")
    }
    
  }, error = function(e) {
    message("Error in CMH test for treatment ", treat, ": ", e$message)
  })
}
