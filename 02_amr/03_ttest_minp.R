
# T-test with min-p
library(tidyverse)
library(doMC)


# Set-up parallel processing
registerDoMC(cores = 12)

# Set seed
set.seed(3012)


# Load Data --------------------------------------------------------------------

counts <- read_rds("01_data/simulated_counts.rds")


## rM -------------------------------------


# Convert counts to reads per million
for (i in 5:ncol(counts)){
  counts[,i] <- (counts[[i]]/counts$total_reads)*1000000
}



# Functions --------------------------------------------------------------------


# start_col: the index of first column with counts
# group_col: grouping column (tx)
# total_reads_col: column with total non-human read numbers for log transforming
# time_col: time (longitudinal) or measurment month/year
  # if not a longitudinal study, use seperate datasets for pre.post intervention min-p calculation
# log2: should the data be log2+e transformed? log2(MLS + 1/total_reads*1000000)
# list: do you want the results to be returned in list format?
ttest_log2rm <- function(dat, group_col, tot_reads_col, start_col, time_col,
                         list = FALSE, log2 = TRUE){
  
  # Store variable (e.g. drug class names)
  rm_names <- colnames(dat)[start_col:length(colnames(dat))]
  
  # Change column names to generic names
  colnames(dat)[tot_reads_col] <- "total_reads"
  colnames(dat)[group_col] <- "group"
  colnames(dat)[time_col] <- "timepoint"
  
  if(list){ # Store in a list
    result <- list()
  }else{ # Store in a single table
    result <- tibble(timepoint = NA, var = NA, pval = NA)
  }
  # Loop for class
  for (i in 1:length(rm_names)){ # t-test for each time for class i
    
    # The index for class i
    dat_col <- which(names(dat) == rm_names[i]) # i
    
    # Subset data and log transform class i to create variable log2
    dat_log2 <- dat[, c(time_col, group_col, tot_reads_col, dat_col)] %>% 
      mutate(log2 = log2(!!sym(rm_names[i]) + 1/total_reads*1000000)) %>% 
      # Also give the original variable a generic name
      dplyr::rename(untrans = !!sym(rm_names[i]))
    
    # Store vector of times for class i
    times <- sort(unique(as.character(dat_log2$timepoint)))
    
    # Create a matrix to hold t-test results for each timepoint
    mttest <- matrix(ncol = 3, nrow = length(times),
                     dimnames = list(NULL, c("timepoint", "var", "pval")))
    # Loop for time
    for (j in 1:length(times)){
      # Filter for just the correct timepoint
      dat_tp <- dat_log2 %>% 
        filter(timepoint == times[j])
      # Store time and class name
      mttest[j,1] <- times[j]
      mttest[j,2] <- rm_names[i]
      # T-test with either the original variable or the log2 transformed
      if (log2){ # log2 
        mttest[j,3] <- t.test(log2 ~ group, data = dat_tp)$p.value
      } else{ # raw data
        mttest[j,3] <- t.test(untrans ~ group, data = dat_tp)$p.value
      }
    } # j
    
    # Store Results
    if(list){ # Store in a list
      result[[i]] <- as_tibble(mttest) %>% 
        mutate(
          timepoint = as.numeric(timepoint),
          pval = as.numeric(pval)
        )
    }else{ # Store in a single table
      # Bind results from class i to the overall result table
      result <- rbind(result,
                      as_tibble(mttest) %>% 
                        mutate(
                          timepoint = as.numeric(timepoint),
                          pval = as.numeric(pval)
                        ))
      # If it's the first iteration, take off the top row with NA used to initialize
      if(i == 1){result <- result[-1,]}
    }
  } # i
  return(result) 
} # end function


# Runs ttest_log2rm() and Min-p for multiple comparisons
# B: number of permutations
ttest_min_p <- function(dat, B = 999, group_col, tot_reads_col, 
                        start_col, time_col, 
                        list = FALSE, log2 = TRUE){
  
  # Make sure group column is named tx
  colnames(dat)[group_col] <- "tx"
  
  # Sort and store timepoints
  times <- sort(as.numeric(unique(as.vector(unlist(dat[,time_col])))))
  
  # Vector of class names
  rM_names <- colnames(dat)[c(start_col:ncol(dat))]
  
  # Call function to compute observed pvalues
  observed_result_table <- ttest_log2rm(
    dat = dat, group_col = group_col, 
    tot_reads_col = tot_reads_col, 
    start_col = start_col, time_col = time_col,
    list = list, log2 = log2
  )
  # Store observed P-values
  observed_p_values <- observed_result_table$pval
  # Extract observed minimum P-value
  observed_min_p <- min(observed_p_values, na.rm = TRUE)
  
  # Storage for permuted minimum P-values (numeric vector)
  min_perm_p_values <- numeric(B)
  
  for (b in 1:B) { # Permutation loop
    
    # Sample IDs without replacement
    perm_dat <- dat %>% 
      mutate(tx = sample(tx, length(tx), replace = FALSE))
    
    # Compute p-values from permuted data
    perm_p_values <- ttest_log2rm(
      dat = perm_dat, group_col = group_col, 
      tot_reads_col = tot_reads_col, 
      start_col = start_col, time_col = time_col,
      list = list, log2 = log2
    )
    # Store minimum P-value from vector of permuted p-values
    min_perm_p_values[b] <- min(perm_p_values$pval, na.rm = TRUE)
    
  } # for b
  
  # Tibble to hold adjusted p-values
  m_padj <- tibble(
    timepoint = rep(times, length(rM_names)),
    var = observed_result_table$var,
    minp_adj = rep(0)
  )
  # Adjust P-values
  for (i in 1:length(observed_p_values)) {
    # Proportion of permuted min p values less than or equal to the observed p-value i
    m_padj[i,3] <- mean(min_perm_p_values <= observed_p_values[i])
    
  } # for i
  
  # Join to raw p-value from observed results
  padj_result_table <- observed_result_table %>% 
    left_join(m_padj, by = c("var", "timepoint"))
  
  return(padj_result_table)
} # end function



# Workflow ----------------------------------------------------------------------


# Function call: want to run 10,000 permutations for actual analysis
result_min_p_ttest <- ttest_min_p(dat = counts, B = 50, group_col = 3, 
                                  tot_reads_col = 4, start_col = 5, 
                                  time_col = 2, list = FALSE, 
                                  log2 = TRUE
                                  )

result_min_p_ttest


# Save -------------------------------------------------------------------------

# If using many permutations, run in parallel in the background
# Save results to analyze in a different script
# save(result_min_p_ttest, file = "data/.../ttest_minp_result.rda")

