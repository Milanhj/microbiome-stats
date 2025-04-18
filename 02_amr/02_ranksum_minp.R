
# Wilcoxon Rank Sum with Min-p


library(tidyverse)
library(doMC)


# Set-up parallel processing
registerDoMC(cores = 12)

set.seed(3012)



# Load Data --------------------------------------------------------------------

counts <- read_rds("01_data/simulated_counts.rds")


## rM -------------------------------------


# Convert counts to reads per million
for (i in 5:ncol(counts)){
  counts[,i] <- (counts[[i]]/counts$total_reads)*1000000
}


# Functions --------------------------------------------------------------------


# dat: dataset with just rM/counts and grouping variable in 1st column
rank_sum <- function(dat){
  
  # Arrange by grouping variable
  dat <- dat %>% arrange(!!sym(colnames(dat)[1]))
  
  # Create a matrix to store results with rows = the number of rM columns
  res_matrix <- matrix(ncol = 2, nrow = ncol(dat)-1,
                       dimnames = list(NULL, c("var", "p_value")))
  
  for (i in 1:nrow(res_matrix)){
    # Store data with tx group (1) and class column i (i+1)
    dat_col <- dat[,c(1,i+1)] 
    
    # Store class name
    res_matrix[i,1] <- colnames(dat_col)[2]
    
    # Skip any columns with all 0s
    if (all(dat_col[,2] == 0)) {
      next
    } 
    # Rank sum test and store in p-value column (2): class_i_counts ~ tx
    res_matrix[i,2] <-
      wilcox.test(dat_col[,2][[1]] ~ dat_col[,1][[1]],
                  data = dat_col)$p.value
  } # for i
  
  # Make output into a tibble
  res_matrix <- as_tibble(res_matrix) %>% 
    mutate(p_value = if_else(!is.nan(p_value), as.numeric(p_value), NA))
  # Return result matrix
  return(res_matrix)
  
} # end function

# dat: dataset with all meta data before count columns 
# start_col: the index of first column with counts
# group_col: grouping column (tx)
# total_reads_col: column with total non-human read numbers for log transforming
# time_col: time (longitudinal) or measurment month/year
  # if not a longitudinal study, use seperate datasets for pre.post intervention min-p calculation
rank_sum_longitudinal <- function(dat, time_col, group_col, start_col){
  
  # Sort and store timepoints
  times <- sort(as.numeric(unique(as.vector(unlist(dat[,time_col])))))
  
  # Names of the classes
  rM_names <- colnames(dat)[start_col:ncol(dat)]
  # Name the grouping column 
  colnames(dat)[group_col] <- "group"
  # Name the time column
  colnames(dat)[time_col] <- "timepoint"
  
  # Tibble for storing results for all time points
  all_result <- tibble(timepoint = NA, var = NA, p_value = NA)
  
  for (t in 1:length(times)){
    # Filter for time t, remove all variables except counts and group column
    # Arrange by grouping variable
    dat_t <- dat %>% 
      filter(timepoint == times[t]) %>% 
      select(group, all_of(rM_names)) %>% 
      arrange(group)
    
    # Create a matrix to store results with rows = the number of rM columns
    res_matrix <- matrix(
      ncol = 3, nrow = length(rM_names), 
      dimnames = list(NULL, c("timepoint", "var", "p_value"))
    )
    
    for (i in 1:length(rM_names)){
      # Store data with tx group (1) and class column i (i+1)
      dat_col <- dat_t[,c(1,i+1)] 
      
      # Store timepoint
      res_matrix[i,1] <- times[t]
      # Store class name
      res_matrix[i,2] <- colnames(dat_col)[2]
      
      # Skip any columns with all 0s
      if (all(dat_col[,2] == 0)) {
        next
      } 
      # Rank sum test and store in p-value column (2): class_i_counts ~ tx
      res_matrix[i,3] <-
        wilcox.test(dat_col[,2][[1]] ~ dat_col[,1][[1]],
                    data = dat_col)$p.value
    } # for i
    
    # Bind to the storing tibble
    all_result <- all_result %>% 
      rbind( # Make output into a tibble
        as_tibble(res_matrix) %>% 
          mutate(
            timepoint = as.numeric(timepoint),
            p_value = if_else(!is.nan(p_value), as.numeric(p_value), NA))
      )
    
    # If first iteration of the loop, remove NA row used to initialize the tibble
    if (t == 1){
      all_result <- all_result[-1,]
    }
    
  } # for t
  
  # Return result matrix
  return(all_result)
  
} # end function


# B: number of permutations
rank_sum_minp <- function(dat, time_col, group_col, start_col, B = 9999){
  
  # Sort and store timepoints
  times <- sort(as.numeric(unique(as.vector(unlist(dat[,time_col])))))
  # Names of the classes
  rM_names <- colnames(dat)[start_col:ncol(dat)]
  # Make sure group column is named tx
  colnames(dat)[group_col] <- "tx"
  
  # Call function to compute observed pvalues
  observed_result_table <- rank_sum_longitudinal(dat, time_col = time_col, 
                                                 group_col = group_col,
                                                 start_col = start_col)
  # Store p-values
  observed_p_values <- observed_result_table$p_value
  
  # Extract observed minimum P-value
  observed_min_p <- min(observed_p_values, na.rm = TRUE)
  
  # Storage for permuted minimum P-values (numeric vector)
  min_perm_p_values <- numeric(B)
  
  #min_perm_p_values <- tibble(nrow = B, ncol = 2)
  
  for (b in 1:B) { # Permutation loop
    
    # Sample IDs without replacement
    perm_dat <- dat %>% 
      mutate(tx = sample(tx, length(tx), replace = FALSE))
    
    # Compute p-values from permuted data
    perm_p_result <- rank_sum_longitudinal(perm_dat, time_col = time_col,
                                           group_col = group_col, 
                                           start_col = start_col
    )
    # Record minimum P-value from vector of permuted p-values
    min_perm_p_values[b] <- min(perm_p_result$p_value, na.rm = TRUE)
    
  } # for b
  
  # Use this look to create the correct time labels
  # The rank sum function does all classes for each time individually, so need 18 of each timepoint sequentially
  time_var <- c()
  # Build a vector of times
  for (t in 1:length(times)){
    time_var <- append(time_var, rep(times[t], length(rM_names)))
  }
  
  # Adjust P-values
  m_padj <- tibble(
    timepoint = time_var,
    var = rep(rM_names, length(times)),
    minp_adj = rep(0)
  )
  for (i in 1:length(observed_p_values)) {
    # Proportion of permuted min p values less than or equal to the observed p-value i
    m_padj[i,3] <- mean(min_perm_p_values <= observed_p_values[i])
    
  } # for i
  
  padj_result_table <- observed_result_table %>% 
    left_join(m_padj, by = c("var", "timepoint"))
  
  # Return results
  return(padj_result_table)
  
} # end function


# Calculate --------------------------------------------------------------------


# Function call: want to run 10,000 permutations for actual analysis
rsum_minp_result <- rank_sum_minp(counts, B = 50, time_col = 2, 
                                   group_col = 3, start_col = 5)



# Save -------------------------------------------------------------------------

# # If using many permutations, run in parallel in the background
# # Save results to analyze in a different script
# save(result_min_p_ranksum, file = "data/.../ranksum_minp_result.rda")


