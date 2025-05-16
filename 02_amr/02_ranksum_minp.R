
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
# select_vars: option to min-p select classes separately, 
  # Vector of class (outcome_var) names to analyze independent of the others
# thresh: the threshold for proportion of non-zero values needed to be included in analysis
rank_sum_minp <- function(dat, time_col, group_col, select_vars = NULL,
                           start_col, thresh = NULL, B = 999){
  
  # Convenience function to remove the too sparse columns
  remove_sparse_cols <- function(dat, start_col, keep_vars, thresh = 0.2){
    remove_indices <- c()
    for (i in start_col:ncol(dat)){
      # Add an NA in position i if the column isnt numeric
      if (!is.numeric(dat[,i])) remove_indices[i] <- NA
      if (!is.null(keep_vars)){
        if(colnames(dat)[i] %in% keep_vars) {
          remove_indices[i] <- NA
        } else{
          # If proportion non-zero rows is below thresh, add its index to remove it
          if (sum(dat[[i]] != 0) / length(dat[[i]]) < thresh){
            remove_indices[i] <- i
          }else{
            remove_indices[i] <- NA
          }
        }
      } else{
        # If proportion non-zero rows is below thresh, add its index to remove it
        if (sum(dat[[i]] != 0) / length(dat[[i]]) < thresh){
          remove_indices[i] <- i
        }else{
          remove_indices[i] <- NA
        }
      }
    } # for i
    # Remove NA for just the column indices
    remove_indices <- remove_indices[!is.na(remove_indices)]
    # Remove columns with too many zero values
    if (length(remove_indices) == 0){
      return(dat)
    }else{
      dat <- dat[,-c(remove_indices)]
      return(dat)
    }
  } # end function
  
  # Clean data by removing columns with too few reads
  if (!is.null(thresh)){
    dat <- remove_sparse_cols(
      dat, start_col = start_col, thresh = thresh, 
      keep_vars = select_vars
    ) #select_vars
  }
  
  # Sort and store timepoints
  times <- sort(as.numeric(unique(as.vector(unlist(dat[,time_col])))))
  # Names of the classes with classes from select_vars removed
  rM_names <- colnames(dat)[start_col:ncol(dat)]
  rM_names <- rM_names[!(rM_names %in% select_vars)]
  # Make sure group column is named tx
  colnames(dat)[group_col] <- "tx"
  
  # List to store all observed minimun p-values from timepoints
  # Length of the number of select variables, plut one for the remaining classes 
  observed_pvals <- vector("list", length = length(select_vars)+1)
  names(observed_pvals) <- append(select_vars, "others")
  
  # 1) Observed
  # Call function to compute observed pvalues
  observed_result_table <- rank_sum_longitudinal(dat, time_col = time_col, 
                                                 group_col = group_col,
                                                 start_col = start_col)
  # Other classes: append observed p-values from time t to the combined vector
  observed_pvals[[length(observed_pvals)]] <- 
    append( # Append result list in final position of observed p
      observed_pvals[[length(observed_pvals)]],  
      # Extract p-values
      as.vector(unlist(
        observed_result_table %>% 
          filter(!(var %in% select_vars)) %>% 
          select(p_value)
      ))
    )
  if (!is.null(select_vars)){
    # Observed p-values for all select variables
    for (v in 1:length(select_vars)){
      observed_pvals[[v]] <- as.vector(unlist(
        observed_result_table %>% 
          filter(var == select_vars[v]) %>% #
          select(p_value)
      ))
    } # end for v
  }
  
  # Storage for permuted minimum P-values (numeric vector)
  min_perm_p_values <-  vector("list", length = length(select_vars)+1)
  names(min_perm_p_values) <- append(select_vars, "others")
  
  for (b in 1:B) { # Permutation loop
    
    # Sample IDs without replacement
    perm_dat <- dat %>% 
      mutate(tx = sample(tx, length(tx), replace = FALSE))
    
    # Compute p-values from permuted data
    perm_p_result <- rank_sum_longitudinal(perm_dat, time_col = time_col,
                                           group_col = group_col, 
                                           start_col = start_col
    )
    # Store minimum p-value for all of rM_names
    min_perm_p_values[[length(min_perm_p_values)]][b] <- 
      min(
        as.vector(unlist(
          perm_p_result %>% 
            filter(!(var %in% select_vars)) %>% 
            select(p_value)
        )),
        na.rm = TRUE
      )
    if (!is.null(select_vars)){
      # Minimum p-values for all select variables
      for (vp in 1:length(select_vars)){
        min_perm_p_values[[vp]][b] <- 
          min(
            as.vector(unlist(
              perm_p_result %>% 
                filter(var == select_vars[vp]) %>%
                select(p_value)
            )),
            na.rm = TRUE
          )
      } # end for v
    }
    
  } # for b
  
  # Create a vector for storing a vector of times/names to build results table
  time_var <- c()
  names_var <- c()
  group_var <- c()
  
  if (!is.null(select_vars)){
    for (v in 1:length(select_vars)){
      # Add in the times for each 
      time_var <- append(time_var, rep(times))
      names_var <- append(names_var, rep(select_vars[v], length(times)))
      group_var <- append(group_var, 
                          rep(str_c("independent", as.character(v)), length(times)))
    }
  }
  # Times, names, and run group for the other classes
  time_var <- append(time_var, sort(rep(times, length(rM_names))))
  names_var <- append(names_var, rep(rM_names, length(times)))
  group_var <- append(group_var, rep("together", length(rM_names)*length(times)))
  
  # Adjust p-values for each group g
  # Vector for storing all adjusted p-values
  padj <- c()
  for (g in 1:length(observed_pvals)){ # for group g
    #print(g)
    for (p in 1:length(observed_pvals[[g]])) { # for observed p-value p
      # Proportion of permuted min p values less than or equal to the observed p-value i
      padj <- append(padj,
                     mean(min_perm_p_values[[g]] <= observed_pvals[[g]][p])
      )
    } # for p
  } # for g
  
  # Table for adjusted P-values
  m_padj <- tibble(
    timepoint = time_var,
    var = names_var,
    group = group_var,
    minp_adj = padj
  )
  # Join to observed results
  padj_result_table <- observed_result_table %>% 
    left_join(m_padj, by = c("var", "timepoint")) %>% 
    dplyr::rename(outcome_var = var) %>% 
    relocate(group, .after = "outcome_var")
  
  # Return results
  return(padj_result_table)
  
} # end function


# Calculate --------------------------------------------------------------------


# Function call: want to run 10,000 permutations for actual analysis
rsum_minp_result <- rank_sum_minp(counts, B = 50, time_col = 2, thresh = 0.1,
                                   group_col = 3, start_col = 5)



# Save -------------------------------------------------------------------------

# # If using many permutations, run in parallel in the background
# # Save results to analyze in a different script
# save(result_min_p_ranksum, file = "data/.../ranksum_minp_result.rda")


