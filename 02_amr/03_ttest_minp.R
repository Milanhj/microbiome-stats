
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

# dat: data with a group variable, total reads variable, and timepoint variable
# group_col: column with grouping variable
# tot_reads_col: column index with total reads number for that sample
# start_col: column index where count data starts
# time_col: column index with timepoint variable
# list: logical, list == TRUE returns results in list format, FALSE for single tibble
# log2: if TRUE, log2 transforms data for the t-test

ttest_log2rm <- function(dat, group_col, tot_reads_col, start_col, time_col,
                         list = FALSE, log2 = TRUE){
  
  get_lm_results <- function(dat, formula) {
    # Model Output
    res <- summary(lm(formula, data = dat))
    
    # Tibble to store lm output for all covariates
    store_result <- tibble(
      var = NA,
      estimate = NA,
      SE = NA, 
      p_value = NA)
    
    # Organize lm output for all covariates, remove intercept
    for (i in 2:nrow(res$coefficients)){
      store_result <- rbind(
        store_result,
        tibble(
          var = rownames(res$coefficients)[i],
          estimate = res$coefficients[i, "Estimate"],
          SE = res$coefficients[i,"Std. Error"],
          p_value = res$coefficients[i,"Pr(>|t|)"]
        )
      )
      # Remove the first row (with the NAs)
      if (i == 2)  store_result <- store_result[-1,]
    } # for i
    return(store_result)
  } # end function
  
  rm_names <- colnames(dat)[start_col:length(colnames(dat))]
  colnames(dat)[tot_reads_col] <- "total_reads"
  colnames(dat)[group_col] <- "tx"
  colnames(dat)[time_col] <- "timepoint"
  
  if(list){ # Store in a list
    result <- list()
  }else{ # Store in a single table
    result <- tibble(
      timepoint = NA, 
      name = NA,
      var = NA, 
      estimate = NA,
      SE = NA,
      p_value = NA
    )
  }
  for (i in 1:length(rm_names)){ # t-test for each time for class i
    
    # The index for class i
    dat_col <- which(names(dat) == rm_names[i]) # i
    
    # Subset data and log transform class i to create variable log2
    dat_log2 <- dat[, c(time_col, group_col, tot_reads_col, dat_col)] %>% 
      mutate(log2 = log2(!!sym(rm_names[i]) + 1/total_reads*1000000)) %>% 
      # Also give the original variable a generic name
      dplyr::rename(untrans = !!sym(rm_names[i]))
    
    # Store vector of times for class i
    times <- sort(unique(dat_log2$timepoint))
    
    # Create a tibble to hold lm t-test results for each timepoint
    lm_tib <- tibble(
      timepoint = NA, 
      name = NA,
      var = NA, 
      estimate = NA,
      SE = NA,
      p_value = NA
    )
    for (t in 1:length(times)){
      # Filter for just the correct timepoint
      dat_tp <- dat_log2 %>% 
        filter(timepoint == times[t])
      
      # T-test with either the original variable or the log2 transformed
      if (log2){
        lm_tib <- rbind(
          lm_tib, 
          as_tibble(cbind(
            tibble(
              timepoint = times[t],
              name = rm_names[i]),
            get_lm_results(dat_tp, as.formula(log2 ~ tx))
          ))
        )
      } else{
        
        lm_tib <- rbind(
          lm_tib, 
          as_tibble(cbind(
            tibble(
              timepoint = times[t],
              name = rm_names[i]),
            get_lm_results(dat_tp, as.formula(untrans ~ tx))
          ))
        )
        
      } # if else
      if (t == 1) lm_tib <- lm_tib[-1,]
    } # t
    
    if(list){ 
      # Store in a list
      result[[i]] <- lm_tib
    }else{ # Store in a single table
      # Bind results from class i to the overall result table
      result <- rbind(result,lm_tib)
      # If it's the first iteration, take off the top row with NA used to initialize
      if(i == 1){result <- result[-1,]
      }
    }
    
  } # i
  if(log2){
    colnames(result)[4] <- "log2estimate"
    colnames(result)[5] <- "log2SE"
  }
  return(result) 
} # end function




# Runs ttest_log2rm() and Min-p for multiple comparisons
# B: number of permutations
ttest_min_p <- function(dat, B = 999, group_col, tot_reads_col, 
                        start_col, time_col,
                        list = FALSE, log2 = TRUE){
  
  # Make sure group column is named tx
  colnames(dat)[group_col] <- "tx"
  
  # Vector of class names
  rM_names <- colnames(dat)[c(start_col:ncol(dat))]
  
  # Call function to compute observed pvalues
  observed_result <- ttest_log2rm(
    dat = dat, group_col = group_col, 
    tot_reads_col = tot_reads_col, 
    start_col = start_col, time_col = time_col,
    list = list, log2 = log2
  )
  observed_p_values <- observed_result$p_value
  
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
    )$p_value
    
    # Record minimum P-value from vector of permuted p-values
    min_perm_p_values[b] <- min(perm_p_values, na.rm = TRUE)
    
  } # for b
  
  # Adjust P-values
  m_padj <- tibble(
    timepoint = observed_result$timepoint,
    name = observed_result$name,
    minp_adj = rep(0)
  )
  for (i in 1:length(observed_p_values)) {
    # Proportion of permuted min p values less than or equal to the observed p-value i
    m_padj[i,3] <- mean(min_perm_p_values <= observed_p_values[i])
    
  } # for i
  
  minp_result <- observed_result %>% 
    left_join(m_padj, by = c("name", "timepoint"))
  
  # Return results
  return(minp_result)
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

