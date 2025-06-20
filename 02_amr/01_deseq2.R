
# Functions and workflow for working with DESeq2

library(tidyverse)
library(DESeq2)

# Install DESeq2 from bioconductor
  #> install.packages("BiocManager")
  #> BiocManager::install("DESeq2")

# Load Data --------------------------------------------------------------------

# Simulated count data
counts <- read_rds("01_data/simulated_counts.rds")

# Store count column names
count_col_names <- colnames(counts[,-c(1:3)])


# Functions --------------------------------------------------------------------

# Raw output from DESeq2 package 
#> dat: dataframe with counts and metadata for design formula
  #> contains atleast IDs, grouping variable, and counts with all integer values
  #> Need id variable in first column
#> formula: design formula with condition of interest last
#> stop_col: index for final column of metadata
#> alpha: significance cutoff used for optimizing the independent filtering
  #> default = 0.05 here, but default is 0.1 for the package function
  #> set to adjusted p-value cutoff
#> test: (two options) Wald for cross-sectional, LRT for longitudinal
#> sf_type: how to handle size factors
  #> sf_type = "median_ratio" uses default deseq2 esimation
  #> sf_type = "custom" uses custom size factors in sf argument
#> Optional arguments:
  #> total_counts: optional vector of total counts for custom size factor estimation
  #> reduced: for LRT, formula without the final interaction term
  #> ordered: is one of the covariates ordered? if yes, set to TRUE

do_deseq <- function(dat, stop_col, formula, alpha = 0.05, test = "Wald",
                     sf_type = "custom", total_counts = NULL,
                     ordered = FALSE, reduced = NULL){
  
  # Rename id column id if it isn't already
  if(colnames(dat)[1] != "id") colnames(dat)[1] <- "id"
  
  # If any IDs are duplicated, ID with unique row numbers
  if (any(duplicated(dat$id) == TRUE)){
    # Replace old ID column with row numbers
    dat <- dat %>% 
      # remove original id
      select(-id) %>% 
      # create new id variable with unique identifiers
      mutate(
        id = row_number(),        # create ID
        id = factor(id)           # make it a factor
      ) %>% 
      # make sure new id is the first column
      relocate(id)
  } # end id if
  
  if (test == "Wald"){
    
    # Transposed Count Data Frame 
    count_data <- t(
      dat[,-c(2:stop_col)] %>% 
        column_to_rownames(var = colnames(dat[,1]))
    ) 
    
    # Column data matrix
    col_data <- dat[, 1:stop_col] %>% 
      column_to_rownames(var = colnames(dat[,1]))
    
    # Make sure col_data and count_data match
    if (unique(rownames(col_data) != colnames(count_data))) break
    
    if (ordered){ 
      # Manually make the design matrix to include ordered factors
      design_matrix <- model.matrix(formula, data = col_data)
      
      # Create DESeq Data Object
      dds <- DESeqDataSetFromMatrix(
        countData = count_data,
        colData = col_data,
        design = design_matrix # group
      )
    } else{
      
      # Create DESeq Data Object
      dds <- DESeqDataSetFromMatrix(
        countData = count_data,
        colData = col_data,
        design = formula # group
      )
    }
    
    if (sf_type == "median_ratio"){
      
      # Store Model object
      deseq_mod <- DESeq(dds, test = test)
      
    } else{
      
      if(sf_type == "custom"){
        
        # Calculate the geometric mean of total counts
        geo_mn <- exp(mean(log(total_counts))) 
        # Custom size factors
        custom_sf <- total_counts/geo_mn
        # Manually set size factors
        sizeFactors(dds) <- custom_sf 
        
        # Store Model object using 
        deseq_mod <- DESeq(dds, 
                           fitType = "parametric",
                           test = test)
      } else{ 
        stop("invalid sf_type input") 
      } # ifelse 2
      
    } # if 1
    
    # Store Results
    mod_results <- results(deseq_mod, alpha = alpha)
    
  } else {
    
    if (test == "LRT"){
      
      # Transposed Count Data Frame 
      count_data <- t(
        dat[,-c(2:stop_col)] %>% 
          column_to_rownames(var = colnames(dat[,1]))
      ) 
      
      # Column data matrix
      col_data <- dat[, c(1:stop_col)] %>% 
        column_to_rownames(var = colnames(dat[,1]))
      
      # Make sure col_data and count_data match
      if (unique(rownames(col_data) != colnames(count_data))) break
      
      if (ordered){ 
        # Manually make the design matrix to include ordered factors
        design_matrix <- model.matrix(formula, data = col_data)
        reduced_matrix <- model.matrix(reduced, data = col_data)
        
        # Create DESeq Data Object
        dds <- DESeqDataSetFromMatrix(
          countData = count_data,
          colData = col_data,
          design = design_matrix # group
        )
      } else{
        
        # Create DESeq Data Object
        dds <- DESeqDataSetFromMatrix(
          countData = count_data,
          colData = col_data,
          design = formula # group
        )
      }
      # 
      # # Create DESeq Data Object
      # dds <- DESeqDataSetFromMatrix(
      #   countData = count_data,
      #   colData = col_data,
      #   # Full model
      #   design = formula # ~ group + timepoint + group:timepoint
      # )
      
      if (sf_type == "median_ratio"){
        
        if (ordered){
          # Store Model object with reduced model for ordinal var
          deseq_mod <- DESeq(dds, test = test, reduced = reduced_matrix) # ~group + timepoint
        } else{
          # Store Model object with reduced model
          deseq_mod <- DESeq(dds, test = test, reduced = reduced) # ~group + timepoint
        }
        
      } else{
        
        if(sf_type == "custom"){
          
          # Calculate the geometric mean of total counts
          geo_mn <- exp(mean(log(total_counts))) 
          # Custom size factors
          custom_sf <- total_counts/geo_mn
          # Manually set size factors
          sizeFactors(dds) <- custom_sf 
          
          if (ordered){
            # Store Model object using 
            deseq_mod <- DESeq(dds, fitType = "parametric",
                               test = test, reduced = reduced_matrix)
          } else{
            # Store Model object using 
            deseq_mod <- DESeq(dds, fitType = "parametric",
                               test = test, reduced = reduced)
          }
        } else{ 
          stop("invalid sf_type input") 
        } # ifelse 2
      }
      
      # Store Results
      mod_results <- results(deseq_mod, alpha = alpha)
      
    } else{ stop("Invalid test input") } # end else 2
    
  } # end else 1
  
  return(mod_results)
  
} # End function


# Display DESeq output in a cleaner, more easy to work with way
  #> mod: model output from my do_deseq() function or package results(DESeq())
  #> vars: the names of which variables to include in the output table
  #> var_label: what to name the column containing vars as rows
  #> digits: how many digits to round to
  #> p.adjust.method: whether or not to correct for multiple comparisons
    #> string for which p.adjust method to use
  #> na.rm: if TRUE, NA rows are removed

display_deseq_results <- function(mod, vars = NULL, var_label = NULL, digits = NULL,
                                  p.adjust.method = NULL, na.rm = FALSE){
  # Start with a generic variable name for levels
  default_label <- "outcome"
  
  # Matrix to store Values
  pvals <- matrix(
    nrow = length(rownames(mod)),
    ncol = 4,
    dimnames = list(NULL, c(default_label, "log2foldchange", "lfcSE", "p_value")))
  
  for (i in 1:length(rownames(mod))){
    # Pull out rowname from deseq output to store dat in matrix
    row_name <- rownames(mod)[i]
    # Store values
    pvals[i, 1] <- row_name
    pvals[i, 2] <- as.numeric(mod$log2FoldChange[i])
    pvals[i, 3] <- as.numeric(mod$lfcSE[i])
    pvals[i, 4] <- as.numeric(mod$pvalue[i])
    
  } # end for i
  
  # Deal with variable level display
  if (!is.null(vars)){ # select just levels included in vars
    
    # Create a temporary tibble for pulling out the variables you want
    # Extra steps so any length of vars can be used
    temp_tib <- as_tibble(pvals)
    
    # Optional padjust for multiple comparisons
    if (!is.null(p.adjust.method)) {
      temp_tib <- temp_tib %>% 
        mutate(p.adjust.method = p.adjust(p_value, method = p.adjust.method) 
        )
    }
    # Filter for just the variables listed in vars
    pvals_tibble <- tibble(
      outcome = temp_tib$outcome[temp_tib$outcome %in% vars]
    ) %>% 
      left_join(temp_tib, by = "outcome")
    
  } else{
    # Make pvals into tibble without filtering
    pvals_tibble <- as_tibble(pvals)
    
    # Optional padjust for multiple comparisons
    if (!is.null(p.adjust.method)) {
      pvals_tibble <- pvals_tibble %>% 
        mutate(p.adjust.method = p.adjust(p_value, method = p.adjust.method) 
        )
    }
  } # end ifelse
  
  # Get variables into the format we want them in
  pvals_tibble <- pvals_tibble %>% 
    mutate(
      p_value = as.numeric(p_value),
      lfcSE = as.numeric(lfcSE),
      log2foldchange = as.numeric(log2foldchange),
      foldchange = 2^log2foldchange
    ) %>% 
    relocate(foldchange, .after = log2foldchange) %>%
    arrange(default_label)
  
  # Rename outcome if a variable label is given
  if (!is.null(var_label)){
    colnames(pvals_tibble)[colnames(pvals_tibble) == "outcome"] <- var_label
  }
  # Round if digits provided
  if (!is.null(digits)){
    pvals_tibble <- pvals_tibble %>% 
      mutate(
        p_value = round(p_value, digits = digits),
        lfcSE = round(lfcSE, digits = digits),
        log2foldchange = round(log2foldchange, digits = digits),
        foldchange = round(foldchange, digits))
    # If padjust in there too, round that as well
    if(!is.null(p.adjust.method)){
      pvals_tibble <- pvals_tibble %>% 
        mutate(p.adjust.method = round(p.adjust.method, digits = digits))
    }
  } # end if
  
  return(pvals_tibble)
  
} # end function


# fit_deseq2(): do_deseq and display_deseq_results wrapped together
# outputs a list of tables for each variable level (gene) with results for each independent variable 
fit_deseq2 <- function(dat, stop_col, formulas, vars = NULL,
                       sf_type = "custom", total_counts = NULL, 
                       p.adjust.method = NULL, list = TRUE,
                       alpha = 0.05, test = "Wald",
                       ordered = FALSE, reduced = NULL, na.rm = FALSE,
                       var_label = NULL, digits = NULL){
  
  formula_outs <- list()
  # Calculate the raw package ouput 
  for (i in 1:length(formulas)){
    # Get the raw DESeq2 output
    raw_out <- do_deseq(dat, stop_col = stop_col, sf_type = sf_type, 
                        total_counts = total_counts, ordered = ordered,
                        formula = formulas[[i]][[1]], test = test) 
    
    # Organize raw output into a easier table
    out <- display_deseq_results(
      raw_out, vars = vars, var_label, digits = digits, 
      p.adjust.method = p.adjust.method, na.rm = na.rm
    )
    formula_outs[[i]] <- out
    # Add the name of the variable to the list
    if (length(formulas) == 1){
      # Pattern for splitting a formula with just 1 covariate
      covs <- str_split_i(
        as.character(formulas[[i]]), pattern = "~", 2
      )
    } else{
      # Pattern for splitting a formula with 2+ covariates
      covs <- str_split_1(
        as.character(formulas[[i]]), pattern = "[\\s]*\\+[\\s]*"
      )
    }
    
    names(formula_outs)[[i]] <- covs[[length(covs)]]
  } # end for i
  
  # List to hold tables with all results for each level
  list_results <- vector(length = nrow(out), "list")
  names(list_results) <- out[,1][[1]]
  
  for (i in 1:nrow(out)){
    # matrix to store the data with variable names in the first column
    store_dat <-  cbind(var = names(formula_outs), 
                        matrix(nrow = length(formula_outs),
                               ncol = ncol(formula_outs[[1]]),
                               dimnames = list(NULL, names(formula_outs[[1]])))
    ) # 5 cols
    for(j in 1:length(formula_outs)){
      # Fill matrix with data for a particular var level
      store_dat[j,2] <- formula_outs[[j]][i,1][[1]] # [i,1]
      store_dat[j,3] <- formula_outs[[j]][i,2][[1]] # [i,1]
      store_dat[j,4] <- formula_outs[[j]][i,3][[1]] # [i,1]
      store_dat[j,5] <- formula_outs[[j]][i,4][[1]] # [i,1]
      store_dat[j,6] <- formula_outs[[j]][i,5][[1]] # [i,1]
      
      if (!is.null(p.adjust.method)){
        store_dat[j,7] <- formula_outs[[j]][i,6][[1]]
      }
      
    } # end for j
    
    # Make outputs numeric
    final_tib <- as_tibble(store_dat) %>% 
      mutate(log2foldchange = as.numeric(log2foldchange),
             foldchange = as.numeric(foldchange),
             lfcSE = as.numeric(lfcSE),
             p_value = as.numeric(p_value)
      )
    # If padjust used, make that column numeric too
    if (!is.null(p.adjust.method)){
      final_tib <- final_tib %>% 
        mutate(p.adjust.method = as.numeric(p.adjust.method))
    } 
    
    # Fill slot for each level
    list_results[[i]] <- final_tib
  } # end for i
  
  # Return either a list or a full table
  if (list){
    return(list_results)
  } else{
    # Also deal with just single value in vars
    if(!is.null(vars) & length(vars) <= 1){
      return(list_results)
    } else{
      # Put into a single table
      tib <- list_results[[1]]
      for (j in 2:length(list_results)){
        tib <- rbind(tib, list_results[[j]])
      } # end j
      return(tib)
    }
  } # end else
  
  
} # end function


# For non-longitudinal data:
  # fit_deseq2() with a permutation loop for min-p
  # B: number of permutations for minP
  # tot_counts_col: column index with total counts
  # select_var: option to run min-p on a single drug class name
  # sep_time: set to true if you're doing parwise comparison at different timepoints
  # Function renames grouping column tx
  # need total counts col

deseq2_min_p <- function(dat, B = 999, stop_col, tot_counts_col, time_col = NULL,
                         formulas, alpha = 0.05, test = "Wald", group_col,
                         sf_type = "custom", sep_time = TRUE, select_var = NULL,
                         ordered = FALSE, reduced = NULL, digits = NULL){
  
  # If just one variable is getting min-p-ed, must be ran at sep times
  if (!is.null(select_var)) sep_time <- TRUE
  # Throw a warning message if time column is null
  if (!is.null(select_var) & is.null(time_col)){
    message("Must input a time column index when running min-P on a single variable.")
  }
  
  # Names of the classes
  rM_names <- colnames(dat)[(stop_col+1):ncol(dat)]
  # Rename group column (permutation column) "tx"
  colnames(dat)[group_col] <- "tx"
  # Deal with the outcome variable, give it a generic name
  var_label <- "outcome_var"
  
  if (sep_time){ # Run timepoints individually
    
    # Sort and store timepoints
    times <- sort(as.numeric(unique(as.vector(unlist(dat[,time_col])))))
    # Make sure time column is named timepoint
    colnames(dat)[time_col] <- "timepoint"
    
    # Initialize a tibble to store all results from DESeq2
    all_observed_result <- tibble(
      timepoint = NA, var = NA, outcome_var = NA, 
      log2foldchange = NA, foldchange = NA, lfcSE = NA, 
      p_value = NA
    )
    # Vector to store all observed minimun p-values from timepoints
    observed_pvals_all <- c()
    # Vector to store all minimun p-values from all permutations for all timepoints
    min_perm_pvals_all <- c()
    for (t in 1:length(times)){
      
      # Filter for that timepoint data
      dat_t <- dat %>% 
        filter(timepoint == times[t]) # t
      
      # Call function to compute observed p-values
      observed_result_table <- fit_deseq2(
        dat_t, stop_col = stop_col, formulas = formulas, alpha = alpha,
        test = test, sf_type = "custom", 
        total_counts = dat_t[[tot_counts_col]],
        var_label = var_label, digits = digits,
        p.adjust.method = NULL, list = FALSE, reduced = reduced,
        ordered = ordered,  na.rm = FALSE
      ) %>% 
        # Time variable
        mutate(timepoint = rep(times[t], length(rM_names))) %>% 
        relocate(timepoint)
      
      # Either filter for a single class or use results for all
      if (!is.null(select_var)){ # Single rM name
        
        # Filter results for just a single class
        limited_observed_result <- observed_result_table %>%
          filter(outcome_var == select_var)
        # Append observed p-values from time t to the combined vectoe
        observed_pvals_all[t] <- limited_observed_result$p_value
        # Bind observed result to combined results table
        all_observed_result <- rbind(all_observed_result, limited_observed_result)
        
      } else{ # All rM names
        
        # Append observed p-values from time t to the combined vectoe
        observed_pvals_all <- append(observed_pvals_all, observed_result_table$p_value)
        # Bind observed result to combined results table
        all_observed_result <- rbind(all_observed_result, observed_result_table)
      }
      
      # Remove the NA row used to create tibble
      if (t == 1){
        all_observed_result <- all_observed_result[-1,]
      } 
      
      # Storage for permuted minimum P-values for time t 
      min_perm_pvals_t <- c()
      for (b in 1:B) { # Permutation loop
        
        # Sample IDs without replacement
        perm_dat_t <- dat_t %>% 
          mutate(tx = sample(tx, length(tx), replace = FALSE))
        
        # Compute p-values from permuted data, same procedure as for observed
        perm_p_result <- fit_deseq2(
          perm_dat_t, stop_col = stop_col, formulas = formulas, alpha = alpha,
          test = test, sf_type = "custom", 
          total_counts = perm_dat_t[[tot_counts_col]],
          var_label = var_label, digits = digits,
          p.adjust.method = NULL, list = FALSE, reduced = reduced,
          ordered = ordered,  na.rm = FALSE
        )
        # Either filter for a single class or use results for all
        if (!is.null(select_var)){ # Single rM name
          
          # Select just a single class
          limited_perm_result <- perm_p_result %>%
            filter(outcome_var == select_var)
          # Store the p-value
          min_perm_pvals_t[b] <- limited_perm_result$p_value
          
        } else{ # All rM names
          # Record minimum P-value from vector of permuted p-values
          min_perm_pvals_t[b] <- min(perm_p_result$p_value, na.rm = TRUE)
        }
      } # for b
      
      # Append the permuted minPs to the vector will minP for all times
      min_perm_pvals_all <- append(min_perm_pvals_all, min_perm_pvals_t)
      
    } # for t
    
    # Build table for holding adjusted p-values
    if (!is.null(select_var)){ # For just a single variable
      m_padj <- tibble(
        timepoint = times,
        outcome_var = rep(select_var, length(times)),
        minp_adj = rep(0)
      )
    } else{ # For all variables
      
      # Build a vector of times for result table
      time_var <- c()
      for (p in 1:length(times)){
        time_var <- append(time_var, rep(times[p], length(rM_names)))
      } # for p
      # Build table for holding adjust P-values for all rM names
      m_padj <- tibble(
        timepoint = time_var,
        outcome_var = rep(rM_names, length(times)),
        minp_adj = rep(0)
      )
    }
    # Adjust p-values
    for (i in 1:length(observed_pvals_all)) {
      # Proportion of permuted min p values less than or equal to the observed p-value i
      m_padj[i,3] <- mean(min_perm_pvals_all <= observed_pvals_all[i])
      
    } # for i
    
    # Format results
    padj_result <- all_observed_result %>% 
      left_join(m_padj, by = c("outcome_var", "timepoint")) %>% 
      mutate(SE = 2^lfcSE) %>% 
      relocate(SE, .after = "lfcSE")
    
    return(padj_result)
    
  } else { # If don't need to calculate for separate timepoints
    
    # Call function to compute observed pvalues
    observed_result_table <- fit_deseq2(
      dat, stop_col = stop_col, formulas = formulas, alpha = alpha,
      test = test, sf_type = "custom", 
      total_counts = dat[[tot_counts_col]],
      var_label = var_label, digits = digits,
      p.adjust.method = NULL, list = FALSE, reduced = reduced,
      ordered = ordered,  na.rm = FALSE
    )
    
    # Store P-values
    observed_p_values <- observed_result_table$p_value
    
    # Storage for permuted minimum P-values (numeric vector)
    min_perm_p_values <- c()
    for (b in 1:B) { # Permutation loop
      
      # Sample IDs without replacement
      perm_dat <- dat %>% 
        mutate(tx = sample(tx, length(tx), replace = FALSE))
      
      # Compute p-values from permuted data, same procedure as for observed
      perm_p_result <- fit_deseq2(
        perm_dat, stop_col = stop_col, formulas = formulas, alpha = alpha,
        test = test, sf_type = "custom", 
        total_counts = perm_dat[[tot_counts_col]],
        var_label = var_label, digits = digits,
        p.adjust.method = NULL, list = FALSE, reduced = reduced,
        ordered = ordered,  na.rm = FALSE
      )
      
      # Record minimum P-value from vector of permuted p-values
      min_perm_p_values[b] <- min(perm_p_result$p_value, na.rm = TRUE)
    } # for b
    
    # Build table to store adjusted P-values
    m_padj <- tibble(
      var = observed_result_table$var,
      outcome_var = observed_result_table$outcome_var,
      minp_adj = rep(0))
    
    for (i in 1:length(observed_p_values)) {
      # Proportion of permuted min p values less than or equal to the observed p-value i
      m_padj[i,3] <- mean(min_perm_p_values <= observed_p_values[i])
      
    } # for i
    
    deseq_padj_res <- observed_result_table %>%
      left_join(m_padj, by = c("var", "outcome_var")) %>% 
      mutate(SE = 2^lfcSE) %>% 
      relocate(SE, .after = "lfcSE")
    
    return(deseq_padj_res)
    
  } # end else (sep_time == FALSE)
  
} # end function




# For longitudinal data:
  # fit_deseq2() with a permutation loop for westfall-young min-p
  # B: number of permutations for minP
  # stop_col: final column index of metadata (before counts start)
  # tot_counts_col: column index with total counts
  # group_col: column index of grouping variable
    # Function renames grouping column tx
  # select_vars: vector of variables to run through min-p procedure separately
    # All other columns will be done together
  # thresh: threshold for proportion of non-zero values needed for removing a count column

#dat <- counts
minp_deseq2_longitudinal <- function(
    dat, B = 999, formulas, stop_col, tot_counts_col, time_col,
    group_col, alpha = 0.05, test = "Wald", select_vars = NULL, thresh = NULL,
    ordered = FALSE, reduced = NULL, digits = NULL){
  
  # Convenience function to remove the too sparse columns
  remove_sparse_cols <- function(dat, start_col, thresh = 0.2, keep_vars = NULL){
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
  
  # Names of the genera, with selected variables removed
  rM_names <- colnames(dat)[(stop_col+1):ncol(dat)]
  rM_names <- rM_names[!(rM_names %in% select_vars)]
  
  # Sort and store timepoints
  times <- sort(as.numeric(unique(as.vector(unlist(dat[,time_col])))))
  
  # Make sure time column is named timepoint
  colnames(dat)[time_col] <- "timepoint"
  # Rename group column (permutation column) "tx"
  colnames(dat)[group_col] <- "tx"
  # Deal with the outcome variable, give it a generic name
  var_label <- "outcome_var"
  
  # Initialize a tibble to store all observed results from DESeq2
  # End in a table with a row for each timepoint
  all_observed_result <- tibble(
    timepoint = NA, var = NA, outcome_var = NA, 
    log2foldchange = NA, foldchange = NA, lfcSE = NA, 
    p_value = NA
  )
  
  ## 1) Observed
  
  # List to store the names of genera at each timepoint for building result table
  togther_group_names <- vector("list", length = length(times))
  # List to store timepoints for building result table
  togther_group_times <- vector("list", length = length(times))
  
  # List to store all observed minimun p-values from timepoints
  # Length of the number of select variables, plut one for the remaining genera 
  observed_pvals <- vector("list", length = length(select_vars)+1)
  names(observed_pvals) <- append(select_vars, "others")
  
  # Time loop for observed result
  for (t in 1:length(times)){
    
    # Filter for that timepoint data
    dat_t <- dat %>% 
      filter(timepoint == times[t]) # t
    
    # Remove sparse columns before running DESeq2
    if (!is.null(thresh)){
      dat_t <- remove_sparse_cols(
        dat_t, start_col = stop_col+1, thresh = thresh, 
        keep_vars = select_vars
      ) #select_vars
    }
    
    # Call function to compute observed p-values
    observed_result_table <- fit_deseq2(
      dat_t, stop_col = stop_col, formulas = formulas, alpha = alpha,
      test = test, sf_type = "custom", 
      total_counts = dat_t[[tot_counts_col]],
      var_label = var_label, digits = digits,
      p.adjust.method = NULL, list = FALSE, reduced = reduced,
      ordered = ordered,  na.rm = FALSE
    ) %>% 
      # Time variable
      mutate(timepoint = times[t]) %>% 
      #mutate(timepoint = rep(times[t], length(rM_names))) %>% 
      relocate(timepoint) #
    
    # Store names and times
    togther_group_names[[t]] <- 
      observed_result_table$outcome_var[!(observed_result_table$outcome_var
                                          %in% select_vars)]
    togther_group_times[[t]] <- rep(times[t], length(togther_group_names[[t]]))
    
    # Other genera: append observed p-values from time t to the combined vector
    observed_pvals[[length(observed_pvals)]] <- append(
      # Append result list in final position of observed p
      observed_pvals[[length(observed_pvals)]],
      # Extract p-values
      as.vector(unlist(
        observed_result_table %>% 
          filter(!(outcome_var %in% select_vars)) %>% 
          select(p_value)
      ))
    )
    
    if (!is.null(select_vars)){
      # Observed p-values for all select variables
      for (v in 1:length(select_vars)){
        observed_pvals[[v]][t] <- as.vector(unlist(
          observed_result_table %>% 
            filter(outcome_var == select_vars[v]) %>% #
            select(p_value)
        ))
      } # end for v
    }
    
    # Bind observed result to combined results table
    all_observed_result <- rbind(all_observed_result, observed_result_table)
    
    # remove the NA row used to create tibble
    if (t == 1){
      all_observed_result <- all_observed_result[-1,]
    } 
  } # for t observed
  
  ## 1) Permutation
  
  # List to store permuted minimum P-values for the groups of variables
  min_permutation_pvals <- vector("list", length = length(select_vars)+1)
  names(min_permutation_pvals) <- append(select_vars, "others")
  
  # Permutation loop
  for (b in 1:B){
    
    # List to store all p-values from permutation b for the groups of variables
    pval_b <- vector("list", length = length(select_vars)+1)
    names(pval_b) <- append(select_vars, "others")
    
    # Time loop for permutation
    for (tp in 1:length(times)){
      
      # The columns removed from the observed data at time t
      to_remove <- rM_names[!(rM_names %in% togther_group_names[[tp]])]
      
      # Filter for that timepoint data
      # Remove columns removed from that observed timepoint
      # Sample IDs without replacement
      perm_dat_t <- dat %>% 
        filter(timepoint == times[tp]) %>% 
        select(-all_of(to_remove)) %>% 
        mutate(tx = sample(tx, length(tx), replace = FALSE)) # t
      
      # Compute p-values from permuted data, same procedure as for observed
      perm_result_t <- fit_deseq2(
        perm_dat_t, stop_col = stop_col, formulas = formulas, alpha = alpha,
        test = test, sf_type = "custom", 
        total_counts = perm_dat_t[[tot_counts_col]],
        var_label = var_label, digits = digits,
        p.adjust.method = NULL, list = FALSE, reduced = reduced,
        ordered = ordered,  na.rm = FALSE
      ) 
      
      # Other genera: append permuted p-values from time t to the combined vector
      pval_b[[length(pval_b)]] <- append(
        # Append result list in final position of observed p
        pval_b[[length(pval_b)]],
        # Extract p-values
        as.vector(unlist(
          perm_result_t %>% 
            filter(!(outcome_var %in% select_vars)) %>% 
            select(p_value)
        ))
      )
      
      if(!is.null(select_vars)) {
        # Permuted p-values for all select variables
        for (vp in 1:length(select_vars)){
          pval_b[[vp]][tp] <- as.vector(unlist(
            perm_result_t %>% 
              filter(outcome_var == select_vars[vp]) %>% #
              select(p_value)
          ))
        } # end for vp
      }
      
    } # for tp 
    
    # Compute and store minimum p-values for each group in permutation b
    ## Each p-value represents results from an individual timepoint
    for(i in 1:length(pval_b)){
      min_permutation_pvals[[i]][b] <- min(pval_b[[i]], na.rm = TRUE)
    } # for i
    
  } # for b
  
  # Make a vector of names and times
  # Need to do select_vars first, replicate the name for each timepoint
  # names must remain in the same order as the indices of the observed pvals list (MLS, betalactam, others)
  var_names <- c()
  timepoints <- c()
  # Vector for variable to indicate if it was run independently
  group <- c()
  
  if(!is.null(select_vars)) {
    # Loop to index for the names run separately 
    for (s in 1:(length(observed_pvals)-1)){
      var_names <- append(var_names, rep(names(observed_pvals)[s], length(times)))
      timepoints <- append(timepoints, times)
      group <- append(group, rep(
        str_c("independent", as.character(s)), length(times))
      )
    } # for s
  }
  
  # Loop to index for the names run together
  for (i in 1:length(togther_group_names)){
    var_names <- append(var_names, togther_group_names[[i]])
    timepoints <- append(timepoints, togther_group_times[[i]])
    # Add on "together" for all genera in the group run together
    group <- append(group, rep("together", length(togther_group_names[[i]])))
  } # for i
  
  # Adjust p-values for each group
  # Vector for storing all adjusted p-values
  padj <- c()
  for (g in 1:length(observed_pvals)){ # for group g
    for (p in 1:length(observed_pvals[[g]])) { # for observed p-value p
      
      # Proportion of permuted min p values less than or equal to the observed p-value i
      padj <- append(padj,
                     mean(min_permutation_pvals[[g]] <= observed_pvals[[g]][p])
      )
    } # for p
  } # for g
  
  # Build table for adjust P-values result
  m_padj <- tibble(
    timepoint = timepoints,
    outcome_var = var_names,
    run_group = group,
    minp_adj = padj
  )
  
  # Join min-p adjusted p-values to the observed results
  padj_result <- all_observed_result %>% 
    left_join(m_padj, by = c("outcome_var", "timepoint")) %>% 
    select(-var) %>% 
    relocate(run_group, .after = "timepoint")
  
  return(padj_result)
  
  
} # end function



# Wald --------------------------------------------------------------


### custom offset --------------------------

# Test differential expression between treatment arms
out_wald_csf <- do_deseq(counts, stop_col = 4, alpha = 0.05, test = "Wald",
                         formula = ~ group, sf_type = "custom",
                         total_counts = counts$total_reads
                         )

# Visualize results for p.adjust.method = bonferroni-holm
result_wald_csf <- display_deseq_results(
  out_wald_csf, vars = count_col_names, var_label = "x",
  p.adjust.method = "holm") 

result_wald_csf



### median ratio ------------------


# Test differential expression between treatment arms
out_wald_mr <- do_deseq(counts, stop_col = 4, alpha = 0.05, test = "Wald",
                     formula = ~ group, sf_type = "median_ratio")

# Visualize results for p.adjust.method = NULL
result_wald_mr <- display_deseq_results(
  out_wald_mr, vars = count_col_names, var_label = "x",
  p.adjust.method = NULL)

result_wald_mr



# Visualize results for p.adjust.method 
result_wald_fdr <- display_deseq_results(
  out_wald_mr, vars = count_col_names, var_label = "x", p.adjust.method = "fdr")

result_wald_fdr


# For just a subset of variables
result_wald_subset <- display_deseq_results(
  out_wald_mr, vars = c("x1", "x5", "x18"), var_label = "x")


result_wald_subset




## fit_deseq2  -------------------------------------------------


# Wrapped function
# fit_deseq2 for a single covariate
out_fit_deseq2 <- fit_deseq2(counts, stop_col = 4, formulas = list(c(~group)),
                             sf_type = "custom", 
                             total_counts = counts$total_reads,
                             vars = count_col_names, var_label = "x",
                             p.adjust.method = "holm", list = FALSE,
                             alpha = 0.05, test = "Wald",
                             ordered = FALSE, reduced = NULL, na.rm = FALSE,
                             digits = NULL)

out_fit_deseq2




## tx and time ------------------------------------------------


# Define both formulas 
formulas <- list(c(~ timepoint + group), c(~ group + timepoint))



# Run wrapped function
result_group_time_csf <- fit_deseq2(dat = counts, stop_col = 4, formulas = formulas, 
                                    alpha = 0.05, test = "Wald", 
                                    sf_type = "custom",
                                    total_counts = counts$total_reads,
                                    na.rm = FALSE, p.adjust.method = "fdr",
                                    vars = count_col_names, 
                                    var_label = "x", digits = 4
)
result_group_time_csf




# LRT -------------------------------------------------------------



# Test treatment_arm:timepoint
out_LRT <- do_deseq(counts, stop_col = 4, 
                    formula = ~ group + timepoint + group:timepoint,
                    reduced = ~ group + timepoint,
                    sf_type = "custom",
                    total_counts = counts$total_reads,
                    alpha = 0.05, test = "LRT")

# Visualize results
result_LRT <- display_deseq_results(
  out_LRT, vars = count_col_names, var_label = "x")


result_LRT




# Min-P ------------------------------------------------------------------------


## Cross-sectional --------------------------------------

# For all variables
out_minp <- deseq2_min_p(counts, B = 5, formulas = list(c(~tx)),
                         stop_col = 4, group_col = 3, tot_counts_col = 4,
                         time_col = NULL, sf_type = "custom", alpha = 0.05,
                         test = "Wald", sep_time = FALSE)


# For a single x variable
out_minp_x1 <- deseq2_min_p(counts, B = 5, formulas = list(c(~tx)), 
                            select_var = "x1", stop_col = 4, group_col = 3,
                            tot_counts_col = 4, time_col = 2, 
                            sf_type = "custom", alpha = 0.05,
                            test = "Wald", sep_time = TRUE
                            )



## Longitudinal --------------------------------------


minp_deseq2_longitudinal(
  counts, B = 5, formulas = list(c(~tx)), select_vars = c("x1", "x5"), 
  stop_col = 4, group_col = 3, tot_counts_col = 4, time_col = 2, alpha = 0.05,
  test = "Wald", thresh = NULL, 
  ordered = FALSE, reduced = NULL, digits = NULL)




minp_deseq2_longitudinal(
  counts, B = 5, formulas = list(c(~tx)), select_vars = c("x1", "x5"), 
  stop_col = 4, group_col = 3, tot_counts_col = 4, time_col = 2, alpha = 0.05,
  test = "Wald", thresh = 0.1, 
  ordered = FALSE, reduced = NULL, digits = NULL)





