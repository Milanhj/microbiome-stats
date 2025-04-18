
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
# Use this function to output results for each variable in a design formula with multiple independent variables
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
    covs <- str_split_1(
      as.character(formulas[[i]]), pattern = "[\\s]*\\+[\\s]*"
    )
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





# Treatment Wald --------------------------------------------------------------


## just tx  -------------------------------------------------


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




## tx and time ------------------------------------------------


# Define both formulas 
formulas <- list(c(~ timepoint + group), c(~ group + timepoint))



### custom offset -----------------------


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



### median ratio ------------------------


# Run wrapped function
result_group_time_mr <- fit_deseq2(dat = counts, stop_col = 4, formulas = formulas, 
                                alpha = 0.05, test = "Wald", 
                                sf_type = "median_ratio",
                                na.rm = FALSE, p.adjust.method = "fdr",
                                vars = count_col_names, 
                                var_label = "x", digits = 4
           )
result_group_time_mr





# Longitudinal LRT -------------------------------------------------------------



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




