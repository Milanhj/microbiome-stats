
# Functions and workflow for working with DESeq2

library(tidyverse)
library(DESeq2)

# Install DESeq2 from bioconductor
  #> install.packages("BiocManager")
  #> BiocManager::install("DESeq2")

# Load Data --------------------------------------------------------------------

# Simulated count data
counts <- read_rds("data/simulated_counts.rds")

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
                     sf_type, total_counts = NULL,
                     ordered = FALSE, reduced = NULL){
  
  # Rename id column id if it isn't already
  if(colnames(dat)[1] != "id") colnames(dat)[1] <- "id"
  
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
      
      # Create DESeq Data Object
      dds <- DESeqDataSetFromMatrix(
        countData = count_data,
        colData = col_data,
        # Full model
        design = formula # ~ group + timepoint + group:timepoint
      )
      
      # Store Model object with reduced model
      deseq_mod <- DESeq(dds, test = test, reduced = reduced) # ~group + timepoint
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
  #> fdr: if TRUE, q-values are calculated and returned
  #> na.rm: if TRUE, NA rows are removed
display_deseq_results <- function(mod, vars, var_label, digits = 4, 
                                  fdr = FALSE, na.rm = FALSE){
  # Matrix to store Values
  pvals <- matrix(
    nrow = length(rownames(mod)),
    ncol = 4,
    dimnames = list(NULL, c(var_label, "log2foldchange", "lfcSE", "p_value")))
  
  for (i in 1:length(rownames(mod))){
    row_name <- rownames(mod)[i]
    if (any(row_name == vars)){
      pvals[i, 1] <- row_name
      pvals[i, 2] <- round(mod$log2FoldChange[i], digits = digits)
      pvals[i, 3] <- round(mod$lfcSE[i], digits = digits)
      pvals[i, 4] <- round(mod$pvalue[i], digits = digits)
    } else { next }
  } # end loop
  
  if (na.rm){ # remove the NA levels
    # Filter NA
    pvals_tibble <- as_tibble(pvals) %>% 
      filter(!is.na(p_value))       
  } else{
    # Make pvals into tibble without filtering
    pvals_tibble <- as_tibble(pvals)
  } # end ifelse
  
  # Get variables into the format we want them in
  pvals_tibble <- pvals_tibble %>% 
    mutate(
      p_value = as.numeric(p_value),
      lfcSE = as.numeric(lfcSE),
      log2foldchange = as.numeric(log2foldchange),
      foldchange = round(2^log2foldchange, digits = digits)
    ) %>% 
    relocate(foldchange, .after = log2foldchange) %>%
    arrange(var_label)
  
  # Run an optional FDR
  if (fdr == TRUE) {
    pvals_tibble <- pvals_tibble %>% 
      mutate(padj = round(p.adjust(p_value, method = "fdr"), 
                          digits = digits))
  } 
  return(pvals_tibble)
} # end function



# fit_deseq2(): do_deseq and display_deseq_results wrapped together
# outputs a list of tables for each variable level (gene) with results for each independent variable 
# Use this function to output results for each variable in a design formula with multiple independent variables
fit_deseq2 <- function(dat, stop_col, formulas, alpha = 0.05, test = "Wald",
                       sf_type, total_counts = NULL,
                       ordered = FALSE, reduced = NULL, 
                       na.rm = FALSE, fdr = FALSE,
                       vars, var_label, digits = 4){
  
  formula_outs <- list()
  # Calculate the raw package ouput 
  for (i in 1:length(formulas)){
    # Get the raw DESeq2 output
    raw_out <- do_deseq(dat, stop_col = stop_col, sf_type = sf_type, 
                        total_counts = total_counts, 
                        formula = formulas[[i]][[1]], test = test) 
    
    # Organize raw output into a easier table
    formula_outs[[i]] <- display_deseq_results( 
      raw_out, vars = vars, var_label, digits, fdr = fdr, na.rm = na.rm
    )
    # Add the name of the variable to the list
    covs <- str_split_1(
      as.character(formulas[[i]]), pattern = "[\\s]*\\+[\\s]*"
    )
    names(formula_outs)[[i]] <- covs[[length(covs)]]
  } # end for i
  
  # List to hold tables with all results for each level
  list_results <- vector(length = length(vars), "list")
  names(list_results) <- vars
  
  for (i in 1:length(vars)){
    # matrix to store the data with variable names in the first column
    store_dat <-  cbind(var = names(formula_outs), 
                        matrix(nrow = length(formula_outs),
                               ncol = ncol(formula_outs[[1]]),
                               dimnames = list(NULL, names(formula_outs[[1]])))
                        
    ) # 5 cols
    for(j in 1:length(formula_outs)){
      # Fill matrix with data for a particular var level
      store_dat[j,2] <- formula_outs[[j]][i,1][[1]]
      store_dat[j,3] <- formula_outs[[j]][i,2][[1]]
      store_dat[j,4] <- formula_outs[[j]][i,3][[1]]
      store_dat[j,5] <- formula_outs[[j]][i,4][[1]] 
      store_dat[j,6] <- formula_outs[[j]][i,5][[1]] 
      
      # If we did an fdr, store it in a 7th column
      if(fdr){store_dat[j,7] <- formula_outs[[j]][i,6][[1]] }
      
    } # end for j
    # Fill slot for each level
    list_results[[i]] <- as_tibble(store_dat)
  } # end for i
  
  return(list_results)
} # end function



# Workflows --------------------------------------------------------------------



## Treatment Wald --------------------------------


### just tx -------------------

# Test differential expression between treatment arms
out_wald <- do_deseq(counts, stop_col = 3, alpha = 0.05, test = "Wald",
                     formula = ~ group, sf_type = "median_ratio")

# Visualize results for fdr = FALSE
result_wald <- display_deseq_results(
  out_wald, vars = count_col_names, var_label = "x", fdr = FALSE)

result_wald


# Visualize results for fdr = TRUE
result_wald_fdr <- display_deseq_results(
  out_wald, vars = count_col_names, var_label = "x", fdr = TRUE)

result_wald_fdr

### tx and time -------------------


# Define both formulas 
formulas <- list(c(~ timepoint + group), c(~ group + timepoint))

# Run wrapped function
result_group_time <- fit_deseq2(dat = counts, stop_col = 3, formulas = formulas, 
                                alpha = 0.05, test = "Wald", 
                                sf_type = "median_ratio",
                                na.rm = FALSE, fdr = TRUE,
                                vars = count_col_names, 
                                var_label = "x", digits = 4
           )
result_group_time

## Longitudinal LRT --------------------------------


# Test treatment_arm:timepoint
out_LRT <- do_deseq(counts, stop_col = 3, 
                    formula = ~ group + timepoint + group:timepoint,
                    reduced = ~ group + timepoint,
                    sf_type = "median-ratio",
                    alpha = 0.05, test = "LRT")

# Visualize results
result_LRT <- display_deseq_results(
  out_LRT, vars = count_col_names, var_label = "x")


result_LRT




