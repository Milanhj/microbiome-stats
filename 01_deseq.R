
# Functions and workflow for working with DESeq2

library(tidyverse)
library(DESeq2)

# Install DESeq2 from bioconductor
  #> install.packages("BiocManager")
  #> BiocManager::install("DESeq2")

# Load Data --------------------------------------------------------------------

# Simulated count data
counts <- read_rds("data/simulated_counts.rds")


# Functions --------------------------------------------------------------------

# do_deseq()
  #> Need ID variable in first spot
  #> df: Dataset with atleast IDs, treatment, and counts with all integer values
  #> tx_col: Treatment column index, must be directly before the count data starts
  #> alpha: significance cutoff used for optimizing the independent filtering
    #> Default = 0.05 here, but default is 0.1 for the package function
    #> Set to adjusted p-value cutoff
  #> test: (two options) Wald for cross-sectional, LRT for longitudinal
  #> time_col: Optional column index of timpoint data (must be between ID and treatment)

do_deseq <- function(df, tx_col, alpha = 0.05, test = "Wald", time_col = NULL){
  
  # Rename group column tx if it isn't already
  if(colnames(df)[tx_col] != "tx") colnames(df)[tx_col] <- "tx"
  
  # Rename id column id if it isn't already
  if(colnames(df)[1] != "id") colnames(df)[1] <- "id"

  if (test == "Wald"){
    
    # Transposed Count Data Frame 
    count_data <- t(
      df[,-c(2:tx_col)] %>% 
        column_to_rownames(var = colnames(df[,1]))
    ) 
    # Column data matrix
    col_data <- df[, c(1,tx_col)] %>% 
      column_to_rownames(var = colnames(df[,1]))
    
    # Make sure col_data and count_data match
    if (unique(rownames(col_data) != colnames(count_data))) break
    
    # Pull out tx variable for creating DESeq Data Object
    # tx <- col_data[,1]
    
    # Create DESeq Data Object
    dds <- DESeqDataSetFromMatrix(
      countData = count_data,
      colData = col_data,
      design = ~ tx
    )
    
    # Store Model object
    deseq_mod <- DESeq(dds, test = test)
    
    # Store Results
    mod_results <- results(deseq_mod, alpha = alpha)
    
  } else {
    
    if (test == "LRT"){
      
      # Rename time column timepoint if it isn't already
      if (colnames(df)[time_col] != "timepoint") colnames(df)[time_col] <- "timepoint"
     
      # If any IDs are duplicated, ID with unique row numbers
      if (any(duplicated(df$id) == TRUE)){
        # Replace old ID column with row numbers
        df <- df %>% 
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
        df[,-c(2:tx_col)] %>% 
          column_to_rownames(var = colnames(df[,1]))
      ) 
      
      # Column data matrix
      col_data <- df[, c(1:tx_col)] %>% 
        column_to_rownames(var = colnames(df[,1]))
      
      # Make sure col_data and count_data match
      if (unique(rownames(col_data) != colnames(count_data))) break
      
      # Pull out tx variable for creating DESeq Data Object
      # timepoint <- col_data[,1]
      # tx <- col_data[,2]
      
      # Create DESeq Data Object
      dds <- DESeqDataSetFromMatrix(
        countData = count_data,
        colData = col_data,
        # Full model
        design = ~ tx + timepoint + tx:timepoint
      )
      
      # Store Model object with reduced model
      deseq_mod <- DESeq(dds, test = test, reduced = ~tx + timepoint)
      
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
  
display_deseq_results <- function(mod, vars, var_label){
  
  # Matrix to store Values
  pvals <- matrix(nrow = length(rownames(mod)),
                  ncol = 3,
                  dimnames = list(NULL, c(var_label, "p_value", "padj")))
  
  for (i in 1:length(rownames(mod))){
    
    row_name <- rownames(mod)[i]
    
    if (any(row_name == vars)){
      pvals[i, 1] <- row_name
      pvals[i, 2] <- mod$pvalue[i]
      pvals[i, 3] <- signif(mod$padj[i], digits = 3)
    } else { next }
    
  } # end loop
  
  # Put Matrix into a more easily usable format
  pvals_tibble <- as_tibble(pvals) %>% 
    # Remove NA
    filter(!is.na(p_value)) %>% 
    # Add a column for padjust q-value
    mutate(
      p_value = as.numeric(p_value),
      padj = as.numeric(padj)
    ) %>% 
    arrange(p_value)
  
  return(pvals_tibble)
  
} # end function



# Workflows --------------------------------------------------------------------



## Treatment Wald --------------------------------


# Test differential expression between treatment arms
out_wald <- do_deseq(counts, tx_col = 3, alpha = 0.05, test = "Wald")

# Visualize results
result_wald <- display_deseq_results(
  out_wald, vars = colnames(counts), var_label = "x")



## Longitudinal LRT --------------------------------


# Test treatment_arm:timepoint
out_LRT <- do_deseq(counts, tx_col = 3, alpha = 0.05, test = "LRT", time_col = 2)

# Visualize results
result_LRT <- display_deseq_results(
  out_LRT, vars = colnames(counts), var_label = "x")


result_LRT



