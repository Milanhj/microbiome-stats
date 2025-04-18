
# Test for correlation/association

library(tidyverse)



# Load Data --------------------------------------------------------------------

# Simulated count data
counts <- read_rds("01_data/simulated_counts.rds")

just_counts <- counts[,-c(1:3)]


# Functions --------------------------------------------------------------------


# Arguments:
  #> dat: dataframe with only count data
  #> method: which correlation coefficient to use (default = spearman)
  #> p.adjust.method: method for p.adjust correction for multiple comparisons
  #> nonzero_thres: threshold of non-zero values to use for filtering out columns with mostly zero counts
  #> var_lab: optional, what column names should be in result table (e.g "gene","class", etc)
    #> defaults to "var"
  #> digits: optional, how many digits to round final outputs to
corr_fdr <- function(dat, method = "spearman", p.adjust.method = "fdr",
                     nonzero_thres = round(0.1*nrow(dat)),
                     var_lab = NULL, digits = NULL){
  
  # if no var_lab provided, use a generic one
  if(is.null(var_lab)) var_lab <- "var"
  
  # First: Filter out columns with too many zeros
  # Vector to store names of columns to remove
  remove_cols <- c()
  # Convenience vector to track the non-zero values in each column
  for (i in 1:ncol(dat)){
    # Create a counter for non-zero values in a column
    non_zero <- 0
    for (j in 1:nrow(dat[,i])){
      if (dat[j,i] > 0){
        non_zero <- non_zero+1
      } else{next}
    } # end for j
    # Store column names with fewer than threshold for non-zero values
    if (non_zero <= nonzero_thres){
      remove_cols[i] <- colnames(dat)[i]
    } else{next}
  } # end for i
  # Take out NA values
  remove_cols <- remove_cols[!is.na(remove_cols)]
  # Remove all zero columns
  dat_clean <- dat %>%
    select(-all_of(remove_cols))
  
  # Second: get spearman correlation coefficients
  corr_coef <- as.data.frame(cor(dat_clean, method = method))
  # Get just the top triangle
  corr_coef_top <- matrix(ncol = ncol(corr_coef), nrow = nrow(corr_coef))
  for (i in 1:ncol(corr_coef)){
    # vector to store values before 1
    vals <- c()
    # counter
    j <- 0
    while (j+1 != i){ # stop loop at diagonal
      vals[j+1] <- corr_coef[j+1,i]
      j <- j + 1
    } # end while
    # Append on NA to create a vector the same length as nrow(corr_coef)
    corr_coef_top[,i] <- append(vals, rep(NA, length(corr_coef[,i]) - j))
  } # end for i
  # Rename rows and columns
  colnames(corr_coef_top) <- rownames(corr_coef_top) <- colnames(corr_coef)
  
  # Validate dimensions of correlation matrix
  if (any(dim(corr_coef) != c(ncol(dat_clean), ncol(dat_clean))) | 
      any(dim(corr_coef_top) != c(ncol(dat_clean), ncol(dat_clean)))){
    break
  } # end if
  
  # Third: correlation p value with cor.test()
  corr_pval <- matrix(nrow = ncol(dat_clean), ncol = ncol(dat_clean))
  for (i in 1:ncol(dat_clean)){
    for (j in 1:ncol(dat_clean)){
      
      corr_pval[i,j] <- cor.test(dat_clean[[i]], dat_clean[[j]],
                                 method = method, exact = FALSE)$p.value
    } # end j
  } # end i
  # Get just the top triangle
  corr_p_top <- matrix(ncol = ncol(corr_pval), 
                       nrow = nrow(corr_pval))
  for (i in 1:ncol(corr_pval)){
    # vector to store values before 1
    vals <- c()
    # counter
    j <- 0
    while (j+1 != i){ # stop loop at diagonal
      vals[j+1] <- corr_pval[j+1,i]
      j <- j + 1
    } # end while
    # append on NA
    corr_p_top[,i] <- append(vals, rep(NA, length(corr_pval[,i]) - j))
  } # end for i
  # Rename rows and columns
  colnames(corr_p_top) <- rownames(corr_p_top) <- colnames(corr_coef)
  # Validate dimensions of p matrix
  if (any(dim(corr_pval) != c(ncol(dat_clean), ncol(dat_clean))) | 
      any(dim(corr_p_top) != c(ncol(dat_clean), ncol(dat_clean)))){
    break
  } # end if
  
  # Fourth: FDR
  # Convert matrix values to vector
  vector_vals <- as.vector(corr_p_top)
  # NA positions
  na_positions <- is.na(vector_vals)
  # Remove NA values
  vector_p_no_na <- vector_vals[!na_positions]
  
  if (p.adjust.method == "fdr"){
    # FDR BH
    vector_padj <- p.adjust(vector_p_no_na, method = "fdr")
  } else{
    if(p.adjust.method == "bonferroni"){
      # Bonferroni correction
      vector_padj <- p.adjust(vector_p_no_na, method = "bonferroni")
    } else{
      if (p.adjust.method == "holm"){
        vector_padj <- p.adjust(vector_p_no_na, method = "holm")
      } else{
        stop("invalid input for p.adjust.method") 
      } # ifelse 3
    } # ifelse 2
  } # ifelse 1
  
  # Re-insert NA at original positions
  vector_vals[!na_positions] <- vector_padj
  # Coerce vector back into matrix
  corr_padjust_top <- matrix(vector_vals, nrow = nrow(corr_coef))
  # Add back in arg_group names
  colnames(corr_padjust_top) <- rownames(corr_padjust_top) <- colnames(corr_coef)
  
  # Validate dimensions of padj matrix
  if (any(dim(corr_padjust_top) != c(ncol(dat_clean), ncol(dat_clean)))) break
  
  # Fifth: Extract name pairs of arg_groups with padj values less than 0.05
  # Matrix to store 5 variables (gene names, cor coef, raw p-value, and q-value)
  pairs <- matrix(nrow = 1, ncol = 5, 
                  dimnames = list(
                    NULL, 
                    c(str_c(var_lab, "1"), str_c(var_lab, "2"), 
                      "cor", "pval", "padj")))
  
    # Skip NA
    for (i in 1:ncol(corr_padjust_top)){
      for(j in 1:nrow(corr_padjust_top)){
        # If it isnt NA
        if(!is.na(corr_padjust_top[j,i])){
          pairs <- rbind(
            pairs, c(rownames(corr_padjust_top)[i], 
                     colnames(corr_padjust_top)[j], 
                     # Add in the coefficient from the other matrix 
                     corr_coef[j,i],  # corr coef
                     # Raw p-value
                     corr_p_top[j,i],      # raw p
                     # padj for selected method
                     corr_padjust_top[j,i] # padj col
            )
          )
          # Convert to tibble and make variables numeric
          pairs <- as_tibble(pairs) %>% 
            mutate(
              cor = as.numeric(cor),
              pval = as.numeric(pval), 
              padj = as.numeric(padj)
            )
        } else{ next } # end ifelse
      } # end j
    } # end i
  
  # Remove first NA row 
  pairs <- pairs[-1,]
  
  if(!is.null(digits)){
    pairs <- as_tibble(pairs) %>% 
      mutate(
        cor = signif(as.numeric(cor), digits),
        pval = signif(as.numeric(pval), digits),
        padj = signif(as.numeric(padj), digits)
      )
  }
  
  # Return a list with the correlation, pvalue, fdr matrices, and a summary table
  return(list("correlation_coefficients" = corr_coef_top, 
              "pvalue_matrix" = corr_p_top,
              "qvalue_matrix" = corr_padjust_top, 
              "pairwise_comparison" = pairs,
              "removed_cols" = remove_cols)
  )
  
  
} # end function




# Workflow ---------------------------------------------------------------------


## Correlation and FDR -----------------------------


# Function call and store
out_fdr <- corr_fdr(dat = just_counts, method = "spearman", 
                p.adjust.method = "fdr") 

out_bf <- corr_fdr(dat = just_counts, method = "spearman", 
               p.adjust.method = "bonferroni") 

out_holm <- corr_fdr(dat = just_counts, method = "spearman", 
                 p.adjust.method = "holm" )



## Inspect -----------------------------------------


# Inspect Result table for each padjust method
result_fdr <- out_fdr[[4]] %>% 
  dplyr::rename(fdr = padj)

result_bf <- out_bf[[4]] %>% 
  dplyr::rename(bonferroni = padj)

result_holm <- out_holm[[4]] %>% 
  dplyr::rename(holm = padj)




