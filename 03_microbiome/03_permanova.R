
# PERMANOVA

library(tidyverse)
library(vegan)


# Load Data --------------------------------------------------------------------


# Simulated count data
counts <- read_rds("01_data/simulated_counts.rds") %>% 
  relocate(group, .after = "total_reads")


# Functions --------------------------------------------------------------------


# Function from 02_pcoa returning Eigenvectors, Eigenvalues, and Centroids from PCoA
  #> dat: count data with time variable and group column preceding count columns
  #> start_col: first column with counts
  #> method: distance measure for vegan distance matrix (takes "l1" or "l2")
pcoa <- function(dat, start_col, method){
  
  # Sort by group
  dat <- arrange(dat, dat[,start_col-1][[1]])
  
  # Store just diversity variables as data
  data <- dat[, start_col:ncol(dat)]
  
  # Store ID variables
  arm <- dat[, 1:start_col-1]
  
  # Compute distance matrix using selected norm as method and store as dis
  if (method == "l1"){
    dist_m <- vegdist(data, method = "manhattan")
  } else if (method == "l2"){
    dist_m <- vegdist(data, method = "euclidean")
  } else {
    stop("Invalid method input: choose 'l1' or 'l2'")
  }
  
  # PCoA with distance matrix
  pcoa_mod <- betadisper(d = dist_m, group = arm[,start_col-1][[1]])
  
  # Site points in multivariate space
  vectors <- data.frame(
    group = pcoa_mod$group,
    data.frame(pcoa_mod$vectors)
  )
  
  # Name each PC
  names(vectors) <- "group"
  for (i in 2:ncol(vectors)){
    names(vectors)[i] <- str_c("v_PCoA", as.numeric(i-1))
  }
  
  # Centroids in multivariate space  
  centroids <- data.frame(
    group = rownames(pcoa_mod$centroids),
    data.frame(pcoa_mod$centroids)
  )
  
  # Eigenvalues
  eigs <- as.data.frame(as.numeric(pcoa_mod$eig))
  names(eigs) <- "value"
  
  out <- list(centroids, vectors, eigs)
  names(out) <- c("centroids", "vectors", "eig")
  
  # Return output eigenvalues and pc1/pc2 eigenvectors 
  return(out)
  
} # end function

# Running PERMANOVA at all timepoints in a dataset
  #> count_start_col: index of first column with counts
  #> time_col: index of column containing timepoints
  #> group_col: column index of the grouping variable
  #> method: distance measure to use
    #> "l1" for manhattan
    #> "l2" for euclidean
  #> pcoa: logical for if permanova should be done after a pcoa or on the observed data
  #> B: number of permutations for permanova
  #> thresh_variance: threshold for determining how many PCs to by percent variance explained

pcoa_permanova <- function(dat, count_start_col, group_col, time_col,
                           run_pcoa = TRUE, B = 9999,
                           method, thresh_variance = 0.8){
  
  # `Vegan` required, load vegan
  library(vegan)
  
  # Sort and store timepoints
  times <- sort(as.numeric(unique(as.vector(unlist(dat[,time_col])))))
  # rename column of group and time variable
  colnames(dat)[group_col] <- "group"
  colnames(dat)[time_col] <- "timepoint"
  
  # 1) PCoA
  
  # list to hold PCoA outputs for each timepoint
  #pcoa_out <- list()
  
  # Tibble to store PERMANOVA results for each timepoint
  permanova_results <- tibble(
    timepoint = NA,
    p_value = NA
  )
  
  # Vector to store the number of PCs explaining the threshold of variance
  PC_thresh_var <- c()
  
  for (t in 1:length(times)){
    print(t)
    # Filter data
    dat_t <- dat %>% 
      filter(timepoint == times[t]) %>% 
      relocate(group, .before = colnames(dat)[count_start_col])
    
    if (run_pcoa) {
      
      # Run PCoA
      pcoa_out <- pcoa(dat_t, count_start_col, method = method)
      
      # Sum of the eigenvalues
      total_var <- sum(pcoa_out$eig[,1])
      # percent variance
      perc_var <- pcoa_out$eig[,1]/total_var
      # Starting place for while loop (percent variance of first PC)
      for_sum_var <- perc_var[1]
      
      # Counter for while loop starting at position 2
      counter <- 2
      while(for_sum_var <= thresh_variance){
        for_sum_var <- sum(for_sum_var, perc_var[counter])
        counter <- counter + 1
      } 
      # Add the count for number of PCs it takes to explain the threshold of variance
      PC_thresh_var[t] <- counter
      
      # Principal Components explaining the threshold for time t
      dat_counts <- as.matrix(pcoa_out$vectors[,2:(counter+1)])
      # Treatment arms (first column in the pcoa vector output)
      arms <- data.frame(group = pcoa_out$vectors[,1]) 
      
      # Rename method for permanova function
      if (method == "l1") distance_measure <- "manhattan"
      if (method == "l2") distance_measure <- "euclidean"
      
      # Calculate permanova results
      permanova_out <- vegan::adonis2(
        dat_counts ~ group, 
        data = arms,
        method = distance_measure,
        permutations = B
      )
      # Format and store results in overall tibble
      permanova_results <- rbind(
        permanova_results,
        tibble(timepoint = times[t],
               p_value = permanova_out$`Pr(>F)`[1])
      )
      
    } else{
      
      # Rename method for permanova function
      if (method == "l1") distance_measure <- "manhattan"
      if (method == "l2") distance_measure <- "euclidean"
      
      # Data for time t
      dat_counts <- dat_t[,-c(1:(count_start_col-1))]
      # Arms
      arms <- dat_t[,c(1:4)]
      
      # Calculate permanova results
      permanova_out <- vegan::adonis2(
        dat_counts ~ group, 
        data = arms,
        method = distance_measure,
        permutations = B
      )
      # Format and store results in overall tibble
      permanova_results <- rbind(
        permanova_results,
        tibble(timepoint = times[t],
               p_value = permanova_out$`Pr(>F)`[1])
      )
      
    } # ifelse pcoa
    
    # Remove first row of NA if t=1
    if (t == 1){
      permanova_results <- permanova_results[-1,]
    }
    
  } # for t
  
  # Return the number of PCs used if pcoa, otherwise return just the permanova results
  if (run_pcoa){
    return(list(PC_thresh_var, permanova_results))
  } else{
    return(permanova_results)
  }
  
} # end function


# Workflows ---------------------------------------------------------------------



## PCoA --------------------------------------


# Manhattan distance
l1_pcoa_permanova_result <- pcoa_permanova(
  counts, count_start_col = 5, group_col = 4, time_col = 2, run_pcoa = TRUE,
  B = 9999, method = "l1", thresh_variance = 0.8)


# Euclidean distance
l2_pcoa_permanova_result <- pcoa_permanova(
  counts, count_start_col = 5, group_col = 4, time_col = 2, run_pcoa = TRUE,
  B = 9999, method = "l2", thresh_variance = 0.8)




# Visualize results
l1_pcoa_permanova_result
l2_pcoa_permanova_result



## Observed Data --------------------------------------


# Manhattan distance
l1_permanova_result <- pcoa_permanova(
  counts, count_start_col = 5, group_col = 4, time_col = 2, run_pcoa = FALSE,
  B = 9999, method = "l1")


# Euclidean distance
l2_permanova_result <- pcoa_permanova(
  counts, count_start_col = 5, group_col = 4, time_col = 2, run_pcoa = FALSE,
  B = 9999, method = "l2")




# Visualize results
l1_permanova_result
l2_permanova_result




