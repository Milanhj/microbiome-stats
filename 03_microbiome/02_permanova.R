
# PERMANOVA

library(tidyverse)
library(vegan)


# Load Data --------------------------------------------------------------------


# Simulated count data
counts <- read_rds("01_data/simulated_counts.rds") %>% 
  relocate(group, .after = "total_reads")


# Functions --------------------------------------------------------------------


# Combines runs PERMANOVA, comparing groups at a single timepoint
  #> dat: count data with time variable and group column preceding count columns
  #> start_col: first column with counts
  #> mthd: distance measure for `vegan::adonis2()`
  #> perm: number of permutations

fit_permanova <- function(dat, start_col, mthd, perm = 9999){
  
  # Filter data
  arm_dat <- dat[,c(1:start_col-1)]
  # Diversity data as matrix
  count_dat <- as.matrix(dat[,-c(1:start_col-1)])
  
  # Treatment arm
  arm <- as.matrix(arm_dat[,start_col-1])
  
  # PERMANOVA model
  permanova_out <- vegan::adonis2(
    count_dat ~ arm,
    #data = arm_dat,
    method = mthd,
    permutations = perm
  )
  return(permanova_out)
  
} # end function





# Workflow ---------------------------------------------------------------------



# Manhattan distance
l1_permanova <- fit_permanova(counts, start_col = 4, 
                              mthd = "manhattan", perm = 999
                              )


# Euclidean distance
l2_permanova <- fit_permanova(counts, start_col = 4, 
                              mthd = "euclidean", perm = 999
                              )

# Table Results
permanova_results <- data.frame(
  time1 = c(l1_permanova$`Pr(>F)`[1], l2_permanova$`Pr(>F)`[1]),
  row.names = c("Manhattan", "Euclidean")
  )


# Visualize results
permanova_results



# # For longitudinal sets: 

# # Vector of possible timepoints
# times <- unique(counts$timepoint)
# 
# 
# # For storing results 
# l1_permanova <- list()             # List for manhattan results
# l2_permanova <- list()             # List for euclidean results
# 
# for (i in 1:length(times)){
#   # Timepoint i
#   time_i <- times[i]
#   # Filter data for just time i
#   dat <- counts %>% 
#     filter(timepoint == time_i)
#   
#   # Manhattan permanova
#   l1_permanova[[i]] <- fit_permanova(dat, start_col = 4,
#                                      mthd = "manhattan", perm = 999)
#   # Euclidean permanova
#   l2_permanova[[i]] <- fit_permanova(dat, start_col = 4,
#                                      mthd = "euclidean", perm = 999)
# } # end loop











