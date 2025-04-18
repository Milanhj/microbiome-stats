
# Calculate alpha diversity 

library(tidyverse)
library(patchwork)

# library(vegan) can return same output with vegan::diversity(vector, index)


# Load Data --------------------------------------------------------------------


# Simulated count data
counts <- read_rds("01_data/simulated_counts.rds") %>% 
  relocate(group, .after = "total_reads")


# Functions --------------------------------------------------------------------


# Shannon: Calculate L1 Diversity for a vector of counts or normalized reads
shannon <- function(v){
  
  prob <- c()
  for (i in 1:length(v)){
    
    # Handle zeros
    if (v[i] != 0){
      # Genus proportion of set
      prob[i] <- v[i]/sum(v)
    } else{ next }
    
  } # end for
  
  # Remove NA
  prob <- prob[!is.na(prob)]
  
  if (is.null(prob)) {
    # Return 1 for shannon diversity if all values are 0 exp(0) = 1
    div <- 1
  } else {
    
    # calculate diversity
    div <- exp(-sum(prob*log(prob)))
    
  } # end else
  
  return(div)
  
} # end function



# Calculate L2 Diversity for a vector of counts or normalized reads
# Returns a vector with [1] simpson's index, and [2] inverse simpson 
simpson <- function(v){
  
  prob <- c()
  for (i in 1:length(v)){
    
    # Deal with 0 sum
    if (sum(v) == 0){
      prob[i] <- 0
    } else{
      
      # Genus proportion of set
      prob[i] <- v[i]/sum(v)
      
    } # end else
    
  } # end for loop
  
  # Simpson's index 
  si <- 1 - sum(prob^2)
  # Calculate inverse Simpson's index
  inv_si <- (sum(prob^2))^(-1)
  
  return(c(si, inv_si))
  
}



# Calculate all diversity indices simultaneously
  #> Output table with treatment arm, shannon, simpson's index, and inverse simpson
  #> dat: Count data with treatment arm in column immediately before counts start
  #> start_col: Index of first count column 

calc_diversity <- function(dat, start_col){
  
  # Store just diversity variables as data
  data <- as.data.frame(dat[, start_col:ncol(dat)])
  
  # Store ID variables
  arm <- as.data.frame(dat[, start_col-1])
  
  # Create empty matrix to store values
  div_out <- matrix(nrow = nrow(arm), ncol = 3,
                    dimnames = list(NULL, c("shannon", "simpson", "inv_simpson"))
  )
  for (i in 1:nrow(arm)){
    
    # Pool data
    vi <- as.numeric(data[i,])
    # Calculate Shannon Diversity for row i
    l1 <- shannon(vi)
    l2 <- simpson(vi)
    
    # Store diversity
    div_out[i,1] <- l1    # Shannon
    div_out[i,2] <- l2[1] # SI
    div_out[i,3] <- l2[2] # inverse SI
    
  } # end for
  
  # Fill row i in 1st column with treatment arm
  div_out <- as_tibble(cbind(arm, as.data.frame(div_out)))
  
  return(div_out)
  
} # end function




# Plot diversity values from calc_diversity() output, requires `patchwork` loaded
  #> dat: Output table from calc_diversity()
  #> xlim: Manually set limits of x-axis (if you don't do this, values get truncated)
  #> timepoint: Character string for plot title (e.g. "baseline", "36 months"..)
  #> groups: Character vector for naming groups in the legend
  #> colors: Vector of colors matching groups

plot_diversity <- function(dat, xlim, timepoint, groups, colors){
  
  # Rename group column tx if it isn't already
  if(colnames(dat)[1] != "tx") colnames(dat)[1] <- "tx"
  
  p_shannon <- ggplot(data = dat) +
    geom_density(
      aes(shannon, fill = tx), 
      alpha = 0.5,
      color = FALSE
    ) + 
    scale_fill_manual(
      labels = groups,
      values = colors
    ) +
    scale_x_continuous(limits = xlim) +
    labs(
      title = str_c("Diversity", timepoint, sep = " "),
      x = "Shannon's Diversity (effective number)",
      y = "Density"
    ) +
    theme_light() +
    theme(
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
      # aspect.ratio = 1
    )
  
  p_inv_si <- ggplot(data = dat) +
    geom_density(
      aes(inv_simpson, fill = tx), 
      alpha = 0.5,
      color = FALSE
    ) + 
    scale_fill_manual(
      labels = groups, 
      values = colors
    ) +
    scale_x_continuous(limits = xlim) +
    labs(
      x = "Inverse Simpson's Diversity (effective number)",
      y = "Density"
    ) +
    theme_light() +
    theme(
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
      # aspect.ratio = 1
    )
  
  p_combined <- p_shannon + p_inv_si + plot_layout(guides = "collect")
  
  return(p_combined)
  
} # end function



# Workflow ---------------------------------------------------------------------

# Calculate diversity for each ID
div_out <- calc_diversity(counts, start_col = 5)



# Now that we have vector data reduced to a scalar, can run a t-test or anova

# Shannon's Diversity test
t_shannon <- t.test(shannon ~ group, data = div_out)

# Simpson's Diversity test
t_inv_simpson <- t.test(inv_simpson ~ group, data = div_out)



# Plots -------------------------------------------------------------------------


# Plot Density distribution of diversity 
plot_diversity(div_out, xlim = c(10, 18), 
               timepoint = "",                 # not including time in analysis
               groups = c("A", "B"), 
               colors = c("#5D3A9B", "#E66100")
               )


