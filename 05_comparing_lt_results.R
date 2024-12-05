
# Direct comparison of results from multiple time points

library(tidyverse)
library(vegan)
library(plotly)


# Table Function ---------------------------------------------------------------


# Combines p-value outputs for different timepoints in a single table
  #> list_results: fit_permanova() output for each timepoint collected in a list
  #> var_col: the column of the result variable to be represented
  #> bind: if no key variable for joining exists, use bind = TRUE to bind columns instead of joining

combine_table <- function(list_results, var_col, var_name, bind = FALSE){
  
  # Loop to rename and store just variable of interest
  new_list <- list()
  for (i in 1:length(list_results)){
    
    tib_i <- list_results[[i]][, c(1, var_col)]
    # Change the name of the variable to reflect timepoint
    colnames(tib_i)[2] <- str_c(var_name, as.character(i-1), sep = "_")
    # Store in the new list
    new_list[[i]] <- tib_i
    
  } # end for
  
  if (bind == FALSE){ # Left Joining when a key variable exists
    
    # Starting place for creating the full table
    tib_full <- new_list[[1]]
    # Store the name of the key variable to join tibbles with
    join_var <- colnames(tib_full)[1]
    
    for(j in 2:length(new_list)){
      
      tib_j <- new_list[[j]]
      
      # Join columns
      tib_full <- tib_full %>% 
        left_join(tib_j, by = join_var)
      
    } # end for
    
  } else{ # If there is no key variable to join, then bind 
    
    # Starting place for creating the full table
    tib_full <- new_list[[1]]
    
    for(j in 2:length(new_list)){
      
      tib_j <- new_list[[j]][,2]
      
      # Bind columns
      tib_full <- as_tibble(tib_full %>% cbind(tib_j))
    } # end for
    
  } # end else
  
  return(tib_full)
  
} # end function



# Plot 3D Function -------------------------------------------------------------



# Plot a 3D scatterplot for visualizing time centroids on 3 principal components
  #> centroids: centroid coordinate output from `pcoa()`
  #> colors: vector of colors matching number of timepoints
  #> labels: text to label reach centroid with timepoint
  #> title_text: add an optional title

plotly_scatterplot <- function(centroids, colors, labels, title_text = NULL){
  
  # use plotly to make a 3d scatterplot
  p <- plot_ly(centroids, #centroids_l2, 
               x = ~PCoA1, y = ~PCoA2, z = ~PCoA3,
               colors = colors
  ) %>% 
    # Add in lines to connect points (needs to go in before points?)
    add_trace(type = "scatter3d",
              mode = "lines",
              line = list(color = "grey", width = 2),
              showlegend = FALSE
    ) %>% 
    # Add timepoints 
    add_trace(type = "scatter3d",
              mode = "markers",
              color = ~ group, 
              marker = list(size = 4, symbol = "arrow-open"),
              colors = colors,
              showlegend = FALSE
    ) %>% 
    # Label timepoints
    add_trace(color = ~ group,
              type = "scatter3d", 
              mode = "text",
              name = labels,
              text = labels,
              showlegend = FALSE
    ) %>% 
    # Add title
    layout(title = title_text)
  
  return(p)
  
} # end function




