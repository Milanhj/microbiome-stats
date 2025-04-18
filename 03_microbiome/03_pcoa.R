
# PCoA: Principal Coordinate Analysis

library(tidyverse)
library(vegan)
library(patchwork)


# Load Data --------------------------------------------------------------------


# Simulated count data
counts <- read_rds("data/simulated_counts.rds")



# Functions --------------------------------------------------------------------


# Function returning Eigenvectors, Eigenvalues, and Centroids from PCoA
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




# Plot Centroids from PCoA1 and PCoA2
  #> vectors: Eigenvectors from pcoa() output
  #> centroids: Centroids from pcoa() output
  #> label: labels for groups in the legend
  #> colors: option to manually set group colors
  #> title_text: optional plot title
  #> subtitle_text: optional plot subtitle

plot_pcoa <- function(vectors, centroids, label, colors = NULL,
                      title_text = NULL, 
                      subtitle_text = NULL){ 
  
  # Eigenvectors and Centroids for PC1 and PC2 
  v1_2 <- as_tibble(vectors[,1:3])
  c1_2 <- as_tibble(centroids[,1:3]) %>% 
    select(group, c_PCoA1 = PCoA1, c_PCoA2 = PCoA2)
  
  # Values to plot lines from segments to points in tibble format
  for_seg <- v1_2 %>% 
    group_by(group) %>% 
    left_join(c1_2, by = "group") %>% 
    # add in a variable labeling row numbers by group
    mutate(point = row_number())
  
  # PCoA 1 and 2
  plot_1_2 <- 
    ggplot() + 
    stat_ellipse(
      data = for_seg,
      aes(x = v_PCoA1, y = v_PCoA2, fill = group),
      geom = "polygon",
      alpha = 0.175
    ) +
    geom_segment(
      data = for_seg,
      aes(x = v_PCoA1, xend = c_PCoA1,
          y = v_PCoA2, yend = c_PCoA2,
          color = group),
      alpha = 0.30
    ) +
    # centroids
    geom_point(
      data = for_seg, # group and first two PC
      aes(x = c_PCoA1, y = c_PCoA2, color = group),
      size = 5,
      # shape = 18
    ) +
    # points
    geom_point(
      data = for_seg,
      aes(x = v_PCoA1, y = v_PCoA2, color = group), 
      size = 2
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "PCoA1",
      y = "PCoA2"
    ) +
    theme_light() +
    theme(
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      aspect.ratio = 1)
  
  if (is.null(colors)){
    
    return(plot_1_2)
    
  } else {
    
    plot_1_2 <- plot_1_2 +     
      scale_color_manual(
        labels = label,
        values = colors
      ) +
      scale_fill_manual(
        labels = label,
        values = colors
      )
    return(plot_1_2)
    
  } # end else
  
} # end function



# Workflow ---------------------------------------------------------------------


# Calculate centroids and eigenvectors
l1_pcoa <- pcoa(counts, start_col = 4, "l1")
l2_pcoa <- pcoa(counts, start_col = 4, "l2")


# Plot centroids L1
p_l1 <- plot_pcoa(vectors = l1_pcoa$vectors, 
          l1_pcoa$centroids,
          label = c("A", "B"),
          colors = c("#5D3A9B", "#E66100"),
          title_text = "Manhattan Distance",
          subtitle_text =  NULL
          )


# Plot centroids L2
p_l2 <- plot_pcoa(vectors = l2_pcoa$vectors, 
                l2_pcoa$centroids,
                label = c("A", "B"),
                colors = c("#5D3A9B", "#E66100"),
                title_text = "Euclidean Distance",
                subtitle_text =  NULL
                )


# Plot separately
p_l1
p_l2

# Plot together
p_l1 + p_l2 + plot_layout(guides = "collect")







