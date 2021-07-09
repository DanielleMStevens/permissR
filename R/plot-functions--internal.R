#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 06/09/2021
# Script Purpose: Find promiscuous sites in bacterial genome
# Inputs: N/A - from input-file.R 
# Outputs: Figure which plots A. total regions with annotaiton, GC content, shannon entropy, and IS/PAI elements and B. inidividual
#           plots for main hits as well as tab-delimited file output of major hit regions ranked by complexity
#-----------------------------------------------------------------------------------------------
#' 
#' 
#' 
#' 
#' @noRd
#' @keywords internal

######################################################################
# functions for subsetiing and plotting data
######################################################################


# function for plotting annotations
.annotation_plot <- function(data_to_set_up_scale, data_to_plot) {
  
  plot <- ggplot() +
    
    geom_rect(data = data_to_set_up_scale, 
              aes(xmin = start, xmax = end, ymin = 0, ymax = 0.5, fill = factor(contig)), 
              stat = "identity", fill = "white") + 
    
    
    geom_rect(data = data_to_plot,  
              aes(xmin = start, xmax = end, ymin = 0, ymax = 0.5, fill = factor(contig)), 
              stat = "identity", fill = "grey35") +
    
    theme_classic() +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme(legend.position = "none", axis.text.y.left = element_blank(), 
          axis.ticks.y.left = element_blank(),
          axis.title.y.left = element_text(angle = 360),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          strip.placement = "outside",
          axis.text.x = element_blank()) +
    facet_grid(. ~ factor(contig, 
                          levels = unique(forcats::fct_inorder(forcats::fct_drop(contigs_file_reformate$contig)))), 
               scales = 'free_x', space = 'free_x', switch = 'x') + 
    ylab("Genes\nPresent") 
  
  
  return(plot)
}



# function for plotting annotations
#.annotation_subplot <- function(data_to_set_up_scale, data_to_plot) {
  
#  plot <-
    
    #geom_rect(data = data_to_set_up_scale, 
    #          aes(xmin = start, xmax = end, ymin = 0, ymax = 0.5, fill = factor(contig)), 
    #          stat = "identity", fill = "white") + 
    
    
  
    #ggplot(data_to_plot, 
    #       aes(xmin = start, xmax = end, label = gene, fill = gene, y = contig, forward = strand)) + 
    #  geom_gene_arrow(arrowhead_height = grid::unit(6, "mm"), 
    #                  arrowhead_width = unit(2, "mm"),
    #                  arrow_body_height = grid::unit(6, "mm")) + 
  #ggfittext::geom_fit_text(min.size = 4, place = "center", vjust = 1)

    
#    theme_classic() +
#    scale_y_continuous(expand = c(0, 0)) +
#    scale_x_continuous(expand = c(0, 0)) +
#    theme(legend.position = "none", axis.text.y.left = element_blank(), 
#          axis.ticks.y.left = element_blank(),
#          axis.title.y.left = element_text(angle = 360),
#          strip.background = element_blank(),
#          strip.text.x = element_blank(),
#          strip.placement = "outside",
#          axis.text.x = element_blank()) +
#    facet_grid(. ~ factor(contig, 
#                          levels = unique(forcats::fct_inorder(forcats::fct_drop(contigs_file_reformate$contig)))), 
#               scales = 'free_x', space = 'free_x', switch = 'x') + 
#    ylab("Genes\nPresent") 
  
  
#  return(plot)
#}


# function for mobile element annotations
.mobile_element_plot <- function(data_to_set_up_scale, data_to_plot){
  
  if (nrow(data_to_plot) > 0){
    plot <- ggplot() +
      geom_rect(data = data_to_set_up_scale, 
                aes(xmin = start, xmax = end, ymin = 0, ymax = 0.5, fill = factor(contig)), 
                stat = "identity", fill = "white") + 
      
      geom_rect(data = data_to_plot,  
                aes(xmin = start, xmax = end, ymin = 0, ymax = 0.5, fill = factor(contig)), 
                stat = "identity", fill = "grey35") +
      
      theme_classic() +
      scale_y_continuous(expand = c(0, 0)) +
      scale_x_continuous(expand = c(0, 0)) +
      theme(legend.position = "none", axis.text.y.left = element_blank(), 
            axis.ticks.y.left = element_blank(),
            axis.title.y.left = element_text(angle = 360),
            axis.text.x = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            strip.placement = "outside") +
      facet_grid(~factor(contig, 
                         levels = unique(forcats::fct_inorder(forcats::fct_drop(contigs_file_reformate$contig)))), 
                 scales = 'free_x', space = 'free_x', switch = 'x') + 
      ylab("IS Elements\nPresent")
  }
  if (nrow(data_to_plot) == 0){
    plot <- ggplot() +
      geom_rect(data = data_to_set_up_scale, 
                aes(xmin = start, xmax = end, ymin = 0, ymax = 0.5, fill = factor(contig)), 
                stat = "identity", fill = "white") + 
      geom_blank() +
      
      theme_classic() +
      scale_y_continuous(expand = c(0, 0)) +
      scale_x_continuous(expand = c(0, 0)) +
      theme(legend.position = "none", axis.text.y.left = element_blank(), 
            axis.ticks.y.left = element_blank(),
            axis.title.y.left = element_text(angle = 360),
            axis.text.x = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            strip.placement = "outside") +
      facet_grid(~factor(contig, 
                         levels = unique(forcats::fct_inorder(forcats::fct_drop(contigs_file_reformate$contig)))), 
                 scales = 'free_x', space = 'free_x', switch = 'x') + 
      ylab("IS Elements\nPresent")
  }
  
  return(plot)
}



# function for plotting shannon entropy
.shannon_entropy_plot <- function(data_to_set_up_scale, data_to_plot) {
  plot <- ggplot() +
    
    geom_rect(data = data_to_set_up_scale, 
              aes(xmin = start, xmax = end, ymin = 0, ymax = 0.5, fill = factor(contig)), 
              stat = "identity", fill = "white") + 
    
    geom_line(data = data_to_plot,
              aes(x = start, y = shannon.entropy, fill = factor(contig)), 
              stat = "identity", size = 0.3) +
    
    facet_grid(. ~ factor(contig, 
                          levels = unique(forcats::fct_inorder(forcats::fct_drop(contigs_file_reformate$contig)))), scales = 'free_x', space = 'free_x', switch = 'x') + 
    theme_classic() +
    scale_x_continuous(expand = c(0, 0)) +
    ylim(1.4,2.0) +
    theme(axis.text.x = element_blank(), 
          axis.title.x = element_blank(),
          axis.text.y = element_text(family = "Arial", size = 11, color = "black"),
          axis.title.y.left = element_text(angle = 360),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          strip.placement = "outside") +
    ylab("Shannon\nEntropy") 
  
  return(plot)
}


#data_to_set_up_scale
# function for plotting gc content with contig info
.gc_content_plot <- function(data_to_set_up_scale, data_to_plot, sites_to_plot) {
  
  plot <- ggplot() +
    
    geom_rect(data = data_to_set_up_scale, 
              aes(xmin = start, xmax = end, ymin = 0, ymax = 1),
              fill = "white") +
    
    geom_line(data = data_to_plot, 
              aes(x = start, y =  gc.content, fill = factor(contig)), 
              stat = "identity", size = 0.3, color = "black") + 
    
    geom_rect(data = sites_to_plot,  
              aes(xmin = position_1, xmax = position_2, ymin = 0, ymax = 1, fill = factor(contig)), 
              stat = "identity", fill = "#0072B2", alpha = 0.8) +
    
    facet_grid(. ~ factor(contig, 
                          levels = unique(forcats::fct_inorder(forcats::fct_drop(data_to_plot$contig)))), scales = 'free_x', space = 'free_x', switch = 'x') + 
    theme_classic() +
    scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme(axis.text.x = element_text(family = "Arial", size = 11, color = "black", 
                                     angle = 90), 
          axis.text.y = element_text(family = "Arial", size = 11, color = "black"),
          axis.title.y.left = element_text(angle = 360),
          strip.background = element_rect(fill = "#B5EB5A"),
          legend.position = "none",
          strip.text.x = element_text(angle = 90, family = "Arial", 
                                      size = 14, color = "black"),
          strip.placement = "outside") +
    ylab("%GC") +
    xlab("\nPosition along Chromosome") 
  
  return(plot)
}


