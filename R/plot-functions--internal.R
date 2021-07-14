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
#' @param data_to_set_up_scale inputs data that helps set the scale of the plot
#' @param data_to_plot the data we are interested in actually plotting
#' @param sites_to_plot provides info on the permissive sites found to be able to plots them
#' 
#' 
#' @keywords internal
#' @noRd
annotation_plot <- function(data_to_set_up_scale, data_to_plot) {
  
  ######################################################################
  # functions for subsetiing and plotting data
  ######################################################################
  
  plot <- ggplot2::ggplot() +
    
    ggplot2::geom_rect(data = data_to_set_up_scale, 
                       ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = 0.5, fill = factor(contig)), 
                       stat = "identity", 
                       fill = "white") + 
    
    ggplot2::geom_rect(data = data_to_plot,  
                       ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = 0.5, fill = factor(contig)), 
                       stat = "identity",
                       fill = "grey35") +
    
    ggplot2::facet_grid(. ~ factor(contig, 
                        levels = unique(forcats::fct_inorder(forcats::fct_drop(contigs_file_reformate$contig)))), 
                        scales = 'free_x', 
                        space = 'free_x', 
                        switch = 'x') + 
    
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme(legend.position = "none", 
                   axis.text.y.left = ggplot2::element_blank(), 
                   axis.ticks.y.left = ggplot2::element_blank(),
                   axis.title.y.left = ggplot2::element_text(angle = 360),
                   strip.background = ggplot2::element_blank(),
                   strip.text.x = ggplot2::element_blank(),
                   strip.placement = "outside",
                   axis.text.x = ggplot2::element_blank()) +

    ggplot2::ylab("Genes\nPresent") 
  
  
  return(plot)
}




# function for mobile element annotations
.mobile_element_plot <- function(data_to_set_up_scale, data_to_plot){
  
  
  # if there are IS elements present
  if (nrow(data_to_plot) > 0){
    
    plot <- ggplot2::ggplot() +
      
      ggplot2::geom_rect(data = data_to_set_up_scale, 
                         ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = 0.5, fill = factor(contig)), 
                         stat = "identity", 
                         fill = "white") + 
      
      ggplot2::geom_rect(data = data_to_plot,  
                         ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = 0.5, fill = factor(contig)), 
                         stat = "identity", 
                         fill = "grey35") +
      
      ggplot2::theme_classic() +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      
      ggplot2::theme(legend.position = "none", 
                     axis.text.y.left = ggplot2::element_blank(), 
                     axis.ticks.y.left = ggplot2::element_blank(),
                     axis.title.y.left = ggplot2::element_text(angle = 360),
                     axis.text.x = ggplot2::element_blank(),
                     strip.background = ggplot2::element_blank(),
                     strip.text.x = ggplot2::element_blank(),
                     strip.placement = "outside") +
      
      ggplot2::facet_grid(~factor(contig, 
                         levels = unique(forcats::fct_inorder(forcats::fct_drop(contigs_file_reformate$contig)))), 
                         scales = 'free_x',
                         space = 'free_x', 
                         switch = 'x') + 
      
      ggplot2::ylab("IS Elements\nPresent")
  }
  
  # if there are no IS elements present
  if (nrow(data_to_plot) == 0){
    
    plot <- ggplot2::ggplot() +
      
      ggplot2::geom_rect(data = data_to_set_up_scale, 
                         ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = 0.5, fill = factor(contig)), 
                         stat = "identity", 
                         fill = "white") + 
      
      ggplot2::geom_blank() +
      ggplot2::theme_classic() +
      
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::theme(legend.position = "none", 
                     axis.text.y.left = ggplot2::element_blank(), 
                     axis.ticks.y.left = ggplot2::element_blank(),
                     axis.title.y.left = ggplot2::element_text(angle = 360),
                     axis.text.x = ggplot2::element_blank(),
                     strip.background = ggplot2::element_blank(),
                     strip.text.x = ggplot2::element_blank(),
                     strip.placement = "outside") +
      
      ggplot2::facet_grid(~factor(contig, levels = unique(forcats::fct_inorder(forcats::fct_drop(contigs_file_reformate$contig)))), 
                          scales = 'free_x', 
                          space = 'free_x', 
                          switch = 'x') + 
      
      ggplot2::ylab("IS Elements\nPresent")
  }
  
  return(plot)
}



# function for plotting shannon entropy
.shannon_entropy_plot <- function(data_to_set_up_scale, data_to_plot) {
  
  plot <- ggplot2::ggplot() +
    
    ggplot2::geom_rect(data = data_to_set_up_scale, 
                       ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = 0.5, fill = factor(contig)), 
                       stat = "identity", 
                       fill = "white") + 
    
    ggplot2::geom_line(data = data_to_plot,
                       ggplot2::aes(x = start, y = shannon.entropy, fill = factor(contig)), 
                       stat = "identity", 
                       size = 0.3) +

    
    ggplot2::theme_classic() +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::ylim(1.4,2.0) +
    
    ggplot2::facet_grid(. ~ factor(contig, 
                          levels = unique(forcats::fct_inorder(forcats::fct_drop(contigs_file_reformate$contig)))), 
               scales = 'free_x', space = 'free_x', switch = 'x') + 
    
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), 
          axis.title.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(family = "Arial", size = 11, color = "black"),
          axis.title.y.left = ggplot2::element_text(angle = 360),
          strip.background = ggplot2::element_blank(),
          strip.text.x = ggplot2::element_blank(),
          strip.placement = "outside") +
    ggplot2::ylab("Shannon\nEntropy") 
  
  return(plot)
}



# function for plotting gc content with contig info
.gc_content_plot <- function(data_to_set_up_scale, data_to_plot, sites_to_plot) {
  
  plot <- ggplot2::ggplot() +
    
    ggplot2::geom_rect(data = data_to_set_up_scale, 
                       ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = factor(contig)),
                       fill = "white") +
    
    ggplot2::geom_line(data = data_to_plot, 
                       ggplot2::aes(x = start, y =  gc.content, fill = factor(contig)), 
                       stat = "identity", 
                       size = 0.3, 
                       color = "black") + 
    
    ggplot2::geom_rect(data = sites_to_plot,  
                       ggplot2::aes(xmin = position_1, xmax = position_2, ymin = 0, ymax = 1, fill = factor(contig)), 
                       stat = "identity", 
                       fill = "#0072B2", 
                       alpha = 0.8) +
    

    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    
    ggplot2::facet_grid(. ~ factor(contig, 
                          levels = unique(forcats::fct_inorder(forcats::fct_drop(data_to_plot$contig)))), 
                          scales = 'free_x', space = 'free_x', switch = 'x') + 
    
    ggplot2::theme(axis.text.x = ggplot2::element_text(family = "Arial", size = 11, color = "black", angle = 90), 
                   axis.text.y = ggplot2::element_text(family = "Arial", size = 11, color = "black"),
                   axis.title.y.left = ggplot2::element_text(angle = 360),
                   strip.background = ggplot2::element_rect(fill = "#B5EB5A"),
                   legend.position = "none",
                   strip.text.x = ggplot2::element_text(angle = 90, family = "Arial", size = 14, color = "black"),
                   strip.placement = "outside") +
    
    ggplot2::ylab("%GC") +
    ggplot2::xlab("\nPosition along Chromosome") 
  
  return(plot)
}






#in development - ignore for now

# function for plotting annotations
#.annotation_subplot <- function(data_to_set_up_scale, data_to_plot) {
#  
#  data_to_set_up_scale <- data_to_set_up_scale[c(1,nrow(data_to_set_up_scale)),]
#  data_to_set_up_scale <- data.frame(data_to_set_up_scale[1,1:2],data_to_set_up_scale[2,3])
#  colnames(data_to_set_up_scale) <- c("contig","position_start","position_end")
#  print(data_to_set_up_scale)
#  
#  data_to_plot <- cbind(data_to_plot, rep(data_to_set_up_scale[,2:3],1))
#  data_to_plot[data_to_plot == "-"] <- -1
#  data_to_plot[data_to_plot == "+"] <- +1
#  
#  print(data_to_plot)
#  
#  plot <- data_to_plot %>%
#    ggplot(aes(xmin = start, xmax = end, y = contig)) +
#    geom_linerange(data = data_to_set_up_scale, aes(xmin = position_start,, xmax = position_end),
#                   color = "black") +
#    geom_gene_arrow(aes(forward = strand),
#                    fill = "grey", color = "black",
#                    arrowhead_height = grid::unit(6, "mm"), 
#                    arrowhead_width = unit(2, "mm"),
#                    arrow_body_height = grid::unit(6, "mm")) +
#    theme_genes()  +
#    theme(axis.title.y = element_blank(),
#          axis.line.x = element_blank(),
#          #axis.ticks.x = element_blank(),
#          #axis.text.x = element_blank(),
#          axis.text.y = element_blank(),
#          legend.position = "none") +
#    ggrepel::geom_text_repel(data = data_to_plot %>% mutate(start = (start + end)/2), aes(x = start, y = contig, label = gene), 
#                             inherit.aes = F, nudge_y = 1)
  #facet_grid(. ~ factor(contig, levels = unique(forcats::fct_inorder(forcats::fct_drop(contigs_file_reformate$contig)))), scales = 'free_x', space = 'free_x', switch = 'x')
  
  #plot <- ggplot(data = data_to_plot, 
  #               aes(xmin = start, xmax = end, label = gene, fill = gene, y = contig, forward = strand)) + 
  #  
  #
  #  geom_linerange(data = data_to_set_up_scale, aes(xmin = position_start,, xmax = position_end)) +
  #  geom_gene_arrow() + 
  #scale_y_discrete(limits = )
  
  #  theme_genes() +
  #  theme(legend.position = "none", 
  #        axis.title.y = element_blank(), 
  #        axis.text.x = element_text(color = "black"), 
  #        
  #ggfittext::geom_fit_text(min.size = 4, place = "center", vjust = 1)
  
  
  
  #geom_rect(data = data_to_set_up_scale, 
  #          aes(xmin = start, xmax = end, ymin = 0, ymax = 0.5, fill = factor(contig)), 
  #          stat = "identity", fill = "white") + 
  
  
  
  #    scale_y_continuous(expand = c(0, 0)) +
  #    scale_x_continuous(expand = c(0, 0)) +
  #    theme(legend.position = "none", axis.text.y.left = element_blank(), 
  #          axis.ticks.y.left = element_blank(),
  #          axis.title.y.left = element_text(angle = 360),
  #          strip.background = element_blank(),
  #          strip.text.x = element_blank(),
  #          strip.placement = "outside",
  #          axis.text.x = element_blank()) +
  #     + 
  #    
  
  
#  return(plot)
#}


