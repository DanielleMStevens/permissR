#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 06/09/2021
# Script Purpose: Find permisive sites in bacterial genome
# Inputs: Genbank and Fasta file
# Outputs: Figure which plots A. total regions with annotaiton, GC content, shannon entropy, and IS/PAI elements and B. inidividual
#           plots for main hits as well as tab-delimited file output of major hit regions ranked by complexity
#-----------------------------------------------------------------------------------------------

#' Load GBFF/GBK and FASTA file
#'
#' This function loads both a DNA contig fasta file and annotation genbank file
#' (gbff/gbk) to search for regions a minimum 1.5 kB in size which carry no
#' encoded information to clone into pSelAct-Express for integrative expression.
#' 
#'
#' @param infile Path to the input file
#' @return Plots and files of potential permissive sites in a provided bacterial genome
#' @export
#' 
#' 
permissR <- function(genbank_file_path, fasta_file_path){

######################################################################
# ascii art for running script
######################################################################
  
  permissR <- r"{

  -------------------------------------------------------------------------
  permissR - An R Script for finding Permissive Sites in Bacterial Genomes 
  -------------------------------------------------------------------------
  
  
  }"
  
  cat(permissR)
  
  strain_name <- readline(prompt = "Before analyzing your genome, what is the strain's name (- or _ only, no spaces): ")
  
  
  
  ######################################################################
  # import files to process
  ######################################################################
  
  # run internal function to load inital files if path is not given
  
  if (missing(genbank_file_path) || missing(fasta_file_path)){
    .input_files()
  }
  else {
    fasta_file <- fasta_file_path
    genbank_file <- genbank_file_path
  }
  
  
  ######################################################################
  # run script for IS element prediction
  ######################################################################
  
    # run internal function to call if is elements need to be predicted/to load output
    .MBE_files()
  
  
  ######################################################################
  # check files and begin analysis
  ######################################################################
  
    cat("Cheking Files....\n")
    #if (str_split(fasta_file, "\\.")[[1]][2] != "fasta"){
    #  return(print("You imputed the files in the wrong order. Please retry in the following order: genbank file path, fasta file path."))
    #}
    #if (any(c("gbff","gbk") %in% str_split(genbank_file, "\\.")[[1]][2]) == FALSE){
    #  return(print("You imputed the files in the wrong order. Please retry in the following order: genbank file path, fasta file path."))
    #}
      
  ######################################################################
  # upload file of gene positions
  ######################################################################
  
  
    cat("Loading Genkbank File...\n") 
    annotation_file <- genbankr::parseGenBank(file = genbank_file, text = readLines(genbank_file), ret.seq = FALSE)
    
  
    # reformate annotation file - determine how many contigs are present
    cat("Determining number of contigs....\n")
    number_of_contigs <- 0
    contigs_file_reformate <- data.frame("contig" = character(0), "start" = numeric(0), "end" = numeric(0))
    for (i in 1:length(annotation_file$FEATURES)){
      if(any(grepl("organism", colnames(annotation_file$FEATURES[[i]]))) == TRUE){
        number_of_contigs <- number_of_contigs + 1
        temp_df <- data.frame("contig" = paste("Contig", number_of_contigs, sep = "_"),
                              "start" = annotation_file$FEATURES[[i]]$start,
                              "end" = annotation_file$FEATURES[[i]]$end)
        contigs_file_reformate <- rbind(contigs_file_reformate, temp_df)
      }
    }
    
    
    # check contig size/genome quality to see if it will work with script
    if (nrow(subset(contigs_file_reformate, end > 15000)) == 0){
      return(print("Your genome is rather fragmented. We do not recommend using this pipeline on your genome. Exiting now..."))
    }
    
  
   ######################################################################
   # extract annotation information from gebank
  ######################################################################
      
    # determine type of genbank file (.gbff vs .gbk)
    cat("Reformating Genbank File....\n")
    
    anno_file_reformate <- data.frame("gene" = character(0), "locus-tag" = character(0), 
                                       'type' = character(0), 'strand' = character(0), 
                                       'start' = numeric(0), 'end' = numeric(0))
    
    anno_file_reformate <- .reformate_genbank()

    
  ######################################################################
  # Scalling gene positions based on chromosome information
  ######################################################################
  
    cat("Adjusting annotation for when multiple contigs are present...\n")
    
    contig_number <- 1
    for (i in 1:nrow(anno_file_reformate)){
      if (i < nrow(anno_file_reformate)){
        anno_file_reformate[i,7] <- contigs_file_reformate[contig_number,1]
        if (anno_file_reformate$locus.tag[i] != anno_file_reformate$locus.tag[i+1]){
          if (anno_file_reformate$start[i+1] < anno_file_reformate$start[i]){
            contig_number <- contig_number + 1
          }
        }
      }
      if (i == nrow(anno_file_reformate)){
        anno_file_reformate[i,7] <- contigs_file_reformate[contig_number,1]
      }
    }
    
    
    # change order of columns and update name
    anno_file_reformate <- anno_file_reformate[,c(7,1,2,3,4,5,6)]
    colnames(anno_file_reformate) <- c("contig","gene","locus-tag",'type','strand','start','end')
    
  
  ######################################################################
  # find differences between coding genes
  ######################################################################
    
    cat("-------------------------------------------------\n")
    cat("Calculating distance Between Protein-coding Genes...\n")
    dist_between_genes <- data.frame("contig" = character(0), "gene-1" = character(0), "gene-2" = character(0), 
                                     "position-1" = numeric(0), "position-2" = numeric(0), "distance" = numeric(0))
    
    
    for (i in 1:nrow(anno_file_reformate)){
      if (i < nrow(anno_file_reformate)){
        if (anno_file_reformate[i,1] != anno_file_reformate[i+1,1]){
          next
        }
          distance_between_genes <- anno_file_reformate[i+1,6] - anno_file_reformate[i,7]
          temp_df <- data.frame("contig" = anno_file_reformate[i,1],
                                "gene-1" = anno_file_reformate[i,2],
                                "gene-2" = anno_file_reformate[i+1,2],
                                "position-1" = anno_file_reformate[i,7],
                                "position-2" = anno_file_reformate[i+1,6],
                                "distance" = distance_between_genes)
          dist_between_genes <- rbind(dist_between_genes, temp_df)
      }
    }
    
    
  ######################################################################
  # filter regions with minimum distance of 1.5 kb
  ######################################################################
    
    
    dist_between_genes <- subset(dist_between_genes, dist_between_genes$distance > 1500)
  
    
  ######################################################################
  # Pull out sequence of regions of interest
  ######################################################################
  
    cat("-------------------------------------------------\n")
    cat("Loading Fasta file...\n")
    # import fasta file
    whole_genome <- seqinr::read.fasta(file = fasta_file, seqtype = "DNA")

    
    # rename names catagory for conistency between fasta files
    reference_names <- data.frame()
    for (i in 1:length(whole_genome)){
      reference_names <- rbind(reference_names, data.frame("Original_name" = names(whole_genome)[i], "contigs" =  paste("Contig", i, sep = "_")))
      names(whole_genome)[i] <- paste("Contig", i, sep = "_")
    }
  
  
    cat("Pulling out sequence of site and calculating GC content...\n")
    seq_between_genes <- data.frame(dist_between_genes, "overall-gc-content" = numeric(nrow(dist_between_genes)), 
                                    "sequence" = character(nrow(dist_between_genes)))
    
    seq_between_genes$sequence <- as.character(seq_between_genes$sequence)
    
    for (i in 1:nrow(seq_between_genes)){
      #temp fix until I figure out how to better filter out plasmids
      if ((any(names(whole_genome) %in% dist_between_genes[i,1])) == FALSE){
        next
      }
      if ((any(names(whole_genome) %in% dist_between_genes[i,1])) == TRUE){
      which_contig <- whole_genome[names(whole_genome) %in% dist_between_genes[i,1]]
      seq_between_genes[i,7] <- seqinr::GC(which_contig[[1]][dist_between_genes$position.1[i]:dist_between_genes$position.2[i]])
      pull_sequence <- as.character(paste(which_contig[[1]][dist_between_genes$position.1[i]:dist_between_genes$position.2[i]], collapse = ""))
      seq_between_genes[i,8] <- as.character(pull_sequence)
      }
    }
    
    
    # rank seq_between_genes and subsequently dist_between_genes by gc-content
    cat("Ranking sites by GC content...\n")
    seq_between_genes <- seq_between_genes[order(seq_between_genes$overall.gc.content),]
    dist_between_genes <- dist_between_genes[match(seq_between_genes$gene.1, dist_between_genes$gene.1),]
    
    
    # remove sequences that con't have contig matches as they are located on plasmids (to be avoided)
    for (i in nrow(seq_between_genes):1){
      if ((any(names(whole_genome) %in% dist_between_genes[i,1])) == FALSE){
        seq_between_genes <- seq_between_genes[-i,]
        dist_between_genes <- dist_between_genes[-i,]
      }
    }
  
    
    
  ######################################################################
  # import ISEScan Output to remove any regions within/near IS Elements
  ######################################################################
    
    if(var1_IS == TRUE){
      cat("-------------------------------------------------\n")
      cat("Loading IS elements file...\n")
      IS_elements <- read.delim(file = is_element_file, header = TRUE, stringsAsFactors = FALSE, sep = "")
      IS_elements <- IS_elements[2:nrow(IS_elements),]
      
      
      cat("Reformating IS Elements....\n")
      IS_elements_reform <- data.frame("contig" = character(0), "start" = numeric(0), "end" = numeric(0))
      for (i in 1:nrow(IS_elements)){
        contig_number <- reference_names[reference_names$Original_name %in% IS_elements$seqID[i],2]
        temp_df <- data.frame("contig" = contig_number,
                              "start" = IS_elements$isBegin[i],
                              "end" = IS_elements$isEnd[i])
        IS_elements_reform <- rbind(IS_elements_reform, temp_df)
      }
    }
    
    
  ######################################################################
  # import IslandViewer4 Output to remove any regions within GIs
  ######################################################################
    
    #if(var2_GI == TRUE){
    #  cat("Loading GI file...\n")
    #  GI_results <- read.delim(file = GI_file, header = TRUE, stringsAsFactors = FALSE)
    #}
    
  ######################################################################
  # run fasta file against iscan elements script
  ######################################################################
    
    
    
    
  
    
    
    
  ######################################################################
  # calculate gc content across the genome
  ######################################################################
  
    cat("-------------------------------------------------\n")
    cat("Calculating GC Content accross the genome...\n")
    gc_content_whole_genome <- data.frame("contig" = character(0), "start" = numeric(0), "end" = numeric(0), 
                                          "gc-content" = numeric(0))
    
    for (j in 1:length(whole_genome)){
      for (i in seq(from = 1, to = length(whole_genome[[j]]), by = 100)){
        gc_value <- seqinr::GC(whole_genome[[j]][i:(i+1000)])
        temp_df <- data.frame("contig" = paste("Contig", j, sep = "_"),
                              "start" = i,
                              "end" = i + 1000,
                              "gc-content" = gc_value)
        gc_content_whole_genome <- rbind(gc_content_whole_genome, temp_df)
      }
    }
  
  
  ######################################################################
  # calculate shannon entropy across the genome
  ######################################################################
  
    cat("Calculating Shannon Entropy accross the genome...\n")
    
    shannon_entropy_whole_genome <- data.frame("contig" = character(0), "start" = numeric(0), "end" = numeric(0), "shannon-entropy" = numeric(0))
    
    for (j in 1:length(whole_genome)){
      for (i in seq(from = 1, to = length(whole_genome[[j]]), by = 100)){
        shannon_value <- DescTools::Entropy(table(whole_genome[[j]][i:(i+1000)]))
        temp_df <- data.frame("contig" = paste("Contig", j, sep = "_"),
                              "start" = i,
                              "end" = i + 1000,
                              "shannon-entropy" = shannon_value)
        shannon_entropy_whole_genome <- rbind(shannon_entropy_whole_genome, temp_df)
      }
    }
  
  
  
  ######################################################################
  # check if that site has high homology to the remaining part of the genome
  ######################################################################
      
    cat("-------------------------------------------------\n")
    cat("Searching to check if sites have homology elsewhere in the genome...\n")
    matches <- data.frame()
    total <- length(whole_genome)*nrow(seq_between_genes)
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    k <- 0
    
    for (j in 1:nrow(seq_between_genes)){  
      for (i in 1:length(whole_genome)){
        output <- Biostrings::vmatchPattern(seq_between_genes[j,8], paste(whole_genome[[i]][1:length(whole_genome[[i]])], collapse = ""), 
                                  max.mismatch = (nchar(seq_between_genes[j,8])*0.5))
  
        output <- as.data.frame(unlist(output))
        if(nrow(output) > 0){
          temp_df <- data.frame(output, names(whole_genome)[[i]], seq_between_genes[j,4], seq_between_genes[j,5])
          matches <- rbind(matches, temp_df)
  
        }
        
        # update progress bar
        k <- k + 1
        setTxtProgressBar(pb, k)
      }
    }
    close(pb)
    
    
    # check if hits are against itself, if so remove
    colnames(matches) <- c("site_start","site_end","site_width","hit_contig","hit_start","hit_end")
    for (i in seq(from = nrow(matches), to = 1, by = -1)) {
      if(matches$site_start[i] == matches$hit_start[i] && matches$site_end[i] == matches$hit_end[i]){
        matches <- matches[-i,]
      }
    }
    
    
  ######################################################################
  # pull out matches to add to final output to user to warn there may be a homologous site
  ######################################################################  
    
    if (nrow(matches) > 0){
      for (i in 1:nrow(matches)){
        which_contig <- whole_genome[names(whole_genome) %in% matches[i,4]]
        pull_sequence <- as.character(paste(which_contig[[1]][matches$hit_start[i]:matches$hit_end[i]], collapse = ""))
        matches[i,7] <- as.character(pull_sequence)
      }
      colnames(matches) <- c("site_start","site_end","site_width","hit_contig","hit_start","hit_end","homologous_site")
    }
    
    
  ######################################################################
  # filter our small-ish contigs which are unlikely to be sites - minimum 10 kB
  ######################################################################  
    
    cat("Filtering out contigs less than 15 kb....\n")
    contigs_to_filter <- subset(contigs_file_reformate, end > 15000)
    
    # remove small contigs
    anno_file_reformate <- anno_file_reformate[anno_file_reformate$contig %in% contigs_to_filter$contig,]
    IS_elements_reform <- IS_elements_reform[IS_elements_reform$contig %in% contigs_to_filter$contig,]
    shannon_entropy_whole_genome <- shannon_entropy_whole_genome[shannon_entropy_whole_genome$contig %in% contigs_to_filter$contig,]
    gc_content_whole_genome <- gc_content_whole_genome[gc_content_whole_genome$contig %in% contigs_to_filter$contig,]
    
    contigs_file_reformate <- contigs_file_reformate[contigs_file_reformate$contig %in% contigs_to_filter$contig,]
  
  
  ######################################################################
  # checking contigs for plasmid components - not a good method atm
  ######################################################################  
    
    contigs_file_reformate <- subset(contigs_file_reformate, contig == unique(gc_content_whole_genome$contig))
    
    #making sure all files match contig number
    anno_file_reformate <- subset(anno_file_reformate, contig == contigs_file_reformate$contig)
    
    
  
  
  ######################################################################
  # plotting and saving data for whole genome
  ######################################################################
  
    cat("-------------------------------------------------\n")
    cat("Plotting data for the whole genome and exporting the figure...\n")
    
    if(var1_IS == FALSE){
      whole_genome_plot <- suppressWarnings(patchwork::wrap_plots(list(
                                              .annotation_plot(contigs_file_reformate, anno_file_reformate) + 
                                                ggtitle(strain_name) +
                                                theme(plot.title = element_text(lineheight = 1.2, face = "bold", family = "Arial", size = 16)), 
                                              .shannon_entropy_plot(contigs_file_reformate, shannon_entropy_whole_genome), 
                                              .gc_content_plot(contigs_file_reformate, gc_content_whole_genome, dist_between_genes)), 
                                              ncol = 1, nrow = 3, heights = c(0.4,0.8,0.8)))
    }
    if(var1_IS == TRUE){
      whole_genome_plot <- suppressWarnings(patchwork::wrap_plots(list(
                                              .annotation_plot(contigs_file_reformate, anno_file_reformate) +
                                                ggtitle(strain_name) +
                                                theme(text = element_text(size = 14),
                                    
                                                  plot.title = element_text(lineheight = 1.2, face = "bold", family = "Arial", size = 18)),
                                              .mobile_element_plot(contigs_file_reformate, IS_elements_reform) + 
                                                theme(text = element_text(size = 14)),
                                              .shannon_entropy_plot(contigs_file_reformate, shannon_entropy_whole_genome) + 
                                                theme(text = element_text(size = 14), axis.text.y = element_text(size = 14)),
                                              .gc_content_plot(contigs_file_reformate, gc_content_whole_genome, dist_between_genes) + 
                                                theme(text = element_text(size = 14), 
                                                      axis.text.x = element_text(size = 14), 
                                                      axis.text.y = element_text(size = 14))),
                                              ncol = 1, nrow = 4, heights = c(0.4,0.8,0.8,0.8)))
    }
    
    
    #print(whole_genome_plot)
    whole_genome_file_name <- paste(strain_name,"whole_genome_plot.pdf", sep = "_")
    
    
    # before saving plots, need to define output directory
    main_path <- getwd()
    output_directory <- paste(strain_name, "outputs", sep = "_")
    dir.create(file.path(main_path, output_directory))
    
    ggsave(filename = file.path(output_directory, whole_genome_file_name), width = 14, height = 6, units = c("in"), device = cairo_pdf)
  
    
  
  ######################################################################
  # regions of intereest
  ######################################################################
  
  
    cat("Plotting data each site predicted...\n")
    # function to filter site data
    sites_filter_for_plotting <- function(data_in_to_filter, filter_info){
      position_1 <- filter_info$position.1[1] - 1000
      position_2 <- filter_info$position.2[1] + 1000
      data_in_to_filter <- subset(data_in_to_filter, contig == filter_info$contig[1])
      data_in_to_filter <- subset(data_in_to_filter, end > filter_info$position.1[1] - 1000)
      data_in_to_filter <- subset(data_in_to_filter, start < filter_info$position.2[1] + 1000)
      
      return(data_in_to_filter)
    }
    
    
  
    # ctreate a list of all subplots then plot them together
    list_of_plots <- list()
    if (var1_IS == FALSE){
      for (i in 1:nrow(dist_between_genes)){
        scale_plot <- sites_filter_for_plotting(gc_content_whole_genome, dist_between_genes[i,]) 
        
        # annotation plot
        plot1 <- .annotation_plot(scale_plot, 
                                 sites_filter_for_plotting(anno_file_reformate, dist_between_genes[i,])) +
          coord_cartesian(xlim = c(dist_between_genes$position.1[i] - 1000, dist_between_genes$position.2[i] + 1000)) +
          annotate("text", x = c(dist_between_genes$position.1[i], dist_between_genes$position.2[i]), y = 0.51, 
                   label = c(dist_between_genes$gene.1[i], dist_between_genes$gene.2[i])) +
          expand_limits(y = 0.52) + 
          xlab("")
        
        
        # shannon entropy plot
        plot2 <- .shannon_entropy_plot(scale_plot, 
                                      sites_filter_for_plotting(shannon_entropy_whole_genome, dist_between_genes[i,]))  +
          coord_cartesian(xlim = c(dist_between_genes$position.1[i] - 1000, dist_between_genes$position.2[i] + 1000))
        
        
        # gc content plot
        plot3 <- .gc_content_plot(scale_plot, 
                                 sites_filter_for_plotting(gc_content_whole_genome, dist_between_genes[i,]), 
                                 dist_between_genes[i,]) + theme(strip.text.x = element_text(angle = 360)) +
          coord_cartesian(xlim = c(dist_between_genes$position.1[i] - 1000, dist_between_genes$position.2[i] + 1000))
        
        
        list_of_plots[[i]] <- patchwork::wrap_plots(list(plot1,plot2,plot3), ncol = 1, nrow = 3, heights = c(0.4,0.8,0.8))
      }
      
      #determine length of the page to plot the subplots onto
      subplots <- patchwork::wrap_plots(list_of_plots, ncol = 2)
    
      length_of_page <- length(subplots$patches$plots)
      if(length_of_page == 1){
        length_of_page <- length_of_page*4.25
      }
      if(length_of_page > 1){
        length_of_page <- (length_of_page/2)*4.25
      }
    }
  
    if (var1_IS == TRUE){
      for (i in 1:nrow(dist_between_genes)){
        scale_plot <- sites_filter_for_plotting(gc_content_whole_genome, dist_between_genes[i,]) 
        
        # annotation plot
        plot1 <- .annotation_plot(scale_plot, 
                                 sites_filter_for_plotting(anno_file_reformate, dist_between_genes[i,])) +
          coord_cartesian(xlim = c(dist_between_genes$position.1[i] - 1000, dist_between_genes$position.2[i] + 1000)) +
          annotate("text", x = c(dist_between_genes$position.1[i], dist_between_genes$position.2[i]), 
                   y = 0.6, 
                   size = 2.5,
                   label = c(dist_between_genes$gene.1[i], dist_between_genes$gene.2[i]), 
                   hjust = 0, 
                   angle = 45) +
          expand_limits(y = 2) +
          xlab("") +
          theme(plot.margin = unit(c(0.4,0.2,-0.7,0.2), "cm"))
        
        
        # is element plot
        
        plot2 <- .mobile_element_plot(scale_plot, 
                                 sites_filter_for_plotting(IS_elements_reform, dist_between_genes[i,])) +
          coord_cartesian(xlim = c(dist_between_genes$position.1[i] - 1000, dist_between_genes$position.2[i] + 1000)) 
  
              
        # shannon entropy plot
        plot3 <- .shannon_entropy_plot(scale_plot, 
                                      sites_filter_for_plotting(shannon_entropy_whole_genome, dist_between_genes[i,]))  +
          coord_cartesian(xlim = c(dist_between_genes$position.1[i] - 1000, dist_between_genes$position.2[i] + 1000)) 
        
        
        # gc content plot
        plot4 <- .gc_content_plot(scale_plot, 
                                 sites_filter_for_plotting(gc_content_whole_genome, dist_between_genes[i,]), 
                                 dist_between_genes[i,]) + theme(strip.text.x = element_text(angle = 360)) +
          coord_cartesian(xlim = c(dist_between_genes$position.1[i] - 1000, dist_between_genes$position.2[i] + 1000)) +
          theme(strip.text = element_text(size = 6))
        
        
        list_of_plots[[i]] <- patchwork::wrap_plots(list(plot1, plot2, plot3, plot4), ncol = 1, nrow = 4, heights = c(0.4, 0.4, 0.4, 0.4)) 
      }
      
      #determine length of the page to plot the subplots onto
      #subplots <- patchwork::wrap_plots(list_of_plots, ncol = 3)
      
      #length_of_page <- length(subplots$patches$plots)
      #if((length_of_page%%2) == 1){
      #  length_of_page <- length_of_page*6
      #}
      #if((length_of_page%%2) == 0){
      #  length_of_page <- (length_of_page/2)*6
      #}
    }
    
    for (i in 1:length(list_of_plots)){
      subplots_file_name <- paste(strain_name,"subplot",i,".pdf", sep = "_")
      
      ggsave(filename = file.path(output_directory, subplots_file_name), plot = list_of_plots[[i]], width = 3.5, height = 6.5, units = c("in"), device = cairo_pdf)
      
    }
    
  
  
  ######################################################################
  # write hits as tab delimited file 
  ######################################################################
  
    top_hits_table_name <- paste(main_path, output_directory, paste(strain_name, "top_hits.txt", sep = "_"), sep = "/")
    write.table(seq_between_genes, file = top_hits_table_name, sep = "\t", row.names = FALSE, quote = FALSE)
    
    # writes out table if any hits have homologous regions elsewhere in the genome
    if(nrow(matches) > 0){
      top_hits_table_name <- paste(main_path, output_directory, paste(strain_name, "homologous_hits.txt", sep = "_"), sep = "/")
      write.table(matches, file = top_hits_table_name, sep = "\t", row.names = FALSE, quote = FALSE)
    }
    
    
}
