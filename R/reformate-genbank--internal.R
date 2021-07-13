#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 02/22/2021
# Script Purpose: Load files to run thorough find-regions.R to search for promiscuous sites
# Inputs: fasta file, genbank file, iSEScan output
# Outputs: N/A
#-----------------------------------------------------------------------------------------------
#' Using GUI to load inital files if path is not provided by user
#' 
#' 
#' 
#' @noRd
#' @keywords internal

######################################################################
# load initial files to start script
######################################################################





#test_genbank <- function(file_in){
  
  
 # if(grepl("\\b.gbff\\b", genbank_file) == TRUE){
#    for (i in 1:nrow(file_in)){
#      
      # skip empty annotations 
 #     if(length(file_in$FEATURES[[i]]) == 0){
#        next
#      }
      
#      if(length(file_in$FEATURES[[i]]) != 0)){
 #       if(any(grepl("organism", colnames(file_in$FEATURES[[i]]))) == TRUE){
  #        next
  #     }
  #      if(file_in$FEATURES[[i]]$type == "CDS"){
  #        gene_name <- file_in$FEATURES[[i]]$locus_tag
  #        locus_tag <- file_in$FEATURES[[i]]$locus_tag
  #      }
  #      if(file_in$FEATURES[[i]]$type == "gene"){
  #        gene_name <- file_in$FEATURES[[i]]$locus_tag
  #        locus_tag <- file_in$FEATURES[[i]]$locus_tag
  #      }
  #      if(file_in$FEATURES[[i]]$type == "regulatory"){
  #        gene_name <- file_in$FEATURES[[i]]$regulatory_class
  #        locus_tag <- file_in$FEATURES[[i]]$db_xref
   #     }
   #   }
        

   #   print(i)
          
  #    temp_df <- data.frame("gene" = gene_name, 
  #                          "locus-tag" = locus_tag,
  #                          'type' = file_in$FEATURES[[i]]$type,
  #                          'strand' = file_in$FEATURES[[i]]$strand,
  #                          'start' = file_in$FEATURES[[i]]$start,
  #                          'end' = file_in$FEATURES[[i]]$end)
  #        
  #    anno_file_reformate <- rbind(anno_file_reformate, temp_df)
          
  #  }
 # }
  
 # return(anno_file_reformate)
  
#}


#Problem_annotations <- c("translation","organism","linkage_evidence")

#for (i in 1:length(annotation_file$FEATURES)){

#if(length(annotation_file$FEATURES[[i]]) != 0){
#  if(any(grepl(paste(Problem_annotations, collapse = "|"), colnames(annotation_file$FEATURES[[i]]))) == FALSE){
#if(any(grepl("locus_tag", colnames(annotation_file$FEATURES[[i]]))) == FALSE){
#  next
#}

# extract gene name when possible, otherwise use locus-tag
#if(any(grepl("gene",colnames(annotation_file$FEATURES[[i]]))) == TRUE){
#  gene_name <- annotation_file$FEATURES[[i]]$gene
#  locus_tag <- annotation_file$FEATURES[[i]]$locus_tag
#}
# if(any(grepl("gene",colnames(annotation_file$FEATURES[[i]]))) == FALSE){
#  gene_name <- annotation_file$FEATURES[[i]]$locus_tag
#  locus_tag <- annotation_file$FEATURES[[i]]$locus_tag
#}

    #}
    
    # occasuinally some CDS sequences get through, this is to remove those redudancies
    #anno_file_reformate <- anno_file_reformate[anno_file_reformate$type %in% c("gene","regulatory"),]
  #}
  
 # if(grepl("\\b.gbk\\b", genbank_file) == TRUE){
  #  Problem_annotations <- c("organism","linkage_evidence")
   # for (i in 1:length(annotation_file$FEATURES)){
  #    if(length(annotation_file$FEATURES[[i]]) == 0){
   #     next
   #   }
   #   
    #  if(length(annotation_file$FEATURES[[i]]) != 0){
    #    if(any(grepl(paste(Problem_annotations, collapse = "|"), colnames(annotation_file$FEATURES[[i]]))) == FALSE){
          #if(any(grepl("locus_tag", colnames(annotation_file$FEATURES[[i]]))) == FALSE){
          #  next
          #}
          
          # extract gene name when possible, otherwise use locus-tag
     #     if(any(grepl("gene",colnames(annotation_file$FEATURES[[i]]))) == FALSE){
    #        gene_name <- annotation_file$FEATURES[[i]]$locus_tag
     #       locus_tag <- annotation_file$FEATURES[[i]]$locus_tag
    #      }
     #     if(annotation_file$FEATURES[[i]]$type == "regulatory"){
    #        gene_name <- annotation_file$FEATURES[[i]]$regulatory_class
    #        locus_tag <- annotation_file$FEATURES[[i]]$db_xref
    #      }
    #      
         # temp_df <- data.frame("gene" = gene_name, 
        #                        "locus-tag" = locus_tag,
        #                        'type' = annotation_file$FEATURES[[i]]$type,
        #                        'strand' = annotation_file$FEATURES[[i]]$strand,
        #                        'start' = annotation_file$FEATURES[[i]]$start,
        #                        'end' = annotation_file$FEATURES[[i]]$end)
        #  
        #  anno_file_reformate <- rbind(anno_file_reformate, temp_df)
          #print(anno_file_reformate)
          
        #}

