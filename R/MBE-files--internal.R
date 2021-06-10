#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 06/09/2021
# Script Purpose: Load files to run thorough find-regions.R to search for promiscuous sites
# Inputs: fasta file, genbank file, iSEScan output
# Outputs: N/A
#-----------------------------------------------------------------------------------------------
#' tries to run ISEScan and asks user for output to load into script
#' 
#' 
#' 
#' @noRd
#' @keywords internal

######################################################################
# library packages need to load
######################################################################

.MBE_files <- function(){
  check_IS_Scan <- readline(prompt = "Do you need to run ISEScan? (y/n): ")
  
  if (check_IS_Scan == "y"){
    number_of_cores <- parallel::detectCores()/2
    update_file_name <- gsub(getwd(), "", fasta_file)
    update_file_name <- paste(".",update_file_name, sep = "")
    
    system(paste('isescan.py', update_file_name, 'proteome hmm --nthread', number_of_cores, sep = ' '))
  }
  
  
  
  
  var1_IS <- readline(prompt = "Is there an output from ISEScan? (y/n): ")
  if (var1_IS == "y"){
    var1_IS <<- TRUE
    cat("Select .fasta.out file from ISEScan Output\n")
    is_element_file <- file.choose()
  }
  if (var1_IS == "n"){
    var1_IS <<- FALSE
  }
}

cat("-------------------------------------------------\n")


#var2_GI <- readline(prompt = "Is there an output from IslandViewer4? (y/n): ")
#if (var2_GI == "y"){
#  var2_GI <- TRUE
#  print("Select IslandViewer4 file Output")
#  GI_file <- file.choose()
#}
#if (var2_GI == "n"){
#  var2_GI <- FALSE
#}

#cat("-------------------------------------------------\n")


