#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 02/22/2021
# Script Purpose: Load files to run thorough find-regions.R to search for promiscuous sites
# Inputs: fasta file, genbank file, iSEScan output
# Outputs: N/A
#-----------------------------------------------------------------------------------------------
# Using GUI to load inital files if path is not provided by user


######################################################################
# load initial files to start script
######################################################################

.input_files <- function(){

  cat("Next Select Genbank - GBK/GBFF File\n")
  genbank_file <<- file.choose()
  
  
  cat("Select Corresponding DNA Fasta File\n")
  fasta_file <<- file.choose()
  
  
  cat("-------------------------------------------------\n")
  #return(genbank_file,fasta_file)
}

