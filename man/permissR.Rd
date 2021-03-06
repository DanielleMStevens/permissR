% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Find-regions.R
\name{permissR}
\alias{permissR}
\title{Load GBFF/GBK and FASTA file}
\usage{
permissR(genbank_file_path, fasta_file_path)
}
\arguments{
\item{genbank_file_path}{takes in path of genbank file}

\item{fasta_file_path}{takes in path of fasta file}
}
\value{
Plots and files of potential permissive sites in a provided bacterial genome
}
\description{
This function loads both a DNA contig fasta file and annotation genbank file
(gbff/gbk) to search for regions a minimum 1.5 kB in size which carry no
encoded information to clone into pSelAct-Express for integrative expression.
}
