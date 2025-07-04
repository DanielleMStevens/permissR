<img align="right" width="140" height="150" src="https://github.com/DanielleMStevens/permissR/blob/main/github_images/permissR_logo.png">

# permissR - an R Package for finding Bacterial Permissive Sites

* __Danielle M. Stevens__ <a itemprop="sameAs" content="https://orcid.org/0000-0001-5630-137X" href="https://orcid.org/0000-0001-5630-137X" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a>   </br>
_Dept. of Plant Pathology, University of California, Davis_, </br>
_Integrative Genetics and Genomics Graduate Group, University of California, Davis_

* __Gitta Coaker__ <a itemprop="sameAs" content="https://orcid.org/0000-0003-0899-2449" href="https://orcid.org/0000-0003-0899-2449" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a> </br>
_Dept. of Plant Pathology, University of California, Davis_

---
__permissR__ is an R script for finding permissive sites, i.e. sites which appear to contain no functionally encoded information, proteins, regulatory regions, and mobile elements.

This package works intandom with pSelAct Express, a vector modified for gene expression in the chromosome of a permissive site. Therefore, the decisions to find such sites are ranked for this purpose.

permissR tries to identify sites which are at least 1.5 kB in length (for relatively easy, specific recombination), that contain no protein encoding genes, no mobile elements which may trigger structural changes, and are ranked by complexity to make cloning a little easier (i.e. reasonable GC-content). The outputs of this pipeline include one folder of figures which inlcudes an asessment of the overall genome, figures for regions of interest, tab-delimited text file with the ranked list of sites, and finally a tab-delimited file of any sites which share some degree of homology with the predicted site (this is only procurred if there are any). This pipeline was designed for bacteria and in its current state, not recommended on eurkaryotes.

![](/github_images/permissR_figure.png)



## 1. Install Package

To install, type the following code into your R console:


```
install.packages("devtools")
devtools::install_github("DanielleMStevens/permissR")
```

## 2. Install and run ISEScan

ISEScan can be installed easily with conda. For info on how to intall conda, [see this Github repo](https://github.com/DanielleMStevens/ROS_production_review/blob/master/process_files.md). We can then run it in one of two ways: 1) on the terminal directly or 2) on the terminal through the R-console. 

On the terminal, you'll still need to use conda to install isescan. 
```shell
conda install isescan
```


### To run it on the terminal follow below:
```
# run the following command, update the number of threads based on your computer

isescan.py --seqfile /path/to/fasta_file.fasta --output fasta_file --nthread 6
```

### If you want to run ISEScan while in R (i.e. the script will input the commands for you), you'll need to run the following commands:

The on the R-console:
```R
usethis::edit_r_profile()

# edit path to conda version of python. Anaconda should be located in your home or user directory.
Sys.setenv(PATH = paste("/home/danimstevens/anaconda3/bin", Sys.getenv("PATH"), sep=":")) 
```

Restart R and check to make sure the right version/path of python is used. 

```R
system('which python')
#/home/danimstevens/anaconda3/bin/python
```

permissR will ask the user if you want to run ISEScacn on the genome selected and will run it, but this will only work if the above steps are caried out.

## 3. Run permissR on the R-console

If the user wants to run the function similar as one would on the command line, you can feed permissR the path (relative or absolute is fine) of each file such that it looks like the following:

```
permmissR::permissR("/path/to/genbank/file.gbff", "/path/to/fasta/file.fasta")
```

Or you can run the function with no arguments and a GUI interface will appear, letting the user select each file which is required to run.
```
permissR::permissR()

# The console will alert the user to select the files to analyze. 
# This will include a fasta file of the DNA contigs, a genbank (.gbff) file, 
# and the output of ISEScan.
```
