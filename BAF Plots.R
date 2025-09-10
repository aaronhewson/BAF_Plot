## This script plots the B Allele Frequency (BAF) value of SNPs at their location in the apple genome HFTH1 v1.0.a1 (https://doi.org/10.1038/s41467-019-09518-x). 
## Each sample is plotted and saved to a PNG separately, and a PDF files with all plots compiled is also saved.

# Load Packages -----------------------------------------------------------

library("Qploidy")
library("dplyr")
library("ggplot2")




# Set Working Directory and Input/Output Paths ----------------------------

#Set working directory
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/BAF Plots/Inputs")

#Set input and output paths
input_dir <- "C:/Users/curly/Desktop/Apple Genotyping/Methods/BAF Plots/CNV_Inputs"
output_dir <- "C:/Users/curly/Desktop/Apple Genotyping/Methods/BAF Plots/Output_Plots"



# Read SNP Locations from BLAST Search ------------------------------------
BLAST <- read.delim("BLAST results.tsv", header = TRUE, sep = "\t")



# List Diploid and Triploid Samples ---------------------------------------
triploids <- read.delim("TriploidSampleNames.txt", header=FALSE)
triploids <- triploids$V1

files <- list.files(input_dir, pattern = ".txt", full.names = TRUE)
files.dip <- files[!(files %in% triploids)]
files.trip <- files[(files %in% triploids)]




# Plot BAF for Diploids ---------------------------------------------------
for (file in  files.dip){
  #Reading in data
  data <- read.table(file, header = TRUE,row.names = NULL)

  #Keep only BLAST-matched SNPs, add locations, remove other columns
  data <- data %>% semi_join(BLAST, by = "ProbeSetName") %>% select(-Sample,-Chromosome,-Position,-Log2Ratio) %>% left_join(BLAST, by = "ProbeSetName") %>% select(-ProbeSetName)
  
  #Rearrange columns
  data <- data[,c("Chromosome", "Position", "BAF")]
  
  #Rename columns
  colnames(data) <- c("Chr", "Position", "sample")
  
  #Removing file extension from filename
  file_base <- tools::file_path_sans_ext(basename(file))
  
  #Plot, and export as PNG
  png(filename = file.path(output_dir, paste0(file_base,".png")), width = 1500, height = 1000)
  print(plot_baf(data, ploidy = 2, area_single = 0, dot.size = 2, font_size = 24, add_expected_peaks = TRUE) + labs(title = paste(file_base)))
  
  dev.off()
  }




# Plot BAF for Triploids --------------------------------------------------


for (file in  files.trip){
  #Reading in data
  data <- read.table(file, header = TRUE,row.names = NULL)
  
  #Keep only BLAST-matched SNPs, add locations, remove other columns
  data <- data %>% semi_join(BLAST, by = "ProbeSetName") %>% select(-Sample,-Chromosome,-Position,-Log2Ratio) %>% left_join(BLAST, by = "ProbeSetName") %>% select(-ProbeSetName)
  
  #Rearrange columns
  data <- data[,c("Chromosome", "Position", "BAF")]
  
  #Rename columns
  colnames(data) <- c("Chr", "Position", "sample")
  
  #Removing file extension from filename
  file_base <- tools::file_path_sans_ext(basename(file))
  
  #Plot, and export as PNG
  png(filename = file.path(output_dir, paste0(file_base,".png")), width = 1500, height = 1000)
  print(plot_baf(data, ploidy = 3, area_single = 0, dot.size = 2, font_size = 24, add_expected_peaks = TRUE) + labs(title = paste(file_base)))
 
   dev.off()
}