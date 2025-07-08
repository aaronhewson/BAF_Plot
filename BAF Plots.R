## This script produces plots of B Allele Frequency (BAF) throughout the apple HFTH1 genome. 
# Each sample is plotted and saved to a PNG separately, with a PDF file saved that combines all


#Load packages
library("Qploidy")
library("dplyr")
library("ggplot2")

#Set working directory
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/BAF Plots")

#Set input and output paths, and PDF output paths
input_dir <- "C:/Users/curly/Desktop/Apple Genotyping/Methods/BAF Plots/CNV_Inputs"
output_dir <- "C:/Users/curly/Desktop/Apple Genotyping/Methods/BAF Plots/Output_Plots"
output_pdf <- "C:/Users/curly/Desktop/Apple Genotyping/Methods/BAF Plots/Output_Plots/All_Plots.pdf"  

#Read SNP BLAST list
BLAST <- read.delim("BLAST results.tsv", header = TRUE, sep = "\t")
  
#List the diploid and triploid CNV files separately
triploids <- read.delim("TriploidSampleNames.txt", header=FALSE)
triploids <- triploids$V1

files <- list.files(input_dir, pattern = ".txt", full.names = TRUE)
files.dip <- files[!(files %in% triploids)]
files.trip <- files[(files %in% triploids)]

#Set PDF parameters - output path, dimensions, plots per page.
pdf(file = output_pdf, width = 6, height = 8)
par(mfrow = c(4,1))

#Setting plot count to zero
plot_count <- 0

#Plot diploid BAF plots
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
  
  #Plotting histogram and exporting as PNG
  png(filename = file.path(output_dir, paste0(file_base,".png")), width = 1000, height = 600)
  par(mar = c(5,5,4,2))
  plot_baf(data, ploidy = 2, area_single = 0.01, add_expected_peaks = TRUE) + labs(title = paste(file_base))
  
  dev.off()
  
  #Plotting histogram for PDF
  par(mar = c(3,4,2,4))
  plot_baf(data, ploidy = 2, area_single = 0.01, add_expected_peaks = TRUE) + labs(title = paste(file_base))
  
  #Counting histograms plotted for PDF  
  plot_count <- plot_count + 1
  
  #Making new page on PDF for every 4 histograms
  if (plot_count %% 4 == 0) {
    par(mfrow = c(4,1))
  }
}

#Plot triploid BAF plots
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
  
  #Plotting histogram and exporting as PNG
  png(filename = file.path(output_dir, paste0(file_base,".png")), width = 1000, height = 600)
  par(mar = c(5,5,4,2))
  plot_baf(data, ploidy = 3, area_single = 0.01, add_expected_peaks = TRUE) + labs(title = paste(file_base))
  
  dev.off()
  
  #Plotting histogram for PDF
  par(mar = c(3,4,2,4))
  plot_baf(data, ploidy = 3, area_single = 0.01, add_expected_peaks = TRUE) + labs(title = paste(file_base))
  
  #Counting histograms plotted for PDF  
  plot_count <- plot_count + 1
  
  #Making new page on PDF for every 4 histograms
  if (plot_count %% 4 == 0) {
    par(mfrow = c(4,1))
  }
}

dev.off()








##WORKING CODE (tested)
  
diploid <- read.delim("1_1_T57L.txt")

triploid <- read.delim("1_1_T57R.txt")

plot_baf(diploid,area_single = 0.01,ploidy=3, add_estimated_peaks = TRUE,add_expected_peaks = TRUE)

plot_baf(triploid,area_single = 0.01,ploidy=3, add_estimated_peaks = TRUE,add_expected_peaks = TRUE)


