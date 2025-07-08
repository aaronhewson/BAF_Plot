setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/BAF Plots")

library("Qploidy")


diploid <- read.delim("1_1_T57L.txt")

triploid <- read.delim("1_1_T57R.txt")

plot_baf(diploid,area_single = 0.01,ploidy=3, add_estimated_peaks = TRUE,add_expected_peaks = TRUE)

plot_baf(triploid,area_single = 0.01,ploidy=3, add_estimated_peaks = TRUE,add_expected_peaks = TRUE)
