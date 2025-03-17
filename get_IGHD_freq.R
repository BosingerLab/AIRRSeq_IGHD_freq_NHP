install.packages("BiocManager")
install.packages("alakazam")

BiocManager::install("alakazam", force = TRUE)

library(stringr)
library(alakazam)

data_path <- "Y:/Runs/Separated_fastq/functional"
files <- list.files(path = data_path,"productive-T.tsv")
files

list_freq <- list()

shortlist <- as.character(unlist(read.table("List.txt")))
shortlist
length(shortlist)

files_short <- files[grepl(paste(shortlist,collapse = '|'),files)]
files_short <- paste0(data_path,"/",files_short)

combined_data <- data.frame()

get_combined <- function(f){
  d <- read.table(f,header= T, sep="\t")
  d$Sample <- gsub(".Rep.*","",d$sample)
  combined_data <<- rbind(combined_data, d)
}
lapply(files_short,get_combined)

colnames(combined_data)
combined_data$d_call[combined_data$d_call == ""] <- "-"

count <- countGenes(combined_data,gene="d_call", mode="gene", groups= "Sample")
write.table(count,"D_gene_counts.txt",quote = F, sep = "\t", row.names = F)

count_d341 <- count[count$gene == "IGHD3-41",]
write.table(count_d341,"D3-41_gene_counts.txt",quote = F, sep = "\t", row.names = F)

