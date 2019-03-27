

library(data.table)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  args[3] = "out.txt"
} 

bimfile <- args[1]
summary_stats <- args[2]
output <- args[3]


adjust_beta <- function(bimfile, summary_stats, output){
  
  bim_file <- fread(bimfile)
  
  summ_stats <- fread(summary_stats)
  col.names <- c("SNP", "A1", "effect", "P")
  summ_stats <- data.frame(summ_stats)
  summ_stats <- summ_stats[,col.names] 
## find common SNPs between bim file and summary stats
  summ_stats <- filter(summ_stats, SNP %in% bim_file$V2)
  bim_file <- filter(bim_file, V2 %in% summ_stats$SNP)
## Order both files and remove duplicates from the bim,
## so to have identical SNPs order  
  summ_stats <- summ_stats[order(summ_stats$SNP),]
  
  bim_file <- bim_file[order(bim_file$V2),]
  bim_file <- bim_file[!duplicated(bim_file$V2), ]
## Change the sign if the effect allele is different  
  summ_stats$OR <- ifelse(summ_stats$A1 == bim_file$V5, 
                          summ_stats$effect, summ_stats$effect*-1)
## Save the summary stats  
  summary_stats_path <- summary_stats
  summary_stats_path <- gsub(".tbl$", "",summary_stats_path)
  
  write.table(summ_stats, 
              output,
              quote = F, 
              col.names = T, 
              row.names = F)
  return(summ_stats)
}


adjust_data <- adjust_beta(bimfile = bimfile, summary_stats = summary_stats, output = output ) 
