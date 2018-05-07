#Normalising genes
Counts <- read.delim('***PATH***/RNA-Seq-counts.txt', header = TRUE, skip =1, row.names = 1)

CountsNorm <- Counts / colSums(Counts)[col(Counts)]

library("xlsx")
write.xlsx(CountsNorm, "***PATH***/countsNorm.xlsx")
