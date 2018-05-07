#Normalising genes
Counts <- read.delim('/home/shane/Documents/Course_11/data/RNA-Seq-counts.txt', header = TRUE, skip =1, row.names = 1)

CountsNorm <- Counts / colSums(Counts)[col(Counts)]

library("xlsx")
write.xlsx(CountsNorm, "/home/shane/Documents/Course_11/data/countsNorm.xlsx")