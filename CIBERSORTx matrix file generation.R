#As input, CIBERSORTx requires a matrix file containing the normalized gene counts, formatted in a data frame.

setwd("C:/Users/user/folder/")

data<- read.csv("C:/Users/user/folder/countsfile.csv")
write.table(as.data.frame(data), "countsfile.tsv", quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)
head(data)
