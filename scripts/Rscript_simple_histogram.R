args = commandArgs(trailingOnly=TRUE)
mydata <- read.table(args[1], sep="\t", header = FALSE)
colnames(mydata) <- c("REF_proportion")

pdf(args[2])
hist(mydata$REF_proportion, breaks = 50)
dev.off()
