args <- commandArgs(TRUE)
#par(mar = c(2.5, 2.5, 1.6, 1.1), mgp = c(1.5, 0.5, 0), mfrow=c(1,1))
library(VennDiagram)
a = read.table(args[1])[,1]
b = read.table(args[2])[,1]

x = list(a,b)
names(x) = args[3:4]
venn.diagram(x=x, file = paste(args[5], "tiff", sep="."), cex = 1.5, cat.cex = 1.5,cat.pos = 0, fill = c("red", "green"))
