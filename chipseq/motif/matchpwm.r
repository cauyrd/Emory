library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg18)	
DIR = "/compbio/data/motif/human-mouse/"
motiflist = read.table(paste(DIR, "matrixlist.txt", sep =""), colClasses = "character")[,1]


induce = read.table("induce.bed", colClasses=c("character", "numeric", "numeric"))
repress = read.table("repress.bed", colClasses=c("character", "numeric", "numeric"))
count.induce = rep(0, length(motiflist))
count.repress = rep(0, length(motiflist))
chrnames= unique(induce[,1])
#chrnames = paste("chr", c(1:22,"X","Y"), sep ="")

#matchpwm <- function (peak, count, chrnames = unique(induce[,1])) {
for (i in 1:length(motiflist)) {
	pwm = t(as.matrix(read.table(paste(DIR, motiflist[i], sep ="/"))))
	rownames(pwm)=c("A","C","G","T")
	for (chr in chrnames) {
		seq.induce = Views(unmasked(Hsapiens[[chr]]), induce[induce[,1] == chr,2], induce[induce[,1] == chr,3])
		count.induce[i] = count.induce[i] + countPWM(pwm, seq.induce, min.score="80%") + countPWM(reverseComplement(pwm), seq.induce, min.score="80%")
		seq.repress = Views(unmasked(Hsapiens[[chr]]), repress[repress[,1] == chr,2], repress[repress[,1] == chr,3])
		count.repress[i] = count.repress[i] + countPWM(pwm, seq.repress, min.score="80%") + countPWM(reverseComplement(pwm), seq.repress, min.score="80%")
    }
}
#}
#   count.induce = rbind(count.induce, count1)
#   count.repress = rbind(count.repress, count2)

#matchpwm(induce, count.induce)
#matchpwm(repress, count.repress)
write.table(data.frame(motiflist, count.induce, count.repress), "count.txt", row.names=F, col.names=c("motif", "induce", "repress"), quote=F, sep="\t")
