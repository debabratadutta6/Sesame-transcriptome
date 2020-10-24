winDialog("ok", "Please select Stringtie Assembled Transcripts file")
options(stringsAsFactors=FALSE)
data1 <- choose.files(default = "", caption = "Select input datafile",
             multi = TRUE, filters = Filters,
             index = nrow(Filters))

data11 <- read.csv(data1, header=F, sep=" ")
data11 <- data11[-1,]
data11 <- data11[-1,]

winDialog("ok", "Please select Stringtie Gene Count file")

data2 <- choose.files(default = "", caption = "Select input datafile",
             multi = TRUE, filters = Filters,
             index = nrow(Filters))

data21 <- read.csv(data2, header=TRUE, sep="\t")

winDialog("ok", "Please select Stringtie Transcript Count file")

data3 <- choose.files(default = "", caption = "Select input datafile",
             multi = TRUE, filters = Filters,
             index = nrow(Filters))

data31 <- read.csv(data3, header=TRUE, sep="\t")

gene <- data11$V2
gene <- gsub(";", "", gene)
transcript <- data11$V4
transcript <- gsub(";", "", transcript)
locus1 <- data11$V6
locus1 <- gsub(";", "", locus1)
locus2 <- data11$V8
locus2 <- gsub(";", "", locus2)

num1 <- as.numeric(locus1)
num2 <- as.numeric(locus2)

bind <- as.data.frame(cbind (gene, transcript, locus1, locus2, num1, num2))

x <- nrow(bind)

bind$seqs <- 1:x
bind$isna1 <- is.na(bind$num1)
bind$isna2 <- is.na(bind$num2)

bind2 <- as.data.frame(bind)
sub1 <- subset(bind2, bind2$isna1==TRUE)
sub2 <- subset(bind2, bind2$isna2==TRUE)

cols1 <- c(1,2,3,7)
sub12 <- sub1[,cols1]
colnames(sub12) <- c("gene", "transcript", "locus", "seqs")

cols2 <- c(1,2,4,7)
sub22 <- sub2[,cols2]
colnames(sub22) <- c("gene", "transcript", "locus", "seqs")

finalish <- rbind(sub12, sub22) 
final <- finalish[order(finalish$seqs),]
final <- final[,-4]

dup <- duplicated(final)
final$dup <- dup
sub3 <- subset(final, final$dup==FALSE)


#############################################################################
sub4 <- sub3
colnames(sub4)[1] <- c("gene_id")

merge1 <- merge(data21, sub4, by.x="gene_id")
merge2 <- merge1
merge2$gene_id <- merge2$locus
new1 <- merge2[,1:2]

#dup1 <- duplicated(new1)
#new1$dup <- dup1
#sub1 <- subset(new1, new1$dup==FALSE)
sub1 <- new1

rows <- nrow(sub1)
sub1$row <- 1:rows
sub1$gene_id <- paste(sub1$row,"_",sub1$gene_id, sep="")
new11 <- sub1[,-c(3:4)]

write.table(new11, "Adjusted Gene Counts.tabular", quote=F, sep="\t", row.names=F)

#############################################################################
sub5 <- sub3
colnames(sub5)[2] <- c("transcript_id")

merge3 <- merge(data31, sub5, by.x="transcript_id")
merge4 <- merge3
merge4$transcript_id <- merge4$locus
new2 <- merge4[,1:2]

#dup2 <- duplicated(new2)
#new2$dup <- dup2
#sub2 <- subset(new2, new2$dup==FALSE)
sub2 <- new2

rows <- nrow(sub2)
sub2$row <- 1:rows
sub2$transcript_id <- paste(sub2$row,"_",sub2$transcript_id, sep="")
new21 <- sub2[,-3]

write.table(new21, "Adjusted Transcript Counts.tabular", quote=F, sep="\t", row.names=F)





