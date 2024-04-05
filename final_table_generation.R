setwd("/Users/minchunchen/lab/data/isoseq/methy")
temp = list.files(pattern="\\.tsv$")
myfiles = lapply(temp, read.delim)
tumour <- myfiles[[2]][,c('seg_id',"seg_name",'X04PM0832.tumour.hs1.aligned_m_methfrac')] 
n <- c(4,6,8,10,12,14,16,18,20,22,24,26)
for (i in n) {
  tumour <- cbind(tumour,myfiles[[i]][,colnames(myfiles[[i]])[10]])
}
colnames(tumour) <- c("seg_id","seg_name","04PM0832.tumour",'06PM0981.tumour','07PM0575.tumour','10PM1163.tumour','10PM1669.tumour','10PM2283.tumour','10PM2356.tumour','12PM0142.tumour','12PM0832.tumour','12PM1049.tumour','13PM0352.tumour','14PM0263.tumour','21DSRB006.tumour')

normal <- myfiles[[1]][,c('seg_id',"seg_name",'X04PM0832.normal.hs1.aligned_m_methfrac')] 
n <- c(3,5,7,9,11,13,15,17,19,21,23,25)
for (i in n) {
  normal <- cbind(normal,myfiles[[i]][,colnames(myfiles[[i]])[10]])
}
colnames(normal) <- c("seg_id","seg_name","04PM0832.normal",'06PM0981.normal','07PM0575.normal','10PM1163.normal','10PM1669.normal','10PM2283.normal','10PM2356.normal','12PM0142.normal','12PM0832.normal','12PM1049.normal','13PM0352.normal','14PM0263.normal','21DSRB006.normal')
tumour$tumour_average_meth <- rowMeans(tumour[,3:15],na.rm = TRUE)
normal$normal_average_meth <- rowMeans(normal[,3:15],na.rm = TRUE)
##read in overlap data
load("~/lab/data/isoseq/isoTe.RData")
setwd("/Users/minchunchen/lab/data/isoseq")
First_ov <- read.delim('First_ov.txt',check.names = FALSE)
##gene name
colnames(First_ov)
First_ov$GENEID <- mcols(se.tvn)$GENEID[match(First_ov$txId,mcols(se.tvn)$TXNAME)]
gene_gtf <- fread("/Users/minchunchen/lab/data/isoseq/Homo_sapiens-GCA_009914755.4-2022_07-genes_name.gtf",sep = '\t',header = FALSE)
setnames(gene_gtf, names(gene_gtf), c("chr","source","type","start","end","score","strand","phase","attributes","gene_id","gene_name") )
refe <- gene_gtf[,c("gene_id","gene_name")]
First_ov$gene <- refe$gene_name[match(First_ov$GENEID,refe$gene_id)]
colnames(df_deg)
First_ov$log2FoldChange <- df_deg$log2FoldChange[match(First_ov$txId,rownames(df_deg))]
First_ov$padj <- df_deg$padj[match(First_ov$txId,rownames(df_deg))]
colnames(normal)
First_ov$tumour_average_meth <- tumour$tumour_average_meth[match(First_ov$rep_id,tumour$seg_name)]
First_ov$normal_average_meth <- normal$normal_average_meth[match(First_ov$rep_id,normal$seg_name)]
colnames(isoclassification)
iso_ratio<- read.delim("Novel_transcripts_ratio.txt",check.names=FALSE)
head(iso_ratio)
First_ov$Iso_ratio <- iso_ratio$Iso_ratio[match(First_ov$txId,rownames(iso_ratio))]
First_ov$length <- isoclassification$length[match(First_ov$txId,isoclassification$isoform)]
First_ov$coding <- isoclassification$coding[match(First_ov$txId,isoclassification$isoform)]
qctranscripts <-read.delim('isofilter_inclusion-list.txt',header = FALSE)
fst_exo <- as.data.frame(first_exByTx)
First_ov$fist_exo_start <- fst_exo$start[match(First_ov$txId,rownames(fst_exo))]
First_ov$fist_exo_end <- fst_exo$end[match(First_ov$txId,rownames(fst_exo))]
head(First_ov)
trancrip <- transcripts(txdb, columns=c("tx_id", "tx_name"), filter=NULL, use.names=FALSE)
as.data.frame(trancrip)$width
tabel <- data.table(tx_name = trancrip$tx_name,
                       range=as.character(trancrip),
                       pre_mRNA = as.data.frame(trancrip)$width)
table <- as.data.frame(tabel)
First_ov <-  as.data.table(First_ov)
trancrip[match(First_ov$txId,trancrip$tx_name)]
First_ov$transcript_coor <- as.character(range(trancrip[match(First_ov$txId,trancrip$tx_name)]))
colnames(data) <- "coor"
First_ov <- read.delim("novel_transcript_first_exon.txt",check.names = FALSE)
First_ov$transcript_coor <- data$coor
export(trancrip, "transcript_first_exo.bed")
transcrp$coo <- apply(trancrip,1, paste0(seqnames,":",ranges(trancrip)))

write.table(First_ov, file="novel_transcript_first_exon.txt", sep = "\t", quote = F, row.names =F, col.names = T)
write.csv(First_ov, file="novel_transcript_first_exon.csv")
First_ov$transcript_coor<- table$range[match(First_ov$txId,table$tx_name)]
First_ov$pre_mRNA<- table$pre_mRNA[match(First_ov$txId,table$tx_name)]
First_ov$fist_exo_start <- NULL
First_ov$fist_exo_end <- NULL
colnames(First_ov)
write.table(First_ov, file="novel_transcript_first_exon.txt", sep = "\t", quote = F, row.names =F, col.names = T)
write.csv(First_ov, file="novel_transcript_first_exon.csv")