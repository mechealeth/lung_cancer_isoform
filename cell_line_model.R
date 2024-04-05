##cancer cell specific transcipts
library(DESeq2)
#set path
setwd("/Users/minchunchen/lab/data/cell_line")
se.multi <- readRDS('/Users/minchunchen/lab/data/cell_line/bambucell_line.rds')
isoclassification <- read.delim("isoseq_classification.txt")
qctranscripts <-read.delim('isofilter_inclusion-list.txt',header = FALSE)
colnames(se.multi) <- sapply(strsplit(colnames(se.multi), "\\."), "[", 1)
head(qctranscripts)
#qc and subset novel transcripts
se.multi.qc <- se.multi[mcols(se.multi)$TXNAME %in% qctranscripts$V1,]
cell <- sapply(strsplit(rownames(colData(se.multi.qc)), "_"), "[", 1)
condition <-  sapply(strsplit(rownames(colData(se.multi.qc)), "_"), "[", 2)
condition <- c("DMSO","DMSO","DNMTi","DNMTi","DNMTi","DNMTi","DNMTi","DNMTi","Ctrl","DMSO","DMSO","DNMTi","DNMTi","DNMTi","DNMTi","DNMTi","Ctrl")
timepoint <- sapply(strsplit(rownames(colData(se.multi.qc)), "_"), "[", 3)
timepoint <- c("D0","D0","D2","D2","D4","D4","D5","D5","D0","D0","D0","D2","D2","D4","D4","D4","D0")
sample <- data.frame(cell,condition,timepoint,row.names=rownames(colData(se.multi.qc)))
head(assays(se.multi.qc)$counts[,c(1,2,9,10,11,17)])
dds <- DESeqDataSetFromMatrix(round(assays(se.multi.qc)$counts[,c(1,2,9,10,11,17)]), colData = sample[c(1,2,9,10,11,17),],design = ~cell)
dds.deseq <- DESeq(dds)
deRes <- DESeq2::results(dds.deseq,contrast = c("cell", "HARA", "NL20"), independentFiltering = FALSE)
head(deRes[order(deRes$padj), ])
summary(deRes)
DESeq2::plotMA(deRes,alpha = 0.01)
deg <- deRes[deRes$padj < 0.05& !is.na(deRes$padj)& deRes$log2FoldChange > 1, ]
df_deg <- as.data.frame(deg)
table(df_deg$log2FoldChange > 0)
novel_dte <- df_deg[rownames(df_deg) %like% "BambuTx" &df_deg$log2FoldChange > 2, ]
##overlaping isoform with TE
comp_gtf <- fread("/Users/minchunchen/lab/data/cell_line/cell_tumorcmp.annotated.gtf",sep = '\t',header = FALSE)
head(comp_gtf)
setnames(comp_gtf, names(comp_gtf), c("chr","source","type","start","end","score","strand","phase","attributes") )
extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- unlist(strsplit(gtf_attributes, " "))
  if(att_of_interest %in% att){
    return(gsub("\"|;","", att[which(att %in% att_of_interest)+1]))
  }else{
    return(NA)}
}
comp_gtf$transcript_id <- unlist(lapply(comp_gtf$attributes, extract_attributes, "transcript_id"))
comp_gtf$gene_id <- unlist(lapply(comp_gtf$attributes, extract_attributes, "gene_id"))
comp_gtf$cmp_ref <- unlist(lapply(comp_gtf$attributes, extract_attributes, "cmp_ref"))
comp_gtf$class_code <- unlist(lapply(comp_gtf$attributes, extract_attributes, "class_code"))
comp_gtf <- comp_gtf[comp_gtf$class_code == '=' & comp_gtf$transcript_id %like% "BambuTx",]
comp_gtf$attributes <- NULL
##cancer cell line specific DTE overlap with Tumor sample
dte_ov <- novel_dte[rownames(novel_dte)%in%comp_gtf$transcript_id,]
###time series analysis
df1 <-assays(se.multi.qc)$counts[,c(-7,-8,-9,-17)]
colnames(df1)
col_oder <- c("NL20_DMSO_R1","NL20_DMSO_R2","NL20_DNMTi_D2_R1","NL20_DNMTi_D2_R2","NL20_DNMTi_D4_R1","NL20_DNMTi_D4_R2","NL20_DNMTi_D4_R3","HARA_DMSO_R1","HARA_DMSO_R2","HARA_DNMTi_D2_R1","HARA_DNMTi_D2_R2",
              "HARA_DNMTi_D4_R1","HARA_DNMTi_D4_R2" )
cell <- sapply(strsplit(colnames(df1), "_"), "[", 1)
replicate <- c("R1","R2","R1","R2","R1","R2","R3","R1","R2","R1","R2","R1","R2")
timepoint <- c("0","0","2","2","4","4","4","0","0","2","2","4","4")
colData <- data.frame(cell,timepoint,replicate,row.names=colnames(df1))
ddsTC <- DESeqDataSetFromMatrix(countData = round(df1), colData = colData ,design = ~ cell + timepoint + cell:timepoint)
cell <- factor(colData$cell)
ddsTC$cell=relevel(cell, "HARA")
ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ cell + timepoint)
resTC <- results(ddsTC)
res <- as.data.frame(ddsTC)
head(resTC[order(resTC$padj),], 7)
p <- plotCounts(ddsTC,"BambuTx4", 
                   intgroup = c("timepoint","cell"), returnData = TRUE)
p$timepoint <- as.numeric(as.character(p$timepoint))
ggplot(p,aes(x = timepoint, y = count, color = cell, group = cell)) + 
  geom_point() + stat_summary(fun.y=mean, geom="line") +
  scale_y_log10()
resultsNames(ddsTC)
res30 <- results(ddsTC, name="cellNL20.timepoint4", test="Wald")
head(res30[order(res30$padj), ])
summary(res30)
DESeq2::plotMA(res30)
sig_n2<- res30[res30$padj < 0.05 & !is.na(res30$padj)& res30$log2FoldChange > 1, ]
sig_n2 <- as.data.frame(sig_n2)
###NL20
colnames(assays(se.multi.qc)$counts)
df1 <-assays(se.multi.qc)$counts[,10:16]
colnames(df1)
col_oder <- c("NL20_DMSO_R1","NL20_DMSO_R2","NL20_DNMTi_D2_R1","NL20_DNMTi_D2_R2","NL20_DNMTi_D4_R1","NL20_DNMTi_D4_R2","NL20_DNMTi_D4_R3")
###"HARA_DMSO_R1","HARA_DMSO_R2","HARA_DNMTi_D2_R1","HARA_DNMTi_D2_R2","HARA_DNMTi_D4_R1","HARA_DNMTi_D4_R2"###
cell <- sapply(strsplit(colnames(df1), "_"), "[", 1)
replicate <- c("R1","R2","R1","R2","R1","R2","R3")
timepoint <- c("0","0","2","2","4","4","4")
colData <- data.frame(cell,timepoint,replicate,row.names=colnames(df1))
dds <- DESeqDataSetFromMatrix(countData = round(df1) , colData, design = ~ timepoint)
dds.deseq <- DESeq(dds)
resultsNames(dds.deseq)
deGeneRes <- DESeq2::results(dds.deseq,name = "timepoint_4_vs_0")
head(deGeneRes[order(deGeneRes$padj), ])
summary(deGeneRes)
DESeq2::plotMA(deGeneRes)
deg2<- deGeneRes[deGeneRes$padj < 0.05 & !is.na(deGeneRes$padj)& deGeneRes$log2FoldChange > 1, ]
deg2 <- as.data.frame(deg2)
NL20_sensitive <- rbind(deg1,deg2)
NL20_sensitive <- NL20_sensitive[rownames(NL20_sensitive) %like%"BambuTx",]
###HARA
colnames(assays(se.multi.qc)$counts)
df1 <-assays(se.multi.qc)$counts[,1:6]
colnames(df1)
col_oder <- c("HARA_DMSO_R1","HARA_DMSO_R2","HARA_DNMTi_D2_R1","HARA_DNMTi_D2_R2","HARA_DNMTi_D4_R1","HARA_DNMTi_D4_R2")
cell <- sapply(strsplit(colnames(df1), "_"), "[", 1)
replicate <- c("R1","R2","R1","R2","R1","R2")
timepoint <- c("0","0","2","2","4","4")
colData <- data.frame(cell,timepoint,replicate,row.names=colnames(df1))
dds <- DESeqDataSetFromMatrix(countData = round(df1) , colData, design = ~ timepoint)
dds.deseq <- DESeq(dds)
resultsNames(dds.deseq)
deGeneRes <- DESeq2::results(dds.deseq,name = "timepoint_4_vs_0")
head(deGeneRes[order(deGeneRes$padj), ])
summary(deGeneRes)
DESeq2::plotMA(deGeneRes)
deg4<- deGeneRes[deGeneRes$padj < 0.05 & !is.na(deGeneRes$padj)& deGeneRes$log2FoldChange > 1, ]
deg4 <- as.data.frame(deg4)
HARA_sensitive <- rbind(deg3,deg4)
HARA_sensitive <- HARA_sensitive[rownames(HARA_sensitive) %like%"BambuTx",]
###cancer_specific_novel
colnames(assays(se.multi.qc)$counts)
df1 <-assays(se.multi.qc)$counts[,c(1,2,9,10,11,17)]
colnames(df1)
col_oder <- c("HARA_DMSO_R1","HARA_DMSO_R2","HARA_DNMTi_D2_R1","HARA_DNMTi_D2_R2","HARA_DNMTi_D4_R1","HARA_DNMTi_D4_R2")
cell <- sapply(strsplit(colnames(df1), "_"), "[", 1)
colData <- data.frame(cell,row.names=colnames(df1))
dds <- DESeqDataSetFromMatrix(countData = round(df1) , colData, design = ~ cell)
cell <- factor(colData$cell)
dds$cell=relevel(cell, "NL20")
dds.deseq <- DESeq(dds)
resultsNames(dds.deseq)
deGeneRes <- DESeq2::results(dds.deseq,name = "cell_HARA_vs_NL20")
head(deGeneRes[order(deGeneRes$padj), ])
summary(deGeneRes)
DESeq2::plotMA(deGeneRes)
deg<- deGeneRes[deGeneRes$padj < 0.001 & !is.na(deGeneRes$padj)& deGeneRes$log2FoldChange > 2, ]
deg <- as.data.frame(deg)
cancer_specific <- deg[rownames(deg) %like%"BambuTx",]
##cancer specific meth sensitive
NL20_sensitive <- NL20_sensitive[!rownames(NL20_sensitive) %in% rownames(HARA_sensitive),]
NL20_sensitive <- NL20_sensitive[rownames(NL20_sensitive) %in% rownames(cancer_specific),]

##overlap with TE
library(data.table)
library(stringr)
gtf <- fread("/Users/minchunchen/lab/data/isoseq/T2T_CHM13_v2_rmsk_TE_chr.gtf",sep = '\t',header = FALSE)
head(gtf)
setnames(gtf, names(gtf), c("chr","source","type","start","end","score","strand","phase","attributes") )
extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- unlist(strsplit(gtf_attributes, " "))
  if(att_of_interest %in% att){
    return(gsub("\"|;","", att[which(att %in% att_of_interest)+1]))
  }else{
    return(NA)}
}
gtf$gene_id <- unlist(lapply(gtf$attributes, extract_attributes, "gene_id"))
gtf$family_id <- unlist(lapply(gtf$attributes, extract_attributes, "family_id"))
gtf$class_id <- unlist(lapply(gtf$attributes, extract_attributes, "class_id"))
rm(gtf)
gr_list <- makeGRangesFromDataFrame(gtf, keep.extra.columns=TRUE,
                                    seqnames.field="chr",
                                    start.field="start",
                                    end.field="end",
                                    strand.field="strand",
                                    starts.in.df.are.0based=TRUE)
names(gr_list) <- paste0('rep',1:length(gr_list))

grangesList <- rowRanges(se.multi.qc[mcols(se.multi.qc)$TXNAME%in%rownames(NL20_sensitive),])
exons_granges <- unlist(grangesList)
ov <- findOverlaps(exons_granges,gr_list )
p <- Pairs(exons_granges, gr_list, hits=ov)
hitIntersect <- pintersect(p)
rm(p)
gc()
qHits <- queryHits(ov)
sHits <- subjectHits(ov)
rm(ov)
gc()
overlapWidth <- width(hitIntersect)
rm(hitIntersect)
gc()

ov <- data.table(txId = names(exons_granges)[qHits],
                 exon_rank = exons_granges[qHits]$exon_rank,
                 anno_status = exons_granges[qHits]$anno_status,
                 nNovelExon = exons_granges[qHits]$nnovel,
                 nAnnotatedExon = exons_granges[qHits]$nannotated,
                 rep_id = names(gr_list)[sHits],
                 strand = as.character(strand(gr_list[sHits])),
                 rep_name = gr_list[sHits]$gene_id,
                 rep_class = gr_list[sHits]$class_id,
                 rep_family = gr_list[sHits]$family_id,
                 tx_width = width(exons_granges)[qHits],
                 rep_width = width(gr_list)[sHits],
                 ov_width = overlapWidth)

##export bed file for methylation
inter <- ov[ov$exon_rank == "1",]
gr_list
names(gr_list)
gr <- gr_list[names(gr_list) %in% inter$rep_id,]
gr
export(gr, "ov.bed")

###read methlation data
meth <- read.table(file = 'ov.NL20.hs1.aligned.mod.segmeth.tsv', sep = '\t', header = TRUE)
rownames(meth) <- meth$seg_id
meth$seg_id <- NULL
rownames(meth) <- c("chr12:121989307-121990791_LTR12C","chr3:33448197-33448761_HERV9NC-int","chr3:33448762-33455009_SVA_F","chr7:74483462-74485480_THE1B","chrX:1985724-1986064_SVA_F","chrY:2049625-2049965_THE1B")
colnames(meth) <- c("NL20","NL20_DNMTi_D4","HARA","HARA_DNMTi_D4")
p <- as.data.frame(t(meth[4,]))
p$day <- c("0","4","0","4")
p$cell <- c("NL20","NL20","HARA","HARA")
colnames(p) <- c("meth_Frac","day","cell")
ggplot(p,aes(x = day, y = meth_Frac, color = cell, group = cell)) + 
  geom_point() + stat_summary(fun.y=mean, geom="line") +
  scale_y_log10() +labs(x="DNMTi_TimePoint", y="meth_Frac", title="chr7:74483462-74485480_THE1B Methylation ")
inter$rep_id <- NULL
inter$coding <- isoclassification[]
inter$coding  <- isoclassification[match(inter$txId,isoclassification$isoform),]$coding
write.csv(inter,"HARA_specific_sen.csv")
