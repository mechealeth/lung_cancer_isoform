#load lib
library(bambu)
require(GenomicAlignments) 
require(AnnotationDbi)
require(data.table)
require(readxl)
require(ggplot2)
require(RColorBrewer)
require(gridExtra)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(DESeq2)
library(Gviz)
library(GenomicRanges)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
#set path and load data
setwd("/Users/minchunchen/lab/data/isoseq")
out_data_dir <- "/Users/minchunchen/lab/data/isoseq/bambu"
se.multi <- readRDS('/Users/minchunchen/lab/data/isoseq/bambu_NDR0.155rebambu.rds')
qctranscripts <-read.delim('isofilter_inclusion-list.txt',header = FALSE)
isoclassification <- read.delim("isoseq_classification.txt")
colnames(se.multi) <- sapply(strsplit(colnames(se.multi), "_"), "[", 1)
#qc and subset novel transcripts
se.multi.qc <- se.multi[mcols(se.multi)$TXNAME %in% qctranscripts$V1,]
se.multi.novel = se.multi.qc[mcols(se.multi.qc)$TXNAME %like% "Bambu",]
mcols(se.multi.novel)
colnames((isoclassification))
mcols(se.multi.novel)$structural_category <- isoclassification[match(mcols(se.multi.novel)$TXNAME,isoclassification$isoform),]$structural_category
unique(mcols(se.multi.novel)$structural_category)
readCount <- read.table("transcript_count_matrix.csv",sep = ',', header = TRUE,row.names = 1,check.names = FALSE)
sampleInfo<- as.data.frame(colnames(readCount))
colnames(sampleInfo) <- "sample_id"
out <- stringr::str_split_fixed(sampleInfo$sample_id,"\\.",2)
colnames(out) <- c("patients",'groups')
sample <- cbind(sampleInfo,out)
rownames(sample) <- sample$sample_id
sample$sample_id <- NULL
sample$groups<- factor(sample$groups, levels = c('tumour','normal'))
###
dds <- DESeqDataSetFromMatrix(countData = readCount, colData = sample, design = ~ groups)
dim(dds)
rowSums(counts(dds))
keep <- rowSums(counts(dds) >= 5) >= 3
table(keep)
dds1 <- dds[keep,]
dds1 <- DESeq(dds1)
normalized_counts <- as.data.frame(counts(dds1,normalized=TRUE))
se.multi.qc.f <- se.multi.qc[mcols(se.multi.qc)$TXNAME%in% rownames(normalized_counts),]
nrows <- 149053
ncols <- 28
rn <- names(rowRanges(se.multi.qc.f))
library(dplyr)
rn2 <- c("04PM0832.normal","08PM0370.normal","10PM2283.normal","12PM0832.normal","14PM0263.normal","06PM0981.normal","10PM1163.normal","10PM2356.normal","12PM1049.normal","21DSRB006.normal","07PM0575.normal","10PM1669.normal","12PM0142.normal","13PM0352.normal",
         "04PM0832.tumour","08PM0370.tumour","10PM2283.tumour","12PM0832.tumour","14PM0263.tumour","06PM0981.tumour","10PM1163.tumour","10PM2356.tumour","12PM1049.tumour","21DSRB006.tumour","07PM0575.tumour","10PM1669.tumour","12PM0142.tumour","13PM0352.tumour" )
counts1 <- as.data.frame(normalized_counts[,match(rn2,colnames(normalized_counts))])
cpm <- counts1 %>% arrange(factor(rownames(counts1), levels = rn)) %>% as.data.frame()
readcount <- readCount[rownames(readCount)%in%rownames(normalized_counts),]
readcount <- readcount[,match(rn2,colnames(readcount))]
readcount <- readcount %>% arrange(factor(rownames(counts1), levels = rn)) %>% as.data.frame()
rowRanges <- rowRanges(se.multi.qc.f)
sampleInfo<- as.data.frame(colnames(readcount))
colnames(sampleInfo) <- "sample_id"
out <- stringr::str_split_fixed(sampleInfo$sample_id,"\\.",2)
colnames(out) <- c("patients",'groups')
sample <- cbind(sampleInfo,out)
colData <- sample
se.tvn <- SummarizedExperiment(assays=list(counts=readcount,CPM=cpm),rowRanges=rowRanges, colData=colData)
rm(se.multi)

isoclassification <- read.delim("isoseq_classification.txt")
colnames((isoclassification))
mcols(se.tvn)$structural_category<- isoclassification[match(mcols(se.tvn)$TXNAME,isoclassification$isoform),]$structural_category
gen_gtf <- fread("/Users/minchunchen/lab/data/isoseq/Homo_sapiens-GCA_009914755.4-2022_07-genes_name.gtf",sep = '\t',header = FALSE)
setnames(gen_gtf, names(gen_gtf), c("chr","source","type","start","end","score","strand","phase","attributes","gene_id","gene_name") )
refe <- gen_gtf[,c("gene_id","gene_name")]
mcols(se.tvn)$gene <- refe$gene_name[match(mcols(se.tvn)$GENEID,refe$gene_id)]
mcols(se.tvn)$gene[mcols(se.tvn)$gene == ''] <- NA 
rm(refe)

###DTE
library(DESeq2)
dds <- DESeqDataSetFromMatrix(round(assays(se.tvn)$counts), colData = colData(se.tvn),
                              design = ~groups)
dds.deseq <- DESeq(dds)
deRes <- DESeq2::results(dds.deseq,contrast = c("groups", "tumour", "normal"), independentFiltering = FALSE)
head(deRes[order(deRes$padj), ])
summary(deRes)
DESeq2::plotMA(deRes,alpha = 0.01)
##PCA
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("groups", "patients"), returnData=TRUE,ntop = 4000)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=patients, shape=groups)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
##
library(FactoMineR)
vsd <-counts(dds.deseq,normalized=TRUE)
vsd_pca <- PCA(t(vsd), ncp = 2, scale.unit = TRUE, graph = FALSE)
pca_sample <- data.frame(vsd_pca $ind$coord[ ,1:2])
pca_eig1 <- round(vsd_pca$eig[1,2], 2)
pca_eig2 <- round(vsd_pca$eig[2,2],2 )
out <- stringr::str_split_fixed(rownames(pca_sample),"\\.",2)
colnames(out) <- c('patients','groups')
pca_sample <- cbind(pca_sample, out)
library(ggplot2)
p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = groups), size = 5) +  
  scale_color_manual(values = c("#F8766D",  "#00C094")) + 
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +  
  labs(x =  paste('PCA1:', pca_eig1, '%'), y = paste('PCA2:', pca_eig2, '%'), color = '') 
library(ggrepel)
p <- p +  geom_text_repel(aes(label = rownames(pca_sample)), size = 3, show.legend = FALSE, 
                          box.padding = unit(0.5, 'lines'))
print(p)
p + stat_ellipse(aes(color = groups), level = 0.95, show.legend = FALSE)

##DTE result 
deg <- deRes[deRes$padj < 0.001& !is.na(deRes$padj)& deRes$log2FoldChange > 2, ]
df_deg <- as.data.frame(deg)
table(df_deg$log2FoldChange > 0)
novel_dte <- df_deg[rownames(df_deg) %like% "BambuTx" &df_deg$log2FoldChange > 2, ]
iso_ratio<- read.delim("Novel_transcripts_ratio.txt",check.names=FALSE)
refe$gene_name[match(dtu$gene_id,refe$gene_id)]
novel_dte$iso_ratio <- iso_ratio[match(rownames(novel_dte),iso_ratio$TXNAME),]$Iso_ratio
## define cancer specifc isoform
novel_dte <- novel_dte[novel_dte$iso_ratio > 0.1,]
novel_expr <- readcount[rownames(readcount)%in%rownames(novel_dte),]
colnames(novel_expr)
keep <- rowSums(novel_expr[,15:28] >= 5) >= 3
table(keep)
novel_expr <- novel_expr[keep,]
table <- as.data.frame(table(geneTxTable_extended$newTxClassAggregated))
colnames(table) <- c('Var1','number')
pct <- round(100*table$number/sum(table$number))
pie(table$number,labels = paste(table$Var1, sep = " ", pct, "%"), col = rainbow(length(table$number)), main = "structural category of 2379 Novel transcripts ")

###cancer specific  isoform characteristic
######TE intersection#######
##prepare TE reference
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

###
trancrip <- transcripts(txdb, columns=c("tx_id", "tx_name"), filter=NULL, use.names=FALSE)

txdb <- makeTxDbFromGFF(file="Homo_sapiens-GCA_009914755.4-2022_07-genes_chr.gtf")
rowRanges(txdb)
anno_exByTx <- exonsBy(txdb, 'tx', use.names = TRUE)
txLengths <- transcriptLengths(txdb)
txLengths.tbldf <- data.table(txLengths)
txLengths.tbldf [, `:=`(nisoform = length(unique(tx_id))),
                 by = gene_id]
geneTxTable <- txLengths.tbldf[,.(tx_name, gene_id, nisoform)]
setnames(geneTxTable, 'gene_id', 'gene_name')
grangesList <- rowRanges(se.tvn[mcols(se.tvn)$TXNAME%in%rownames(novel_expr),])
exons_granges <- unlist(grangesList)
anno_exons <- unlist(anno_exByTx)
gr_list <- makeGRangesFromDataFrame(gtf, keep.extra.columns=TRUE,
                                    seqnames.field="chr",
                                    start.field="start",
                                    end.field="end",
                                    strand.field="strand",
                                    starts.in.df.are.0based=TRUE)
names(gr_list) <- paste0('rep',1:length(gr_list))


###
ov <- findOverlaps(exons_granges, anno_exons, type = "any")
exons_granges$anno_status <- "novel"
exons_granges[queryHits(ov)]$anno_status <- "annotated"
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


nnovel <- data.table(table(names(exons_granges[exons_granges$anno_status=="novel"])))
nannotated <- data.table(table(names(exons_granges[exons_granges$anno_status=="annotated"])))
exons_granges$nnovel <- nnovel[match(names(exons_granges),V1)]$N
exons_granges$nannotated <- nannotated[match(names(exons_granges),V1)]$N  
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

ov[, repRatio := ov_width/tx_width]
ov[is.na(nNovelExon), nNovelExon := 0]
ov[is.na(nAnnotatedExon), nAnnotatedExon := 0]
isoTEratio_anno <- unique(ov[anno_status == "annotated", list(averageRepRatio = sum(repRatio)/nAnnotatedExon,
                                                              strand = strand,
                                                              rep_class = rep_class,
                                                              n_rep = length(rep_id)), by = list(txId, rep_name)])
isoTEratio_novel <- unique(ov[anno_status == "novel", list(averageRepRatio = sum(repRatio)/nNovelExon,
                                                           strand = strand,
                                                           rep_class = rep_class,
                                                           n_rep = length(rep_id)), by = list(txId, rep_name)])
isoTEratio_anno[, repRatioIso_all := sum(averageRepRatio), by = txId]
isoTEratio_novel[, repRatioIso_all := sum(averageRepRatio), by = txId]
isoTEratio_anno[, anno_status := "annotated"]
isoTEratio_novel[, anno_status := "novel"]
isoTEratio <- rbind(isoTEratio_anno, isoTEratio_novel)
isoTEratio[, repRatio_corrected:=ifelse(repRatioIso_all>1, 1, repRatioIso_all)]
isoTEratio_final <- isoTEratio[!(grepl("tx.",txId)&(anno_status == "annotated"))]
isoTEratio_final_wide <- dcast(isoTEratio_final, txId ~ rep_name, value.var = "averageRepRatio")
isoTEratio_final_wide[is.na(isoTEratio_final_wide)] <- 0
write.table(isoTEratio, file="cancer_associate_isoTEratio.txt",sep="\t", quote=F)
se.tvn_novel_dte <- se.tvn[mcols(se.tvn)$TXNAME%in%rownames(novel_expr),]
unique(mcols(se.tvn_novel_dte)$txClassDescription)
geneTxTable_extended <- as.data.table(rowData(se.tvn_novel_dte))
setnames(geneTxTable_extended, c("GENEID","TXNAME"),c("gene_name","tx_name"))
geneTxTable_extended[,nisoform:=length(unique(tx_name)), by = gene_name]
geneTxTable <- txLengths.tbldf[,.(tx_name, gene_id, nisoform)]
setnames(geneTxTable, 'gene_id', 'gene_name')
geneTxTable_extended <- geneTxTable[geneTxTable_extended, on = c("gene_name","tx_name")]
geneTxTable_extended[, newTxClassAggregated:=ifelse(grepl("newFirstExon",txClassDescription)&(grepl("newLastExon",txClassDescription)),"newFirstLastExon",
                                                    ifelse(grepl("newFirstExon",txClassDescription), "first exon",
                                                           ifelse(grepl("newLastExon",txClassDescription), "last exon", 
                                                                  ifelse(grepl("Junction|allNew",txClassDescription),"internal exon",
                                                                         ifelse(grepl("unspliced",txClassDescription),"unspliced",
                                                                                ifelse(grepl("newGene-",txClassDescription),"new gene",txClassDescription))))))]



novel_dte_ov <- novel_dte[rownames(novel_dte)%in%isoTEratio$txId,]

####
setnames(isoTEratio_final_wide, "txId","tx_name")
isoTEratio_final_wide <- geneTxTable_extended[isoTEratio_final_wide, on = "tx_name"]
setnames(isoTEratio,"txId","tx_name")
isoTEratio <- geneTxTable_extended[isoTEratio, on = "tx_name"]
unique(isoTEratio$newTxClassAggregated)
isoTEratio_anno
cutoff_value <- 0.5
canonical.labs <- c("Annotated", "Novel")
names(canonical.labs) <- c("annotated", "novel")
plot_tmp <- unique(isoTEratio[anno_status=="annotated"][,.(tx_name, repRatio_corrected,anno_status, newTxClassAggregated)])
p_anno <- ggplot(plot_tmp, aes(repRatio_corrected))+
  geom_histogram(aes(fill = newTxClassAggregated),binwidth = 0.02,col = "white",boundary = 0, position ="stack")+
  geom_vline(xintercept=cutoff_value, color = 'steelblue', linetype = "dashed")+
  scale_fill_manual(values = c("#1B9E77B3","lightblue","steelblue","darkseagreen1","forestgreen"), 
                    breaks = c("annotation","first exon","last exon","internal exon","newWithin"),
                    # label = c("Novel","Canonical"),
                    name = "Isoform type")+
  ggtitle("")+
  #coord_trans(y = "log10")+
  #facet_wrap(~anno_status,scales = "free_y", labeller = labeller(anno_status = canonical.labs))+
  #scale_y_log10()+
  xlab("Overlap percentage with repeats")+
  ylab("Number of transcripts")+
  theme_classic()

plot_tmp <- unique(isoTEratio[anno_status=="novel"][,.(tx_name, repRatio_corrected,anno_status, newTxClassAggregated)])
table(plot_tmp$newTxClassAggregated)
plot_tmp[, newTxClassAggregated_factor := factor(newTxClassAggregated, 
                                                 levels = rev(c("first exon","last exon","internal exon","new gene","newFirstLastExon")))]
p_novel <- ggplot(plot_tmp, aes(repRatio_corrected))+
  geom_histogram(aes(fill = newTxClassAggregated_factor),binwidth = 0.02,col = "white",boundary = 0, position ="stack")+
  geom_vline(xintercept=cutoff_value, color = 'steelblue', linetype = "dashed")+
  scale_fill_manual(values = c("#1B9E77B3","lightblue","steelblue","darkseagreen1","forestgreen"), 
                    breaks = c("newFirstLastExon","first exon","last exon","internal exon","new gene"),
                    # label = c("Novel","Canonical"),
                    name = "Isoform type")+
  ggtitle("")+
  #coord_trans(y = "log10")+
  #facet_wrap(~anno_status,scales = "free_y", labeller = labeller(anno_status = canonical.labs))+
  #scale_y_log10()+
  xlab("Overlap percentage with repeats")+
  ylab("Number of transcripts")+
  theme_classic()
##plotting
chr3:121634644-121634930
chr3:121633530-121633834
library(Gviz)
library(GenomicRanges)
mcols(se.tvn[mcols(se.tvn)$TXNAME == 'BambuTx1502',])
rowRanges(se.tvn[mcols(se.tvn)$GENEID == 'ENSG05220015048',])
plotBambu(se.tvn, type = "annotation", gene_id = "ENSG05220053473")
mcols(se.tvn[mcols(se.tvn)$GENEID == 'ENSG05220015048',])
as.data.frame(gr_list[names(gr_list) == "rep3683766",])
polt <- as.data.frame(rowRanges(se.tvn[mcols(se.tvn)$GENEID == 'ENSG05220053473',]))
colnames(polt) <- c('group','transcript','chromosome','start','end','width','strand','exon','exon_endRank')
erv <- data.frame(group = 1,
                  transcript = "AluJb",
                  chromosome = 'chr3',
                  start = 121634644,
                  end =121634930,
                  width =287,
                  strand ='-',
                  exon = 1,
                  exon_endRank =1)

novel_dte[rownames(novel_dte) %in% isoTEratio[isoTEratio$newTxClassAggregated=='new gene' &isoTEratio$coding =='coding',]$tx_name,]
erv <- data.frame(group = c(4,2,3,1),
                  transcript = c("AluSp","MER58A",'MER5A','MIR3'),
                  chromosome = c("chr6","chr6","chr6","chr6"),
                  start =c(14424027,14426247,14426479,14430158),
                  end =c(14424308,14426476,14426534,14430216),
                  width =c(282,230,56,59),
                  strand =c("+","+","+","+"),
                  exon = c(1,2,3,4),
                  exon_endRank =c(4,3,2,1))
data <- rbind(polt,erv)
grtrack <- GeneRegionTrack(polt[1:4,], transcriptAnnotation = "transcript", background.title = "brown",name ="Transcript")
grtrack1 <- GeneRegionTrack(erv,transcriptAnnotation = "transcript", background.title = "brown",name ="TE")
plotTracks(grtrack)
grtrack2<- GeneRegionTrack(erv, transcriptAnnotation = "transcript", background.title = "brown",name ="TE")
plotTracks(list(grtrack,grtrack1),collapseTranscripts = "longest")

ideoTrack <- IdeogramTrack(genome = "hs1", chromosome = "chr3")
plotTracks(ideoTrack, from = 85e6, to = 129e6)
gen <- genome(cpgIslands)
atrack <- AnnotationTrack(cpgIslands, name = "CpG")
plotTracks(list(ideoTrack,grtrack,grtrack1),collapseTranscripts = "longest")

switch <-as.data.frame(rowRanges(se.multi.qc[mcols(se.multi.qc)$TXNAME == 'ENST05220074466',]))[1,]
ov <- findOverlaps(gr_list,rowRanges(se.multi.qc[mcols(se.multi.qc)$TXNAME == 'ENST05220074466',]))
rowRanges(se.multi.qc[mcols(se.multi.qc)$GENEID == 'ENSG05220020076',])



se.tvn.dtu <-se.tvn[mcols(se.tvn)$TXNAME %in% rownames(novel_dte),]
mcols(se.tvn.dtu)
table(isoclassification[isoclassification$isoform %like% "BambuTx",]$coding)
table(isoclassification[isoclassification$isoform %in% rownames(novel_dte),]$coding)
normalized_counts <- counts(dds.deseq,normalized=TRUE)
colnames(isoclassification)
rowRanges(se.tvn.dtu)
tumor_novel_gene <- mcols(se.tvn.dtu)[mcols(se.tvn.dtu)$GENEID %like% "BambuGene",]$GENEID
mcols(se.tvn.dtu)[mcols(se.tvn.dtu)$GENEID%in%tumor_novel_gene,]
table(mcols(se.tvn.dtu)[mcols(se.tvn.dtu)$GENEID%in%tumor_novel_gene,]$structural_category)
plotBambu(se.tvn, type = "annotation", gene_id = "ENSG05220012074")

##cancer specific novel coding transcript structural_category distribution
dis <- as.data.frame(table(isoclassification[isoclassification$isoform %in% rownames(novel_dte),]$coding))
rownames(dis) <- dis$Var1
dis$Var1 <- NULL
boxplot(t(dis))
cancer_specific_novel <- isoclassification[isoclassification$isoform %in% rownames(novel_dte)&isoclassification$coding == 'coding',]
table <- as.data.frame(table(cancer_specific_novel$structural_category))
colnames(table) <- c('Var1','number')
pct <- round(100*table$number/sum(table$number))
pie(table$number,labels = paste(table$Var1, sep = " ", pct, "%"), col = rainbow(length(table$number)), main = "Cancer overexpression coding Novel transcrpts")


##DTU
library(DEXSeq)
colData(se.tvn)$condition = factor(colData(se.tvn)$groups)
dxd <- DEXSeqDataSet(countData = round(assays(se.tvn)$counts), sampleData = as.data.frame(colData(se.tvn)),
                     design = ~sample + exon + condition:exon, featureID = rowData(se.tvn)$TXNAME,
                     groupID = rowData(se.tvn)$GENEID)

rowRanges(dxd)
dxr <- DEXSeq(dxd)
dxr
mcols(dxr)$description
table ( dxr$padj < 0.05 )
table ( tapply( dxr$padj < 0.05, dxr$groupID, any ) )
df<- dxr[(!is.na(dxr$padj) & (dxr$padj < 0.010000)), ]
df <- as.data.frame(df)
df1 <- as.data.frame(df)
novel_df <- df1[df1$featureID %like% "BambuTx", ]
rownames(novel_df)
out <- stringr::str_split_fixed(rownames(novel_df),"\\:",2)
colnames(out) <- c("gene_id",'transcripts_id')
novel_df <- cbind(novel_df,out)
table ( novel_df$padj < 0.05 & novel_df$log2fold_tumour_normal > 0 )
#characterisitc of dtu
df<- as.data.frame(dxr[(!is.na(dxr$padj) & (dxr$padj < 0.050000)) & !is.na(dxr$log2fold_tumour_normal) & dxr$log2fold_tumour_normal > 0, ])
out <- stringr::str_split_fixed(rownames(df),"\\:",2)
colnames(out) <- c("gene_id",'transcripts_id')
df <- cbind(df,out)
table(mcols(se.multi.qc)[mcols(se.multi.qc)$TXNAME %in% df$transcripts_id,]$txClassDescription)
mcols(se.multi.qc)[mcols(se.multi.qc)$TXNAME %in% df$transcripts_id & mcols(se.multi.qc)$txClassDescription =='newFirstJunction:newFirstExon',]
plotBambu(se.tvn, type = "annotation", gene_id = "ENSG05220034315")
#dtu gene set enrichment
dtu <- as.data.frame(dxr[(!is.na(dxr$padj) & (dxr$padj < 0.050000)) & abs(dxr$log2fold_tumour_normal) > 0.5 & !is.na(dxr$log2fold_tumour_normal), ])
out <- stringr::str_split_fixed(rownames(dtu),"\\:",2)
colnames(out) <- c("gene_id",'transcripts_id')
dtu <- cbind(dtu,out)
gene_gtf <- fread("/Users/minchunchen/lab/data/isoseq/Homo_sapiens-GCA_009914755.4-2022_07-genes_name.gtf",sep = '\t',header = FALSE)
setnames(gene_gtf, names(gene_gtf), c("chr","source","type","start","end","score","strand","phase","attributes","gene_id","gene_name") )
refe <- gene_gtf[,c("gene_id","gene_name")]
dtu$gene <- refe$gene_name[match(dtu$gene_id,refe$gene_id)]

ranks <- dtu$stat
names(ranks)<-dtu$gene
library(fgsea)
pathways.hallmark <- gmtPathways("/Users/minchunchen/lab/h.all.v7.2.symbols.gmt")
head(pathways.hallmark)
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) 

# To see what genes are in each of these pathways:
fgseaResTidy$adjPvalue <- ifelse(fgseaResTidy$padj <= 0.1, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways Enrichment Score from GSEA")

#dtu_se
dtu_se <- se.tvn[mcols(se.tvn)$TXNAME %in% novel_df$transcripts_id,]
isoclassification <- read.delim("isoseq_classification.txt")
colnames((isoclassification))
mcols(dtu_se)$ORF_seq <- isoclassification[match(mcols(dtu_se)$TXNAME,isoclassification$isoform),]$ORF_seq

###Gviz visulization
library(Gviz)
library(GenomicRanges)
library(rtracklayer)
library(GenomicRanges)
gr_obj =  import("chm13v2.0_rmsk.bb")
gr_list = split(gr_obj, gr_obj$name)
chr3_TE <-gr_list[seqnames(gr_list) == 'chr3',]
BambuTx286_granges <- unlist(rowRanges(se.tvn[mcols(se.tvn)$TXNAME == 'BambuTx286',]))

## export cds for mhcferry
isoclassification[isoclassification$isoform == "BambuTx1519", ]$ORF_seq
aa <- isoclassification[isoclassification$isoform %in%rownames(novel_expr)&isoclassification$coding == "coding", ][,c("isoform","ORF_seq")]
write.table(aa, file="cancer_associate_isoTEratio_aa.txt",sep="\t", quote=F)
seqnames(gr_list) == 'chr3'

