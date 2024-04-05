setwd("/lustre/tmchen4/ref")
library(data.table)
library(stringr)
gtf <- fread("/lustre/tmchen4/ref/Homo_sapiens-GCA_009914755.4-2022_07-genes_chr.gtf",sep = '\t',header = FALSE)
setnames(gtf, names(gtf), c("chr","source","type","start","end","score","strand","phase","attributes") )
# which will drastically reduce the file size
genes <- genes[type == "gene"]

# the problem is the attributes column that tends to be a collection
# of the bits of information you're actually interested in
# in order to pull out just the information I want based on the 
# tag name, e.g. "gene_id", I have the following function:
extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- unlist(strsplit(gtf_attributes, " "))
  if(att_of_interest %in% att){
    return(gsub("\"|;","", att[which(att %in% att_of_interest)+1]))
  }else{
    return(NA)}
}
gtf$gene_id <- unlist(lapply(gtf$attributes, extract_attributes, "gene_id"))
df$transcript_id <- unlist(lapply(gtf$attributes, extract_attributes, "transcript_id"))
df$transcript_id <- unlist(lapply(gtf$attributes, extract_attributes, "transcript_id"))
#write gtf
fwrite(gtf,'/Users/minchunchen/lab/data/modified.sofilter.filtered.gtf',sep = '\t',quote = FALSE,col.names = FALSE)

/lustre/tmchen4/iso-seq/bambu/bambu_NDR0.155/SQANTI3_output_hs1SJ.out.tab/isofilter.filtered.gtf
/lustre/tmchen4/ref/Homo_sapiens-GCA_009914755.4-2022_07-genes_chr.gtf

projection_parent_gene <- unlist(lapply(gtf$attributes, extract_attributes, "projection_parent_gene"))
df <- unlist(lapply(gtf$attributes, extract_attributes, "gene_id"))
gtf[x %in% month.name[c(2)], 
gene_gtf <- gtf[type == "gene", ]
gene_gtf$gene_name <- unlist(lapply(gene_gtf$attributes, extract_attributes, "gene_name"))
gene_gtf$gene_id <- unlist(lapply(gene_gtf$attributes, extract_attributes, "gene_id"))
gene_gtf$projection_parent_gene  <- unlist(lapply(gene_gtf$attributes, extract_attributes, "projection_parent_gene")) 
fwrite(gene_gtf,'/lustre/tmchen4/ref/Homo_sapiens-GCA_009914755.4-2022_07-genes_name.gtf',sep = '\t',quote =
FALSE,col.names = FALSE)

#prepare annotation
inDir = '/Users/minchunchen/lab/data/isoseq'
flattenedFile = list.files(inDir, pattern="gtf$", full.names=TRUE)
basename(flattenedFile) 
colnames(df)
df <- as.data.frame(round(assays(se.multi.qc)$counts),row.names = paste0(rowData(se.multi.qc)$GENEID,":",rowData(se.multi.qc)$TXNAME))
lapply(colnames(df)[1:8],function(i){
  a <- as.data.frame(df[,colnames(df)==i],row.names = rownames(df))
  rownames(a) <- rownames(df)
  write.table(a, file = paste0(i,".txt"),sep = "\t",quote = F,row.names = T,col.names = F)
})
countFiles = list.files("/Users/minchunchen/lab/data/isoseq/sample",pattern = 'txt',full.names=TRUE)
setwd("/Users/minchunchen/lab/data/isoseq/sample")

basename(countFiles)
dxd <- DEXSeqDataSetFromHTSeq(countFiles,
  sampleData=as.data.frame(colData(se.multi.qc)),
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile)