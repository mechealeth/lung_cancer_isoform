if (!requireNamespace("IOBR", quietly = TRUE))
  devtools::install_github("IOBR/IOBR")
library(IOBR) 
data("anno_grch38", package = "IOBR")
eset <- anno_eset(eset = tumor_expr, annotation = anno_grch38, probe = "symbol")
eset <- count2tpm(countMat = eset, source = "local", idType = "symbol")
head(eset)
##transform id to symbol
gen_gtf <- fread("/Users/minchunchen/lab/data/isoseq/Homo_sapiens-GCA_009914755.4-2022_07-genes_name.gtf",sep = '\t',header = FALSE)
setnames(gen_gtf, names(gen_gtf), c("chr","source","type","start","end","score","strand","phase","attributes","gene_id","gene_name") )
refe <- gen_gtf[,c("gene_id","gene_name")]
rm(gen_gtf)
logTPM <- function(x) {return(log2(x+1))}
##transform gene count to tpm (HARA/NL20)
cell_line_expr <- read.delim("counts_gene.txt")
head(cell_line_expr)
colnames(cell_line_expr) <- sapply(strsplit(colnames(cell_line_expr), "\\."), "[", 1)
cell_line_expr <- cell_line_expr[,c(1,10,18)]
cell_line_expr$GENEID <- refe$gene_name[match(cell_line_expr$GENEID,refe$gene_id)]
cell_line_expr$GENEID[cell_line_expr$GENEID == ''] <- NA 
cell_line_expr <- na.omit(cell_line_expr)
k = !duplicated(cell_line_expr$GENEID);table(k)
cell_line_expr <- cell_line_expr[k,]
rownames(cell_line_expr) <- cell_line_expr$GENEID
cell_line_expr$GENEID <- NULL
rm(eset)
# count to log2(TPM+1)
cell_line_expr <- count2tpm(countMat = cell_line_expr, source = "local", idType = "symbol")
cell_line_expr <- cell_line_expr %>% mutate_if(is.numeric, logTPM)
head(cell_line_expr)
tumor_expr <- read.csv("/Users/minchunchen/lab/data/isoseq/gene_count_matrix.csv",check.names = FALSE)
head(tumor_expr)
rownames(tumor_expr) <- tumor_expr$gene_id
tumor_expr$gene_id <- NULL
tumor_expr <- tumor_expr[colnames(tumor_expr)%like%".tumour"]
tumor_expr$gene <-refe$gene_name[match(rownames(tumor_expr),refe$gene_id)]
tumor_expr$gene[tumor_expr$gene == ''] <- NA 
tumor_expr <- na.omit(tumor_expr)
head(tumor_expr)
k = !duplicated(tumor_expr$gene);table(k)
tumor_expr <- tumor_expr[k,]
rownames(tumor_expr) <- tumor_expr$gene
tumor_expr$gene <- NULL
head(tumor_expr)
# count to log2(TPM+1)
tumor_expr <- count2tpm(countMat = tumor_expr, source = "local", idType = "symbol")
tumor_expr <- tumor_expr %>% mutate_if(is.numeric, logTPM)
##common gene
list <- intersect(rownames(tumor_expr),rownames(cell_line_expr))
)
common_gene <- as.data.frame(intersect(list,rownames(lusc_expression)))
colnames(common_gene) <- "gene"
head(common_gene)
##
lusc_expression <- lusc_expression[rownames(lusc_expression)%in%common_gene$gene,]
cell_line_expr <- cell_line_expr[rownames(cell_line_expr)%in%common_gene$gene,]
colnames(cell_line_expr) <- c("HARA_lab","NL20")
tumor_expr <- tumor_expr[rownames(tumor_expr)%in%common_gene$gene,]
cell_line_expr$gene <- rownames(cell_line_expr)
lusc_expression$gene <-rownames(lusc_expression)
cell_line_expr <- merge(cell_line_expr,lusc_expression,by="gene")
cell_line_expr <- as.data.frame(cell_line_expr)
rownames(cell_line_expr) <- cell_line_expr$gene
cell_line_expr$gene <- NULL
write.table(cell_line_expr,"cell_line_expr.txt",sep = "\t")
write.table(tumor_expr,"tumo_expr.txt",sep = "\t")
###mutation data
cell_line_mut <- read.csv("OmicsSomaticMutations .csv")
