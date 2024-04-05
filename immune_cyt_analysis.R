setwd("/Users/minchunchen/lab/data/isoseq")
#data loading and sample information prepare
readCount <- read.table("gene_count_matrix.csv",sep = ',', header = TRUE,row.names = 1,check.names = FALSE)
sampleInfo<- as.data.frame(colnames(readCount))
colnames(sampleInfo) <- "sample_id"
out <- stringr::str_split_fixed(sampleInfo$sample_id,"\\.",2)
colnames(out) <- c("patients",'groups')
sample <- cbind(sampleInfo,out)
rownames(sample) <- sample$sample_id
sample$sample_id <- NULL
sample$groups<- factor(sample$groups, levels = c('tumour','normal'))
#deseq2
#Create deseqData
dds <- DESeqDataSetFromMatrix(countData = readCount, colData = sample, design = ~ groups)
#gene qc and normalization
dim(dds)
rowSums(counts(dds))
keep <- rowSums(counts(dds) >= 10) >= round(ncol(dds)*0.05)
table(keep)
dds1 <- dds[keep,]
dds1 <- DESeq(dds1)
vsd <-as.data.frame(counts(dds1,normalized=TRUE))
rn2 <- c("04PM0832.normal","08PM0370.normal","10PM2283.normal","12PM0832.normal","14PM0263.normal","06PM0981.normal","10PM1163.normal","10PM2356.normal","12PM1049.normal","21DSRB006.normal","07PM0575.normal","10PM1669.normal","12PM0142.normal","13PM0352.normal",
         "04PM0832.tumour","08PM0370.tumour","10PM2283.tumour","12PM0832.tumour","14PM0263.tumour","06PM0981.tumour","10PM1163.tumour","10PM2356.tumour","12PM1049.tumour","21DSRB006.tumour","07PM0575.tumour","10PM1669.tumour","12PM0142.tumour","13PM0352.tumour" )
vsd <- as.data.frame(vsd[,match(rn2,colnames(vsd))])
gene_gtf <- fread("/Users/minchunchen/lab/data/isoseq/Homo_sapiens-GCA_009914755.4-2022_07-genes_name.gtf",sep = '\t',header = FALSE)
setnames(gene_gtf, names(gene_gtf), c("chr","source","type","start","end","score","strand","phase","attributes","gene_id","gene_name") )
refe <- gene_gtf[,c("gene_id","gene_name")]
vsd$gene <- refe$gene_name[match(rownames(vsd),refe$gene_id)]
vsd[vsd == ''] <- NA 
vsd <- as.data.frame(na.omit(vsd))
k = !duplicated(vsd$gene)
table(k)
vsd <- vsd[k,]
rownames(vsd) <- vsd$gene
vsd$gene <- NULL  
write.table(vsd[,15:28], file="gene_count.txt",sep="\t", quote=F)
normalized_counts1 = rownames_to_column(vsd)
write.table(normalized_counts1 ,file = "exp.txt",row.names = F,quote = F,sep = "\t")
f = "ciber_GSE201050.Rdata"
if(!file.exists(f)){
  #devtools:: install_github ("Moonerss/CIBERSORT")
  library(CIBERSORT)
  lm22f = system.file("extdata", "LM22.txt", package = "CIBERSORT")
  TME.results = cibersort(lm22f, 
                          "exp.txt" , 
                          perm = 1000, 
                          QN = T)
  save(TME.results,file = f)
}
load(f)
TME.results[1:4,1:4]
T.results <- TME.results[TME.results[,23] < 0.05,]
re <- T.results[,-(23:25)]
write.table(re2[,13:25], file="gene_count.txt",sep="\t", quote=F)
library(pheatmap)
k <- apply(re,2,function(x) {sum(x == 0) < nrow(TME.results)/2})
table(k)
re2 <- as.data.frame(t(re[,k]))
ano=data.frame(ici$cluster[match(rownames(re), ici$SampleID)])
sample <-ici$cluster[match(rownames(re), ici$SampleID)]
Group<- ici$Group[match(rownames(re), ici$SampleID)]
respond <- ici$group[match(rownames(re), ici$SampleID)]
ano=cbind(ano,Group,respond)
colnames(ano) <- c('Cluster','Group','Respond')
rownames(ano)=colnames(re2)
ano <- ano %>% arrange(Cluster)
###
re2_ordered <- re2[, rownames(ano)]
dim(re2)
pheatmap(re2[,13:25],scale = "row",
         show_colnames = T,
         cluster_cols = F,
         cluster_rows  = F,
         drop_levels = TRUE,
         color = c(rep("blue",11),colorRampPalette(colors = c("blue","white","red"))(50),
                   rep("red",11)),
         #legend_breaks = c(-5,-2,0,2,5)
         filename = 'ciber_hmap.png'
)

#Hist
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))

dat <- re %>% 
  as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)
dat$group <-ici$cluster[match(dat$Sample, ici$SampleID)]
dat <- dat %>% arrange(group)
group <- as.data.frame(ano$Cluster)
dat$Sample = factor(dat$Sample,ordered = T,levels = unique(dat$Sample)) 
dat2 = data.frame(a = 1:ncol(re2),
                  b = 1,
                  group = unique(dat$Sample) )
p1 = ggplot(dat2,aes(x = a, y = b)) + 
  geom_tile(aes(fill = group)) +
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        axis.title = element_blank()) + 
  scale_x_continuous(expand = c(0, 0)) +
  labs(fill = "group")

p2 = ggplot(dat,aes(Sample, Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x=a,y = "Estiamted Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  ) + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(22))
 a <- '04PM0832.tumour 08PM0370.tumour 10PM2283.tumour 12PM0832.tumour 14PM0263.tumour 14PM0263.tumour 06PM0981.tumour 10PM2356.tumour 12PM1049.tumour'
library(patchwork)
p1/ p2 + plot_layout(heights = c(1,12),guides = "collect" ) &
  theme(legend.position = "bottom") 

sample.order <- dat$sample
ggbarplot(dat,  x = "Sample",  y = "Proportion",  size = 0,  fill = "Cell_type",  color = "Cell_type",order = sample.order)  +  theme(  axis.text.x = element_text(      angle = 90,      hjust = 1,      vjust = 1,      size = 1    ),    legend.position = "bottom",    legend.text = element_text(size= 8),    legend.title= element_text(size= 8)  )
#compare
##immune analysis
geneCount <- read.table("gene_count_matrix.csv",sep = ',', header = TRUE,row.names = 1,check.names = FALSE)
gene_gtf <- fread("/Users/minchunchen/lab/data/isoseq/Homo_sapiens-GCA_009914755.4-2022_07-genes_name.gtf",sep = '\t',header = FALSE)
setnames(gene_gtf, names(gene_gtf), c("chr","source","type","start","end","score","strand","phase","attributes","gene_id","gene_name") )
refe <- gene_gtf[,c("gene_id","gene_name")]
geneCount$gene <- refe$gene_name[match(rownames(geneCount),refe$gene_id)]
geneCount[geneCount == ''] <- NA
geneCount <- na.omit(geneCount)
k = !duplicated(geneCount$gene)
table(k)
geneCount <- geneCount[k,]
rownames(geneCount) <- geneCount$gene
geneCount$gene <- NULL
sampleInfo<- as.data.frame(colnames(geneCount))
colnames(sampleInfo) <- "sample_id"
out <- stringr::str_split_fixed(sampleInfo$sample_id,"\\.",2)
colnames(out) <- c("patients",'groups')
sample <- cbind(sampleInfo,out)
rownames(sample) <- sample$sample_id
sample$sample_id <- NULL
sample$groups<- factor(sample$groups, levels = c('tumour','normal'))
#deseq2
#Create deseqData
dds <- DESeqDataSetFromMatrix(countData = geneCount, colData = sample, design = ~ groups)
dds.deseq <- DESeq(dds)
normalize_count <- counts(dds.deseq,normalized=TRUE)
write.table(as.data.frame(normalize_count), file="gene_symbol_tumor_vs_normal_iso.txt",sep="\t", quote=F)
###CTY calculate
compute.CYT <- function(RNA.tpm){
  CYT.read <- c("GZMA", "PRF1")
  match_CYT.genes <- match(CYT.read, rownames(RNA.tpm))
  if (anyNA(match_CYT.genes)){
    warning(paste0("differenty named or missing signature genes : \n", paste(CYT.read[!CYT.read %in% rownames(RNA.tpm)], collapse = "\n")))
    match_CYT.genes <- stats::na.omit(match_CYT.genes)
  }
  subset_RNA.tpm <- RNA.tpm[match_CYT.genes, ]
  score <- as.matrix(apply(subset_RNA.tpm + 0.01, 2, function(X) exp(mean(log(X)))))
  return(data.frame(CYT = score, check.names = FALSE))
}
CYT_score <- compute.CYT(normalize_count)
rownames(sample[sample$groups == 'tumour',])
CYT_score <- CYT_score[rownames(CYT_score) %in% rownames(sample[sample$groups == 'tumour',]),]
CYT_score <- as.data.frame(CYT_score,row.names = rownames(sample[sample$groups == 'tumour',]))
immune_scoe <- read.table("/Users/minchunchen/Downloads/xcell_result/xCell_gene_count_xCell_0734022824.txt",check.names=FALSE ,sep = '\t',header = TRUE,row.names = 1)
colnames(immune_scoe) <- colnames(readcount)
rownames(immune_scoe) <-  immune_scoe[,1]
immune_scoe[,1] <- NULL
immune_score <- t(immune_scoe[rownames(immune_scoe)=='ImmuneScore',colnames(immune_scoe) %in% rownames(sample[sample$groups == 'tumour',])])
sig_tumor <- cbind(CYT_score,immune_score)
#expr data transform
c_expr <- cbind(expr,immune_score)
colnames(c_expr)
library(ggstatsplot)
y <- as.numeric(c_expr[,"CYT_score"])
colnames <- colnames(c_expr)
cor_data_df5<- data.frame(colnames)
for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(c_expr[,i]),y,type="spearman")
  cor_data_df5[i,2] <- test$estimate
  cor_data_df5[i,3] <- test$p.value 
}
names(cor_data_df5) <- c("symbol","correlation","pvalue")
###

  hist(as.data.frame(c_expr)$CYT_score, main = "Distribution of CYT score",
     xlab = "CYT_score", ylab = "sample numbers", col = "skyblue")  # 直方图
curve(dnorm(x, mean = mean(as.data.frame(c_expr)$CYT_score), sd = sd(as.data.frame(c_expr)$CYT_score)), add = TRUE, col = "red")  # 正态分布的概率密度函数
# shapiro-wilk检验
shapiro.test(as.data.frame(c_expr)$CYT_score)  # 适用于小样本数据，N≤50，若p值<0.05，则数据不符合正态分布

time <- dataL[,"month"]status <- dataL[,"Status"]x <- as.matrix(dataL[,-c(1,2)])  #x为输入特征,是矩阵格式y <- as.matrix(dataL$Status)lasso <- glmnet(x = x, y = y,                family = "binomial",                alpha = 1,  # alpha = 1为LASSO回归，= 0为岭回归，0和1之间则为弹性网络                  nlambda = 100)  # nlambda表示正则化路径中的个数，这个参数就可以起到一个阈值的作用，决定有多少基因的系数可以留下来。默认值为100。print(lasso)plot(lasso, xvar = "lambda", label = TRUE)  # 系数分布图，由“log-lambda”与变量系数作图，展示根据lambda的变化情况每一个特征的系数变化，展示了Lasso回归筛选变量的动态过程plot(lasso, xvar = "dev", label = TRUE)  # 也可以对%dev绘图plot(lasso, xvar = "norm", label = TRUE)  # “L1范数”与变量系数作图

set.seed(1234)
lasso_cv <- cv.glmnet(x = as.matrix(c_expr), y = as.numeric(c_expr[,"CYT_score"]), family = "poisson", alpha = 1, nlambda = 100)  # 交叉验证，如果结果不理想，可以重新单独运行这一行代码，或者换一下种子数
plot(lasso_cv)
lambda <- lasso_cv$lambda.min
coef_lasso_cv <- coef(lasso_cv, s = lambda)
coef_lasso_cv[,1][coef_lasso_cv[,1]!=0]

##immune score
set.seed(1234)
lasso_cv <- cv.glmnet(x = as.matrix(c_expr), y = as.numeric(c_expr[,"ImmuneScore"]), family = "gaussian", alpha = 1, nlambda = 100)  # 交叉验证，如果结果不理想，可以重新单独运行这一行代码，或者换一下种子数
plot(lasso_cv)
lambda <- lasso_cv$lambda.min
coef_lasso_cv <- coef(lasso_cv, s = lambda)
coef_lasso_cv[,1][coef_lasso_cv[,1]!=0]

#Estimate tumor purity
install.packages("tidyestimate")
library(tidyestimate) 
rownames(normalized_counts1) <- normalized_counts1$rowname
normalized_counts1$rowname <- NULL
scores <- vsd[,15:28] |> 
  filter_common_genes(id = "hgnc_symbol", tell_missing = FALSE, find_alias = TRUE) |> 
  estimate_score(is_affymetrix = TRUE)
scores |> 
  plot_purity(is_affymetrix = TRUE)
scores$immune_score <- c(0.9403,0.2052, 0.4036,0.5095, 0.9257,0.8507,0.0973,0.1648,0.8361, 0.8319,0.4799,0.3923,0.6911,0.4715)
immune_scoe[rownames(immune_scoe)=='ImmuneScore',]

write_csv(scores ,file = "tumour.immune_score_purity.csv")
