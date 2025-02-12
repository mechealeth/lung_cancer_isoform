#####for SVA enrichment hist plot
# Load necessary library
library(ggplot2)

# Create dataframe
df <- data.frame(
  Category = c("SVA", "SVA_F"),
  Enrichment_Score = c(5.395144, 7.019241),
  P_Value = c(0.564850253, 7.23776255456158e-07)
)

# Generate the bar plot
ggplot(df, aes(x = Category, y = Enrichment_Score, fill = Category)) +
  geom_bar(stat = "identity", width = 0.6) + 
  geom_text(aes(label = sprintf("%.2e", P_Value)), vjust = -0.5, size = 5) + 
  labs( x = "Category", y = "Enrichment Score") +
  theme_classic() +  # Removes grid but keeps x and y axis
  theme(legend.position = "none")  # Hide legend
##### for gene transcript composition hist plot
mcols(se.multi.qc[mcols(se.multi.qc)$TXNAME == 'BambuTx548',])
rowRanges(se.multi.qc[mcols(se.multi.qc)$GENEID == 'ENSG05220040131',])
plotBambu(se.multi.qc, type = "annotation", gene_id = "ENSG05220040131")
mcols(se.multi.qc[mcols(se.multi.qc)$GENEID == 'ENSG05220040131',])

df <- round(assays(se.multi.qc)$counts)
colnames(df)
df <- df[,c(-7,-8,-9,-17)]
colnames(df)
erv_transcrtipt <- c("BambuTx647","BambuTx648")
df <- df[erv_transcrtipt,]
df <- prop.table(df, margin = 2)
df <- as.data.frame(t(df))
df$sample <- rownames(df)
library(data.table)
long <- melt(setDT(df), id.vars = "sample",variable.name = "isoform")
long <- as.data.frame(long)

allcolour=c("#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

name <- c("HARA_D0_R1","HARA_D0_R2","HARA_D2_R1","HARA_D2_R2","HARA_D4_R1","HARA_D4_R2","NL20_D0_R1","NL20_D0_R2","NL20_D2_R1","NL20_D2_R2","NL20_D4_R1","NL20_D4_R2","NL20_D4_R3")
long$sample <-rep(c("HARA_D0_R1","HARA_D0_R2","HARA_D2_R1","HARA_D2_R2","HARA_D4_R1","HARA_D4_R2","NL20_D0_R1","NL20_D0_R2","NL20_D2_R1","NL20_D2_R2","NL20_D4_R1","NL20_D4_R2","NL20_D4_R3"),times=2) 
library(ggplot2)
ggplot(long) + 
  geom_bar(aes(x =sample, y= value, fill = isoform),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Expression Ratio')+
  scale_fill_manual(values = allcolour)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid")) + theme(axis.text.x = element_text(size = 6))

####sunburst plot(two levels)
library(plotly)

# Define total values
total_value <- 1231
te_derived_value <- 588
other_novels_value <- total_value - te_derived_value

# Compute first-level percentages (full circle)
te_derived_percent <- (te_derived_value / total_value) * 100
other_novels_percent <- (other_novels_value / total_value) * 100

# Define Cancer-Specific proportions within TE-derived Novels
cancer_specific_value <- 290
none_value <- 1062
total_te_derived <- cancer_specific_value + none_value

# Compute second-level percentages within "TE-derived Novels" (relative percentages)
cancer_specific_percent <- (cancer_specific_value / total_te_derived) * te_derived_percent
none_percent <- (none_value / total_te_derived) * te_derived_percent

# Define hierarchical structure with only two levels
data <- data.frame(
  labels = c(
    "Total", "Other Novels", "TE-derived Novels",
    "Cancer-specific", "None"
  ),
  parents = c(
    "", "Total", "Total",
    "TE-derived Novels", "TE-derived Novels"
  ),
  values = c(
    100, other_novels_percent, te_derived_percent,
    cancer_specific_percent, none_percent
  ),
  colors = c(
    "#ffffff", "#4682B4", "#FF8C00",
    "#FFD700", "#B0C4DE"
  ) # Different colors for clarity
)

# Create sunburst plot ensuring correct hierarchical structure
fig <- plot_ly(
  data,
  labels = ~labels,
  parents = ~parents,
  values = ~values,
  type = 'sunburst',
  branchvalues = 'total',  # Ensures correct scaling at each level
  textinfo = 'label+percent entry',
  marker = list(colors = data$colors)  # Apply custom colors
)

# Show plot
fig
