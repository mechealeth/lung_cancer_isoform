library(readr)
library(dplyr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(pheatmap)

df <- read.csv("/scratch/user/uqmche33/utilis/CCLE_meta.csv", header = FALSE) 
df$V11 <- sub("^Non-Small Cell Lung Cancer.*?,\\s*", "", df$V11)
read_counts <- read.table("/scratch/project/te_lung/CCLE_lung/bamfile/sample_total_reads.txt", header = FALSE, sep = "\t",
                          col.names = c("sample", "sample", "total_mapped_reads"))

head(read_counts)
junction_counts <- read.delim("/scratch/project/te_lung/CCLE_lung/junction_counts_matrix.tsv", header = TRUE, sep = "\t")
# Join junction counts with total mapped reads
norm_counts <- junction_counts %>%
  rename(sample = Sample) %>%
  left_join(read_counts %>% select(sample, total_mapped_reads), by = "sample") %>%
  mutate(across(starts_with("BambuTx"), ~ .x / total_mapped_reads * 1e6, .names = "norm_{.col}")) # CPM scaling
plot_df <- norm_counts %>%
  select(sample, starts_with("norm_")) %>%
  pivot_longer(cols = starts_with("norm_"), names_to = "Transcript", values_to = "CPM") %>%
  left_join(df %>% select(V1, V11, V15), by = c("sample" = "V1"))

####cancer subtype#### Filter to just one transcript, e.g., BambuTx548
single_tx <- plot_df %>%
  filter(Transcript == "norm_BambuTx548_121990050_121995520")

p1 <- ggplot(single_tx, aes(x = V11, y = CPM, fill = V11)) +
  geom_violin(trim = FALSE, scale = "width", color = "black", alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, 
               color = "black", fill = "yellow") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"  # hide legend
  ) +
  labs(
    x = "Cancer subtype (V11)",
    y = "Normalized junction count (CPM)"
  )

ggsave("/scratch/project/te_lung/CCLE_lung/cancer_type_violin_plot.pdf", plot = p1, width = 6, height = 4, units = "in", device = "pdf")

####primary vs metastss #### Filter to one transcript
single_tx <- plot_df %>%
  filter(Transcript == "norm_BambuTx548_121990050_121995520")

p <- ggplot(single_tx, aes(x = V15, y = CPM, fill = V15)) +
  geom_violin(trim = FALSE, scale = "width", color = "black", alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, 
               color = "black", fill = "yellow") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    x = "Tissue type",
    y = "Normalized junction count (CPM)"
  )

ggsave("/scratch/project/te_lung/CCLE_lung/cancer_stage_violin_plot.pdf", plot = p, width = 6, height = 4, units = "in", device = "pdf")

###merge all CCLE TEcount data###
library(data.table)

# Folder with your .cntTable files
cnt_dir <- "/scratch/project/te_lung/CCLE_lung/tecount"   # <-- change if needed

files <- list.files(cnt_dir, pattern = "\\.cntTable$", full.names = TRUE)
stopifnot(length(files) > 0)

# read one file -> 2 columns: gene + sample
read_cnt <- function(f) {
  dt <- fread(f)
  setnames(dt, 1, "gene")  # first col = gene name
  sample_id <- sub("\\.cntTable$", "", basename(f))  # e.g. SRR8616192
  setnames(dt, 2, sample_id)
  dt
}

dt_list <- lapply(files, read_cnt)

# merge by gene, keep all genes
merged <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), dt_list)

# replace NAs with 0 and coerce to integer
for (j in 2:ncol(merged)) set(merged, which(is.na(merged[[j]])), j, 0L)

# make rownames = gene, drop gene column
mat <- as.data.frame(merged)
rownames(mat) <- mat$gene
mat$gene <- NULL

# optional: order columns by sample id
mat <- mat[, order(colnames(mat)), drop = FALSE] 
mat <- mat[ !grepl("ENSG", rownames(mat), ignore.case = FALSE), , drop = FALSE ]    
write.table(mat,"/scratch/project/te_lung/CCLE_lung/merged_tecount_matrix.tsv", sep = "\t", quote = FALSE, col.names = NA)

# ---------- Align samples ----------
if (!"sample" %in% names(read_counts) && "sample.1" %in% names(read_counts)) {
  names(read_counts)[names(read_counts) == "sample.1"] <- "sample"
}
samps <- Reduce(intersect, list(colnames(mat), as.character(read_counts$sample), as.character(df$V1)))
stopifnot(length(samps) > 0)

mat   <- data.matrix(mat[, samps, drop = FALSE])
tot   <- setNames(read_counts$total_mapped_reads, read_counts$sample)[samps]
meta  <- df[df$V1 %in% samps, c("V1","V11")]
colnames(meta) <- c("sample","group")
meta  <- meta[match(samps, meta$sample), , drop = FALSE]

# ---------- Collapse rows to TE subfamily (token before first ':') ----------
subfam <- sub(":.*$", "", rownames(mat))   # "L1HS:L1:LINE" -> "L1HS"
counts_by_sf <- rowsum(mat, group = subfam, reorder = FALSE)

# (Optional) restrict to a fixed panel (comment this block to keep all)
te_set <- c("L1HS","L1PA2","L1PA3","L1PA5","L1PA8","LTR12C","LTR12A",
            "MLT2D","MLT1B","SVA_A","SVA_B","SVA_C","SVA_D","SVA_F")
keep <- intersect(te_set, rownames(counts_by_sf))
if (length(keep) > 0) counts_by_sf <- counts_by_sf[keep, , drop = FALSE]

# ---------- CPM normalize ----------
cpm <- t( t(counts_by_sf) / tot ) * 1e6

# ---------- Order columns by cancer group (df$V11) ----------
ord <- order(meta$group, meta$sample)
cpm <- cpm[, ord, drop = FALSE]
meta <- meta[ord, , drop = FALSE]

# ---------- Group means by cancer type ----------
grp_lvls  <- unique(meta$group)
grp_means <- sapply(grp_lvls, function(g) rowMeans(cpm[, meta$group == g, drop = FALSE], na.rm = TRUE))
colnames(grp_means) <- grp_lvls

# ---------- Transform for heatmaps ----------
row_z      <- function(m) t(scale(t(m)))           # z-score across columns per row
mat_grp_z  <- row_z(log1p(grp_means))              # group-mean heatmap
mat_samp_z <- row_z(log1p(cpm))                    # per-sample heatmap

# (Optional) cap extremes for better contrast
cap <- function(x, lo=-2.5, hi=2.5) pmax(pmin(x, hi), lo)
mat_grp_plot  <- cap(mat_grp_z)
mat_samp_plot <- cap(mat_samp_z)

# ---------- Annotations & colors ----------
hm_cols  <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
ann_col  <- data.frame(CancerType = meta$group); rownames(ann_col) <- meta$sample
gaps     <- cumsum(table(meta$group)); if (length(gaps) > 0) gaps <- gaps[-length(gaps)]

# ---------- Output dir ----------
out_dir <- "/scratch/project/te_lung/CCLE_lung"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Auto sizes so labels aren’t cramped
w_samp <- max(7, ncol(mat_samp_plot) * 0.08)
h_all  <- max(4, nrow(mat_samp_plot) * 0.35)

# ---------- NO-CLUSTER heatmaps ----------
# Group means (columns = cancer types)
pheatmap(
  mat_grp_plot,
  color = hm_cols,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  angle_col = 45,
  main = "TE subfamilies — row-scaled log1p(CPM) (group means)",
  filename = file.path(out_dir, "TE_family_heatmap_group_mean.pdf"),
  width = 6, height = 4
)

# Per-sample (one column per sample), with group separators
pheatmap(
  mat_samp_plot,
  color = hm_cols,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  show_colnames = FALSE,
  annotation_col = ann_col,
  gaps_col = gaps,
  main = "TE subfamilies — row-scaled log1p(CPM) (per sample)",
  filename = file.path(out_dir, "TE_family_heatmap_per_sample_no_cluster.pdf"),
  width = w_samp, height = h_all
)

message("Wrote:\n  ",
        file.path(out_dir, "TE_family_heatmap_group_mean.pdf"), "\n  ",
        file.path(out_dir, "TE_family_heatmap_per_sample_no_cluster.pdf"))                 
