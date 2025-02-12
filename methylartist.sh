#source ~/.bashrc
#tissue sample chr12:121987,694-121996313 SVA_F methylation plot
methylartist locus \
--bams \
/lustre/tgfaulkn/Drop/Barun/ONT/07PM0575.normal.hs1.aligned.mod.bam,\
/lustre/tgfaulkn/Drop/Barun/ONT/10PM1163.normal.hs1.aligned.mod.bam,\
/lustre/tgfaulkn/Drop/Barun/ONT/10PM2283.normal.hs1.aligned.mod.bam,\
/lustre/tgfaulkn/Drop/Barun/ONT/10PM2356.normal.hs1.aligned.mod.bam,\
/lustre/tgfaulkn/Drop/Barun/ONT/12PM1049.normal.hs1.aligned.mod.bam,\
/lustre/tgfaulkn/Drop/Barun/ONT/13PM0352.normal.hs1.aligned.mod.bam,\
/lustre/tgfaulkn/Drop/Barun/ONT/14PM0263.normal.hs1.aligned.mod.bam,\
/lustre/tgfaulkn/Drop/Barun/ONT/21DSRB006.normal.hs1.aligned.mod.bam,\
/lustre/tgfaulkn/Drop/Barun/ONT/07PM0575.tumour.hs1.aligned.mod.bam,\
/lustre/tgfaulkn/Drop/Barun/ONT/10PM1163.tumour.hs1.aligned.mod.bam,\
/lustre/tgfaulkn/Drop/Barun/ONT/10PM2283.tumour.hs1.aligned.mod.bam,\
/lustre/tgfaulkn/Drop/Barun/ONT/10PM2356.tumour.hs1.aligned.mod.bam,\
/lustre/tgfaulkn/Drop/Barun/ONT/12PM1049.tumour.hs1.aligned.mod.bam,\
/lustre/tgfaulkn/Drop/Barun/ONT/13PM0352.tumour.hs1.aligned.mod.bam,\
/lustre/tgfaulkn/Drop/Barun/ONT/14PM0263.tumour.hs1.aligned.mod.bam,\
/lustre/tgfaulkn/Drop/Barun/ONT/21DSRB006.tumour.hs1.aligned.mod.bam \
--ref /lustre/tmchen4/ref/hs1.fa \
--motif CG \
-i chr12:121987,694-121996313 \
-g /lustre/tmchen4/ref/Homo_sapiens_chr-GCA_009914755.4-2022_07-genes.gtf.gz \
-l 121989307-121990791 \
-p 1,6,1,3,4 \
-b \
/lustre/tgfaulkn/Drop/Barun/ONT/07PM0575.normal.hs1.aligned.mod.bam:#1f77b4,\
/lustre/tgfaulkn/Drop/Barun/ONT/10PM1163.normal.hs1.aligned.mod.bam:#377eb8,\
/lustre/tgfaulkn/Drop/Barun/ONT/10PM2283.normal.hs1.aligned.mod.bam:#6a3d9a,\
/lustre/tgfaulkn/Drop/Barun/ONT/10PM2356.normal.hs1.aligned.mod.bam:#4c72b0,\
/lustre/tgfaulkn/Drop/Barun/ONT/12PM1049.normal.hs1.aligned.mod.bam:#005ea8,\
/lustre/tgfaulkn/Drop/Barun/ONT/13PM0352.normal.hs1.aligned.mod.bam:#2986cc,\
/lustre/tgfaulkn/Drop/Barun/ONT/14PM0263.normal.hs1.aligned.mod.bam:#2077b4,\
/lustre/tgfaulkn/Drop/Barun/ONT/21DSRB006.normal.hs1.aligned.mod.bam:#165a9e,\
/lustre/tgfaulkn/Drop/Barun/ONT/07PM0575.tumour.hs1.aligned.mod.bam:#d62728,\
/lustre/tgfaulkn/Drop/Barun/ONT/10PM1163.tumour.hs1.aligned.mod.bam:#e41a1c,\
/lustre/tgfaulkn/Drop/Barun/ONT/10PM2283.tumour.hs1.aligned.mod.bam:#b2182b,\
/lustre/tgfaulkn/Drop/Barun/ONT/10PM2356.tumour.hs1.aligned.mod.bam:#a50f15,\
/lustre/tgfaulkn/Drop/Barun/ONT/12PM1049.tumour.hs1.aligned.mod.bam:#e31a1c,\
/lustre/tgfaulkn/Drop/Barun/ONT/13PM0352.tumour.hs1.aligned.mod.bam:#fc9272,\
/lustre/tgfaulkn/Drop/Barun/ONT/14PM0263.tumour.hs1.aligned.mod.bam:#de2d26,\
/lustre/tgfaulkn/Drop/Barun/ONT/21DSRB006.tumour.hs1.aligned.mod.bam:#f03b20 \
--labelgenes \
--svg

# HARA and NL20 chr12:121987,694-121996313 SVA_F methylation plot

methylartist locus \
--bams \
/lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_ont/NL20.hs1.aligned.mod.bam,\
/lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_ont/NL20_DNMTi_D2.hs1.aligned.mod.bam,\
/lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_ont/NL20_DNMTi_D4.hs1.aligned.mod.bam,\
/lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_ont/HARA.hs1.aligned.mod.bam,\
/lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_ont/HARA_DNMTi_D2.hs1.aligned.mod.bam,\
/lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_ont/HARA_DNMTi_D4.hs1.aligned.mod.bam \
--ref /lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_ont/hs1.fa \
--motif CG \
-i chr12:121987,694-121996313 \
-g /lustre/tmchen4/ref/Homo_sapiens_chr-GCA_009914755.4-2022_07-genes.gtf.gz \
-l 121989307-121990791 \
-p 1,6,1,3,4 \
-b \
/lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_ont/NL20.hs1.aligned.mod.bam:#1f77b4,\
/lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_ont/NL20_DNMTi_D2.hs1.aligned.mod.bam:#377eb8,\
/lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_ont/NL20_DNMTi_D4.hs1.aligned.mod.bam:#6a3d9a,\
/lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_ont/HARA.hs1.aligned.mod.bam:#d62728,\
/lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_ont/HARA_DNMTi_D2.hs1.aligned.mod.bam:#e41a1c,\
/lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_ont/HARA_DNMTi_D4.hs1.aligned.mod.bam:#b2182b \
--labelgenes \
--svg 


