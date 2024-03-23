#!/usr/bin/Rscript
#setting path
setwd("/lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_pacbio")
out_data_dir <- "/lustre/tmchen4/iso-seq/cell_lines/bambu"
#load library
library(Rsamtools)
library(bambu)
library(BSgenome)
library(data.table)
#load data
dir() -> files
filename <- files[grep("*hs1.aligned.bam$",files)]
bamFiles <- Rsamtools::BamFileList(filename)
#create fai file with samtool
fa.file <- "/lustre/tmchen4/ref/hs1.fa"
bambuAnnotations <- prepareAnnotations("/lustre/tmchen4/ref/Homo_sapiens_chr-GCA_009914755.4-2022_07-genes.gtf")
setDTthreads(4)
se.multi <- bambu(reads = bamFiles, annotations = bambuAnnotations, genome = fa.file, NDR = 0.1 ,verbose=TRUE)
saveRDS(se.multi, file=paste0(out_data_dir, 'cell_line.rds'))
writeBambuOutput(se.multi, path = "/lustre/tmchen4/iso-seq/cell_line_bambu_out")