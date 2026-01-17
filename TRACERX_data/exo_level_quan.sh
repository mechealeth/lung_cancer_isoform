cd /scratch/user/uqmche33/Ref
grep -E 'gene_name "CFAP251"|gene_id "[^"]*CFAP251[^"]*"' Homo_sapiens_chr-GCA_009914755.4-2022_07-genes.gtf > CFAP251.gtf

featureCounts \
  -a /scratch/user/uqmche33/Ref/CFAP251.gtf \
  -o tracerx_exon_counts.txt \
  -f -O -p \
  -R BAM \
  /QRISdata/Q7816/TRACERX_bam/*.bam
