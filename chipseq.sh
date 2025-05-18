### hs1 TSS
awk '$3=="transcript"' /lustre/tmchen4/ref/Homo_sapiens_chr-GCA_009914755.4-2022_07-genes.gtf | \
awk 'BEGIN{OFS="\t"} {
    match($0, /gene_id "[^"]+"/, a); gene_id=a[0];
    if ($7=="+") { start=$4-3000; end=$4+3000; }
    else { start=$5-3000; end=$5+3000; }
    if (start<0) start=0;
    print $1, start, end, gene_id, ".", $7;
}' > promoters_TSS_3kb.bed

computeMatrix reference-point \
  -S 1HD-H3K4M3-LHC26213_L2_1.srt.nodup.bw 2HD-H3K4M3-LHC26231_L8_1.srt.nodup.bw 3HD-H3K4M3-LHC26248_L8_1.srt.nodup.bw \
  -R promoters_TSS_3kb.bed \
  --referencePoint TSS -b 3000 -a 3000 -bs 50 --skipZeros \
  -out matrix_HD_H3K4me3_TSS.gz

plotHeatmap \
  -m matrix_H3K4me3_TSS.gz \
  -out H3K4me3_promoter_heatmap.pdf \
  --colorMap Reds 

#!/bin/bash
module load anaconda3 
source activate /scratch/project/te_lung/esm/envs
cd /scratch/project/te_lung

# Config
CREDENTIALS="/scratch/project/te_lung/EGA_credential_file.json"
OUTPUT_DIR="/QRISdata/Q7816/TRACERX"
EGAF_LIST="/scratch/project/te_lung/egaf_list.txt"  # Text file with EGAF IDs, one per line

# Loop through each EGAF ID and fetch it
while read -r EGAF_ID; do
    echo "Fetching $EGAF_ID..."
    pyega3 -cf "$CREDENTIALS" fetch "$EGAF_ID" --output-dir "$OUTPUT_DIR"
done < "$EGAF_LIST"
