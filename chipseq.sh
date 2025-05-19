### hs1 TSS

awk 'BEGIN{OFS="\t"} 
$3 == "transcript" {
    match($0, /gene_id "[^"]+"/, gid);
    gsub(/gene_id "|"/, "", gid[0]);
    if ($7 == "+") 
        print $1, $4-1, $4, gid[0], ".", $7; 
    else if ($7 == "-") 
        print $1, $5-1, $5, gid[0], ".", $7;
}' /lustre/tmchen4/ref/Homo_sapiens_chr-GCA_009914755.4-2022_07-genes.gtf > TSS_transcripts.bed


computeMatrix reference-point \
  -S 1HI-H3K4M3-LHC26216_L2_1.srt.nodup.bw 3HI-H3K4M3-LHC26253_L8_1.srt.nodup.bw \
  -R TSS_transcripts.bed \
  --referencePoint TSS -b 3000 -a 3000 -bs 50 --skipZeros \
  -out matrix_HI_H3K4me3_TSS.gz


plotHeatmap \
  -m matrix_HD_H3K4me3_TSS.gz \
  -out matrix_HI_H3K4me3_TSS.gz \
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
