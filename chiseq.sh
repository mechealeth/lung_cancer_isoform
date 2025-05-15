cd /lustre/tmchen4/ref
grep -w 'GAPDH' Homo_sapiens-GCA_009914755.4-2022_07-genes.gtf | grep 'transcript' | head -1
###TSS -1000 to TSS +100
echo -e "12\t6543868\t6544968\tGAPDH_promoter\t0\t+" > /lustre/tmchen4/ref/canonical_GAPDH_promoter_hs1.bed
#bam to bw
#!/bin/bash
#SBATCH --job-name=bash
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --partition=quick
#SBATCH --ntasks-per-node=30
#SBATCH --error=%j.err
#SBATCH --output=%j.out

source ~/.bashrc
conda activate FLAMES
# Input and output directories
INPUT_DIR="/lustre/tmchen4/chipseq_cell/30-1153071837/results"
OUTPUT_DIR="/lustre/tmchen4/chipseq_cell/bw_file"

# Create output dir if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all nodup.bam files
find "$INPUT_DIR" -type f -name "*nodup.bam" | while read bam; do
    base=$(basename "$bam" .bam)
    out_bw="${OUTPUT_DIR}/${base}.bw"

    echo "Processing $bam → $out_bw"

    bamCoverage \
      -b "$bam" \
      -o "$out_bw" \
      --normalizeUsing CPM \
      --binSize 25 \
      --extendReads \
      -p 30
done


