#!/bin/bash
#SBATCH --job-name=make_bw
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --partition=quick
#SBATCH --ntasks-per-node=60
#SBATCH --error=%j.err
#SBATCH --output=%j.out

source ~/.bashrc
conda activate FLAMES

# === Input/Output Directories ===
INPUT_PATTERN="/lustre/tmchen4//chip_new/N2530988_30-1235730417_ChIP_2025-12-03/result/*/chip/*/call-filter_ctl/*/execution/*nodup.bam"
OUTPUT_DIR="/lustre/tmchen4/chipseq_cell/bw_file/new_H3K27ac"

# === Find BAMs, deduplicate by unique basename ===
find $INPUT_PATTERN -type f | while read -r bam; do
    base=$(basename "$bam")
    echo "$base|$bam"
done | sort -t'|' -k1,1 -u | cut -d'|' -f2 | while read -r bam; do

    base=$(basename "$bam" .bam)
    out_bw="${OUTPUT_DIR}/${base}.bw"

    echo "Processing $bam -> $out_bw"

    bamCoverage \
        -b "$bam" \
        -o "$out_bw" \
        --normalizeUsing CPM \
        --binSize 25 \
        --extendReads \
        -p 60
done
