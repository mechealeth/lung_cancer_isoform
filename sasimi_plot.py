#!/usr/bin/bash
set -euo pipefail

# Recursively find all matching BAMs (change "." to a specific dir if you want)
find . -type f -name "*_Pacbio.hs1.aligned.bam" -print0 |
while IFS= read -r -d '' bam; do
  base="${bam%.bam}"
  junc="${base}.junc.bed"
  bedpe="${base}.junc.bedpe"

  echo "[*] Processing: $bam"

  # Ensure BAM is indexed (regtools expects coordinate-sorted, indexed BAM)
  if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
    samtools index -@ 8 "$bam"
  fi

  # Junctions (BED12-like); tweak -a/-m/-M if needed
   regtools junctions extract -s XS -t ts -a 3 -m 20 -M 1000000 "$bam" > "$junc"

  # Convert to BEDPE for pyGenomeTracks links (score = junction count)
  awk 'BEGIN{OFS="\t"}{
         # regtools: chrom start end name score strand ...
         print $1, $2, $2+1, $1, $3-1, $3, $4, $5, $6, $6
       }' "$junc" > "$bedpe"

  # (optional) filter weak junctions for a cleaner plot, e.g. score >= 5
  # awk '$8>=5' "$bedpe" > "${base}.junc.min5.bedpe"
done


regtools junctions extract -a 8 -m 50 -M 500000 -s XS -t ts 07PM0575.normal_Pacbio.hs1.aligned.bam > 07PM0575.normal.junc.bed


###inifile ###

make_tracks_file --trackFiles \
    /lustre/tmchen4/iso-seq/tissue_sample/filter.filtered.gtf \
    /lustre/tmchen4/ref/T2T_CHM13_v2_rmsk_TE_chr.gtf \
    /lustre/tmchen4/iso-seq/allin/align/NL20_DMSO_Rep1.bw \
    /lustre/tmchen4/iso-seq/allin/align/NL20_DMSO_Rep3.bw\
    /lustre/tmchen4/iso-seq/allin/align/NL20_DNMTi_D4_Rep1.bw\
    /lustre/tmchen4/iso-seq/allin/align/NL20_DNMTi_D4_Rep2.bw \
    /lustre/tmchen4/iso-seq/allin/align/HARA_DMSO_Rep1.bw \
    /lustre/tmchen4/iso-seq/allin/align/HARA_DMSO_Rep2.bw\
    /lustre/tmchen4/iso-seq/allin/align/HARA_DNMTi_D4_Rep1.bw \
    /lustre/tmchen4/iso-seq/allin/align/HARA_DNMTi_D4_Rep2.bw \
    /lustre/tmchen4/iso-seq/allin/align/12PM1049.normal.bw \
    /lustre/tmchen4/iso-seq/allin/align/12PM1049.tumour.bw \
    /lustre/tmchen4/iso-seq/allin/align/13PM0352.normal.bw \
    /lustre/tmchen4/iso-seq/allin/align/13PM0352.tumour.bw \
    /lustre/tmchen4/iso-seq/allin/align/14PM0263.normal.bw \
    /lustre/tmchen4/iso-seq/allin/align/14PM0263.tumour.bw \
    /lustre/tmchen4/iso-seq/allin/align/21DSRB006.normal.bw \
    /lustre/tmchen4/iso-seq/allin/align/21DSRB006.tumour.bw \
    --out cl_tiss_genome_tracks.ini



pyGenomeTracks --tracks patient_genome_tracks.ini --region chr1:154758828-154792136 --outFileName patient-erk.pdf

normal :#3246a8
tumor:  #ed1a1a
