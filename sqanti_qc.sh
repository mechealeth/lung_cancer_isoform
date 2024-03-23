# checek gtf strand Before input the bambu out gtf into Sqanti4 qc (Sqanti don't support unstarnded isoform)
cut -f7 extend_annotations_corrected.gtf | sort | uniq -c
#If there are too many unstrand isoform need to check carefully to delete or not
#if only few isoform unstarnded and the expression is quite loa ang only present in very few sample, I chose to delete this unstarnd isoform
awk '$7 != "." || NR < 3' extended_annotations.gtf  >  extend_annotations_corrected.gtf
### the following is the Sqanti3 QC sh file
source ~/.bashrc
source activate SQANTI3.env

cd /lustre/tmchen4/iso-seq/SQANTI3-5.2

export PYTHONPATH=$PYTHONPATH:/lustre/tmchen4/iso-seq/SQANTI3-5.2/cDNA_Cupcake/sequence/
export PYTHONPATH=$PYTHONPATH:/lustre/tmchen4/iso-seq/SQANTI3-5.2/cDNA_Cupcake/

python sqanti3_qc.py  /lustre/tmchen4/iso-seq/cell_line/extend_annotations_corrected.gtf  /lustre/tmchen4/ref/Homo_sapiens_chr-GCA_009914755.4-2022_07-genes.gtf  /lustre/tmchen4/ref/hs1.fa \
                     --polyA_motif_list /lustre/tmchen4/iso-seq/SQANTI3-5.2/data/polyA_motifs/mouse_and_human.polyA_motif.txt    \
                     -o isoseq -d /lustre/tmchen4/iso-seq/cell_line/SQANTI3_output_hs1SJ.out.tab     \
		                 --force_id_ignore \
                     --CAGE_peak /lustre/tmchen4/iso-seq/SQANTI3-5.2/data/ref_TSS_annotation/hglft_refTSS_v3.1_human_coordinate.hg1.bed \
                     --fl_count /lustre/tmchen4/iso-seq/cell_line/fullLengthCounts_transcript.csv \
                     --cpus 30 --report both