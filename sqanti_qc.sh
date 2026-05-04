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



module load sqanti3/5.2
python /sw/local/rocky8/noarch/qcif/software/SQANTI3-5.2/sqanti3_qc.py  /QRISdata/Q7820/i3N/bambu/extended_annotations.gtf /scratch/user/uqmche33/Ref/Homo_sapiens_chr-GCA_009914755.4-2022_07-genes.gtf  /scratch/user/uqmche33/Ref/hs1.fa \
                     --polyA_motif_list /scratch/user/uqmche33/Ref/data/polyA_motifs/mouse_and_human.polyA_motif.txt    \
                     -o isoseq -d /QRISdata/Q7820/i3N/bambu/SQANTI3_output_hs1SJ.out.tab     \
		             --force_id_ignore \
                     --CAGE_peak /scratch/user/uqmche33/Ref/data/ref_TSS_annotation/hglft_refTSS_v3.1_human_coordinate.hg1.bed \
                     --fl_count /QRISdata/Q7820/i3N/bambu/fullLengthCounts_transcript.csv \
                     --cpus 40 --report both



python /sw/local/rocky8/noarch/qcif/software/SQANTI3-5.2hea/sqanti3_filter.py rules -j /lustre/tmchen4/iso-seq/bambu/bambu_NDR0.155/filter_rule.json tissue_isoseq_classification.txt \
	                           --gtf /lustre/tmchen4/iso-seq/tissue_sample/SQANTI3_output_hs1SJ.out.tab/tissue_isoseq_corrected.gtf \
                                -o filter

