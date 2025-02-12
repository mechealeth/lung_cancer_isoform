#TRI HPC
source ~/.bashrc
conda actiuvate FLAMES
# care about spacing
#bam2Bingwing
#bamfile must index before bam2bw
samtools index NL20_DNMTi_D4_Rep2_Pacbio.hs1.aligned.bam
bamCoverage --bam NL20_DNMTi_D4_Rep2_Pacbio.hs1.aligned.bam -o NL20_DNMTi_D4_Rep2.bw --normalizeUsing CPM
#after finish preprocess of the input filen generate configuration file for pyGenomeTracks
make_tracks_file --trackFiles \
    /lustre/tmchen4/iso-seq/tissue_sample/SQANTI3_output_hs1SJ.out.tab/filter.filtered.gtf \
    /lustre/tmchen4/ref/T2T_CHM13_v2_rmsk_TE_chr.gtf \
    /lustre/tmchen4/H1299/Transcripts/LR_CAGE/H1299_DMSO_LRCAG.bw \
    /lustre/tmchen4/H1299/Transcripts/LR_CAGE/H1299_DACSB_LRCAGE.bw \
    /lustre/tmchen4/H1299/ng_wang/chip/f01774f2-30e3-4d2d-aa65-c69ca73f9b52/call-macs2_signal_track/shard-0/execution/DMSO_H3K4me1.srt.nodup_x_DMSO_input.srt.nodup.pval.signal.bigwig \
    /lustre/tmchen4/H1299/ng_wang/chip/3c7ae6e0-6582-4a1b-ab63-a84e52fdff31/call-macs2_signal_track/shard-0/execution/DACSB_H3K4me1.srt.nodup_x_DACSB_Input.srt.nodup.pval.signal.bigwig \
    /lustre/tmchen4/H1299/ng_wang/chip/69ec2cde-c239-4be0-b930-054fada05260/call-macs2_signal_track/shard-0/execution/DMSO_H3K4me3.srt.nodup_x_DMSO_input.srt.nodup.pval.signal.bigwig \
    /lustre/tmchen4/H1299/ng_wang/chip/6f79b64e-198c-4c8f-b2b2-421eeba0bccc/call-macs2_signal_track/shard-0/execution/DACSB_H3K4me3.srt.nodup_x_DACSB_Input.srt.nodup.pval.signal.bigwig \
    /lustre/tmchen4/H1299/ng_wang/chip/d151c0d0-48b6-4bc1-bca5-92b21a225f0d/call-macs2_signal_track/shard-0/execution/DMSO_H3K27ac.srt.nodup_x_DMSO_input.srt.nodup.pval.signal.bigwig \
    /lustre/tmchen4/H1299/ng_wang/chip/8166d340-1fcb-40de-b7a1-85bf8aae2d8e/call-macs2_signal_track/shard-0/execution/DACSB_H3K27ac.srt.nodup_x_DACSB_Input.srt.nodup.pval.signal.bigwig \
    /lustre/tmchen4/H1299/ng_wang/chip/2301b736-9c05-4605-8dc1-236148ac58d6/call-macs2_signal_track/shard-0/execution/DMSO_H3K9me3.srt.nodup_x_DMSO_input.srt.nodup.pval.signal.bigwig \
    /lustre/tmchen4/H1299/ng_wang/chip/41eee49d-7cb3-409c-a3a1-694c2608454e/call-macs2_signal_track/shard-0/execution/DACSB_H3K9me3.srt.nodup_x_DACSB_Input.srt.nodup.pval.signal.bigwig \
    --out histone_genome_tracks.ini
#custimized ini file 
vi histone_genome_tracks.ini  # can change track name & track y axis maxium value & track color & track spacing
# after define the plot style and generate track fig
pyGenomeTracks --tracks histone_genome_tracks.ini --region chr12:121984033-122001272 --outFileName histone_image.pdf





