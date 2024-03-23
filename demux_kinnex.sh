HARA	bc01	bcM0001										
NL20	bc02	bcM0001										
HARA_DMSO	bc03	bcM0001										
HARA_DNMTi_D2	bc04	bcM0001										
HARA_DNMTi_D4	bc05	bcM0001										
HARA_DNMTi_D5	bc06	bcM0001										
NL20_DMSO	bc07	bcM0001										
NL20_DNMTi_D2	bc08	bcM0001										
NL20_DNMTi_D4	bc09	bcM0001										
NL20_DNMTi_D4	bc10	bcM0002										
HARA_DMSO	bc11	bcM0002										
HARA_DNMTi_D2	bc12	bcM0002										
HARA_DNMTi_D4	bc01	bcM0002										
HARA_DNMTi_D5	bc02	bcM0002										
NL20_DMSO	bc03	bcM0002										
NL20_DNMTi_D2	bc04	bcM0002										
NL20_DNMTi_D4	bc05	bcM0002					



#!/bin/bash
#$ -S /bin/bash
#$ -terse
#$ -cwd
#$ -N masseq
#$ -l mem_requested=4G
#$ -l h_vmem=4G
#$ -l tmp_requested=4G
#$ -pe smp 16

THREADS=16

## sample name
SAMPLE=${1}
BARCODE=${2}

## names and locations for input and output files
MASSEQ_PRIMERS=/directflow/KCCGGenometechTemp/projects/iradev/masseq/ref_files/MAS-Seq_Adapter_v3__MAS8.fa
ISOSEQ_BARCODES=/directflow/KCCGGenometechTemp/projects/iradev/masseq/ref_files/IsoSeqX_bc01_5p.fa

## load Tim's conda environment
module add centos7.8/timamo/miniconda3/22.11.1 2> /dev/null; conda init bash > /dev/null; source /home/timamo/.bashrc; eval "$(conda shell.bash hook)"; conda activate snakemake
module load centos7.8/phuluu/samtools/1.11

mkdir ${SAMPLE}_${BARCODE}_kinnex
cd ${SAMPLE}_${BARCODE}_kinnex

## de-concatenate HiFi reads generated with Kinnex protocol into individual read segments, using Skera
skera split --num-threads ${THREADS} ../${SAMPLE}.hifi_reads.${BARCODE}.bam ${MASSEQ_PRIMERS} ${SAMPLE}.${BARCODE}.segmented.bam

## demux sample barcodes and arrange all transcripts into 5'->3' direction, using Lima with --isoseq preset
lima --num-threads ${THREADS} ${SAMPLE}.${BARCODE}.segmented.bam ${ISOSEQ_BARCODES} ${SAMPLE}.${BARCODE}.segmented.lima.bam --isoseq

### process
ls movie*.flnc.bam movie*.flnc.bam movie*.flnc.bam > flnc.fofn
 
#check
find /lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_pacbio -name \*.bcM0001.segmented.lima.IsoSeqX_bc01_5p--IsoSeqX_3p.bam -print


HARA								
NL20								
HARA_DMSO_R1									
HARA_DNMTi_D2_R1									
HARA_DNMTi_D4_R1								
HARA_DNMTi_D5_R1								
NL20_DMSO_R1								
NL20_DNMTi_D2_R1								
NL20_DNMTi_D4_R1									
NL20_DNMTi_D4_R2							
HARA_DMSO_R2								
HARA_DNMTi_D2_R2								
HARA_DNMTi_D4_R2								
HARA_DNMTi_D5_R2								
NL20_DMSO_R3								
NL20_DNMTi_D2_R3								
NL20_DNMTi_D4_R3	

#merge 
find /lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_pacbio -name \*.bcM0002.segmented.lima.IsoSeqX_bc03_5p--IsoSeqX_3p.bam > NL20_DMSO_R3.flnc.fofn
pbmerge HARA.flnc.fofn > HARA.flnc.bam
find /lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_pacbio -name \*.bcM0001.segmented.lima.IsoSeqX_bc02_5p--IsoSeqX_3p.bam > NL20.flnc.fofn
pbmerge NL20.flnc.fofn > NL20.flnc.bam
find /lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_pacbio -name \*.bcM0001.segmented.lima.IsoSeqX_bc03_5p--IsoSeqX_3p.bam > HARA_DMSO.flnc.fofn
find /lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_pacbio -name \*.bcM0002.segmented.lima.IsoSeqX_bc11_5p--IsoSeqX_3p.bam >> HARA_DMSO.flnc.fofn
pbmerge HARA_DMSO.flnc.fofn > HARA_DMSO.flnc.bam
find /lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_pacbio -name \*.bcM0001.segmented.lima.IsoSeqX_bc04_5p--IsoSeqX_3p.bam > HARA_DNMTi_D2.flnc.fofn
find /lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_pacbio -name \*.bcM0002.segmented.lima.IsoSeqX_bc12_5p--IsoSeqX_3p.bam >> HARA_DNMTi_D2.flnc.fofn
pbmerge HARA_DNMTi_D2.flnc.fofn > HARA_DNMTi_D2.flnc.bam
find /lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_pacbio -name \*.bcM0001.segmented.lima.IsoSeqX_bc05_5p--IsoSeqX_3p.bam > HARA_DNMTi_D4.flnc.fofn
find /lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_pacbio -name \*.bcM0002.segmented.lima.IsoSeqX_bc01_5p--IsoSeqX_3p.bam >> HARA_DNMTi_D4.flnc.fofn
pbmerge HARA_DNMTi_D4.flnc.fofn > HARA_DNMTi_D4.flnc.bam
find /lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_pacbio -name \*.bcM0001.segmented.lima.IsoSeqX_bc06_5p--IsoSeqX_3p.bam > HARA_DNMTi_D5.flnc.fofn
find /lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_pacbio -name \*.bcM0002.segmented.lima.IsoSeqX_bc02_5p--IsoSeqX_3p.bam >> HARA_DNMTi_D5.flnc.fofn
pbmerge HARA_DNMTi_D5.flnc.fofn > HARA_DNMTi_D5.flnc.bam
find /lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_pacbio -name \*.bcM0001.segmented.lima.IsoSeqX_bc07_5p--IsoSeqX_3p.bam > NL20_DMSO.flnc.fofn
find /lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_pacbio -name \*.bcM0002.segmented.lima.IsoSeqX_bc03_5p--IsoSeqX_3p.bam >> NL20_DMSO.flnc.fofn
pbmerge NL20_DMSO.flnc.fofn > NL20_DMSO.flnc.bam
find /lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_pacbio -name \*.bcM0001.segmented.lima.IsoSeqX_bc08_5p--IsoSeqX_3p.bam > NL20_DNMTi_D2.flnc.fofn
find /lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_pacbio -name \*.bcM0002.segmented.lima.IsoSeqX_bc04_5p--IsoSeqX_3p.bam >> NL20_DNMTi_D2.flnc.fofn
pbmerge NL20_DNMTi_D2.flnc.fofn > NL20_DNMTi_D2.flnc.bam
find /lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_pacbio -name \*.bcM0001.segmented.lima.IsoSeqX_bc09_5p--IsoSeqX_3p.bam > NL20_DNMTi_D4.flnc.fofn
find /lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_pacbio -name \*.bcM0002.segmented.lima.IsoSeqX_bc10_5p--IsoSeqX_3p.bam >> NL20_DNMTi_D4.flnc.fofn
find /lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_pacbio -name \*.bcM0002.segmented.lima.IsoSeqX_bc05_5p--IsoSeqX_3p.bam >> NL20_DNMTi_D4.flnc.fofn
pbmerge NL20_DNMTi_D4.flnc.fofn >> NL20_DNMTi_D4.flnc.bam





#list sample bamfile
find /lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_pacbio -name \*.${sadapter_barcode}.segmented.lima.IsoSeqX_${sample_barcode}_5p--IsoSeqX_3p.bam > ${sample}.flnc.fofn
#merger bamfile from different flowcell
pbmerge -j 30 -o ${sample}.bam ${sample}.flnc.fofn
#pacbio formatted bamfile to fastq
bam2fastq -o ${sample}  -j 30 ${sample}.fastq.gz
#minimap2 mapping / sam to bam /samtools sort
minimap2 -ax splice:hq -uf -t 30 /lustre/tgfaulkn/Drop/Minchun/hs1.fa  ${sample}.fastq.gz | samtools view -bh - | samtools sort -@ 30 - -o  ${sample}.hs1.aligned.bam
#quick qc
cramino ${sample}.hs1.aligned.bam

##subsample
samtools view -s 0.2 -b HARA_DMSO_R2.hs1.aligned.bam > HARA_DMSO_R2_0.2.hs1.aligned.bam
samtools view -s 0.4 -b HARA_DMSO_R2.hs1.aligned.bam > HARA_DMSO_R2_0.4.hs1.aligned.bam  
samtools view -s 0.6 -b HARA_DMSO_R2.hs1.aligned.bam > HARA_DMSO_R2_0.6.hs1.aligned.bam
samtools view -s 0.8 -b HARA_DMSO_R2.hs1.aligned.bam > HARA_DMSO_R2_0.8.hs1.aligned.bam
samtools view -s 0.3 -b HARA_DMSO_R2.hs1.aligned.bam > HARA_DMSO_R2_0.3.hs1.aligned.bam
samtools view -s 0.5 -b HARA_DMSO_R2.hs1.aligned.bam > HARA_DMSO_R2_0.5.hs1.aligned.bam
samtools view -s 0.9 -b HARA_DMSO_R2.hs1.aligned.bam > HARA_DMSO_R2_0.9.hs1.aligned.bam
samtools view -s 0.15 -b HARA_DMSO_R2.hs1.aligned.bam > HARA_DMSO_R2_0.15.hs1.aligned.bam
samtools view -s 0.55 -b HARA_DMSO_R2.hs1.aligned.bam > HARA_DMSO_R2_0.55.hs1.aligned.bam
#!/bin/bash
#SBATCH --job-name=bash
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --partition=medium
#SBATCH --ntasks-per-node=30
#SBATCH --error=%j.err
#SBATCH --output=%j.out

CURDIR=`pwd`
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID
NODES=`scontrol show hostnames $SLURM_JOB_NODELIST`
for i in $NODES
do
echo "$i:$SLURM_NTASKS_PER_NODE" >> $CURDIR/nodelist.$SLURM_JOB_ID
done
echo $SLURM_NPROCS

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"

echo "start .."
#setting environment
source ~/.bashrc
source activate FLAMES

cd /lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_pacbio

cat /lustre/tgfaulkn/Drop/Minchun/lung_cancer_cell_line_pacbio/sample.txt | while read id
do
name=$(echo $id )
echo $name

minimap2 -ax splice:hq -uf -t 30 /lustre/tgfaulkn/Drop/Minchun/hs1.fa ${name}.fastq.gz | samtools view -bh - | samtools sort -@ 30 - -o ${name}.hs1.aligned.bam	 

done


echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID