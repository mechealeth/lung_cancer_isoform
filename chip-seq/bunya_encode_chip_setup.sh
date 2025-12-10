module load anaconda3
conda create --prefix /scratch/user/uqmche33/conda/encode_chip
source activate /scratch/user/uqmche33/conda/encode_chip
conda install bioconda::caper
mkdir -p ~/.caper
#this step will generate default Caper configuration file ~/.caper/default.conf
caper init slurm
#edit Caper configuration file ~/.caper/default.conf
vim ~/.caper/default.conf
--partition=general
--account=a_faulkner
cd /scratch/project/te_lung/
mkdir chip-seq
cd chip-seq

######slurm script######multiple processing
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=60
#SBATCH --mem=200G
#SBATCH --job-name=chipseq
#SBATCH --time=14:00:00
#SBATCH --partition=general
#SBATCH --account=a_faulkner
#SBATCH -o slurm-%j.output
#SBATCH -e slurm-%j.error

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
conda activate encd-chip
module load singularity/3.8.3


for fq in /lustre/tmchen4/H1299/ChIP-Seq/*_R1_1.fastq.gz;do
    pre=${fq%_R1_1.fastq.gz}
    caper run /scratch/project/te_lung/chip-seq/chip-seq-pipeline2/chip.wdl -i ${pre}.json --singularity
done



echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 30


