#!/bin/bash
#SBATCH --job-name=bash
#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH --partition=quick
#SBATCH --ntasks-per-node=8
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
source activate SQANTI3.env

cd /lustre/tmchen4/iso-seq/tissue_sample/SQANTI3_output_hs1SJ.out.tab

export PYTHONPATH=$PYTHONPATH:/lustre/tmchen4/iso-seq/SQANTI3-5.2/cDNA_Cupcake/sequence/
export PYTHONPATH=$PYTHONPATH:/lustre/tmchen4/iso-seq/SQANTI3-5.2/cDNA_Cupcake/


python /lustre/tmchen4/iso-seq/SQANTI3-5.2/sqanti3_filter.py rules -j /lustre/tmchen4/iso-seq/bambu/bambu_NDR0.155/filter_rule.json tissue_isoseq_classification.txt \
	                           --gtf /lustre/tmchen4/iso-seq/tissue_sample/SQANTI3_output_hs1SJ.out.tab/tissue_isoseq_corrected.gtf \
                                -o filter

echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 30
