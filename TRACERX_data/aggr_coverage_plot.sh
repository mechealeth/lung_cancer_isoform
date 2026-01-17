module load anaconda3
source activate /scratch/user/uqmche33/myenv/bambu_env

wiggletools mean $(cat LUAD.txt)   > LUAD.mean.wig
wiggletools mean $(cat LUSC.txt)   > LUSC.mean.wig
wiggletools mean $(cat Normal.txt) > Normal.mean.wig

wigToBigWig LUAD.mean.wig   /path/to/chrom.sizes LUAD.mean.bw
wigToBigWig LUSC.mean.wig   /path/to/chrom.sizes LUSC.mean.bw
wigToBigWig Normal.mean.wig /path/to/chrom.sizes Normal.mean.bw
