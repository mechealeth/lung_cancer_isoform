#L1HS 
L1HS.gtf

L1HS	custom	gene	        1	6037	.	+	.	gene_id "L1HS"; gene_name "L1HS";
L1HS	custom	transcript	    1	6037	.	+	.	gene_id "L1HS"; transcript_id "L1HS_tx";
L1HS	custom	exon	        1	910	    .	+	.	gene_id "L1HS"; transcript_id "L1HS_tx"; exon_number "1"; feature "5UTR";
L1HS	custom	exon	        911	2450	.	+	.	gene_id "L1HS"; transcript_id "L1HS_tx"; exon_number "2"; feature "ORF1";
L1HS	custom	exon	        2451	2700	.	+	.	gene_id "L1HS"; transcript_id "L1HS_tx"; exon_number "3"; feature "interORF";
L1HS	custom	exon	        2701	5950	.	+	.	gene_id "L1HS"; transcript_id "L1HS_tx"; exon_number "4"; feature "ORF2";
L1HS	custom	exon	        5951	6037	.	+	.	gene_id "L1HS"; transcript_id "L1HS_tx"; exon_number "5"; feature "3UTR";



##align to L1HS

minimap2 -ax map-hifi /lustre/tmchen4/ref/L1consensus.mmi NL20_Pacbio.fastq.gz | samtools sort -o NL20_Pacbio_L1_aligned.bam
samtools index NL20_Pacbio_L1_aligned.bam

minimap2 -ax map-hifi /lustre/tmchen4/ref/L1consensus.mmi HARA_Pacbio.fastq.gz | samtools sort -o HARA_Pacbio_L1_aligned.bam
samtools index HARA_Pacbio_L1_aligned.bam

minimap2 -ax map-hifi /lustre/tmchen4/ref/L1consensus.mmi A549_D.fastq.gz | samtools sort -o A549_Pacbio_L1_aligned.bam
samtools index A549_Pacbio_L1_aligned.bam

##Coverage plot

bedtools genomecov -d -ibam NL20_Pacbio_L1_aligned.bam -g <(echo -e "L1consensus\t6208") > NL20_Pacbio_L1_coverage.txt
samtools view -c -F 260 NL20_Pacbio_L1_aligned.bam
# Sample 2
bedtools genomecov -d -ibam HARA_Pacbio_L1_aligned.bam -g <(echo -e "L1consensus\t6208") > HARA_Pacbio_L1_coverage.txt
samtools view -c -F 260 HARA_Pacbio_L1_aligned.bam

# Sample 3
bedtools genomecov -d -ibam A549_Pacbio_L1_aligned.bam -g <(echo -e "L1consensus\t6208") > A549_Pacbio_L1_coverage.txt
samtools view -c -F 260 A549_Pacbio_L1_aligned.bam

###plot (CPM)



bedtools genomecov -d -ibam your.bam -strand + > sample_sense.txt
bedtools genomecov -d -ibam your.bam -strand - > sample_antisense.txt


mapped_reads = {
    "NL20": 5097,
    "HARA": 16101,
    "A549": ls L1HS
}

import pandas as pd
import matplotlib.pyplot as plt

samples = ["NL20", "HARA", "A549"]
colors = ["gray", "gray", "gray"]

# Total mapped reads per sample (fill with actual counts)
mapped_reads = {
    "sample1": 5097,
    "sample2": 16101,
    "sample3": 1300000
}

fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(7, 4), sharex=True)

for i, sample in enumerate(samples):
    df = pd.read_csv(f"{sample}_Pacbio_L1_coverage.txt", sep="\t", header=None, names=["chr", "pos", "cov"])
    
    # Normalize to CPM
    cpm_scale = mapped_reads[sample] / 1e6
    df["cpm"] = df["cov"] / cpm_scale

    ax = axes[i]
    ax.fill_between(df["pos"], df["cpm"], color=colors[i])
    ax.set_xlim([0, 6208])
    ax.set_ylim([0, df["cpm"].max() * 1.1])
    ax.set_ylabel(f"{sample}", fontsize=9, rotation=0, labelpad=40)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

axes[-1].set_xticks(range(0, 7000, 1000))
axes[-1].set_xlabel("Human L1 consensus sequence (bp)")

# Arrow + label
axes[-1].annotate('', xy=(0, -5), xytext=(6208, -5),
                  arrowprops=dict(arrowstyle='<->', lw=1.5), annotation_clip=False)
axes[-1].text(3104, -10, '6208 bp', ha='center', va='center', fontsize=10)

# Panel label
fig.text(0.01, 0.95, 'H', fontsize=14, fontweight='bold')

plt.tight_layout()
plt.savefig("L1_coverage_CPM_three_samples.pdf", dpi=300)
plt.show()


# Get reads that fully map within 1–910 on L1consensus

samtools view NL20_Pacbio_L1_aligned.bam L1HS:1-910 | \
awk '{ 
  if ($4 >= 1 && ($4 + length($10) - 1) <= 910) 
    print ">"$1"\n"$10 
}' > NL20_reads_within_1_910.fasta

samtools view HARA_Pacbio_L1_aligned.bam L1HS:1-910 | \
awk '{ 
  if ($4 >= 1 && ($4 + length($10) - 1) <= 910) 
    print ">"$1"\n"$10 
}' > HARA_reads_within_1_910.fasta

###get bamfile
samtools view -h NL20_Pacbio_L1_aligned.bam L1HS:1-910 | \
awk 'BEGIN {OFS="\t"} 
     /^@/ {print; next} 
     {
       start = $4
       len = length($10)
       if (start >= 1 && (start + len - 1) <= 910)
         print
     }' | \
samtools view -b -o NL20_within_1_910.bam
samtools sort -o NL20_within_1_910.sorted.bam NL20_within_1_910.bam
samtools index NL20_within_1_910.sorted.bam


samtools view -h HARA_Pacbio_L1_aligned.bam L1HS:1-910 | \
awk 'BEGIN {OFS="\t"} 
     /^@/ {print; next} 
     {
       start = $4
       len = length($10)
       if (start >= 1 && (start + len - 1) <= 910)
         print
     }' | \
samtools view -b -o HARA_within_1_910.bam
samtools sort -o HARA_within_1_910.sorted.bam HARA_within_1_910.bam
samtools index HARA_within_1_910.sorted.bam


samtools view -h A549_Pacbio_L1_aligned.bam L1HS:1-910 | \
awk 'BEGIN {OFS="\t"} 
     /^@/ {print; next} 
     {
       start = $4
       len = length($10)
       if (start >= 1 && (start + len - 1) <= 910)
         print
     }' | \
samtools view -b -o A549_within_1_910.bam
samtools sort -o A549_within_1_910.sorted.bam A549_within_1_910.bam
samtools index A549_within_1_910.sorted.bam
