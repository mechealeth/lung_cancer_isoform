#this script better run in ipython or jupyter
#source ~/.bashrc
#conda activate FLAMEs
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode
def extract_sequence(genome_file, chrom, start, end):
    # Load the genome from a FASTA file
    genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
    
    # Extract the sequence from the specified chromosome and coordinates
    sequence = genome[chrom].seq[start:end]
    return str(sequence)

genome_file = "/lustre/tmchen4/ref/hs1.fa"  # Replace with the path to your genome file
chromosome = "chr12"  # Replace with the chromosome name

# Coordinates of the larger fragment (replace with your actual coordinates)
large_start = 121989307
large_end = 121990791
longer_seq = extract_sequence(genome_file, chromosome, large_start, large_end)
context = 10000
paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
models = [load_model(resource_filename('spliceai', x)) for x in paths]
x = one_hot_encode('N'*(context//2) + longer_seq + 'N'*(context//2))[None, :]
y = np.mean([models[m].predict(x) for m in range(5)], axis=0)
acceptor_prob = y[0, :, 1]
donor_prob = y[0, :, 2]
#sequence sel alignment to find SVA highly repeat region which may be the VNTR region
def plot_self_alignment_dotplot_with_splicing(sequence, window=10):
    # Generate dot matrix for forward vs forward
    forward_dot_matrix = generate_dotplot(sequence, sequence, window)
    
    # Generate dot matrix for forward vs reverse
    reverse_complement = str(Seq(sequence).reverse_complement())
    reverse_dot_matrix = generate_dotplot(sequence, reverse_complement, window)
    
    fig, ax1 = plt.subplots(figsize=(12, 9))
    
    # Plot the dot plot (forward vs forward and forward vs reverse)
    coords_forward = np.argwhere(forward_dot_matrix == 1)
    ax1.scatter(coords_forward[:, 1], coords_forward[:, 0], s=1, color='red', label='Forward')
    
    coords_reverse = np.argwhere(reverse_dot_matrix == 1)
    ax1.scatter(coords_reverse[:, 1], coords_reverse[:, 0], s=1, color='blue', label='Reverse')
    
    ax1.set_xlabel('Sequence Position', fontsize=16)
    ax1.set_ylabel('Sequence Position',  fontsize=16)
    ax1.tick_params(axis='y')
    ax1.legend(loc='upper left')
    
    # Create a second y-axis for splicing probability
    ax2 = ax1.twinx()
    ax2.plot(acceptor_prob, color='black', label='Acceptor Probability')
    ax2.plot(donor_prob, color='blue', label='Donor Probability')
    
    ax2.set_ylabel('Splicing Probability', color='black', fontsize=16)
    ax2.tick_params(axis='y', labelcolor='black')
    
    # Add title and format the plot
    plt.title('Self Align and Splicing Probability plot', fontsize=20)
    
    # Adding legend
    ax2.legend(loc='upper right')
    
    plt.show()

# Run the plot function with the sequence and splicing probability data
plot_self_alignment_dotplot_with_splicing(longer_seq, window=10)

# acording to the SVA annotation result and  spicing AI spilcing site prediction plot SVA derived novel transcripts TSS site in SVA
# Plotting the data
fig, ax1 = plt.subplots(figsize=(20, 10))
# Plotting the first dataset
ax1.plot(acceptor_prob, color='black', label='acceptor_prob')

# Plotting the second dataset
ax1.plot(donor_prob, color='blue', label='donor_prob')

# Adding labels and title
ax1.set_title('Splicing probability',fontsize=20)
ax1.set_xlabel('Base Position',fontsize=16)
ax1.set_ylabel('probability', color='black',fontsize=16)

# Highlighting and annotating specific points in the first dataset
highlight_points1 = [98,140,137,107,115,688,671,746,747]  # example indices to highlight in the first dataset
annotations1 = ["TSS", "TSS","TSS","TSS","TSS","SS","SS","SS","SS"]

for i, point in enumerate(highlight_points1):
    ax1.annotate(annotations1[i], 
                 xy=(point, acceptor_prob[point]), 
                 xytext=(point + 10, acceptor_prob[point] + 0.02),  # Position offset
                 textcoords='data', 
                 arrowprops=dict(arrowstyle="->", color='red' if i == 0 else 'green'),
                 color='red' if i == 0 else 'green',
                 fontsize=14)

# Annotating x-axis regions
ax1.axvspan(1528, 2000, color='yellow', alpha=0.3, label='Alu-like')
ax1.axvspan(482, 1527, color='lightgreen', alpha=0.3, label='VNTR')
ax1.axvspan(1, 481, color='red', alpha=0.3, label='SINE-R')
# Adding legend
ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=5,fontsize=16)
#ax1.legend(fontsize=12, loc='best', frameon=True, fancybox=True)
# Removing the frame
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.spines['left'].set_visible(False)


# Removing ticks
ax1.xaxis.set_ticks_position('none') 
ax1.yaxis.set_ticks_position('none') 
ax1.xaxis.set_tick_params(length=0,labelsize=12)
ax1.yaxis.set_tick_params(length=0,labelsize=12)

# Displaying the plot
plt.show()

