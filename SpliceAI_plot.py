import matplotlib.pyplot as plt
import numpy as np
from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode
import numpy as np
input_sequence = "CACGGTCTCCCTCTCATGCGGAGCCGAAGCTGGACTGTACTGCTGCCATCTCGGCTCACTGCAACCTCCCTGCCTGATTCTCCTGCCTCAGCCTGCCGAGTGCCTGCGATTGCAGGCACGCGCCGCCACGCCTGACTGGTTTTGGTGGAGACCGGGTTTCGCTGTGTTGGCCGGGCCGGTCTCCAGCCCCTAACCGCGAGTGATCCGCCAACCTCGGCCTCCCGAGGTGCCGGGATTGCAGACGGAGTCTCGTTCACTCAGTGCTCAATGGTGCCCAGGCTGGAGTGCAGTGGCGTGATCTCGGCTCACTACAACCTACACCTCCCAGCCGCCTGCCTTGGCCTCCCAAAGTGCCGAGATTGCAGCCTCTGCCCGGCCGCCACCCCGTCTGGGAAGTGAGGAGTGTCTCTGCCTGGCCGCCCATCGTCGGGGATGTGAGGAGCCCCTCTGCCTGGCTGCCCAGTCTGGAAAGTGAGGAGCGTCTCCACCCGGCCGCCATCCCATCTAGGAAGTGAGGAGCGCCTCTTCCCAGCCGCCATCACATCTAGGAAGTGAGGAGCGTCTCTGCCCGGCCGCCCATCGTCTGAGATGTGGGGAGCGCCTCTGCCCCGCCGCCCCATCTGGGATGTGAGGAGTGCCTCTGCCCGGCCGAGACCCCGTCTGGGAGGTGAGGAGCGTCTCTGCCCGGCCGCCCCGTCTGAGAAGTGAGGAGACCCTCTGCCTGGCAACCACCCCGTCTGAGAAGTGAGGAGCCCCTCTGCCCGGCCAGCACCCCGTCCGGGAGGGAGGTGGGGGGGTCAGCCCCCCGCCCGGCCAGCCGCCCCGTCCGGGAGGGAGGTGGGGGGGGTCAGCCCCCCCGCCCGGCCAGCCGCCCCGTCCGGGAGGTGAGGGGCGCCTCTGCCCGGCCGCCCCTACTGGGAAGTGAGGAGCCCCTCTGCCCGGCCACCACCCCGTCTGGGACGTGTGCCCAACAGCTCATTGAGAACGGGCCAGGATGACAATGGCGGCTTTGTGGAATAGAAAGGCGGGAAAGGTGGGGAAAAGATTGAGAAATCGGATGGTTGCCGTGTCTGTGTAGAAAGAAGTAGACATGGGAGACTTTTCATTTTGTTCTGCACTAAGAAAAATTCCTCTGCCTTGGGATCCTGTTGATCTGTGACCTTACCCCCAACCCTGTGCTCTCTGAAACATGTGCTGTGTCCACTCAGGGTTAAATGGATTAAGGGCGGTGCAAGATGTGCTTTGTTAAACAGATGCTTGAAGGCAGCATGCTCGTTAAGAGTCATCACCAATCCCTAATCTCAAGTAATCAGGGACACAAACACTGCGGAAGGCCGCAGGGTCCTCTGCCTAGGAAAACCAGAGACCTTTGTTCACTTGTTTATCTGCTGACCTTCCCTCCACTATTGTCCCATGACCCTGCCAAATCCCCCTCTGTGAGAAACACCCAAGAATTATCAATAAAAAAATAAATTAAAAAAAAAA"
context = 10000
paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
models = [load_model(resource_filename('spliceai', x)) for x in paths]
x = one_hot_encode('N'*(context//2) + input_sequence + 'N'*(context//2))[None, :]
y = np.mean([models[m].predict(x) for m in range(5)], axis=0)
acceptor_prob = y[0, :, 1]
donor_prob = y[0, :, 2]
# Plotting the data
fig, ax1 = plt.subplots(figsize=(20, 15))
# Plotting the first dataset
ax1.plot(acceptor_prob, color='black', label='acceptor_prob')

# Plotting the second dataset
ax1.plot(donor_prob, color='blue', label='donor_prob')

# Adding labels and title
ax1.set_title('Splicing probability',fontsize=20)
ax1.set_xlabel('Base Position',fontsize=16)
ax1.set_ylabel('probability', color='black',fontsize=16)

# Highlighting and annotating specific points in the first dataset
highlight_points1 = [744, 667]  # example indices to highlight in the first dataset
annotations1 = ["BambuTx619ss", "BambuTx620ss"]
for i, point in enumerate(highlight_points1):
    ax1.annotate(annotations1[i], 
                 xy=(point, acceptor_prob[point]), 
                 xytext=(point, acceptor_prob[point] + 320), 
                 textcoords='offset points', 
                 arrowprops=dict(arrowstyle="->", color='red' if i == 0 else 'green'),
                 color='red' if i == 0 else 'green')

highlight_points1 = [107, 114]  # example indices to highlight in the first dataset
annotations1 = ["BambuTx619TSS", "BambuTx620TSS"]
for i, point in enumerate(highlight_points1):
    ax1.annotate(annotations1[i], 
                 xy=(point, donor_prob[point]), 
                 xytext=(point, donor_prob[point] + 320), 
                 textcoords='offset points', 
                 arrowprops=dict(arrowstyle="->", color='red' if i == 0 else 'green'),
                 color='red' if i == 0 else 'green')

# Annotating x-axis regions
ax1.axvspan(16, 365, color='yellow', alpha=0.3, label='Alu-like')
ax1.axvspan(366, 966, color='lightgreen', alpha=0.3, label='VNTR')
ax1.axvspan(967, 1475, color='red', alpha=0.3, label='SINE-R')
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
