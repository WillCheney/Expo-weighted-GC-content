

#Will Cheney Sept 2019

#import Libraries
import matplotlib.pyplot as plt
import numpy as np



# Sequence of transcript to analyze
example = 'ATGGACAGCAGCGCTGCCCCCACGAACGCCAGCAATTGCACTGATGCCTTGGCGTACTCAAGTTGCTCCCCAGCACCCAGCCCCGGTTCCTGGGTCAACTTGTCCCACTTAGATGGCAACCTGTCCGACCCATGCGGTCCGAACCGCACCGACCTGGGCGGGAGAGACAGCCTGTGCCCTCCGACCGGCAGTCCCTCCATGATCACGGCCATCACGATCATGGCCCTCTACTCCATCGTGTGCGTGGTGGGGCTCTTCGGAAACTTCCTGGTCATGTATGTGATTGTCAGATACACCAAGATGAAGACTGCCACCAACATCTACATTTTCAACCTTGCTCTGGCAGATGCCTTAGCCACCAGTACCCTGCCCTTCCAGAGTGTGAATTACCTAATGGGAACATGGCCATTTGGAACCATCCTTTGCAAGATAGTGATCTCCATAGATTACTATAACATGTTCACCAGCATATTCACCCTCTGCACCATGAGTGTTGATCGATACATTGCAGTCTGCCACCCTGTCAAGGCCTTAGATTTCCGTACTCCCCGAAATGCCAAAATTATCAATGTCTGCAACTGGATCCTCTCTTCAGCCATTGGTCTTCCTGTAATGTTCATGGCTACAACAAAATACAGGCAAGGTTCCATAGATTGTACACTAACATTCTCTCATCCAACCTGGTACTGGGAAAACCTGCTGAAGATCTGTGTTTTCATCTTCGCCTTCATTATGCCAGTGCTCATCATTACCGTGTGCTATGGACTGATGATCTTGCGCCTCAAGAGTGTCCGCATGCTCTCTGGCTCCAAAGAAAAGGACAGGAATCTTCGAAGGATCACCAGGATGGTGCTGGTGGTGGTGGCTGTGTTCATCGTCTGCTGGACTCCCATTCACATTTACGTCATCATTAAAGCCTTGGTTACAATCCCAGAAACTACGTTCCAGACTGTTTCTTGGCACTTCTGCATTGCTCTAGGTTACACAAACAGCTGCCTCAACCCAGTCCTTTATGCATTTCTGGATGAAAACTTCAAACGATGCTTCAGAGAGTTCTGTATCCCAACCTCTTCCAACATTGAGCAACAAAACTCCACTCGAATTCGTCAGAACACTAGAGACCACCCCTCCACGGCCAATACAGTGGATAGAACTAATCATCAGCTAGAAAATCTGGAAGCAGAAACTGCTCCGTTGCCCTAA'
sequence = ''

# Returns GC-content percentage of a sequence
def get_gc(sequence):
    sequence = sequence.upper()
    gc = 0
    for i in sequence:
        if i in ['G','C']:
            gc += 1

    return (gc/len(sequence))*100


#Returns lists x and y which contain basepair index and weighted average GC-content at corresponding x index
#Uses equation y[i] = theta * y[i - 1] + (1 - theta) * {1 if C/G else 0 if A/T}
#Theta value approximates the number of bp to average over, ~= 1/(1 - theta) 
def gc_weighted_average(sequence, theta = 0.96):
    sequence = sequence.lower()
    sequence = sequence.replace('u','t')
    x = []
    y = []
    for i in range(len(sequence)):
        if len(y)>0:
            gc = (1-theta)*get_gc(sequence[i]) + theta*y[-1]
            y.append(gc)
            x.append(i)
        else:
            y.append(get_gc(sequence[i]))
            x.append(i)
    return x,y
    




#Generate GC content Graph using matplotlib
x,y = gc_weighted_average(example)
plt.plot(x,y, c = 'royalblue', linewidth = 1)
plt.ylabel('GC Content(%)')
plt.xlabel('Base Pair Index')


     
plt.ylim(0,100)


fig = plt.gcf()
plt.show()



# Save Figure 
#fig.savefig('2020-06-15 ste2',dpi = 300, bbox_inches='tight')








