
from itertools import izip

chrom = 19

#convert output from rabbit and happy to a bed file for scoring

# RABBIT
#  read output
ancestors = [] #ancestors in order of their output number
with open('./rabbit/rabbit_viterbi_Summary.csv') as f:
    line = f.readline()
    #find ancestor output order
    while line[:9] != 'Haplotype':
        line = f.readline()
    line = f.readline()
    while line[:9] == 'haplotype':
        ancestors.append(line.strip().split()[-1])
        line = f.readline()

    #read desc output
    while line[:4] != 'DESC':
        line = f.readline()
    assert(len(ancestors) > 0)
    desc_path = line.strip().split()[1]

#  write bed file
with open('./rabbit/rabbit_haplotypes.bed', 'w') as f:
    #path string is in the form pos-diplocode-pos-diplocode-...-diplocode-pos
    # diplocodes increment every pair of ancestors, eg. for 8 ancestors:
    #  1 - anc1:anc1
    #  2 - anc1:anc2
    #   ...
    #  8 - anc1:anc8
    #  9 - anc2:anc1
    #  10 - anc2:anc2
    #   ...
    # essentially it follows the formula: ((#-1)/8):((#-1)%8) for the indexes of anc1:anc2
    entries = desc_path.split('-')
    print len(entries)
    for i in range(0, len(entries)-1, 2):
        print ''
        print i, entries[i]
        diplocode = int(entries[i+1])
        print diplocode
        diplo_ancs = set((ancestors[(diplocode-1) / 8], ancestors[(diplocode-1) % 8]))
        print diplo_ancs
        f.write('chr%i\t%s\t%s\t%s\n' % (chrom, entries[i], entries[i+2], '_'.join(sorted(diplo_ancs))))





