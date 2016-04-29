
from collections import defaultdict

#convert desc_segments into a list the size of the desc genome with the ancestor at each position
desc_ancs = []
print 'reading descendant segments'
with open('./data/desc_segments_chr1.bed', 'r') as f:
    for line in f:
        _, anc, start, end = line.strip().split()
        desc_ancs.extend([anc] * (int(end) - int(start)))

#start at first SNP and end at last SNP, no informative markers past either end
print 'comparing hmm output'
for filname in ('./data/sim_SNPs_hmm-out.bed',
                './rabbit/rabbit_haplotypes.bed',
                './happy/happy_haplotypes.bed'):
    print filname
    correct_SNPs = 0
    incorrect_SNPs = 0
    with open(filname, 'r') as f:
        for line in f:
            _, start, end, ancs_str = line.strip().split()[:4]
            ancs = ancs_str.split('_')
            for i in range(int(start), int(end)):
                if desc_ancs[i] in ancs:
                    correct_SNPs += 1
                else:
                    incorrect_SNPs += 1

    print 'correct_SNPs =', correct_SNPs
    print 'incorrect_SNPs =', incorrect_SNPs
    print ''

