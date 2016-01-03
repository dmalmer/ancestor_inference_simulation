
from collections import defaultdict


strain_ids = {'330443391': 'S288c',
              '455768006': 'W303',
              '834774811': 'NCIM3186',
              '759090333': 'YJM450',
              '759134535': 'YJM451',
              '759334484': 'YJM555',
              '759349870': 'YJM981',
              '874346690': 'ySR127'}

ref = '330443391'
group_size = 8

bed_SNPs = [] #list of lists: [[chr, pos, pos+1, ancestors],...]

counts = defaultdict(int)

with open('data/kalign-yeast.clustalw', 'r') as f:
#with open('data/tmp_head.txt', 'r') as f:
    f.readline()
    f.readline()
    f.readline()

    p = 0
    line = f.readline().strip()
    while line != '':
        group_lines = {} 
        for _ in range(group_size):
            strain, seq = line.split()
            line = f.readline().strip()

            group_lines[strain] = seq

        curr_SNPs = defaultdict(list)
        for strain, seq in group_lines.items():
            if strain == ref:
                continue
            for i, (nuc_s, nuc_r) in enumerate(zip(group_lines[strain], group_lines[ref])):
                if nuc_s != nuc_r:
                    if nuc_r == '-': #insertion
                        pass
                    elif nuc_s == '-': #deletion
                        pass
                    else: #SNP
                        if strain not in ('455768006', '834774811'):
                            curr_SNPs[i].append(strain)
                            counts[strain] += 1

        for i in sorted(curr_SNPs.keys()):
            print 'CURR SNPS'
            print curr_SNPs[i]
            SNP_str = ''
            for s in sorted(curr_SNPs[i]):
                SNP_str += strain_ids[s] + '_'
            SNP_str = SNP_str[:-1]
            bed_SNPs.append(['chr1', p+i, p+i+1, SNP_str])
        
        p += len(group_lines[ref]) 
        print p

        print group_lines
        print curr_SNPs
        print '\n\n'
        f.readline()
        line = f.readline().strip()


with open('data/ancestor_SNPs.bed', 'w') as f:
    for curr_chr, pos_s, pos_e, ancs in bed_SNPs:
        f.write('%s\t%i\t%i\t%s\n' % (curr_chr, pos_s, pos_e, ancs))

print counts
