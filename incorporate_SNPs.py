
import sys
from collections import defaultdict

wkdir = '/scratch/Users/dama9282/sim_SNPs/'

strain_names = ('AKR/J', 'A/J', 'BALB/cJ', 'C3H/HeJ', 'CAST/EiJ', 'CBA/J', 'DBA/2J', 'LP/J')

#read chrom parameter
chrom = sys.argv[1] if len(sys.argv) > 1 else '19'
print 'chromosome ' + chrom

#genotypes organized as dictionary of lists for each chrom
# chrom lists are a sequence of tuples with (position, [AKR-nuc, A-nuc, BALBc-nuc, ..., LP-nuc])
print 'reading genotypes'
genotypes = []
with open(wkdir + 'SC_MOUSE_GENOMES.genotype.vcf') as f_geno:
    #skip header
    line = f_geno.readline()
    while line[:2] == '##':
        line = f_geno.readline()

    #get strain indexes
    cols = line.strip().split()
    strain_inds = [cols.index(n) for n in strain_names]

    #read genotypes
    for j, line in enumerate(f_geno):
        cols = line.strip().split()
        if cols[0] != chrom:
            continue
        
        #first get genotypes (0/0, 1/1, etc.)
        genos = [cols[i].replace('.','0') for i in strain_inds]
        
        #now determine nucleotides of each strain
        alleles = [cols[3]] + cols[4].split(',')
        nucs = [alleles[int(g[0])] for g in genos]
        
        #add position list to chrom list
        genotypes.append((int(cols[1]), nucs))

        if j % 100000 == 0:
            print j

#create a new fasta for each strain for each chromosome
print 'writing fastas'

#read in reference fasta
ref_fasta = ''
with open(wkdir + 'chr%s.fa' % chrom) as f_in:
    f_in.readline()
    for line in f_in:
        ref_fasta += line.strip()

#write strain-specific fasta by replacing reference nucleotides with strain-specific genotypes
for i, name in enumerate([n.replace('/','') for n in strain_names]):
    print name
    strain_genome = list(ref_fasta)
    for j, (pos, nucs) in enumerate(genotypes):
        strain_genome[pos] = nucs[i]
        if j % 10000 == 0:
            print 'j = ' + str(j)

    with open(wkdir + '%s_chr%s.fa' % (name.replace('/',''), chrom), 'w') as f_out:
        f_out.write('>' + name + '\n')
        f_out.write(''.join(strain_genome) + '\n')



