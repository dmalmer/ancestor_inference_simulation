
from random import choice

import simuPOP as sim
from simuPOP import *


def call_SNPs(alleles):
    max_a = []
    max_ct = -1
    for a in range(4):
        ct = alleles.count(str(a))
        if ct > max_ct:
            max_a = [str(a)]
            max_ct = ct
        elif ct == max_ct:
            max_a.append(str(a))
    ref_allele = choice(max_a)

    SNPs = ['1' if a != ref_allele and a != '5' else '0' for a in alleles]

    return ref_allele, SNPs


def process_seq(seq):
    seq = seq.replace('A','0')
    seq = seq.replace('T','1')
    seq = seq.replace('C','2')
    seq = seq.replace('G','3')
    # not sure what to do about these...
    seq = seq.replace('-','5')
    seq = seq.replace('N','5')
    seq = seq.replace('W','5')
    seq = seq.replace('M','5')
    seq = seq.replace('R','5')
    seq = seq.replace('Y','5')
    seq = seq.replace('K','5')
    seq = seq.replace('S','5')
    return seq


strains = ('W303', 'YJM555', 'NCIM3186', 'YJM450', 'YJM981', 'S288c', 'YJM451', 'ySR127')

# read alignment file
tot_alignments = 8
strain_seqs = {i: '' for i in range(tot_alignments)}

with open('kalign-yeast.clustalw', 'r') as f:
    f.readline()
    f.readline()
    f.readline()
    line = f.readline().strip()
    j = 0
    while line != '':
        for i in range(tot_alignments):
            strain_seqs[i] += process_seq(line.split()[1])
            line = f.readline().strip()
        f.readline()
        line = f.readline().strip()
        j += 1
        if j > 3:
            break
         
# call SNPs on alignment file and construct reference sequence
ref_seq = ''
SNPs_by_strain = {i: '' for i in range(tot_alignments)}
for i in range(len(strain_seqs[0])):
    alleles = [s[i] for s in [strain_seqs[j] for j in range(tot_alignments)]]
    ref_allele, SNPs = call_SNPs(alleles)

    ref_seq += ref_allele
    for j, s in enumerate(SNPs):
        SNPs_by_strain[j] += s

print ''
print ref_seq
for SNPs in SNPs_by_strain.values():
    print SNPs
    


exit(0)

seq_alleles = {}
for strain in strains:
    seq = ''
    with open('S_cerevisiae_chr1/' + strain + '.chr1.fasta', 'r') as f:
        f.readline()
        for line in f:
            seq += line.strip()
    seq = seq.replace('A','0')
    seq = seq.replace('T','1')
    seq = seq.replace('C','2')
    seq = seq.replace('G','3')
    # not sure what to do about these...
    seq = seq.replace('N','5')
    seq = seq.replace('W','5')
    seq = seq.replace('M','5')
    seq = seq.replace('R','5')
    seq = seq.replace('Y','5')
    seq = seq.replace('K','5')
    seq = seq.replace('S','5')
    seq_alleles[strain] = [int(a) for a in seq]

# need to align coordinates to reference (so most positions are homozygous)

# create population
pop = Population(size=2, ploidy=2, loci=100, alleleNames=['A','T','C','G','N'], infoFields=['father_idx', 'mother_idx'])
#for i in range(0, pop.popSize(), 2):
#    pop.individual(i).setSex(MALE)
#    pop.individual(i+1).setSex(FEMALE)

pop.individual(0).setGenotype(seq_alleles['S288c'][:100])
pop.individual(1).setGenotype(seq_alleles['YJM451'][:100])

print '\n\npre evolve:'
dump(pop)

#pop.evolve(initOps=[InitSex(), InitGenotype(freq=[.25,.25,.25,.25])], gen=2)
#pop.evolve(initOps=[InitGenotype(freq=[.25,.25,.25,.25])], gen=2)
pop.evolve(
        initOps=sim.InitSex(sex=(sim.MALE, sim.FEMALE)),
        matingScheme=sim.MonogamousMating(
            numOffspring=2,
            sexMode=(sim.NUM_OF_MALES, 1),
            ops=[
                    sim.MendelianGenoTransmitter(),
                    sim.ParentsTagger(),
                ],
        ),
        gen=1
    )



#pop.evolve(gen=5)


print pop.subPopSize()

print '\n\npost evolve:'
dump(pop)

for indv in pop.individuals():
    print indv.sex()


