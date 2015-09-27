
from random import choice

import simuPOP as sim
from simuPOP import *


strain_ids = {'834774811': 'NCIM3186',
              '330443391': 'S288c',
              '455768006': 'W303',
              '759090333': 'YJM450',
              '759134535': 'YJM451',
              '759334484': 'YJM555',
              '759349870': 'YJM981',
              '874346690': 'ySR127'}

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

    SNPs = ['1' if a != ref_allele and a != '4' else '0' for a in alleles]

    return ref_allele, SNPs


def process_seq(seq):
    seq = seq.replace('A','0')
    seq = seq.replace('T','1')
    seq = seq.replace('C','2')
    seq = seq.replace('G','3')
    # not sure what to do about these...
    seq = seq.replace('-','4')
    seq = seq.replace('N','4')
    seq = seq.replace('W','4')
    seq = seq.replace('M','4')
    seq = seq.replace('R','4')
    seq = seq.replace('Y','4')
    seq = seq.replace('K','4')
    seq = seq.replace('S','4')
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
    


# create population
num_loci = 100
pop = Population(size=2, ploidy=2, loci=num_loci, alleleNames=['A','T','C','G','N'], infoFields=['father_idx', 'mother_idx'])
for i in range(0, pop.popSize(), 2):
    pop.individual(i).setSex(MALE)
    pop.individual(i+1).setSex(FEMALE)

strain_alleles = {}
for k, v in strain_seqs.items():
    strain_alleles[k] = [int(a) for a in v]

print strain_alleles[0]

pop.individual(0).setGenotype(strain_alleles[0][:num_loci])
pop.individual(1).setGenotype(strain_alleles[2][:num_loci])

print '\n\npre evolve:'
dump(pop)

pop.evolve(
        initOps=sim.InitSex(sex=(sim.MALE, sim.FEMALE)),
        matingScheme=sim.RandomMating(
            numOffspring=2,
            sexMode=(sim.NUM_OF_MALES, 1),
            ops=[
                    #sim.Recombinator(rates=.1, convMode=(sim.NUM_MARKERS, 1, 10)),
                    sim.Recombinator(rates=.01),
                    #sim.ParentsTagger(),
                ]
        ),
        gen=2
    )

'''
pop.evolve(
        initOps=sim.InitSex(sex=(sim.MALE, sim.FEMALE)),
        matingScheme=sim.HomoMating(
        #matingScheme=sim.RandomMating(
            #sim.SequentialParentChooser(),
            sim.RandomParentChooser(),
            sim.OffspringGenerator(
                ops=[
                        sim.Recombinator(rates=1, convMode=(sim.NUM_MARKERS, 1, 10)),
                        #sim.ParentsTagger(),
                    ],
                numOffspring=2,
                sexMode=(sim.NUM_OF_MALES, 1),
            )
        ),
        gen=1
    )
'''

#pop.evolve(gen=5)


print pop.subPopSize()

print '\n\npost evolve:'
dump(pop)

for indv in pop.individuals():
    print indv.sex()


