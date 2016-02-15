
from random import choice

import simuPOP as sim
from simuPOP import *

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


if __name__ == '__main__':
    strain_names = {'834774811': 'NCIM3186',
                    '330443391': 'S288c',
                    '455768006': 'W303',
                    '759090333': 'YJM450',
                    '759134535': 'YJM451',
                    '759334484': 'YJM555',
                    '759349870': 'YJM981',
                    '874346690': 'ySR127'}

    #alignment_file = 'data/kalign-yeast.clustalw'
    #strain_ids = ('834774811', '330443391', '455768006', '759090333', '759134535', '759334484', '759349870', '874346690')
    #used_strains = ('759090333', '759134535', '759334484', '759349870', '874346690')
    strain_ids = ('gi|329138864|tp', 'gi|874346693|gb', 'gi|768752667|gb', 'gi|768739925|gb', 'gi|768744865|gb', 'gi|768740643|gb')
    used_strains = ('gi|874346693|gb', 'gi|768752667|gb', 'gi|768739925|gb', 'gi|768744865|gb', 'gi|768740643|gb')
    alignment_file = 'data/aligned.aln'

    recomb_rate = .00001
    num_gens = 20

    # read alignment file
    tot_alignments = len(strain_ids)
    strain_seqs = {s_id: '' for s_id in strain_ids}

    with open(alignment_file, 'r') as f:
        f.readline()
        f.readline()
        f.readline()
        line = f.readline().strip()
        while line != '':
            for _ in range(tot_alignments):
                cols = line.split()
                strain_seqs[cols[0]] += process_seq(cols[1])
                line = f.readline().strip()
            f.readline()
            line = f.readline().strip()
            
    # create population
    num_loci = len(strain_seqs[strain_ids[0]])
    #num_loci = 10000
    pop = Population(size=len(used_strains), ploidy=2, loci=num_loci, alleleNames=['A','T','C','G','N'], 
                    #infoFields=['father_idx', 'mother_idx', 'child_idx'])
                    infoFields=['ind_id'])
    for i in range(0, pop.popSize(), 2):
        pop.individual(i).setSex(MALE)
        if i != pop.popSize()-1:
            pop.individual(i+1).setSex(FEMALE)

    strain_alleles = {}
    for k, v in strain_seqs.items():
        strain_alleles[k] = [int(a) for a in v]

    for i, s_id in enumerate(used_strains):
        pop.individual(i).setGenotype(strain_alleles[s_id][:num_loci])
    
    #for i, indv in enumerate(pop.individuals()):
    #    print str(i) + ': ' + str(indv.genotype())
    #    print str(i) + ': ' + str(len(indv.genotype()))
    #    print indv

    #print '\n\npre evolve:'
    #dump(pop)

    pop.evolve(
            initOps=[
                sim.InitSex(sex=(sim.MALE, sim.FEMALE)),
                sim.IdTagger(),
                ],
            matingScheme=sim.RandomMating(
                numOffspring=2,
                sexMode=(sim.NUM_OF_MALES, 1),
                ops=[
                        sim.IdTagger(),
                        sim.Recombinator(rates=recomb_rate, output='>>data/recombs_r%s_g%i.log' % (str(recomb_rate), num_gens), 
                                         infoFields=['ind_id']),

                    ]
            ),
            gen=num_gens
        )

    #print pop.subPopSize()

    #print '\n\npost evolve:'
    #dump(pop)

    #for i, indv in enumerate(pop.individuals()):
    #    print str(i) + ': ' + str(indv.genotype(chroms=0))
    #    print str(i) + ': ' + str(len(indv.genotype(chroms=0)))
    #    print indv


