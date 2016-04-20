
import simuPOP as sim
from simuPOP import *


def read_recomb_rates(filename, num_loci):
    print 'reading'
    #return list of recomb rate at each position
    per_site_recombs = []
    with open(filename) as f:
        f.readline()
        _, prev_kb_pos, _ = f.readline().strip().split(',')
        for line in f:
            #values in file are in 4*Ne*r/kb, where r = (physical distance (l) * per site rate of recomb (p))
            #4*Ne*r = x, 4*Ne*l*p = x, p = x/(4*Ne*l)
            curr_chr, curr_kb_pos, curr_rate = line.strip().split(',') 
            region_kb_length = float(curr_kb_pos) - float(prev_kb_pos)
            p = float(curr_rate) / (4. * (float(curr_kb_pos) - float(prev_kb_pos)))

            #add rate for every position in region
            per_site_recombs.extend([p for _ in range(int(1000. * region_kb_length))])

            prev_kb_pos = curr_kb_pos
            if curr_chr != 'chr1':
                break

    return per_site_recombs[:num_loci]
             

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
    num_gens = 20
    recomb_file_or_val = 'data/mouse_recomb_rates.csv'
    #recomb_file_or_val = '.00001'

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
    try:
        recomb_rates = read_recomb_rates(recomb_file_or_val, num_loci)
        recomb_str = 'mouse'
    except IOError:
        recomb_rates = float(recomb_file_or_val)
        recomb_str = str(recomb_rates)

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
                        sim.Recombinator(rates=recomb_rates, output='>>data/recombs_r%s_g%i.log' % (recomb_str, num_gens), 
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


