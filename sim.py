
import simuPOP as sim
from simuPOP import *
from shutil import copyfile

def read_recomb_rates(filename, num_loci, chrom=None):
    #return list of recomb rate at each position
    read_chrom = 'chr1' if not chrom else chrom

    per_site_recombs = []
    with open(filename) as f:
        f.readline()
        curr_chr = ''

        #read until up to chromosome we care about
        while curr_chr != read_chrom:
            curr_chr, prev_kb_pos, _ = f.readline().strip().split(',')

        for line in f:
            curr_chr, curr_kb_pos, curr_rate = line.strip().split(',') 
            if curr_chr != read_chrom:
                break

            #values in file are in 4*Ne*r/kb, where r = (physical distance (l) * per site rate of recomb (p))
            #4*Ne*r = x, 4*Ne*l*p = x, p = x/(4*Ne*l)
            region_kb_length = float(curr_kb_pos) - float(prev_kb_pos)
            p = (float(curr_rate) / 100.) / (4. * region_kb_length * 1000.)

            #add rate for every position in region
            per_site_recombs.extend([p for _ in range(int(1000. * region_kb_length))])

            prev_kb_pos = curr_kb_pos

    #fill in the last positions
    if len(per_site_recombs) < num_loci:
        per_site_recombs.extend([per_site_recombs[-1] for _ in range(num_loci - len(per_site_recombs))])

    assert(len(per_site_recombs) > 0)
    return per_site_recombs[:num_loci]
             

def nuc_int(nuc):
    if nuc == 'A':
        return 0
    if nuc == 'T':
        return 1
    if nuc == 'C':
        return 2
    if nuc == 'G':
        return 3
    return 4


def compare_seqs(seqA, seqB):
    assert(len(seqA) == len(seqB))
    matches = 0
    for a, b in zip(seqA, seqB):
        if a == b:
            matches += 1
    return float(matches)/len(seqA)


if __name__ == '__main__':
    # intialize
    debug = True

    #wkdir = '/Users/dama9282/simulation/'
    wkdir = './'

    recomb_file_or_val = wkdir + 'data/mouse_recomb_rates.csv'
    if debug:
        recomb_file_or_val = '.01'
    pop_size = 96
    num_rand_gens = 20
    num_sib_gens = 30
    chrom = 19
    num_loci = -1
    if debug:
        num_loci = 100

    strain_seqs = {}
    strain_names = ('AKRJ', 'AJ', 'BALBcJ', 'C3HHeJ', 'CASTEiJ', 'CBAJ', 'DBA2J', 'LPJ')
    #strain_names = ('AKRJ', 'AJ', 'BALBcJ', 'C3HHeJ')

    # read fasta files
    print 'reading fastas'
    for name in strain_names:
        with open(wkdir + './data/mouse_fastas/%s_chr%i.fa' % (name, chrom)) as f:
            f.readline()
            strain_seqs[name] = f.readline().strip()
            
    # read recomb rate file
    print 'reading recomb file'
    num_loci = len(strain_seqs[strain_names[0]]) if num_loci < 1 else num_loci
    try:
        recomb_rates = read_recomb_rates(recomb_file_or_val, num_loci, 'chr%i' % chrom)
        recomb_str = 'mouse'
    except IOError:
        recomb_rates = float(recomb_file_or_val)
        recomb_str = str(recomb_rates)

    # create population
    print 'creating population' 
    pop = Population(size=pop_size, ploidy=2, loci=num_loci, alleleNames=['A','T','C','G','N'], 
                    infoFields=['ind_id'])
    for i in range(0, pop.popSize(), 2):
        pop.individual(i).setSex(MALE)
        if i != pop.popSize()-1:
            pop.individual(i+1).setSex(FEMALE)

    # set numeric genotypes
    for i in range(pop.popSize()):
        name = strain_names[i % len(strain_names)]
        if debug:
            tmp = 3500000 + i*100
            pop.individual(i).setGenotype([nuc_int(n) for n in strain_seqs[name][tmp:tmp+num_loci]])
        else:
            pop.individual(i).setGenotype([nuc_int(n) for n in strain_seqs[name][:num_loci]])

    # random mating
    print 'random mating'
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
                        sim.Recombinator(rates=recomb_rates, output='>>' + wkdir + 'data/recombs_rand-mating_r%s_g%i.log' % (recomb_str, num_rand_gens), 
                                         infoFields=['ind_id']),
                    ]
            ),
            gen=num_rand_gens
        )

    print 'output random mating descendant sequences'
    #grab last individual from population
    indv = pop.individual(pop.popSize()-1)

    with open(wkdir + './data/sim-desc-seq_rand-mating_0.nuc', 'w') as f:
        f.write(''.join(map(str, indv.genotype(ploidy=0))).replace('0','A').replace('1','T').replace('2','C').replace('3','G').replace('4','N'))
        f.write('\n')

    with open(wkdir + './data/sim-desc-seq_rand-mating_1.nuc', 'w') as f:
        f.write(''.join(map(str, indv.genotype(ploidy=1))).replace('0','A').replace('1','T').replace('2','C').replace('3','G').replace('4','N'))
        f.write('\n')

    print 'sibling inbreeding'

    #create a new population with two copies of that individual (can be considered siblings)
    del pop
    pop = Population(size=2, ploidy=2, loci=num_loci, alleleNames=['A','T','C','G','N'], 
                    infoFields=['ind_id'])
    print 'pop created'
    #brother
    pop.individual(0).setSex(MALE)
    pop.individual(0).setGenotype(indv.genotype())
    #sister
    pop.individual(1).setSex(FEMALE)
    pop.individual(1).setGenotype(indv.genotype())
    print 'genotypes set'

    pop.evolve(
            initOps=[
                sim.InitSex(sex=(sim.MALE, sim.FEMALE)),
                sim.IdTagger(),
                ],
            matingScheme=sim.MonogamousMating(
                numOffspring=2,
                sexMode=(sim.NUM_OF_MALES, 1),
                ops=[
                        sim.IdTagger(),
                        sim.Recombinator(rates=recomb_rates, output='>>' + wkdir + 'data/recombs_sib-mating_r%s_g%i.log' % (recomb_str, num_sib_gens), 
                                         infoFields=['ind_id']),

                    ]
            ),
            gen=num_sib_gens
        )

    print 'output sibling mating descendant sequences'
    #grab last individual from population
    indv = pop.individual(pop.popSize()-1)
    with open(wkdir + './data/sim-desc-seq_sib-mating_0.nuc', 'w') as f:
        f.write(''.join(map(str, indv.genotype(ploidy=0))).replace('0','A').replace('1','T').replace('2','C').replace('3','G').replace('4','N'))
        f.write('\n')
    with open(wkdir + './data/sim-desc-seq_sib-mating_1.nuc', 'w') as f:
        f.write(''.join(map(str, indv.genotype(ploidy=1))).replace('0','A').replace('1','T').replace('2','C').replace('3','G').replace('4','N'))
        f.write('\n')

    print 'combine recomb files'
    #combine random mating and sibling mating recombination files with correct ancestor numbers
    copyfile(wkdir + './data/recombs_rand-mating_r%s_g%i.log' % (recomb_str, num_rand_gens), 
             wkdir + './data/recombs_full_r%s_g%i.log' % (recomb_str, num_rand_gens))

    with open(wkdir + './data/recombs_full_r%s_g%s.log' % (recomb_str, num_rand_gens), 'a') as f_full:
        #add the first siblings as clones of the last descendant (inherit both chromosomes with no crossovers)
        last_desc = pop_size * (num_rand_gens + 1)
        f_full.write('%i %i 0\n' % (last_desc+1, last_desc))
        f_full.write('%i %i 1\n' % (last_desc+1, last_desc))
        f_full.write('%i %i 0\n' % (last_desc+2, last_desc))
        f_full.write('%i %i 1\n' % (last_desc+2, last_desc))
        #very conviently, simuPOP increments recomb file generations across multiple populations/evolutions, so all we 
        # have to do is append the sibling mating recomb file to the full file
        with open(wkdir + './data/recombs_sib-mating_r%s_g%i.log' % (recomb_str, num_sib_gens)) as f_sib:
            for line in f_sib:
                f_full.write(line)




