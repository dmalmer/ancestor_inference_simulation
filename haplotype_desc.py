
'''
each indv is a dictionary of the left and right chromosomes
each chrom is a pair of lists, the first showing the starting position of a new ancestor, the other showing which ancestor is at that position

ancestors 1 and 2 (homozygous at all positions):
1 = {'0': [[0],
          [1]]
     '1': [[0],
          [1]]}
2 = {'0': [[0],
          [2]]
     '1': [[0],
          [2]]}

first generation (special case, each chrom is the exact sequence from an ancestor, so ignore recombs):
3 = {'0': [[0],
          [1]]
     '1': [[0],
          [2]]}
4 = {'0': [[0],
          [2]]
     '1': [[0],
          [1]]}

next generations will have recombination events:
5 = {'0': [[0, 20, 80],
          [1, 2, 1]]
     '1': [[0, 60],
          [2, 1]]}
6 = {'0': [[0, 50],
          [2, 1]]
     '1': [[0, 30],
          [1, 2]]}
'''

from itertools import tee, izip


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)


if __name__ == '__main__':

    #RECOMB_FILE = 'data/recombs_r1e-05_g20.log'
    RECOMB_FILE = 'data/recombs_rand-mating_rmouse_g20.log'
    #RECOMB_FILE = 'data/recombs_r0.01_g20.log'
    CHROM = 19

    # re-trace what parts of each genome came from which ancestor using the recombination file output by simuPOP
    print 'reading recomb file'
    genomes = [] #index+1 = indv ID
    genomes = open(RECOMB_FILE,'r').readlines()
    #for g in genomes:
    #    print g
    #print '\n\n\n'

    anc_pop = int(genomes[0].split()[0]) - 1

    def parent_index(parent, ploidy, anc_pop):
        #return index of genomes[] list that corresponds to parent number in recombs.log file
        return (parent - anc_pop - 1) * 2 + ploidy


    #only care about last genome, start from there and work backwords
    #only looking at ploidy 1 for now
    #haplotypes = [(pos1, parent1, ploidy1), (pos2, parent2, ploidy2), ..., (posN, parentN, ploidyN)]]
    # last line
    init_line = genomes[-1].strip().split()
    parent = int(init_line[1])
    ploidy = int(init_line[2])
    curr_haplotypes = [(0, parent, ploidy)]
    for recomb in init_line[3:]:
        curr_haplotypes.append((int(recomb), parent, abs(ploidy-1)))
    
    print 'starting haplotypes:'
    print curr_haplotypes

    # replace segments with their parent segments until we're only left with original ancestors
    while not all([hap[1] <= anc_pop for hap in curr_haplotypes]):
        print '\n'
        print 'truth check'
        print [hap[1] <= anc_pop for hap in curr_haplotypes]
        print curr_haplotypes
        new_haplotypes = list(curr_haplotypes)

        #keep track of i offset as we insert new parental haplotypes
        offset = 0 
        for i, (curr_pos, curr_parent, curr_ploidy) in enumerate(curr_haplotypes):
            print ''
            print 'hap', i
            #get ending position of current haplotype we are replacing with parent haplotypes
            hap_end_pos = curr_haplotypes[i+1][0] if i+1 < len(curr_haplotypes) else -1
            #get parent line of current haplotype
            print 'parent index', parent_index(curr_parent, curr_ploidy, anc_pop)
            parent_line = genomes[parent_index(curr_parent, curr_ploidy, anc_pop)].strip().split()
            print parent_line
            
            #replace the current haplotype with its parent haplotypes
            new_parent = int(parent_line[1])
            new_ploidy = int(parent_line[2])
            new_recombs = [int(r) for r in parent_line[3:]]
            # first, loop through parent recombinations until we get to the current position
           # this gets us to the correct ploidy of the parent at the current position
            for recomb in [r for r in new_recombs if r < curr_pos]:
                new_ploidy = abs(new_ploidy-1)
            # now, replace the current haplotype with the first parent haplotype
            new_haplotypes[i+offset] = (curr_pos, new_parent, new_ploidy)
            # next, loop through the rest of the parent recombinations until we reach the end of the current haplotype
            for recomb in [r for r in new_recombs if r > curr_pos and (hap_end_pos == -1 or r < hap_end_pos)]:
                offset += 1
                new_ploidy = abs(new_ploidy-1)
                new_haplotypes.insert(i+offset, (recomb, new_parent, new_ploidy))
            print 'end of round haplotypes:'
            print new_haplotypes
            
        curr_haplotypes = new_haplotypes
        #break
        
    print ''
    print 'final haplotypes'
    print curr_haplotypes
    print ''

    for h in curr_haplotypes:
        print h

    # output descendant sequence and ancestral haplotypes
    strain_seqs = {}
    #strain_names = ('AKRJ', 'AJ', 'BALBcJ', 'C3HHeJ', 'CASTEiJ', 'CBAJ', 'DBA2J', 'LPJ')
    strain_names = ('AKRJ', 'AJ', 'BALBcJ', 'C3HHeJ')

    #  read fasta files
    print 'reading fastas'
    for name in strain_names:
        with open('./data/mouse_fastas/%s_chr%i.fa' % (name, CHROM)) as f:
            f.readline()
            strain_seqs[name] = f.readline().strip()

    #  write descendant sequence and each segment's ancestral origin
    print 'writing output'
    num_loci = len(strain_seqs[strain_names[0]])
    with open('data/desc_segments.bed', 'w') as f_bed:
        with open('data/desc_seq.nuc', 'w') as f_nuc:
            for (start, anc_num, _), (end, _, _) in pairwise(curr_haplotypes + [(num_loci, -1, -1)]):
                print anc_num, start, end
                strain = strain_names[(anc_num-1) % len(strain_names)] #ancs are indexed starting at 1 in recomb file
                f_bed.write('chr%i\t%s\t%i\t%i\n' % (CHROM, strain, start, end))
                f_nuc.write(strain_seqs[strain][start:end]) 
            f_nuc.write('\n')
'''
        curr_start, curr_end, curr_anc = haplotypes[0]
        for start, end, anc in haplotypes[1:]:
            if COLLAPSE_FINAL_HAPLOTYPES:
                #collapse adjacent haplotypes from the same ancestor strain
                # different numbered ancestors can be the same strain (larger initial ancestor population adds diversity),
                # so may want to collapse them together for the final bed file. these separations are still recomb
                # events though, so may want to keep them in the final bed file.
                if anc != curr_anc:
                    f.write('chr%i\t%s\t%s\t%s\n' % (CHROM, curr_start, curr_end, curr_anc))
                    curr_start, curr_end, curr_anc = start, end, anc
                else:
                    curr_end = end
            else:
                f.write('chr%i\t%s\t%s\t%s\n' % (CHROM, curr_start, curr_end, curr_anc))
                curr_start, curr_end, curr_anc = start, end, anc
        f.write('chr%i\t%s\t%s\t%s\n' % (CHROM, curr_start, curr_end, curr_anc))
'''
