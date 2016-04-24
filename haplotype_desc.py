
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
from glob import glob
from os import path


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)


def parent_index(parent, ploidy, anc_pop):
    #return index of genomes[] list that corresponds to parent number in recombs.log file
    return (parent - anc_pop - 1) * 2 + ploidy


if __name__ == '__main__':
    #get most recently created recombs file
    recomb_file = sorted(filter(path.isfile, glob('./data/recombs_*.log')), key=lambda x: path.getmtime(x))[-1]
    strain_names = ('AKRJ', 'AJ', 'BALBcJ', 'C3HHeJ', 'CASTEiJ', 'CBAJ', 'DBA2J', 'LPJ')
    chrom = 19

    # re-trace what parts of each genome came from which ancestor using the recombination file output by simuPOP
    print 'reading recomb file: ' + recomb_file
    genomes = [] #index+1 = indv ID
    genomes = open(recomb_file,'r').readlines()
    anc_pop = int(genomes[0].split()[0]) - 1

    #only care about last genome, start from there and work backwords
    #loop over both ploidy's to check homozygosity
    final_segments = {-2: [], -1: []}
    for desc_ind in (-2, -1):
        #haplotypes = [(pos1, parent1, ploidy1), (pos2, parent2, ploidy2), ..., (posN, parentN, ploidyN)]]
        # add parental haplotypes of descendant genome
        init_line = genomes[desc_ind].strip().split()
        parent = int(init_line[1])
        ploidy = int(init_line[2])
        curr_haplotypes = [(0, parent, ploidy)]
        for recomb in init_line[3:]:
            curr_haplotypes.append((int(recomb), parent, abs(ploidy-1)))

        # replace parental haplotypes with their parents' haplotypes until we're only left with original ancestors
        while not all([hap[1] <= anc_pop for hap in curr_haplotypes]):
            #create new replacement
            new_haplotypes = list(curr_haplotypes)

            #keep track of i offset as we insert new parental haplotypes
            offset = 0 
            for i, (curr_pos, curr_parent, curr_ploidy) in enumerate(curr_haplotypes):
                #get ending position of current haplotype we are replacing with parent haplotypes
                hap_end_pos = curr_haplotypes[i+1][0] if i+1 < len(curr_haplotypes) else -1
                #get parent line of current haplotype
                parent_line = genomes[parent_index(curr_parent, curr_ploidy, anc_pop)].strip().split()
                
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
                
            #replace haplotypes
            curr_haplotypes = new_haplotypes

        # output descendant sequence and ancestral haplotypes
        strain_seqs = {}

        #  read fasta files
        for name in strain_names:
            with open('./data/mouse_fastas/%s_chr%i.fa' % (name, chrom)) as f:
                f.readline()
                strain_seqs[name] = f.readline().strip()

        #  write descendant sequence and each segment's ancestral origin
        num_loci = len(strain_seqs[strain_names[0]])
        with open('data/desc_segments.bed', 'w') as f_bed:
            with open('data/desc_seq.nuc', 'w') as f_nuc:
                for (start, anc_num, _), (end, _, _) in pairwise(curr_haplotypes + [(num_loci, -1, -1)]):
                    strain = strain_names[(anc_num-1) % len(strain_names)] #ancs are indexed starting at 1 in recomb file
                    f_bed.write('chr%i\t%s\t%i\t%i\n' % (chrom, strain, start, end))
                    f_nuc.write(strain_seqs[strain][start:end]) 

                    final_segments[desc_ind].append((strain, anc_num, start, end))
                f_nuc.write('\n')

    #make final segments same length
    if len(final_segments[-2]) > len(final_segments[-1]):
        final_segments[-1].extend([('_', -1, -1, -1)] * (len(final_segments[-2]) - len(final_segments[-1])))
    else:
        final_segments[-2].extend([('_', -1, -1, -1)] * (len(final_segments[-1]) - len(final_segments[-2])))

    #output final segments
    # determine collapsed final segments as we're outputing uncollapsed segments
    collapsed_final_segments = {-2: [], -1: []}
    curr_strain_0, _, curr_start_0, curr_end_0 = final_segments[-2][0]
    curr_strain_1, _, curr_start_1, curr_end_1 = final_segments[-1][0]

    print '\nuncollapsed:'
    print 'ploidy 0\t\t\t\t\tploidy 1'
    for i, ((strain_0, anc_num_0, start_0, end_0), (strain_1, anc_num_1, start_1, end_1)) in enumerate(zip(final_segments[-2], final_segments[-1])):
        print '%i - %i: %i_%s' % (start_0, end_0, anc_num_0, strain_0),
        print '\t\t\t',
        print '%i - %i: %i_%s' % (start_1, end_1, anc_num_1, strain_1)

        if i == 0:
            continue

        if strain_0 != '_' and curr_strain_0 != strain_0:
            collapsed_final_segments[-2].append((curr_strain_0, curr_start_0, curr_end_0))
            curr_strain_0, curr_start_0, curr_end_0 = strain_0, start_0, end_0
        elif end_0 != -1:
            curr_end_0 = end_0

        if strain_1 != '_' and curr_strain_1 != strain_1:
            collapsed_final_segments[-1].append((curr_strain_1, curr_start_1, curr_end_1))
            curr_strain_1, curr_start_1, curr_end_1 = strain_1, start_1, end_1
        elif end_1 != -1:
            curr_end_1 = end_1

    #make collapsed final segments same length
    collapsed_final_segments[-2].append((curr_strain_0, curr_start_0, curr_end_0))
    collapsed_final_segments[-1].append((curr_strain_1, curr_start_1, curr_end_1))
    if len(collapsed_final_segments[-2]) > len(collapsed_final_segments[-1]):
        collapsed_final_segments[-1].extend([('_', -1, -1)] * (len(collapsed_final_segments[-2]) - len(collapsed_final_segments[-1])))
    else:
        collapsed_final_segments[-2].extend([('_', -1, -1)] * (len(collapsed_final_segments[-1]) - len(collapsed_final_segments[-2])))

    print '\ncollapsed:'
    print 'ploidy 0\t\t\t\t\tploidy 1'
    for ((strain_0, start_0, end_0), (strain_1, start_1, end_1)) in zip(collapsed_final_segments[-2], collapsed_final_segments[-1]):
        print '%i - %i: %s' % (start_0, end_0, strain_0),
        print '\t\t\t',
        print '%i - %i: %s' % (start_1, end_1, strain_1)

