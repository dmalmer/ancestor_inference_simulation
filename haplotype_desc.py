
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
    RECOMB_FILE = 'data/recombs_rmouse_g20.log'
    #RECOMB_FILE = 'data/recombs_r0.01_g20.log'

    CHROM = 19
    
    COLLAPSE_FINAL_HAPLOTYPES = False

    # re-trace what parts of each genome came from which ancestor using the recombination file output by simuPOP
    print 'reading recomb file'
    genomes = [] #index+1 = indv ID
    with open(RECOMB_FILE,'r') as f:
        cols = f.readline().strip().split()
        first_gen = int(cols[0])

        #set ancestor genomes
        for a in range(1, first_gen):
            genomes.append({c: [[0],[a]] for c in (0, 1)})  

        #set first gen genomes (ancestors are homozygous, so ignore recombs)
        while int(cols[1]) < first_gen:
            genomes.append({})
            genomes[-1][0] = [[0],[int(cols[1])]]  
            cols = f.readline().strip().split()
            genomes[-1][1] = [[0],[int(cols[1])]]  
            
            cols = f.readline().strip().split()

        #set rest of genomes
        while len(cols) > 0:
            genomes.append({})
            for g_chrom in (0, 1):
                parent = genomes[int(cols[1])-1] 
                curr_c = int(cols[2])
                new_pos = [parent[curr_c][0][0]] #always 0
                new_anc = [parent[curr_c][1][0]]
                p_i = [0, 0]
                p_i[curr_c] += 1
                for r in cols[3:]:
                    r = int(r)
                    #add anc info from before recomb
                    while p_i[curr_c] < len(parent[curr_c][0]) and parent[curr_c][0][p_i[curr_c]] < r:
                        if parent[curr_c][0][p_i[curr_c]] > new_pos[-1]:
                            new_pos.append(parent[curr_c][0][p_i[curr_c]])
                            new_anc.append(parent[curr_c][1][p_i[curr_c]])
                        p_i[curr_c] += 1
                    #add recomb
                    new_pos.append(r)
                    curr_c = abs(curr_c - 1)
                    #  catch other chrom parent index to current position 
                    while p_i[curr_c] < len(parent[curr_c][0])-1 and parent[curr_c][0][p_i[curr_c]+1] < r:
                        p_i[curr_c] += 1
                    if p_i[curr_c] < len(parent[curr_c][1]):
                        new_anc.append(parent[curr_c][1][p_i[curr_c]])
                    else:
                        #if the recombination switches to a chrom that is already completely iterated over, we need to simply
                        # add the last ancestor
                        new_anc.append(parent[curr_c][1][-1])

                #add rest of parent after all recombinations
                new_pos.extend(parent[curr_c][0][p_i[curr_c]+1:])
                new_anc.extend(parent[curr_c][1][p_i[curr_c]+1:])
                
                #collapse adjacent positions with the same ancestor and add to genomes list
                genomes[-1][g_chrom] = [[new_pos[0]], [new_anc[0]]]
                for n_p, n_a in zip(new_pos[1:], new_anc[1:]):
                    if n_a != genomes[-1][g_chrom][1][-1]:
                        genomes[-1][g_chrom][0].append(n_p)
                        genomes[-1][g_chrom][1].append(n_a)

                cols = f.readline().strip().split()

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
    haplotypes = [] #list of tuples, (start, end, strain_name) of each segment
    with open('data/desc_seq.nuc', 'w') as f:
        for (start, end), anc_num in zip(pairwise(genomes[-1][0][0] + [num_loci]), genomes[-1][0][1]):
            anc_index = (anc_num-1) % len(strain_names) #ancs are indexed starting at 1 in recomb file
            #write nucleotide sequence
            f.write(strain_seqs[strain_names[anc_index]][start:end]) 
            #record ancestral haplotype
            haplotypes.append((start, end, strain_names[anc_index]))
        #f.write('-' * (len(strain_seqs[strain_names[0]]) - len(desc_seq)) + '\n')

    with open('data/desc_segments.bed', 'w') as f:
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

