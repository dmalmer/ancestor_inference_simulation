
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


genomes = [] #index+1 = indv ID
with open('data/recombs.log','r') as f:
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
    
    print '6:'
    print genomes[5][0]
    print genomes[5][1]
    print '7:'
    print genomes[6][0]
    print genomes[6][1]

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

        print '\ngenome ' + str(len(genomes))
        print genomes[-1][0]
        print genomes[-1][1]

    
exit(0)
for i, g in enumerate(genomes):
    print str(i+1) + ':'
    print '   ' + str(g['0'])
    print '   ' + str(g['1'])


