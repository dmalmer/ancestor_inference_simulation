
from collections import defaultdict


def read_recomb_rates(filename):
    recomb_rates_by_chr = defaultdict(list)
    with open(filename, 'r') as f:
        f.readline()
        line = f.readline()
        while line != '':
            chrom, pos_kb, rate = line.strip().split(',')
            recomb_rates_by_chr[chrom].append([int(float(pos_kb)*1000), float(rate)])
            line = f.readline()

    return recomb_rates_by_chr


def read_bed(filename):
    SNPs_by_chr = defaultdict(list)
    unique_ancs = set()
    with open(filename) as f:
        for line in f:
            chrom, pos, _, ancs_str = line.strip().split()
            ancs = ancs_str.split('_')
            unique_ancs = unique_ancs.union(ancs)
            SNPs_by_chr[chrom].append((int(pos), ancs))

    return SNPs_by_chr, unique_ancs


def genetic_distance(SNP_start, SNP_end, recomb_main_i, recomb_map):
    # In recomb file, each position has an associated rate
    #   chr1,300500,.03
    #   chr1,300510,.05
    # The rate of regions between markers is the value next to the starting marker, eg. above: 300500-300510 = .03

    # Find recomb_map starting position
    if recomb_main_i is None:
        recomb_main_i = 0
        while recomb_main_i < len(recomb_map) and int(recomb_map[recomb_main_i][0]) < SNP_start:
            recomb_main_i += 1
        recomb_start_i = max(recomb_main_i-1, 0)
    else:
        recomb_start_i = recomb_main_i - 1

    # Quick check to make sure recomb_start_i >= 0 (should be, as recomb_index should always be >=1 in else statement)
    if recomb_start_i < 0:
        raise Exception('recomb_start_i should never be less than 0')

    # Find recomb_map ending position
    while (recomb_main_i < len(recomb_map) and int(recomb_map[recomb_main_i][0]) < SNP_end) or recomb_main_i == 0:
        recomb_main_i += 1
    recomb_end_i = recomb_main_i

    # Calc recomb rates
    #  First, special case for SNPs between adjacent genetic markers
    if recomb_end_i - recomb_start_i == 1:
        genetic_dist = ((SNP_end - SNP_start) / 1000.) * recomb_map[recomb_start_i][1]

    #  Otherwise, calc with all recomb rates between SNPs
    else:
        # Proportional rate of first SNP
        genetic_dist = ((recomb_map[recomb_start_i+1][0] - SNP_start) / 1000.) * recomb_map[recomb_start_i][1]

        # Rates in the middle
        for i in range(recomb_start_i+1, recomb_end_i-1):
            genetic_dist += ((recomb_map[i+1][0] - recomb_map[i][0]) / 1000.) * recomb_map[i][1]

        # Proportional rate of second SNP
        genetic_dist += ((SNP_end - recomb_map[recomb_end_i-1][0]) / 1000.) * recomb_map[recomb_end_i-1][1]

    return genetic_dist, recomb_main_i


if __name__ == '__main__':
    chrom = 'chr19'

    recomb_rates_by_chr = read_recomb_rates('./data/mouse_recomb_rates.csv')
    SNPs_by_chr, unique_ancs = read_bed('./data/sim_SNPs.bed')
    unique_ancs.discard('DESC')
    sorted_ancs = sorted(unique_ancs)

    #determine genetic distance of each SNP position
    genetic_distances = [] #genetic distance from previous marker
    recomb_i = None
    prev_pos = recomb_rates_by_chr[chrom][0][0]
    skipped = 0
    for curr_pos, ancs in SNPs_by_chr[chrom][1:]:
        if prev_pos >= curr_pos:
            skipped += 1
            continue

        gen_dist, recomb_i = genetic_distance(prev_pos, curr_pos, recomb_i, recomb_rates_by_chr[chrom])
        genetic_distances.append(gen_dist)

        prev_pos = curr_pos
    # since first position doesn't have a previous marker, give it the genetic distance of the second position
    genetic_distances.insert(0, genetic_distances[0])

    #write RABBIT input files
    # each position is a column, so need to gather all columns before writing each row
    SNP_names = []
    cM_dists = []
    anc_alleles = {anc: [] for anc in unique_ancs}
    desc_alleles = []
    curr_cM_pos = 0.
    for (pos, ancs), dist in zip(SNPs_by_chr[chrom], genetic_distances):
        SNP_names.append('pos%i' % pos)
        curr_cM_pos += dist
        cM_dists.append('%.8f' % curr_cM_pos)

        for anc in unique_ancs:
            if anc in ancs:
                anc_alleles[anc].append('2')
            else:
                anc_alleles[anc].append('1')
        if 'DESC' in ancs:
            desc_alleles.append('2')
        else:
            desc_alleles.append('1')

    with open('./data/RABBIT_markers.csv', 'w') as f:
        f.write('#founders,%i\n' % len(unique_ancs))
        f.write('SNP,%s\n' % ','.join(SNP_names))
        f.write('Chromosome,%s\n' % ','.join([chrom.replace('chr','')] * len(SNP_names)))
        f.write('cM,%s\n' % ','.join(cM_dists))
        for anc in sorted_ancs:
            f.write('%s,%s\n' % (anc, ','.join(anc_alleles[anc])))
        f.write('DESC,%s\n' % ','.join(desc_alleles))

    #write HAPPY input files
    # write alleles file
    with open('./data/HAPPY_markers.csv', 'w') as f:
        f.write('markers %i strains %i\n' % (len(SNP_names), len(unique_ancs)))
        f.write('strain_names %s\n' % '\t'.join(sorted_ancs))
        curr_cM_pos = 0.
        for i, ((pos, ancs), dist) in enumerate(zip(SNPs_by_chr[chrom], genetic_distances)):
            curr_cM_pos += dist
            f.write('marker pos%i 2 %.8f\n' % (pos, curr_cM_pos))
            #create arrays of 1. and 0. corresponding to if the ancestor contains or doesn't contain the current SNP
            missing_SNP = [float(a not in ancs) for a in sorted_ancs]
            have_SNP = [float(a in ancs) for a in sorted_ancs]
            #divide the 1's by the total number of 1's to get the probability of that allele (eg. five 1's turn into .2's)
            tot_missing = len([val for val in missing_SNP if val == 1.])
            tot_have = len(sorted_ancs) - tot_missing
            f.write('allele\t0\t%s\n' % '\t'.join(['%.3f' % (val/tot_missing) for val in missing_SNP]))
            f.write('allele\t1\t%s\n' % '\t'.join(['%.3f' % (val/tot_have) for val in have_SNP]))
    
    # write data file (ped format)

         







