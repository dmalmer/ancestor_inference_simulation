
from glob import glob
from matplotlib import pyplot
from collections import defaultdict
from itertools import izip

#@profile
def compare_seqs(seqA, seqB, filename):
    positions = []
    misses_ct = defaultdict(int)
    try:
        assert(len(seqA) == len(seqB))
        matches = 0
        for i, (a, b) in enumerate(izip(seqA, seqB)):
            if a == b:
                matches += 1
                positions.append(1)
            else:
                positions.append(0)
                misses_ct[a + '-' + b] += 1
                pass#print i, 'miss'

        print misses_ct

        pyplot.figure()
        pyplot.plot(range(len(seqA)), positions)
        pyplot.title(filename)
        pyplot.ylim([-.1,1.1])
        pyplot.savefig('./data/plot_' + filename.rsplit('/',1)[1].split('.')[0] + '.png')

        del positions
        del misses_ct

        return float(matches)/len(seqA)
    except AssertionError:
        return -1


#@profile
if __name__ == '__main__':
    with open('./data/desc_seq.nuc') as f:
        seqB = f.readline().strip()

    for filename in glob('./data/last-indv-seq_1.nuc'):
        with open(filename) as f:
            seqA = f.readline().strip()
        print filename.rsplit('/',1)[1][14:]
        print compare_seqs(seqA, seqB, filename)
    


