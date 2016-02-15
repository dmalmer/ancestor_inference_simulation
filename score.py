
from intervaltree import IntervalTree

true_segments = IntervalTree()

correct = 0
incorrect = 0

with open('data/desc_segments.bed', 'r') as f:
    for line in f:
        cols = line.strip().split()
        true_segments[int(cols[1]):int(cols[2])] = cols[3]

with open('data/sim_SNPs_hmm-out.bed', 'r') as f:
    for line in f:
        cols = line.strip().split()
        ancs = cols[3].split('_')
        # ugly and slow, but good enough
        for i in range(int(cols[1]), int(cols[2])+1):
            intv_set = true_segments.search(i)
            if len(intv_set) > 1:
                raise Exception('big set')
            for interval in intv_set:
                if interval.data in ancs:
                    correct += 1
                else:
                    incorrect += 1


print correct
print incorrect
print float(correct) / (correct + incorrect)

