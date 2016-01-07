

ancs = ('YJM450', 'YJM451', 'YJM981', 'ySR127', 'DESC')
unk = 'YJM981'

with open('data/sim_SNPs.bed', 'w') as f_out:
    with open('data/all_SNPs.bed', 'r') as f_in:
        for line in f_in:
            cols = line.strip().split()
            new_SNPs = ''
            for strain in cols[-1].split('_'):
                if strain in ancs and strain != unk:
                    new_SNPs += strain + '_'
            if new_SNPs != '':
                f_out.write('\t'.join(cols[:3] + [new_SNPs[:-1]]) + '\n')



