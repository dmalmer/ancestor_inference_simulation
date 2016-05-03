
import sys

if __name__ == '__main__':
    wkdir = '/scratch/Users/dama9282/simulation/'
    #wkdir = './'

    append_str = '_%s' % sys.argv[1] if len(sys.argv) > 1 else ''

    strain_names = ['AKRJ', 'AJ', 'BALBcJ', 'C3HHeJ', 'CASTEiJ', 'CBAJ', 'DBA2J', 'LPJ']
    chrom = 1

    strain_seqs = {}

    # read fasta files
    print 'reading fastas'
    print ' ref'
    strain_seqs['ref'] = ''
    with open(wkdir + './data/mouse_fastas/chr%i.fa' % chrom) as f:
        f.readline()
        strain_seqs['ref'] = ''.join(f.readlines()).replace('\n','')

    for name in strain_names:
        print ' ' + name
        with open(wkdir + './data/mouse_fastas/%s_chr%i.fa' % (name, chrom)) as f:
            f.readline()
            strain_seqs[name] = f.readline().strip()

    # call SNPs on files with and without genotyping error applied
    for err_str in ('', '_err'):
        # read descendant sequence
        with open(wkdir + './data/desc_seq%s%s.nuc' % (err_str, append_str)) as f:
            strain_seqs['DESC'] = f.readline().strip()

        # create bed file
        print 'creating bed file'
        with open(wkdir + './data/sim_SNPs%s%s.bed' % (err_str, append_str), 'w') as f:
            for i, ref_nuc in enumerate(strain_seqs['ref']):
                ancs = [name for name in ['DESC'] + sorted(strain_names) if strain_seqs[name][i] != ref_nuc]
                #ignore lines where all ancestors and desc have same SNPs
                if len(ancs) > 0 and len(ancs) != len(strain_names) + 1:
                    f.write('chr%i\t%i\t%i\t%s\n' % (chrom, i, i+1, '_'.join(ancs)))
                if i % 1000000 == 0:
                    print i

