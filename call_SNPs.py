
if __name__ == '__main__':
    #wkdir = '/scratch/Users/dama9282/simulation/'
    wkdir = './'

    strain_names = ['AKRJ', 'AJ', 'BALBcJ', 'C3HHeJ', 'CASTEiJ', 'CBAJ', 'DBA2J', 'LPJ']
    chrom = 19

    strain_seqs = {}

    # read descendant sequence
    with open(wkdir + './data/desc_seq.nuc') as f:
        strain_seqs['DESC'] = f.readline().strip()
    
    # read fasta files
    print 'reading fastas'
    print ' ref'
    strain_seqs['ref'] = ''
    with open(wkdir + './data/mouse_fastas/chr%i.fa' % chrom) as f:
        f.readline()
        strain_seqs['ref'] = ''.join(f.readlines()).replace('\n','')
        print len(strain_seqs['ref'])

    for name in strain_names:
        print ' ' + name
        with open(wkdir + './data/mouse_fastas/%s_chr%i.fa' % (name, chrom)) as f:
            f.readline()
            strain_seqs[name] = f.readline().strip()
            print len(strain_seqs[name])

    # create bed file
    print 'creating bed file'
    with open(wkdir + './data/sim_SNPs.bed', 'w') as f:
        for i, ref_nuc in enumerate(strain_seqs['ref']):
            ancs = [name for name in ['DESC'] + sorted(strain_names) if strain_seqs[name][i] != ref_nuc]
            if len(ancs) > 0:
                f.write('chr%i\t%i\t%i\t%s\n' % (chrom, i, i+1, '_'.join(ancs)))
            if i % 100000 == 0:
                print i

