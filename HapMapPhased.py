#!/usr/bin/env python

# parse HapMap genotypes data
# infer the haplotype from trios data
# input files are HapMap genotypes data (CEU or YRI) and pedigree information

# used python packages
import os
import HapMap

op = dict()
op['parse'] = True
op['sample'] = True

sample_dir = 'HapMap/'
sample_sizes = [1000]

chromosomes = list()
for num in range(1, 23):
    chromosomes.append(str(num))
chromosomes.append('x')
chromosomes.append('y')

populations = ['jpt+chb']

filter = [0, 0.95, 0.05]

for chr in chromosomes:
    for pop in populations:
        print 'Chromosome ' + chr + ' in ' + pop.upper() + ' population'
        pop_dir = pop.upper() + '/'
        genotypes = pop_dir + 'genotypes_chr' + chr + '_' + pop.upper() + '_r22_nr.b36_fwd.phase.gz'
        legend = pop_dir + 'genotypes_chr' + chr + '_' + pop.upper() + '_r22_nr.b36_fwd_legend.txt.gz'
        haplotypes = pop_dir + 'hapmap_chr' + chr + '_' + pop + '.txt.gz'
        if os.access(genotypes, os.F_OK) == False or os.access(legend, os.F_OK) == False:
            continue
        if op['parse']:
            HapMap.parse_phased(genotypes, legend, haplotypes, 'chr' + chr)
        if op['sample']:
            for size in sample_sizes:
                dest_dir = sample_dir + str(size) + '/' + pop_dir + 'chr' + chr + '/'
                if os.access(dest_dir, os.F_OK) == False:
                    os.makedirs(dest_dir)
                prefix = dest_dir + 'hapmap_chr' + chr + '_' + pop
                HapMap.trios_samples(haplotypes, prefix, size)
