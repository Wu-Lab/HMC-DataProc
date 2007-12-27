#!/usr/bin/env python

# parse HapMap genotypes data
# infer the haplotype from trios data
# input files are HapMap genotypes data (CEU or YRI) and pedigree information
# output files are PHASE format files

# used python packages
import os
import HapMap

op = dict()
op['phase'] = True
op['filter'] = True
op['sort'] = True
op['sample'] = True

sample_dir = 'HapMap/'
sample_sizes = [1000]

chromosomes = list()
for num in range(1, 23):
    chromosomes.append(str(num))
chromosomes.append('x')
chromosomes.append('y')

chromosomes = ['1']

populations = ['ceu', 'yri']

filter = [0, 0.95, 0.05]

for chr in chromosomes:
    for pop in populations:
        print 'Chromosome ' + chr + ' in ' + pop.upper() + ' population'
        pop_dir = pop.upper() + '/'
        genotypes = pop_dir + 'genotypes_chr' + chr + '_' + pop + '_r22_nr.b36_fwd.txt.gz'
        pedinfo = 'pedinfo/pedinfo2sample_' + pop.upper() + '.txt.gz'
        haplotypes = pop_dir + 'hapmap_unfiltered_chr' + chr + '_' + pop + '.txt.gz'
        filter_str = '[' + str(filter[0]) + '_' + str(filter[1]) + '_' + str(filter[2]) + ']'
        filtered = pop_dir + 'hapmap_chr' + chr + '_' + pop + '_' + filter_str + '.txt.gz'
        sorted = pop_dir + 'hapmap_chr' + chr + '_' + pop + '_' + filter_str + '_sorted.txt.gz'
        if op['phase']:
            HapMap.trios_phase(genotypes, pedinfo, haplotypes)
        if op['filter']:
            HapMap.trios_filter(haplotypes, filtered, filter)
        if op['sort']:
            HapMap.trios_sort(filtered, sorted)
        if op['sample']:
            for size in sample_sizes:
                dest_dir = sample_dir + str(size) + '/' + pop_dir + 'chr' + chr + '/'
                if os.access(dest_dir, os.F_OK) == False:
                    os.makedirs(dest_dir)
                prefix = dest_dir + 'hapmap_chr' + chr + '_' + pop
                HapMap.trios_samples(sorted, prefix, size)
