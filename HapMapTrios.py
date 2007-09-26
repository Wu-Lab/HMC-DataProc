#!/usr/bin/env python

# parse HapMap genotypes data
# infer the haplotype from trios data
# input files are HapMap genotypes data (CEU or YRI) and pedigree information
# output files are PHASE format files

# used python packages
import HapMap

chromosomes = list()
for num in range(1, 23):
    chromosomes.append(str(num))
chromosomes.append('x')
chromosomes.append('y')

populations = ['ceu', 'yri']

for chr in chromosomes:
    for pop in populations:
        print 'Chromosome ' + chr + ' in ' + pop.upper() + ' population'
        genotypes = pop.upper() + '/genotypes_chr' + chr + '_' + pop + '_r22_nr.b36_fwd.txt.gz'
        pedinfo = 'pedinfo/pedinfo2sample_' + pop.upper() + '.txt.gz'
        haplotypes = pop.upper() + '/hapmap_unfiltered_chr' + chr + '_' + pop + '.txt.gz'
        filtered = pop.upper() + '/hapmap_chr' + chr + '_' + pop + '_[0.8_0.8_0.05].txt.gz'
        sorted = pop.upper() + '/hapmap_chr' + chr + '_' + pop + '_[0.8_0.8_0.05]_sorted.txt.gz'
        HapMap.trios_phase(genotypes, pedinfo, haplotypes)
        HapMap.trios_filter(haplotypes, filtered, [0.8, 0.8, 0.05])
        HapMap.trios_sort(filtered, sorted)
        samples = 'HapMap/100/hapmap_chr' + chr + '_' + pop
        start = 0
        while(start >= 0):
            start = HapMap.select_samples(sorted, samples, start, 100)
