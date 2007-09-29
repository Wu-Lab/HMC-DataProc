#!/usr/bin/env python

# parse HapMap genotypes data
# infer the haplotype from trios data
# input files are HapMap genotypes data (CEU or YRI) and pedigree information
# output files are PHASE format files

# used python packages
import HapMap

op = dict()
op['phase'] = True
op['filter'] = True
op['sort'] = True
op['sample'] = True

sample_dir = 'HapMap/'
sample_sizes = [100, 200, 500, 1000]

chromosomes = list()
for num in range(1, 23):
    chromosomes.append(str(num))
chromosomes.append('x')
chromosomes.append('y')

populations = ['ceu', 'yri']

for chr in chromosomes:
    for pop in populations:
        print 'Chromosome ' + chr + ' in ' + pop.upper() + ' population'
        pop_dir = pop.upper() + '/'
        genotypes = pop_dir + 'genotypes_chr' + chr + '_' + pop + '_r22_nr.b36_fwd.txt.gz'
        pedinfo = 'pedinfo/pedinfo2sample_' + pop.upper() + '.txt.gz'
        haplotypes = pop_dir + 'hapmap_unfiltered_chr' + chr + '_' + pop + '.txt.gz'
        filtered = pop_dir + 'hapmap_chr' + chr + '_' + pop + '_[0.8_0.8_0.05].txt.gz'
        sorted = pop_dir + 'hapmap_chr' + chr + '_' + pop + '_[0.8_0.8_0.05]_sorted.txt.gz'
        if op['phase']:
            HapMap.trios_phase(genotypes, pedinfo, haplotypes)
        if op['filter']:
            HapMap.trios_filter(haplotypes, filtered, [0.8, 0.8, 0.05])
        if op['sort']:
            HapMap.trios_sort(filtered, sorted)
        if op['sample']:
            for size in sample_sizes:
                dest_dir = sample_dir + str(size) + '/'
                if os.access(dest_dir, os.F_OK) == False:
                    os.mkdirs(dest_dir)
                samples = dest_dir + 'hapmap_chr' + chr + '_' + pop
                start = 0
                while(start >= 0):
                    start = HapMap.select_samples(sorted, samples, start, 100)
