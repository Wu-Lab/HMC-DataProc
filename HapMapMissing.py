#!/usr/bin/env python

import os
import shutil
import glob
import HapMap

sample_dir = 'HapMap/'
sample_names = ['100', '200', '500', '1000']
missing_factors = [0.01, 0.05, 0.1]

chromosomes = list()
for num in range(1, 23):
    chromosomes.append(str(num))
chromosomes.append('x')
chromosomes.append('y')

populations = ['ceu', 'yri']

for name in sample_names:
    for factor in missing_factors:
        missing_name = name + '_' + str(factor)
        for chr in chromosomes:
            for pop in populations:
                work_dir = pop.upper() + '/chr' + chr + '/'
                source_dir = sample_dir + name + '/' + work_dir
                missing_dir = sample_dir + missing_name + '/' + work_dir
                phase_dir = sample_dir + missing_name + '.inp/' + work_dir
                hpm2_dir = sample_dir + missing_name + '.hpm2/' + work_dir
                if os.access(missing_dir, os.F_OK) == False:
                    os.makedirs(missing_dir)
                if os.access(phase_dir, os.F_OK) == False:
                    os.makedirs(phase_dir)
                if os.access(hpm2_dir, os.F_OK) == False:
                    os.makedirs(hpm2_dir)
                for source in glob.iglob(source_dir + '*.txt'):
                    basename = source[len(source_dir):len(source)-4]
                    missing_file = missing_dir + basename + '.txt'
                    missing_mask = missing_dir + basename + '.mask'
                    shutil.copy(source, missing_file)
                    genos = HapMap.read_hapmap(missing_file)
                    mask = Genotype.MissingMask(len(genos), genos[0].len(), factor)
                    mask.write(missing_mask)
                for source in glob.iglob(missing_dir + '*.txt'):
                    basename = source[len(missing_dir):len(source)-4]
                    mask_file = missing_dir + basename + '.mask'
                    phase_file = phase_dir + basename + '.inp'
                    hpm2_file = hpm2_dir + basename + '.hpm2'
                    mask.read(mask_file)
                    HapMap.convert_format_to_phase(source, phase_file, True, mask)
                    HapMap.convert_format_to_hpm2(source, hpm2_file, True, mask)
