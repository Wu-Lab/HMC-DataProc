#!/usr/bin/env python

import os
import shutil
import glob
import HapMap

sample_dir = 'HapMap/'
sample_names = ['1000']
split_factors = [2, 4, 6, 8, 10]

chromosomes = list()
for num in range(1, 23):
    chromosomes.append(str(num))
chromosomes.append('x')
chromosomes.append('y')

populations = ['ceu', 'yri']

for name in sample_names:
    for factor in split_factors:
        split_name = name + '_' + str(factor)
        for chr in chromosomes:
            for pop in populations:
                work_dir = pop.upper() + '/chr' + chr + '/'
                source_dir = sample_dir + name + '/' + work_dir
                split_dir = sample_dir + split_name + '/' + work_dir
                phase_dir = sample_dir + split_name + '.inp/' + work_dir
                hpm2_dir = sample_dir + split_name + '.hpm2/' + work_dir
                if os.access(split_dir, os.F_OK) == False:
                    os.makedirs(split_dir)
                if os.access(phase_dir, os.F_OK) == False:
                    os.makedirs(phase_dir)
                if os.access(hpm2_dir, os.F_OK) == False:
                    os.makedirs(hpm2_dir)
                for source in glob.iglob(source_dir + '*.txt'):
                    prefix = split_dir + 'hapmap_chr' + chr + '_' + pop
                    HapMap.split_samples(source, prefix, factor)
                for source in glob.iglob(split_dir + '*.txt'):
                    basename = source[len(split_dir):len(source)-4]
                    phase_file = phase_dir + basename + '.inp'
                    hpm2_file = hpm2_dir + basename + '.hpm2'
                    if HapMap.filter_samples(source, 0.1):
                        print source
                        msg_file = phase_dir + 'HMC_1/' + basename + '.message'
                        if os.access(msg_file, os.F_OK):
                            os.remove(msg_file)
                        msg_file = phase_dir + 'HMC_2/' + basename + '.message'
                        if os.access(msg_file, os.F_OK):
                            os.remove(msg_file)
                        msg_file = hpm2_dir + 'haplorec_1/' + basename + '.message'
                        if os.access(msg_file, os.F_OK):
                            os.remove(msg_file)
                        msg_file = hpm2_dir + 'haplorec_2/' + basename + '.message'
                        if os.access(msg_file, os.F_OK):
                            os.remove(msg_file)
                    HapMap.convert_format_to_phase(source, phase_file, True)
                    HapMap.convert_format_to_hpm2(source, hpm2_file, True)
