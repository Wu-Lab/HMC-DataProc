#!/usr/bin/env python

import os
import shutil
import subprocess
import glob
import HapMap
import Genotype

size = '1000'
source_dir = 'HapMap/' + size + '/'
phase_dir = 'HapMap/' + size + '.inp/'
hpm2_dir = 'HapMap/' + size + '.hpm2/'

chromosomes = list()
for num in range(1, 23):
    chromosomes.append(str(num))
chromosomes.append('x')
chromosomes.append('y')

populations = ['ceu', 'yri']

for chr in chromosomes:
    for pop in populations:
        work_dir = pop.upper() + '/' + 'chr' + chr + '/'
        real_dir = source_dir + work_dir
        infer_dir = phase_dir + work_dir + 'HMC_1/'
        if os.access(infer_dir, os.F_OK) == False:
            raise RuntimeError, 'Can not find inferred genotype data!'
        compare_file = open(infer_dir + 'compare.txt', 'w')
        pattern = '*_chr' + chr + '_' + pop + '_*.txt'
        for source in glob.iglob(real_dir + pattern):
            print source
            data_name = source[len(real_dir):len(source)-4]
            infer_file = infer_dir + data_name + '.out.1.inp'
            genos_real = HapMap.read_hapmap(source)
            genos_infer = HapMap.read_phase(infer_file)
            comp = Genotype.compare(genos_real, genos_infer)
            result = data_name + ' ' + 'chr' + chr + ' ' + pop
            result += ' ' + data_name.split('_')[3]
            result += ' ' + data_name.split('_')[4]
            result += ' ' + str(comp['SE'])
            result += ' ' + str(comp['IGP'])
            result += ' ' + str(comp['IHP'])
            compare_file.write(result + '\n')

for chr in chromosomes:
    for pop in populations:
        work_dir = pop.upper() + '/' + 'chr' + chr + '/'
        real_dir = source_dir + work_dir
        infer_dir = hpm2_dir + work_dir + 'haplorec_1/'
        if os.access(infer_dir, os.F_OK) == False:
            raise RuntimeError, 'Can not find inferred genotype data!'
        compare_file = open(infer_dir + 'compare.txt', 'w')
        pattern = '*_chr' + chr + '_' + pop + '_*.txt'
        for source in glob.iglob(real_dir + pattern):
            print source
            data_name = source[len(real_dir):len(source)-4]
            infer_file = infer_dir + data_name + '.out.1.hpm2'
            genos_real = HapMap.read_hapmap(source)
            genos_infer = HapMap.read_hpm2(infer_file)
            comp = Genotype.compare(genos_real, genos_infer)
            result = data_name + ' ' + 'chr' + chr + ' ' + pop
            result += ' ' + data_name.split('_')[3]
            result += ' ' + data_name.split('_')[4]
            result += ' ' + str(comp['SE'])
            result += ' ' + str(comp['IGP'])
            result += ' ' + str(comp['IHP'])
            compare_file.write(result + '\n')
