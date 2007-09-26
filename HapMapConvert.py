#!/usr/bin/env python

import os
import shutil
import glob
import HapMap

size = '100'
source_dir = 'HapMap/' + size + '/'
phase_dir = 'HapMap/' + size + '.inp/'
hpm2_dir = 'HapMap/' + size + '.hpm2/'

if os.access(phase_dir, os.F_OK) == False:
    os.mkdir(phase_dir)
if os.access(hpm2_dir, os.F_OK) == False:
    os.mkdir(hpm2_dir)

for source in glob.iglob(source_dir + '*.txt'):
    basename = source[len(source_dir):len(source)-4]
    phase = phase_dir + basename + '.inp'
    hpm2 = hpm2_dir + basename + '.hpm2'
    HapMap.convert_format_to_phase(source, phase)
    HapMap.convert_format_to_hpm2(source, hpm2)

chromosomes = list()
for num in range(1, 23):
    chromosomes.append(str(num))
chromosomes.append('x')
chromosomes.append('y')

populations = ['ceu', 'yri']

file_types = ['.txt', '.inp', '.hpm2']
file_dirs = [source_dir, phase_dir, hpm2_dir]

for chr in chromosomes:
    for pop in populations:
        for type in [0, 1, 2]:
            dest_dir = file_dirs[type]
            dest_dir += pop.upper() + '/'
            if os.access(dest_dir, os.F_OK) == False:
                os.mkdir(dest_dir)
            dest_dir += 'chr' + chr + '/'
            if os.access(dest_dir, os.F_OK) == False:
                os.mkdir(dest_dir)
            pattern = '*_chr' + chr + '_' + pop + '_*' + file_types[type]
            for source in glob.iglob(file_dirs[type] + pattern):
                dest = dest_dir + source[len(file_dirs[type]):len(source)]
                shutil.move(source, dest)

