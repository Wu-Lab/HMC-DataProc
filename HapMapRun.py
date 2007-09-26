#!/usr/bin/env python

import os
import shutil
import subprocess
import glob
import HapMap

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

chromosomes = chromosomes[6:24]

for chr in chromosomes:
    for pop in populations:
        work_dir = phase_dir + pop.upper() + '/' + 'chr' + chr + '/'
        result_dir = work_dir + 'HMC_1/'
        if os.access(result_dir, os.F_OK) == False:
            os.mkdir(result_dir)
        msg_file = open(result_dir + 'message.txt', 'w')
        pattern = '*_chr' + chr + '_' + pop + '_*.inp'
        for source in glob.iglob(work_dir + pattern):
            print source
            p = subprocess.Popen('HMC8.exe --nologo -i 1 -a 0.5 ' + source,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            msg_file.write(source + '\n')
            msg_file.writelines(p.stderr)
            msg_file.writelines(p.stdout)
            output_file = source + '.reconstructed'
            dest_file = result_dir + source[len(work_dir):len(source)-3] + 'out.1.inp'
            if os.access(output_file, os.F_OK):
                shutil.move(output_file, dest_file)

for chr in chromosomes:
    for pop in populations:
        work_dir = hpm2_dir + pop.upper() + '/' + 'chr' + chr + '/'
        result_dir = work_dir + 'haplorec_1/'
        if os.access(result_dir, os.F_OK) == False:
            os.mkdir(result_dir)
        msg_file = open(result_dir + 'message.txt', 'w')
        pattern = '*_chr' + chr + '_' + pop + '_*.hpm2'
        for source in glob.iglob(work_dir + pattern):
            print source
            p = subprocess.Popen('haplorec.bat -n 1 ' + source,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            msg_file.write(source + '\n')
            msg_file.writelines(p.stderr)
            msg_file.writelines(p.stdout)
            output_file = source[len(work_dir):len(source)] + '.reconstructed'
            dest_file = result_dir + source[len(work_dir):len(source)-4] + 'out.1.hpm2'
            if os.access(output_file, os.F_OK):
                shutil.move(output_file, dest_file)
