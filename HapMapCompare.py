#!/usr/bin/env python

import os
import shutil
import glob
import HapMap
import Genotype

sample_dir = 'HapMap/'
sample_names = ['1000']

params = {'HMC':{}, 'haplorec':{}, 'PHASE':{}}

params['HMC']['enable'] = True
params['HMC']['output'] = 'HMC_1'
params['HMC']['suffix'] = '.inp'
params['HMC']['readfile'] = HapMap.read_phase

params['haplorec']['enable'] = True
params['haplorec']['output'] = 'haplorec_1'
params['haplorec']['suffix'] = '.hpm2'
params['haplorec']['readfile'] = HapMap.read_hpm2

params['PHASE']['enable'] = False
params['PHASE']['output'] = 'PHASE_1'
params['PHASE']['suffix'] = '.inp'
params['PHASE']['readfile'] = None

chromosomes = list()
for num in range(1, 23):
    chromosomes.append(str(num))
chromosomes.append('x')
chromosomes.append('y')

populations = ['ceu', 'yri']

for method in params.keys():
    p = params[method]
    if p['enable']:
        for name in sample_names:
            source_dir = sample_dir + name + '/'
            work_dir = sample_dir + name + p['suffix'] + '/'
            output_dir = work_dir + p['output'] + '/'
            if os.access(output_dir, os.F_OK) == False:
                os.makedirs(output_dir)
            for chr in chromosomes:
                for pop in populations:
                    sub_dir = pop.upper() + '/' + 'chr' + chr + '/'
                    real_dir = source_dir + sub_dir
                    infer_dir = work_dir + sub_dir + p['output'] + '/'
                    if os.access(infer_dir, os.F_OK) == False:
                        raise RuntimeError, 'Can not find inferred genotype data!'
                    output_file = 'error_' + pop + '_chr' + chr + '_' + p['output'] + '.txt'
                    output_file = output_dir + output_file
                    output = open(output_file, 'w')
                    pattern = '*_chr' + chr + '_' + pop + '_*.txt'
                    for source in glob.glob(real_dir + pattern):
                        print source
                        data_name = source[len(real_dir):len(source)-4]
                        infer_file = infer_dir + data_name + '.out'
                        if os.access(infer_file, os.F_OK) == False:
                            continue
                        genos_real = HapMap.read_hapmap(source)
                        genos_infer = p['readfile'](infer_file)
                        comp = Genotype.compare(genos_real, genos_infer)
                        result = data_name + ' ' + 'chr' + chr + ' ' + pop
                        result += ' ' + data_name.split('_')[3]
                        result += ' ' + data_name.split('_')[4]
                        for k in ['SE', 'IGP', 'IHP']:
                            result += ' ' + str(comp[k])
                        for k in ['SE', 'IGP', 'IHP']:
                            result += ' ' + str(comp[k + 'raw'][0])
                            result += ' ' + str(comp[k + 'raw'][1])
                        output.write(result + '\n')
