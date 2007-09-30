#!/usr/bin/env python

import os
import time
import shutil
import glob
import HapMap

sample_dir = 'HapMap/'
sample_names = ['1000']

def get_HMC_output(source, dir):
    output = source + '.reconstructed'
    return output

def get_haplorec_output(source, dir):
    output = source[len(dir):len(source)] + '.reconstructed'
    return output

params = {'HMC':{}, 'haplorec':{}, 'PHASE':{}}

params['HMC']['enable'] = True
params['HMC']['command'] = 'HMC8.exe --nologo -a 0.5 -i 1'
params['HMC']['output'] = 'HMC_1'
params['HMC']['suffix'] = '.inp'
params['HMC']['temp'] = get_HMC_output

params['haplorec']['enable'] = True
params['haplorec']['command'] = 'java -Xmx1024m -Xms128m -jar HaploRec.jar -n 1'
params['haplorec']['output'] = 'haplorec_1'
params['haplorec']['suffix'] = '.hpm2'
params['haplorec']['temp'] = get_haplorec_output

params['PHASE']['enable'] = False
params['PHASE']['command'] = 'PHASE.exe'
params['PHASE']['output'] = 'PHASE_1'
params['PHASE']['suffix'] = '.inp'
params['PHASE']['temp'] = None

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
            for chr in chromosomes:
                for pop in populations:
                    work_dir = sample_dir + name + p['suffix'] + '/'
                    input_dir = work_dir + pop.upper() + '/' + 'chr' + chr + '/'
                    output_dir = input_dir + p['output'] + '/'
                    if os.access(output_dir, os.F_OK) == False:
                        os.makedirs(output_dir)
                    pattern = '*_chr' + chr + '_' + pop + '_*' + p['suffix']
                    for source in glob.glob(input_dir + pattern):
                        print source
                        basename = source[len(input_dir):len(source)-len(p['suffix'])]
                        msg_file = output_dir + basename + '.message'
                        out_file = output_dir + basename + '.out'
                        if os.access(msg_file, os.F_OK):
                            continue
                        if p['temp'] is None:
                            cmdline = p['command'] + ' ' + source + ' ' + out_file
                        else:
                            cmdline = p['command'] + ' ' + source
                        message = open(msg_file, 'w')
                        begin_time = time.clock()
                        try:
                            p_in, p_out, p_err = os.popen3(cmdline)
                            msg_err = p_err.readlines()
                            msg_out = p_out.readlines()
                            p_out.close()
                            p_err.close()
                            p_in.close()
                        finally:
                            end_time = time.clock()
                            message.write(str(end_time-begin_time) + '\n')
                            message.write('\nstdout:\n\n')
                            message.writelines(msg_out)
                            message.write('\nstderr:\n\n')
                            message.writelines(msg_err)
                            message.close()
                        if not(p['temp'] is None):
                            temp = p['temp'](source, input_dir)
                            if os.access(temp, os.F_OK):
                                shutil.move(temp, out_file)
