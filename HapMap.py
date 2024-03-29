#!/usr/bin/env python

# this package contains the common functions for handling HapMap data

import gzip
import random
import Genotype

# extract and infer gentoype phase from HapMap trios data
def trios_phase(genotypes, pedinfo, haplotypes):
    # parse pedigree information
    pedinfo_file = gzip.open(pedinfo, 'rb')
    samples = dict()
    for line in pedinfo_file:
        line = line.split()
        for i in (0, 1, 2, 3, 4):
            line[i] = int(line[i])
        for i in (5, 6):
            line[i] = line[i].split(':')[4]
        samples[(line[0], line[1])] = line
    trios = list()
    for child in samples.values():
        if (child[2] != 0) and (child[3] != 0):
            father = samples[(child[0], child[2])]
            mother = samples[(child[0], child[3])]
            trios.append([child[6], father[6], mother[6]])
    pedinfo_file.close()
    
    # parse genotypes data
    haplotype_name = haplotypes.split('/')
    haplotype_name = haplotype_name[len(haplotype_name)-1]
    haplotype_obj = open(haplotypes, 'wb')
    haplotype_file = gzip.GzipFile(haplotype_name, 'wb', fileobj = haplotype_obj)
    genotype_file = gzip.open(genotypes, 'rb')
    inconsistent_count = dict()
    for line in genotype_file:
        line = line.split()
        if line[0] == 'rs#':
            headline = line
            map_id = dict([(headline[i], i) for i in range(11, 101)])
            trios_id = [[map_id[child], map_id[father], map_id[mother]]
                for child, father, mother in trios]
            continue
        buffer = ''
        phased_num = 0
        for cid, fid, mid in trios_id:
            child, father, mother = line[cid], line[fid], line[mid]
            try:
                haplos, phased = Genotype.trios_phase(child, father, mother)
            except RuntimeWarning, w:
                #print w.message
                #print cid, fid, mid, headline[cid], headline[fid], headline[mid]
                haplos, phased = [child, 'NN'], [0, 0]
                if inconsistent_count.has_key((headline[cid], headline[fid], headline[mid])):
                    inconsistent_count[(headline[cid], headline[fid], headline[mid])] += 1
                else:
                    inconsistent_count[(headline[cid], headline[fid], headline[mid])] = 1
            buffer += ' ' + haplos[0] + ' ' + str(phased[0])
            buffer += ' ' + haplos[1] + ' ' + str(phased[1])
            phased_num += phased[0] + phased[1]
        buffer = line[2] + ' ' + line[3] + ' ' + str(phased_num) + ' ' + buffer
        haplotype_file.write(buffer + '\n')
    genotype_file.close()
    haplotype_file.close()
    haplotype_obj.close()

    for child, father, mother in trios:
        if inconsistent_count.has_key((child, father, mother)):
            print child, father, mother, inconsistent_count[(child, father, mother)]

# filter the unphased and missing haplotypes
def trios_filter(haplotypes, filtered, filter = [0.8, 0.8, 0.05]):
    for i in range(0, 3):
        if filter[i] > 1:
            filter[i] = 1
        elif filter[i] < 0:
            filter[i] = 0
    filtered_name = filtered.split('/')
    filtered_name = filtered_name[len(filtered_name)-1]
    filtered_obj = open(filtered, 'wb')
    filtered_file = gzip.GzipFile(filtered_name, 'wb', fileobj = filtered_obj)
    haplotype_file = gzip.open(haplotypes, 'rb')
    for line in haplotype_file:
        buffer = line
        line = line.split()
        for i in [1, 2] + range(4, len(line), 2):
            line[i] = int(line[i])
        num = (len(line) - 3) / 2
        phased = line[2]
        non_missing = num*2
        alleles = dict()
        for i in range(3, len(line), 2):
            for j in [0, 1]:
                if line[i][j] == 'N':
                    non_missing -= 1
                elif alleles.has_key(line[i][j]):
                    alleles[line[i][j]] += 1
                else:
                    alleles[line[i][j]] = 1
        threshold = [int(filter[0]*num), int(filter[1]*num*2), int(filter[2]*non_missing)]
        if phased >= threshold[0] and non_missing >= threshold[1] and \
                len(alleles) >= 2 and min(alleles.values()) >= threshold[2]:
            filtered_file.write(buffer)
    haplotype_file.close()
    filtered_file.close()
    filtered_obj.close()

# sort trios data by position
def trios_sort(haplotypes, sorted):
    sorted_name = sorted.split('/')
    sorted_name = sorted_name[len(sorted_name)-1]
    sorted_obj = open(sorted, 'wb')
    sorted_file = gzip.GzipFile(sorted_name, 'wb', fileobj = sorted_obj)
    haplotype_file = gzip.open(haplotypes, 'rb')
    lines = haplotype_file.readlines()
    positions = list()
    map_id = dict()
    for i in range(0, len(lines)):
        line = lines[i].split()
        line[1] = int(line[1])
        positions.append(line[1])
        map_id[line[1]] = i
    positions.sort()
    for i in range(0, len(lines)):
        sorted_file.write(lines[map_id[positions[i]]])
    haplotype_file.close()
    sorted_file.close()
    sorted_obj.close()

# parse HapMap phased genotypes into samples format
def parse_phased(genotypes, legend, haplotypes, chr):
    legend_file = gzip.open(legend, 'rb')
    snp = list()
    for line in legend_file:
        line = line.split()
        if line[0] == 'rs' or line[0][0:2] != 'rs':
            continue
        snp.append([line[2], line[3], line[1]])
    legend_file.close()
    genotypes_file = gzip.open(genotypes, 'rb')
    genos = genotypes_file.readlines()
    genotypes_file.close()
    for i in range(0, len(genos)):
        genos[i] = genos[i].split()
    geno_num = len(genos) / 2
    samples_file = gzip.open(haplotypes, 'wb')
    for i in range(0, len(snp)):
        line = chr + ' ' + snp[i][2] + ' ' + str(geno_num)
        for j in range(0, geno_num):
            alleles = [int(genos[2*j][i]), int(genos[2*j+1][i])]
            line += ' ' + snp[i][alleles[0]] + snp[i][alleles[1]] + ' 1'
        samples_file.write(line + '\n')
    samples_file.close()

# split trios data to small samples
def trios_samples(haplotypes, prefix, length):
    if length < 1:
        raise RuntimeError, 'Error length argument in trios_samples!'
    haplotype_file = gzip.open(haplotypes, 'rb')
    samples = list()
    count = 0
    for line in haplotype_file:
        samples.append(line)
        count += 1
        if count >= length:
            first_pos = int(samples[0].split()[1])
            last_pos = int(samples[len(samples)-1].split()[1])
            spacing = int((last_pos - first_pos) / len(samples))
            filename = prefix + '_' + str(len(samples)) + '_' + str(spacing) + '_' + \
                    str(first_pos) + '_' + str(last_pos) + '.txt'
            print filename
            sample_file = open(filename, 'w')
            sample_file.writelines(samples)
            sample_file.close()
            samples = list()
            count = 0
    haplotype_file.close()
    return len(samples)

# select samples
def select_samples(haplotypes, prefix, start, length):
    haplotype_file = gzip.open(haplotypes, 'rb')
    samples = list()
    line_num = 0
    for line in haplotype_file:
        if line_num >= (start + length):
            break
        elif line_num >= start:
            samples.append(line)
        line_num += 1
    if len(samples) <= 0:
        return -1
    first_pos = int(samples[0].split()[1])
    last_pos = int(samples[len(samples)-1].split()[1])
    spacing = int((last_pos - first_pos) / len(samples))
    filename = prefix + '_' + str(len(samples)) + '_' + str(spacing) + '_' + \
            str(first_pos) + '_' + str(last_pos) + '.txt'
    print filename
    sample_file = open(filename, 'w')
    sample_file.writelines(samples)
    sample_file.close()
    haplotype_file.close()
    if len(samples) == length:
        return (start + length)
    else:
        return -1

# split samples to smaller data with larger spacing
def split_samples(haplotypes, prefix, factor):
    if factor < 2:
        raise RuntimeError, 'Error factor argument in split_samples!'
    haplotype_file = open(haplotypes, 'r')
    samples = [[] for i in range(factor)]
    line_num = 0
    for line in haplotype_file:
        samples[line_num % factor].append(line)
        line_num += 1
    for i in range(factor):
        if len(samples[i]) <= 0:
            continue
        first_pos = int(samples[i][0].split()[1])
        last_pos = int(samples[i][len(samples[i])-1].split()[1])
        spacing = int((last_pos - first_pos) / len(samples[i]))
        filename = prefix + '_' + str(len(samples[i])) + '_' + str(spacing) + '_' + \
                str(first_pos) + '_' + str(last_pos) + '.txt'
        print filename
        sample_file = open(filename, 'w')
        sample_file.writelines(samples[i])
        sample_file.close()
    haplotype_file.close()

# filter samples by individual missing rate
def filter_samples(filename, filter = 0.1):
    snp = read_snp_info(filename)
    genos = read_hapmap(filename)
    max_missing = genos[0].len() * 2.0 * filter
    genos_filtered = list()
    for g in genos:
        if (g.haplos[0].count('N') + g.haplos[1].count('N')) <= max_missing:
            genos_filtered.append(g)
        else:
            for i in range(0, len(g.status)):
                if g.status[i] == 1:
                    snp[i][2] -= 1
    if len(genos) != len(genos_filtered):
        dest_file = open(filename, 'w')
        for i in range(0, len(snp)):
            line = snp[i][0] + ' ' + snp[i][1] + ' ' + str(snp[i][2]) + ' '
            for g in genos_filtered:
                line += ' ' + g.haplos[0][i] + g.haplos[1][i] + ' ' + str(g.status[i])
            dest_file.write(line + '\n')
        dest_file.close()
        return True
    else:
        return False

# convert samples file to PHASE format
def convert_format_to_phase(samples, output, randomize = False, mask = None):
    sample_file = open(samples, 'r')
    output_file = open(output, 'w')
    sample_num = snp_num = 0
    for line in sample_file:
        line = line.replace('N', '?').split()
        if line[0][0] != '#':
            snp_num += 1
            if sample_num == 0:
                sample_num = (len(line) - 3) / 2
                haplotypes = ['' for i in range(0, sample_num * 2)]
                phase_status = ['' for i in range(0, sample_num)]
                positions = 'P '
                allele_types = ''
            elif sample_num != (len(line) - 3) / 2:
                raise RuntimeWarning, 'Inconsistent sample number!' + ' ' + sample_num
            positions += line[1] + ' '
            alleles = dict()
            for i in range(0, sample_num):
                if randomize:
                    r = random.randint(0, 1)
                else:
                    r = 0
                if mask is None or mask.masks[i][snp_num-1] == 1:
                    haplotypes[2*i]   += line[2*i+3][r]
                    haplotypes[2*i+1] += line[2*i+3][1-r]
                else:
                    haplotypes[2*i]   += '?'
                    haplotypes[2*i+1] += '?'
                for j in [0, 1]:
                    a = line[2*i+3][j]
                    if a != '?':
                        if alleles.has_key(a):
                            alleles[a] += 1
                        else:
                            alleles[a] = 1
                phase_status[i] += line[2*i+4]
            if len(alleles) <= 2:
                allele_types += 'S'
            else:
                print samples
                print snp_num, alleles
                print line
                raise RuntimeError, 'Not bi-allelic SNP site!'
    output_file.write(str(sample_num) + '\n')
    output_file.write(str(snp_num) + '\n')
    output_file.write(positions + '\n')
    output_file.write(allele_types + '\n')
    for i in range(0, sample_num):
        output_file.write('#' + str(i+1) + '\n')
        output_file.write(' ' + haplotypes[2*i] + '\n')
        output_file.write(' ' + haplotypes[2*i+1] + '\n')
    output_file.close()
    sample_file.close()

# convert samples file to HPM2 format
def convert_format_to_hpm2(samples, output, randomize = False, mask = None):
    sample_file = open(samples, 'r')
    output_file = open(output, 'w')
    sample_num = snp_num = 0
    for line in sample_file:
        line = line.replace('N', '0')
        line = line.replace('A', '1')
        line = line.replace('C', '2')
        line = line.replace('G', '3')
        line = line.replace('T', '4')
        line = line.split()
        if line[0][0] != '#':
            snp_num += 1
            if sample_num == 0:
                sample_num = (len(line) - 3) / 2
                haplotypes = ['#' + str(int(i/2)) + ' ' for i in range(0, sample_num * 2)]
                allele_names = 'Id '
            elif sample_num != (len(line) - 3) / 2:
                raise RuntimeWarning, 'Inconsistent sample number!' + ' ' + sample_num
            allele_names += 'M' + str(snp_num) + ' '
            for i in range(0, sample_num):
                if randomize:
                    r = random.randint(0, 1)
                else:
                    r = 0
                if mask is None or mask.masks[i][snp_num-1]:
                    haplotypes[2*i]   += line[2*i+3][r] + ' '
                    haplotypes[2*i+1] += line[2*i+3][1-r] + ' '
                else:
                    haplotypes[2*i]   += '0' + ' '
                    haplotypes[2*i+1] += '0' + ' '
    output_file.write(allele_names + '\n')
    for i in range(0, sample_num):
        output_file.write(haplotypes[2*i] + '\n')
        output_file.write(haplotypes[2*i+1] + '\n')
    output_file.close()
    sample_file.close()

# read HapMap SNP information
def read_snp_info(filename):
    snp_info = list()
    source_file = open(filename, 'r')
    for line in source_file:
        line = line.split()
        if line[0][0] != '#':
            snp_info.append([line[0], line[1], int(line[2])])
    source_file.close()
    return snp_info

# read HapMap data to Genotype class
def read_hapmap(filename):
    genos = list()
    source_file = open(filename, 'r')
    sample_num = snp_num = 0
    for line in source_file:
        line = line.split()
        if line[0][0] != '#':
            snp_num += 1
            if sample_num == 0:
                sample_num = (len(line) - 3) / 2
                genos = [Genotype.Genotype() for i in range(0, sample_num)]
            elif sample_num != (len(line) - 3) / 2:
                raise RuntimeWarning, 'Inconsistent sample number!' + ' ' + sample_num
            for i in range(0, sample_num):
                genos[i].append([line[2*i+3][0], line[2*i+3][1]], int(line[2*i+4]))
    source_file.close()
    return genos

# read PHASE format data to Genotype class
def read_phase(filename):
    genos = list()
    source_file = open(filename, 'r')
    sample_num = int(source_file.readline())
    snp_num = int(source_file.readline())
    line = source_file.readline()
    if line.lstrip()[0] == 'P':
        line = source_file.readline()
    genos = [Genotype.Genotype() for i in range(0, sample_num)]
    lines = source_file.readlines()
    for i in range(0, sample_num):
        genos[i].haplos[0] = list(lines[3*i+1].strip().replace(' ', '').replace('?', 'N'))
        genos[i].haplos[1] = list(lines[3*i+2].strip().replace(' ', '').replace('?', 'N'))
        genos[i].status = [1 for j in range(0, snp_num)]
    source_file.close()
    return genos

# read HPM2 format data to Genotype class
def read_hpm2(filename):
    genos = list()
    source_file = open(filename, 'r')
    line = source_file.readline()
    line = line.split()
    if line[0] != 'Id':
        raise RuntimeError, 'Not HPM2 format!'
    snp_num = len(line) - 1
    lines = source_file.readlines()
    for i in range(0, len(lines)):
        lines[i] = lines[i].replace('1', 'A')
        lines[i] = lines[i].replace('2', 'C')
        lines[i] = lines[i].replace('3', 'G')
        lines[i] = lines[i].replace('4', 'T')
        lines[i] = lines[i].split()
        lines[i] = lines[i][1:len(lines[i])]
    sample_num = len(lines) / 2
    genos = [Genotype.Genotype() for i in range(0, sample_num)]
    for i in range(0, sample_num):
        genos[i].haplos[0] = lines[2*i]
        genos[i].haplos[1] = lines[2*i+1]
        genos[i].status = [1 for j in range(0, snp_num)]
    source_file.close()
    return genos
