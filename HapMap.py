#!/usr/bin/env python

# this package contains the common functions for handling HapMap data

import gzip
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
    
    # parse genotypes data
    genotype_file = gzip.open(genotypes, 'rb')
    haplotype_file = gzip.open(haplotypes, 'wb')
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
            haplos = [child, 'NN']
            phased = [0, 0]
            try:
                haplos, phased = Genotype.trios_phase(child, father, mother)
            except RuntimeWarning, w:
                #print w.message
                #print cid, fid, mid, headline[cid], headline[fid], headline[mid]
                if inconsistent_count.has_key((headline[cid], headline[fid], headline[mid])):
                    inconsistent_count[(headline[cid], headline[fid], headline[mid])] += 1
                else:
                    inconsistent_count[(headline[cid], headline[fid], headline[mid])] = 1
            buffer += ' ' + haplos[0] + ' ' + str(phased[0])
            buffer += ' ' + haplos[1] + ' ' + str(phased[1])
            phased_num += phased[0] + phased[1]
        buffer = line[2] + ' ' + line[3] + ' ' + str(phased_num) + ' ' + buffer
        haplotype_file.write(buffer + '\n')

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
    haplotype_file = gzip.open(haplotypes, 'rb')
    filtered_file = gzip.open(filtered, 'wb')
    for line in haplotype_file:
        buffer = line
        line = line.split()
        for i in [1, 2] + range(4, len(line), 2):
            line[i] = int(line[i])
        alleles = dict()
        missing = 0
        for i in range(3, len(line), 2):
            for j in [0, 1]:
                if line[i][j] == 'N':
                    missing += 1
                elif alleles.has_key(line[i][j]):
                    alleles[line[i][j]] += 1
                else:
                    alleles[line[i][j]] = 1
        num = (len(line) - 3) / 2
        threshold = [int(filter[0]*num), int(filter[1]*num*2), int(filter[2]*num*2)]
        if line[2] >= threshold[0] and missing >= filter[1] and \
                len(alleles) >= 2 and min(alleles.values()) >= filter[2]:
            filtered_file.write(buffer)

# sort trios data by position
def trios_sort(haplotypes, sorted):
    haplotype_file = gzip.open(haplotypes, 'rb')
    sorted_file = gzip.open(sorted, 'wb')
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

# select samples
def select_samples(haplotypes, prefix, start, length):
    haplotype_file = gzip.open(haplotypes, 'rb')
    samples = list()
    line_num = -1
    for line in haplotype_file:
        line_num += 1
        if line_num < start:
            continue
        elif line_num >= (start + length):
            break
        samples.append(line)
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
    if len(samples) == length:
        return (start + length)
    else:
        return -1

# convert samples file to PHASE format
def convert_format_to_phase(samples, output):
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
            allele_types += 'S'
            for i in range(0, sample_num):
                haplotypes[2*i] += line[2*i+3][0]
                haplotypes[2*i+1] += line[2*i+3][1]
                phase_status[i] += line[2*i+4]
    output_file.write(str(sample_num) + '\n')
    output_file.write(str(snp_num) + '\n')
    output_file.write(positions + '\n')
    output_file.write(allele_types + '\n')
    for i in range(0, sample_num):
        output_file.write('#' + str(i+1) + '\n')
        output_file.write(' ' + haplotypes[2*i] + '\n')
        output_file.write(' ' + haplotypes[2*i+1] + '\n')

# convert samples file to HPM2 format
def convert_format_to_hpm2(samples, output):
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
                haplotypes[2*i] += line[2*i+3][0] + ' '
                haplotypes[2*i+1] += line[2*i+3][1] + ' '
    output_file.write(allele_names + '\n')
    for i in range(0, sample_num):
        output_file.write(haplotypes[2*i] + '\n')
        output_file.write(haplotypes[2*i+1] + '\n')

#read HapMap data to Genotype class
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
    return genos

#read PHASE format data to Genotype class
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
    return genos

#read HPM2 format data to Genotype class
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
    return genos
