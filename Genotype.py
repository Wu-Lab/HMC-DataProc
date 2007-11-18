#!/usr/bin/env python

# this package contains the common functions for handling genotype data

# trios_phase, phase the child's genotypes from trios

def trios_phase(child, father, mother):
    if child[0] == 'N':
        child = child[1] + child[0]
    if father[0] == 'N':
        father = father[1] + father[0]
    if mother[0] == 'N':
        mother = mother[1] + mother[0]
    inconsistent = False
    haplos = [child, 'NN']
    phased = [0, 0]
    if father != 'NN' and father[0] == father[1]:
        phased[0] = phased[1] = 1
        if child[0] == father[0]:
            if child[1] == 'N':
                if mother[0] == mother[1]:
                    haplos[0] = haplos[1] = father[0] + mother[0]
                else:
                    haplos[0] = haplos[1] = father[0] + 'N'
            elif child[1] == mother[0]:
                haplos[0] = father[0] + child[1]
                haplos[1] = father[0] + mother[1]
            elif child[1] == mother[1] or mother[1] == 'N':
                haplos[0] = father[0] + child[1]
                haplos[1] = father[0] + mother[0]
            else:
                inconsistent = True
        elif child[1] == father[0] or child[1] == 'N':
            if child[0] == 'N':
                if mother[0] == mother[1]:
                    haplos[0] = haplos[1] = father[0] + mother[0]
                else:
                    haplos[0] = haplos[1] = father[0] + 'N'
            elif child[0] == mother[0]:
                haplos[0] = father[0] + child[0]
                haplos[1] = father[0] + mother[1]
            elif child[0] == mother[1] or mother[1] == 'N':
                haplos[0] = father[0] + child[0]
                haplos[1] = father[0] + mother[0]
            else:
                inconsistent = True
        else:
            inconsistent = True
    elif mother != 'NN' and mother[0] == mother[1]:
        phased[0] = phased[1] = 1
        if child[0] == mother[0]:
            if child[1] == 'N':
                if father[0] == father[1]:
                    haplos[0] = haplos[1] = father[0] + mother[0]
                else:
                    haplos[0] = haplos[1] = 'N' + mother[0]
            elif child[1] == father[0]:
                haplos[0] = child[1] + mother[0]
                haplos[1] = father[1] + mother[0]
            elif child[1] == father[1] or father[1] == 'N':
                haplos[0] = child[1] + mother[0]
                haplos[1] = father[0] + mother[0]
            else:
                inconsistent = True
        elif child[1] == mother[0] or child[1] == 'N':
            if child[0] == 'N':
                if father[0] == father[1]:
                    haplos[0] = haplos[1] = father[0] + mother[0]
                else:
                    haplos[0] = haplos[1] = 'N' + mother[0]
            elif child[0] == father[0]:
                haplos[0] = child[0] + mother[0]
                haplos[1] = father[1] + mother[0]
            elif child[0] == father[1] or father[1] == 'N':
                haplos[0] = child[0] + mother[0]
                haplos[1] = father[0] + mother[0]
            else:
                inconsistent = True
        else:
            inconsistent = True
    elif child == 'NN':
        if father != 'NN' and father[0] == father[1]:
            haplos[0] = father[0]
            phased[0] = 1
        else:
            haplos[0] = 'N'
        if mother != 'NN' and mother[0] == mother[1]:
            haplos[0] += mother[0]
            phased[0] = 1
        else:
            haplos[0] += 'N'
        haplos[1] = haplos[0]
        phased[1] = phased[0]
    elif child[0] == child[1]:
        haplos[0] = child
        phased[0] = 1
        if father[0] == child[0]:
            haplos[1] = father[1]
        else:
            haplos[1] = father[0]
        if mother[0] == child[1]:
            haplos[1] += mother[1]
        else:
            haplos[1] += mother[0]
        if haplos[1] == 'NN':
            phased[1] = 0
        else:
            phased[1] = 1
    if inconsistent:
        haplos[0] = child
        haplos[1] = 'NN'
        phased[0] = phased[1] = 0
        raise RuntimeWarning, 'Inconsistent trios!' + ' ' + child + ' ' + father + ' ' + mother
    return haplos, phased

def test_trios_phase():
    genotypes = ['NN', 'AA', 'BB', 'AB', 'AN', 'BN']
    for child in genotypes:
        for father in genotypes:
            for mother in genotypes:
                haplos, phased = trios_phase(child, father, mother)
                print child, father, mother, ' ==> ', haplos[0], phased[0], haplos[1], phased[1]

class Genotype:
    def __init__(self, other = None):
        if other is None:
            self.haplos = [list(), list()]
            self.status = list()
        else:
            self.haplos = [list(other.haplos[0]), list(other.haplos[1])]
            self.status = list(other.status)
    def append(self, genos, status = None):
        if status is None or status == 1:
            status = [1 for i in range(0, len(genos[0]))]
        elif status == 0:
            status = [0 for i in range(0, len(genos[0]))]
        if len(genos[0]) != len(status) or len(genos[1]) != len(status):
            raise RuntimeError, 'Incorrect argument number!'
        self.haplos[0] += genos[0]
        self.haplos[1] += genos[1]
        self.status += status
    def len(self):
        return len(self.status)
    def geno(self, locus):
        return [self.haplos[0][locus], self.haplos[1][locus]]
    def geno_reverse(self, locus):
        return [self.haplos[1][locus], self.haplos[0][locus]]
    def reverse(self):
        g = Genotype(self)
        g.haplos.reverse()
        return g
    def getSwitchDistance(self, other, mask = None):
        if mask is None:
            mask = [1 for i in range(0, self.len())]
        switch_distance = [0, 0]
        phase_status = 0
        for i in range(0, self.len()):
            if self.status[i] == 1 and mask[i] == 1:
                if self.haplos[0][i] != 'N' and self.haplos[1][i] != 'N':
                    if self.haplos[0][i] != self.haplos[1][i]:
                        switch_distance[1] += 1
                        if phase_status == 1:
                            if self.geno(i) == other.geno(i):
                                continue
                            elif self.geno(i) == other.geno_reverse(i):
                                phase_status = 2
                                switch_distance[0] += 1
                            else:
                                raise RuntimeWarning, 'Inconsistent genotype at locus ' + str(i)
                        elif phase_status == 2:
                            if self.geno(i) == other.geno_reverse(i):
                                continue
                            elif self.geno(i) == other.geno(i):
                                phase_status = 1
                                switch_distance[0] += 1
                            else:
                                raise RuntimeWarning, 'Inconsistent genotype at locus ' + str(i)
                        elif phase_status == 0:
                            if self.geno(i) == other.geno(i):
                                phase_status = 1
                            elif self.geno(i) == other.geno_reverse(i):
                                phase_status = 2
                            else:
                                raise RuntimeWarning, 'Inconsistent genotype at locus ' + str(i)
                        else:
                            raise RuntimeError, 'Unknown phase_status value!'
        if switch_distance[1] >= 1:
            switch_distance[1] -= 1
        return switch_distance
    def getDiffNumber(self, other, mask = None):
        if mask is None:
            mask = [1 for i in range(0, self.len())]
        diff_number = [0, 0]
        for i in range(0, self.len()):
            if self.status[i] == 1 and mask[i] == 1:
                if self.haplos[0][i] != 'N' and self.haplos[1][i] != 'N':
                    if self.haplos[0][i] != self.haplos[1][i]:
                        diff_number[1] += 1
                        if self.geno(i) == other.geno(i):
                            continue
                        elif self.geno(i) == other.geno_reverse(i):
                            diff_number[0] += 1
                        else:
                            raise RuntimeWarning, 'Inconsistent genotype at locus ' + str(i)
        if diff_number[1] >= 1:
            diff_number[1] -= 1
        return diff_number
    def compare(self, other, mask = None):
        comparison = dict()
        comparison['SE'] = self.getSwitchDistance(other, mask)
        diff_num_f = self.getDiffNumber(other, mask)
        diff_num_r = self.getDiffNumber(other.reverse(), mask)
        if diff_num_f[1] != diff_num_r[1]:
            raise RuntimeWarning, 'Different possible error numbers!'
        if diff_num_f[0] < diff_num_r[0]:
            comparison['IGP'] = diff_num_f
        else:
            comparison['IGP'] = diff_num_r
        comparison['IHP'] = [int(comparison['IGP'][0]>0), int(comparison['IGP'][1]>0)]
        return comparison

def compare(genos_real, genos_infer, masks = None):
    if len(genos_real) != len(genos_infer):
        raise RuntimeError, 'Inconsistent genotype data!'
    if masks is None:
        masks = [None for i in range(0, len(genos_real))]
    comp_sum = dict()
    for k in ['SE', 'IGP', 'IHP']:
        comp_sum[k] = [0, 0]
    for i in range(0, len(genos_real)):
        comp = genos_real[i].compare(genos_infer[i], masks[i])
        for k in ['SE', 'IGP', 'IHP']:
            for j in [0, 1]:
                comp_sum[k][j] += comp[k][j]
    comparison = dict()
    for k in ['SE', 'IGP', 'IHP']:
        comparison[k + 'raw'] = comp_sum[k]
        if comp_sum[k][1] > 0:
            comparison[k] = float(comp_sum[k][0]) / comp_sum[k][1]
        else:
            comparison[k] = 0
    return comparison
