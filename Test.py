#!/usr/bin/env python

import HapMap
import Genotype

g1 = HapMap.read_phase('HapMap/hapmap_chr1_ceu_1000_1078_85180574_86259382.inp')
g2 = HapMap.read_phase('HapMap/hapmap_chr1_ceu_1000_1078_85180574_86259382.out.1.inp')

c = Genotype.compare(g1, g2)
print c
