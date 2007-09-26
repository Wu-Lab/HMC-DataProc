#!/usr/bin/env python

import Genotype
import HapMap

HapMap.trios_phase('CEU/genotypes_chr1_ceu_r22_nr.b36_fwd.txt.gz',
                   'pedinfo/pedinfo2sample_CEU.txt.gz',
                   'test.chr1.txt.gz')
