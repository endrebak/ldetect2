#!/usr/bin/env python3
#cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True, 

import sys, math, gzip
import numpy as np
# cimport numpy as cnp
import pandas as pd

from time import time

from libc.math cimport exp, fabs
from libc.stdint cimport int8_t

# # calculate Wen/Stephens shrinkage LD estimate
# gmapfile = sys.argv[1] # genetic map
# indfile = sys.argv[2] #list of individuals
# # NE = 11418.0
# # CUTOFF = 1e-7
# outfile = sys.argv[5] # outfile file

# NE = float(sys.argv[3])
# CUTOFF = float(sys.argv[4])
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
cpdef mat2vec(const int64_t [::1] ipos, const int64_t [::1] jpos, const double [::1] vals, partitions):

    # note that partitions need to be precomputed with new ends halfway between i.end and i+1.start
    # last end is smallest < partitions.last_end
    # but also keep last end, since it is used
    # start, end, newend

    # should only take in ipos >= partitions.first_start

    cdef:
       int pstart, pend, curr_locus, curr_locus_index, end_locus, x, y
       double val

    curr_locus_index = 0
    curr_locus = ipos[curr_locus_index]

    # create vector to hold corr_coeff in

    for p in partitions.itertuples(index=False):
        pstart, pend = p[0], p[1]
        end_locus = int((pstart + pend) / 2)
        x = curr_locus
        y = curr_locus
        while curr_locus <= end_locus:

            curr_locus = curr_locus

            while x >= pstart and y <= pend:

                pass
                # check that x and y have a value - how?
                


        curr_locus += 1
