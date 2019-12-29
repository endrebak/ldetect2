#!/usr/bin/env python3
#cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True, 

import sys, math, gzip
import numpy as np
# cimport numpy as cnp
import pandas as pd

from time import time

from libc.math cimport exp, fabs, sqrt
from libc.stdint cimport int8_t
from libc.stdio cimport printf


# # calculate Wen/Stephens shrinkage LD estimate
# gmapfile = sys.argv[1] # genetic map
# indfile = sys.argv[2] #list of individuals
# # NE = 11418.0
# # CUTOFF = 1e-7
# outfile = sys.argv[5] # outfile file

# NE = float(sys.argv[3])
# CUTOFF = float(sys.argv[4])
cimport cython

from libc.stdint cimport int32_t, int64_t

from cython.operator cimport dereference as deref

from libcpp.map cimport map as cpp_map
from libcpp.vector cimport vector as cpp_vector
from libcpp.utility cimport pair as cpp_pair

ctypedef cpp_pair[int32_t, int32_t] loc_pair
ctypedef cpp_map[loc_pair, double] locs_to_val

cdef locs_to_val array_to_map(const	int32_t [::1] ipos, const int32_t [::1] jpos, const double [::1] vals):

    cdef:
        int i, length = len(vals)
        loc_pair locs
        locs_to_val loc_map


    for i in range(length):
        locs = loc_pair(ipos[i], jpos[i])
        loc_map[locs] = vals[i]

    return loc_map



@cython.cdivision(True)
cdef inline double corrcoeff(locs_to_val loc_map, int32_t x, int32_t y):

    cdef:
        loc_pair locs
        double val, double_x, double_y
        cpp_map[loc_pair, double].iterator it

    locs = loc_pair(x, y)

    print("  ", x, y)

    it = loc_map.find(locs)
    if it != loc_map.end(): # found
        val = deref(it).second
        print("  val", val)
        locs = loc_pair(x, x)
        double_x = deref(loc_map.find(locs)).second
        print("  double_x", double_x)
        locs = loc_pair(y, y)
        double_y = deref(loc_map.find(locs)).second
        print("  double_y", double_y)

        result = (val / sqrt(double_x * double_y)) ** 2
        print("  ", result)
        return result
    else:
        return 0.0



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
cpdef mat2vec(const int32_t [::1] ipos, const int32_t [::1] jpos, const double [::1] vals, partitions):

    # should only take in ipos >= partitions.first_start
    # should only take in ipos <= partitions.last_start
    ipos_unique_arr = np.unique(np.array(ipos))
    cdef int32_t [::1] ipos_unique
    ipos_unique = ipos_unique_arr

    cdef:
        int32_t pstart, pend, curr_locus, curr_locus_index, end_locus, x, y, delta, length = len(ipos_unique)
        double val, double_x, double_y, corr
        loc_pair locs
        locs_to_val loc_map = array_to_map(ipos, jpos, vals)

    curr_locus_index = 0

    # out_arr = np.zeros(length)
    # # cdef double [::1] out
    # # out = out_arr

    # # create vector to hold corr_coeff in

    for p in partitions.itertuples(index=False):
        pstart, pend = p[0], p[1]
        end_locus = pend # int((pstart + pend) / 2)
        curr_locus = ipos_unique[curr_locus_index]
        x = curr_locus
        y = curr_locus
        while curr_locus <= end_locus:

            print("----" * 5)
            print("curr_locus", curr_locus)
            corr = 0
            delta = 0

            curr_locus = curr_locus
            print("pstart", pstart)
            print("pend", pend)

            print("x", x)
            print("y", y)

            while x >= pstart and y <= pend:

                print("x", x, "y", y)
                corr += corrcoeff(loc_map, x, y)

                if delta != 0:

                    x = ipos_unique[curr_locus_index - delta + 1]
                    corr += corrcoeff(loc_map, x, y)

                delta += 1

                if ((curr_locus_index - delta) >= 0) and ((curr_locus_index + delta) < length):
                    print("    are in if" * 5)
                    x = ipos_unique[(curr_locus_index - delta)]
                    y = ipos_unique[(curr_locus_index + delta)]
                else:
                    break

            # printf(r"%d\t%f\n", curr_locus, corr)
            print(curr_locus, corr, sep="\t")

            if curr_locus_index + 1 < length:
                curr_locus_index += 1
                curr_locus = ipos_unique[curr_locus_index]
            else:
                break

