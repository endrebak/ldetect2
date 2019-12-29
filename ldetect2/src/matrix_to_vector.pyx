#!/usr/bin/env python3
#cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True, 

from time import time

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

# from cython.operator cimport dereference as deref

# from libcpp.map cimport map as cpp_map
# from libcpp.vector cimport vector as cpp_vector
# from libcpp.utility cimport pair as cpp_pair

# ctypedef cpp_pair[int32_t, int32_t] loc_pair
# ctypedef cpp_map[loc_pair, double] locs_to_val

# cdef locs_to_val array_to_map(const	int32_t [::1] ipos, const int32_t [::1] jpos, const double [::1] vals):

#     cdef:
#         int i, length = len(vals)
#         loc_pair locs
#         locs_to_val loc_map


#     for i in range(length):
#         locs = loc_pair(ipos[i], jpos[i])
#         loc_map[locs] = vals[i]

#     return loc_map



@cython.cdivision(True)
cdef inline double corrcoeff(loc_map, int32_t x, int32_t y):

    cdef:
        double val, double_x, double_y, result

    # print("x, y", x, y)
    val = loc_map.get((x, y), 0)
    double_x = loc_map.get((x, x), 0)
    double_y = loc_map.get((y, y), 0)
    if val and double_x and double_y:

        result = (val / sqrt(double_x * double_y)) ** 2
        # print("adding corrcoeff", result)
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
        int32_t pstart, pend, curr_locus, curr_locus_index, end_locus, x, y, delta, length = len(ipos_unique), i
        double val, double_x, double_y, corr

    loc_map = {}
    for i in range(len(ipos)):
        loc_map[ipos[i], jpos[i]] = vals[i]

    curr_locus_index = 0

    for p in partitions.itertuples(index=False):
        pstart, pend = p[0], p[1]
        end_locus = pend # int((pstart + pend) / 2)
        curr_locus = ipos_unique[curr_locus_index]
        while curr_locus <= end_locus:
            x = curr_locus
            y = curr_locus

            # start = time()
            # print("----" * 5)
            # print("curr_locus", curr_locus)
            corr = 0
            delta = 0

            curr_locus = curr_locus
            # print("pstart", pstart)
            # print("pend", pend)

            # print("x", x)
            # print("y", y)

            while x >= pstart and y <= pend:

                # print("x", x, "y", y)
                if (x, y) in loc_map:
                    # print("here", file=sys.stderr)
                    corr += corrcoeff(loc_map, x, y)
                    # print("there", file=sys.stderr)

                if delta != 0:

                    x = ipos_unique[curr_locus_index - delta + 1]
                    corr += corrcoeff(loc_map, x, y)

                delta += 1

                if ((curr_locus_index - delta) >= 0) and ((curr_locus_index + delta) < length):
                    # print("    are in if" * 5)
                    x = ipos_unique[(curr_locus_index - delta)]
                    y = ipos_unique[(curr_locus_index + delta)]
                else:
                    break

            printf("%d\t%f\n", curr_locus, corr)

            if curr_locus_index + 1 < length:
                curr_locus_index += 1
                curr_locus = ipos_unique[curr_locus_index]
            else:
                break
