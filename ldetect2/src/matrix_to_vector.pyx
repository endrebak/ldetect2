#!/usr/bin/env python3

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
#         loc_map[locs] = vals[return]

#     i loc_map



@cython.cdivision(True)
cdef inline double corrcoeff(loc_map, int32_t x, int32_t y, double theta2):

    cdef:
        double val, double_x, double_y, result

    # print("x, y", x, y)
    val = loc_map.get((x, y), 0)
    double_x = loc_map.get((x, x), theta2)
    double_y = loc_map.get((y, y), theta2)
    if val and double_x and double_y:

        result = (val / sqrt(double_x * double_y)) ** 2
        return result
    else:
        return 0.0


def fill_loc_map(df, loc_map):

    print("Inserting keys", file=sys.stderr)
    s = df.copy().set_index(["i", "j"]).val
    print("Converting to dict", file=sys.stderr)
    # print("filling loc_map", s, file=sys.stderr)
    d = s.to_dict()
    print("Updating dict", file=sys.stderr)
    loc_map.update(d)
    print("Done inserting keys", file=sys.stderr)


def remove_from_loc_map(delete_before, loc_map):

    print("Removing keys", file=sys.stderr)
    keys_not_needed = []
    for (i, j) in loc_map.keys():
        if i < delete_before:
            keys_not_needed.append((i, j)) 

    for k in keys_not_needed:
        del loc_map[k]

    print("Done removing keys", file=sys.stderr)



def read_covariances(f):
    print("Reading file", file=sys.stderr)
    df = pd.read_parquet(f)
    print("Done reading file", file=sys.stderr)
    return df

# def update_ipos(ipos, ipos_add):
#     return np.concatenate([ipos, ipos_add])

@cython.boundscheck(True)
@cython.wraparound(False)
@cython.initializedcheck(True)
cpdef mat2vec(covariances, partitions, double theta2):


    # should only take in ipos >= partitions.first_start
    # should only take in ipos <= partitions.last_start

    cdef:
        int32_t pstart, pend, curr_locus, curr_locus_index, end_locus, x, y, delta, p_num, delete_locus
        double val, double_x, double_y, corr
        int32_t [::1] ipos
        int32_t [::1] jpos
        double [::1] vals

    curr_locus_index = 0

    """When finding end_locus use (pend[k] + pstart[k+1])/2

# end_locus = int((self.partitions[p_num][1] + self.partitions[p_num+1][0]) / 2) # diag - specific

When checking whether x and y within range, use regular partition ends.

# while x >= self.partitions[p_num][0] and y <= self.partitions[p_num][1]:"""

    """
snp_first and snp_last is just the first and last entry of the partitions/intervals
    """

    cdef int32_t [::1] ipos_unique

    loc_map = {}
    df = read_covariances(covariances[0])
    fill_loc_map(df, loc_map)
    ipos_unique_arr = np.unique(df.i.values)
    ipos_unique = ipos_unique_arr
    # n_to_remove = 0

    for p_num, p in enumerate(partitions.itertuples(index=False), 0):

        # print("-----" * 5, file=sys.stderr)

        if p_num + 1 < len(partitions):
            df = read_covariances(covariances[p_num + 1])
            fill_loc_map(df, loc_map)
            ipos_unique_arr = np.concatenate([ ipos_unique_arr, np.unique(df.i.values) ])
            ipos_unique = ipos_unique_arr

        # curr_locus_index -= n_to_remove
        length = len(ipos_unique)

        print("At partiton {}, so we are {}% done".format(str(p), p_num/len(partitions)), file=sys.stderr)
        pstart, pend, end_locus = p[0], p[1], p[2]
        # print("At start curr_locus_index is", curr_locus_index, file=sys.stderr)
        curr_locus = ipos_unique[curr_locus_index]
        # print("At start curr_locus is", curr_locus, file=sys.stderr)
        while curr_locus <= end_locus:
            x = curr_locus
            y = curr_locus

            # start = time()
            # print("----" * 5)
            # print("curr_locus", curr_locus)
            corr = 0
            delta = 0

            # print("pstart", pstart)
            # print("pend", pend)

            # print("x", x)
            # print("y", y)

            while x >= pstart and y <= pend:

                # print("x", x, "y", y)
                if (x, y) in loc_map:
                    # print("here", file=sys.stderr)
                    corr += corrcoeff(loc_map, x, y, theta2)
                    # print("there", file=sys.stderr)
                if delta != 0:

                    x = ipos_unique[curr_locus_index - delta + 1]
                    corr += corrcoeff(loc_map, x, y, theta2)

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


        if p_num + 1 < len(partitions):
            delete_locus = partitions.iloc[p_num + 1][0]
        else:
            delete_locus = end_locus

        remove_from_loc_map(delete_locus, loc_map)
        # to_remove = ipos_unique_arr < end_locus
        # n_to_remove = (ipos_unique_arr < end_locus).sum()

        # print("At end curr_locus_index is", curr_locus_index, file=sys.stderr)
        # curr_locus_index -= n_to_remove
        # print("After removing {} elements it is".format(n_to_remove), curr_locus_index, file=sys.stderr)
        # ipos_unique_arr = ipos_unique_arr[~to_remove]
        # ipos_unique = ipos_unique_arr
        # curr_locus_indexes = np.where(ipos_unique_arr == curr_locus)[0]
        # assert len(curr_locus_indexes) == 1
        # curr_locus_index = [0] + 1
        # print("At end curr_locus is", curr_locus, file=sys.stderr)
