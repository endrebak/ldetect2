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
    if val:
        double_x = loc_map.get((x, x), 0)
        double_y = loc_map.get((y, y), 0)

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
                corr += corrcoeff(loc_map, x, y)

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
            # printf("%d\t%f\t%d\n", curr_locus, corr, delta)

            # end = time()
            # print(curr_locus, corr, end - start, delta, sep="\t")

            if curr_locus_index + 1 < length:
                curr_locus_index += 1
                curr_locus = ipos_unique[curr_locus_index]
            else:
                break


# regular map
# 16050607        1.0     0.07001876831054688
# 16050840        1.0001658474807615      0.3061106204986572
# 16050847        0.004952597023447374    0.47174859046936035
# 16050958        0.00916724354529276     0.6771810054779053
# 16051249        0.3351658815679538      0.886406421661377
# 16051453        0.10693964554132061     1.0469534397125244
# 16051477        0.12357279221171176     1.2537975311279297
# 16051796        0.07130721301728887     1.423884391784668
# 16052080        1.349161739654159       1.6078472137451172
# 16052167        0.1267288147088067      1.7935221195220947
# 16052684        0.2210650949215177      2.010793924331665
# 16052730        0.7016378081744832      2.2011640071868896
# 16052962        1.1169942616863808      2.4065589904785156
# 16052986        1.076809744500464       2.7378647327423096
# 16053139        1.5727499733173513      2.786099672317505
# 16053249        2.505072074956979       2.987943410873413
# 16053254        1.631801123304757       3.204648494720459
# 16053444        2.854333791321056       3.4111592769622803
# 16053659        2.724271088576541       3.6162593364715576
# 16053730        0.8224091678720482      3.835418939590454
# 16053791        0.2107946981785794      4.073134660720825
# 16053862        2.3781599349078593      4.1966118812561035
# 16053863        1.397634333762032       4.828949689865112
# 16054311        1.8968998793773642      4.723801374435425
# 16054454        0.4806830149148064      5.048082590103149
# 16054713        2.1966416739946957      5.02083683013916
# 16054740        0.253149499719289       5.283098220825195
# 16054751        3.1898460869692546      5.4561402797698975
# 16055070        0.3847804318314368      5.766130208969116
# 16055105        0.3392071939444714      5.921470642089844
# 16055122        0.45308988956812457     6.096620321273804
# 16055171        1.4554468986153046      6.338205814361572
# 16055294        0.4393503333691098      6.502943992614746
