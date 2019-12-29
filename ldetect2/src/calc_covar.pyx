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
cpdef calc_covar(haplos, gmapfile, indfile, double NE, double CUTOFF):

    inds = pd.read_table(indfile, header=None, squeeze=True).to_list()

    haps= list()
    nind_int = len(inds)
    s = 0

    for _i in range(1, 2*nind_int):
            s = s+ 1.0/float(_i)

    s = 1/s
    #print "s", s


    cdef:
        int i, j, k, pos1, pos2, len_allpos, toofar, n11, n10, n01, len_g1_int
        double gpos1, gpos2, ee_const, len_g1, ee, df, theta, f11, f1, f2, b, Ds, Ds2, D, nind, thetas, thetas2

    nind = float(nind_int)
    pos2gpos = pd.read_table(gmapfile, index_col=0, header=None, sep=" ", squeeze=True)
    pos2gpos = pos2gpos

    allpos = haplos.pop(0) #.tolist()
    allrs = haplos.pop(1) #.tolist()
    haps = haplos.values

    # for i, (k, v) in enumerate(pos2gpos.items()):
    #     print(k, v)
    #     if i == 5:
    #         break

    pos2gpos = pos2gpos[allpos]
    # print(pos2gpos)
    # raise

    records = []

    len_g1_int = haps.shape[1]
    len_g1 = float(haps.shape[1])
    ee_const = NE*4.0 / (2.0*nind)
    len_allpos = len(allpos)
    theta = s/(2.0*float(nind)+s)
    thetas = (1-theta)*(1-theta)
    thetas2 = (theta/2.0)*(1-theta/2.0)

    cdef double[:] pos2gpos_view
    cdef long[:] allpos_view
    cdef long[:] g1_view, g2_view
    cdef long[:, :] haps_view

    pos2gpos_view = pos2gpos.values
    allpos_view = allpos.values
    haps_view = haps

    for i in range(len_allpos):

            pos1 = allpos_view[i]

            gpos1 = pos2gpos_view[i]

            toofar = 0
            j = i

            while j < len_allpos and not toofar:
                    pos2 = allpos_view[j]
                    gpos2 = pos2gpos_view[j]

                    df = gpos2 - gpos1
                    ee = exp( - df * ee_const)

                    if ee < CUTOFF:
                            toofar = 1
                            j = j+1
                            continue


                    n11, n01, n10 = 0, 0, 0
                    for k in range(len_g1_int):
                        if haps_view[i][k] == 1 and haps_view[j][k] == 1:
                            n11 += 1
                        elif haps_view[i][k] == 0 and haps_view[j][k] == 1:
                            n01 += 1
                        elif haps_view[i][k] == 1 and haps_view[j][k] == 0:
                            n10 += 1

                    f11 = n11/len_g1
                    f1 = (n11+n10)/len_g1
                    f2 = (n11+n01)/len_g1
                    D = f11 - f1*f2
                    Ds = D*ee
                    Ds2 = thetas*Ds

                    if fabs(Ds2) < CUTOFF:
                            j = j+1
                            continue
                    if i == j:
                            Ds2 = Ds2 + thetas2

                    print(allrs[i], allrs[j], pos1, pos2, gpos1, gpos2, D, Ds2)
                    # records.append(result)
                    j = j+1

    # outdf = pd.DataFrame.from_records(records)
    # return outdf
