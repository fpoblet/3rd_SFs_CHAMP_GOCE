# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 16:48:48 2023

1d plots for the combined CHAMP and GOCE observations.

@author: ma042
"""



import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
from scipy.optimize import curve_fit


################################################################################
### FUNCTIONS
################################################################################
def enstrophy(x, u0):
        return (u0/8)*x**(3)
    
def enstrophy_linear(log10_x, log10_u0):
        return 3*log10_x + log10_u0 - np.log10(4)


################################################################################
### PARAMETERS
################################################################################

#### Altitude ranges
# hmin = 220
# hmax = 380
# step = 1
# hvalues = np.arange(hmin, hmax, step)
fix_coef_s2 = 0.1


#### s ranges
s_min = 50
s_max = 1e4
num_bins = 50

## FontSize for every label
fs = 12

ylimm_nc = [0.01, 1e7]

xlimm = [1, 3e4]
ylimm = [1e-2, 1e9]
msize = 6

################################################################################
### CHAMP
################################################################################
path = 'C:/Poblet/IAP/JAPAN2023/research_topic/codes/read_outs/'

t0_t1_lon_lt_hmin_hmaxN = np.load(path+'t0_t1_lon_lt_hmin_hmax_u_lat-0-60N_07080910.npy', allow_pickle=True)

t0N = t0_t1_lon_lt_hmin_hmaxN[0]
t1N = t0_t1_lon_lt_hmin_hmaxN[1]
lonmN = t0_t1_lon_lt_hmin_hmaxN[2]
ltmN = t0_t1_lon_lt_hmin_hmaxN[3]
hminsN = t0_t1_lon_lt_hmin_hmaxN[4]
hmaxsN = t0_t1_lon_lt_hmin_hmaxN[5]

N_allN = np.load(path+'N_u_lat-0-60N_07080910.npy', allow_pickle=True)
s_allN = np.load(path+'s_u_lat-0-60N_07080910.npy', allow_pickle=True)
D_q_allN = np.load(path+'D_3_u_lat-0-60N_LT_07080910.npy', allow_pickle=True)
D_q_allN_median = np.load(path+'D_3_umedian_lat-0-60N_LT_07080910.npy', allow_pickle=True)

t0_t1_lon_lt_hmin_hmaxNs = np.load(path+'t0_t1_lon_lt_hmin_hmax_u_lat-30-60N_07080910.npy', allow_pickle=True)

t0Ns = t0_t1_lon_lt_hmin_hmaxNs[0]
t1Ns = t0_t1_lon_lt_hmin_hmaxNs[1]
lonmNs = t0_t1_lon_lt_hmin_hmaxNs[2]
ltmNs = t0_t1_lon_lt_hmin_hmaxNs[3]
hminsNs = t0_t1_lon_lt_hmin_hmaxNs[4]
hmaxsNs = t0_t1_lon_lt_hmin_hmaxNs[5]

N_allNs = np.load(path+'N_u_lat-30-60N_07080910.npy', allow_pickle=True)
s_allNs = np.load(path+'s_u_lat-30-60N_07080910.npy', allow_pickle=True)
D_q_allNs = np.load(path+'D_3_u_lat-30-60N_LT_07080910.npy', allow_pickle=True)
D_q_allNs_median = np.load(path+'D_3_umedian_lat-30-60N_LT_07080910.npy', allow_pickle=True)
# 

t0_t1_lon_lt_hmin_hmaxNseq = np.load(path+'t0_t1_lon_lt_hmin_hmax_u_lat-0-30N_07080910.npy', allow_pickle=True)

t0Nseq = t0_t1_lon_lt_hmin_hmaxNseq[0]
t1Nseq = t0_t1_lon_lt_hmin_hmaxNseq[1]
lonmNseq = t0_t1_lon_lt_hmin_hmaxNseq[2]
ltmNseq = t0_t1_lon_lt_hmin_hmaxNseq[3]
hminsNseq = t0_t1_lon_lt_hmin_hmaxNseq[4]
hmaxsNseq = t0_t1_lon_lt_hmin_hmaxNseq[5]

N_allNseq = np.load(path+'N_u_lat-0-30N_07080910.npy', allow_pickle=True)
s_allNseq = np.load(path+'s_u_lat-0-30N_07080910.npy', allow_pickle=True)
D_q_allNseq = np.load(path+'D_3_u_lat-0-30N_LT_07080910_cyc.npy', allow_pickle=True)
D_q_allNseq_median = np.load(path+'D_3_umedian_lat-0-30N_LT_07080910_cyc.npy', allow_pickle=True)
# 


t0_t1_lon_lt_hmin_hmaxS = np.load(path+'t0_t1_lon_lt_hmin_hmax_u_lat-0-60S_07080910.npy', allow_pickle=True)

t0S = t0_t1_lon_lt_hmin_hmaxS[0]
t1S = t0_t1_lon_lt_hmin_hmaxS[1]
lonmS = t0_t1_lon_lt_hmin_hmaxS[2]
ltmS = t0_t1_lon_lt_hmin_hmaxS[3]
hminsS = t0_t1_lon_lt_hmin_hmaxS[4]
hmaxsS = t0_t1_lon_lt_hmin_hmaxS[5]

N_allS = np.load(path+'N_u_lat-0-60S_07080910.npy', allow_pickle=True)
s_allS = np.load(path+'s_u_lat-0-60S_07080910.npy', allow_pickle=True)
D_q_allS = np.load(path+'D_3_u_lat-0-60S_LT_07080910.npy', allow_pickle=True)
D_q_allS_median = np.load(path+'D_3_umedian_lat-0-60S_LT_07080910.npy', allow_pickle=True)

t0_t1_lon_lt_hmin_hmaxSs = np.load(path+'t0_t1_lon_lt_hmin_hmax_u_lat-30-60S_07080910.npy', allow_pickle=True)

t0Ss = t0_t1_lon_lt_hmin_hmaxSs[0]
t1Ss = t0_t1_lon_lt_hmin_hmaxSs[1]
lonmSs = t0_t1_lon_lt_hmin_hmaxSs[2]
ltmSs = t0_t1_lon_lt_hmin_hmaxSs[3]
hminsSs = t0_t1_lon_lt_hmin_hmaxSs[4]
hmaxsSs = t0_t1_lon_lt_hmin_hmaxSs[5]

N_allSs = np.load(path+'N_u_lat-30-60S_07080910.npy', allow_pickle=True)
s_allSs = np.load(path+'s_u_lat-30-60S_07080910.npy', allow_pickle=True)
D_q_allSs = np.load(path+'D_3_u_lat-30-60S_LT_07080910.npy', allow_pickle=True)
# D_q_allSs = np.load(path+'D_3_u_lat-30-60S_LT_07080910_cyc.npy', allow_pickle=True)
# D_q_allSs = np.load(path+'D_3_u_lat-30-60S_LT_07080910_HWM07.npy', allow_pickle=True)
D_q_allSs_median = np.load(path+'D_3_umedian_lat-30-60S_LT_07080910.npy', allow_pickle=True)


t0_t1_lon_lt_hmin_hmaxSseq = np.load(path+'t0_t1_lon_lt_hmin_hmax_u_lat-0-30S_07080910.npy', allow_pickle=True)

t0Sseq = t0_t1_lon_lt_hmin_hmaxSseq[0]
t1Sseq = t0_t1_lon_lt_hmin_hmaxSseq[1]
lonmSseq = t0_t1_lon_lt_hmin_hmaxSseq[2]
ltmSseq = t0_t1_lon_lt_hmin_hmaxSseq[3]
hminsSseq = t0_t1_lon_lt_hmin_hmaxSseq[4]
hmaxsSseq = t0_t1_lon_lt_hmin_hmaxSseq[5]

N_allSseq = np.load(path+'N_u_lat-0-30S_07080910.npy', allow_pickle=True)
s_allSseq = np.load(path+'s_u_lat-0-30S_07080910.npy', allow_pickle=True)
D_q_allSseq = np.load(path+'D_3_u_lat-0-30S_LT_07080910_cyc.npy', allow_pickle=True)
D_q_allSseq_median = np.load(path+'D_3_umedian_lat-0-30S_LT_07080910_cyc.npy', allow_pickle=True)




### USE only 2007-2009 ############################################################
#####
t0N_datetime = pd.to_datetime(t0N, unit='s')
t0N_year = t0N_datetime.year

iyN = np.where((t0N_year == 2007) | (t0N_year == 2008) | (t0N_year == 2009))

t0N = t0N[iyN]
t1N = t1N[iyN]
lonmN = lonmN[iyN]
ltmN = ltmN[iyN]
hminsN = hminsN[iyN]
hmaxsN = hmaxsN[iyN]

s_allN = s_allN[iyN]
N_allN = N_allN[iyN]


#####
t0Ns_datetime = pd.to_datetime(t0Ns, unit='s')
t0Ns_year = t0Ns_datetime.year

iyNs = np.where((t0Ns_year == 2007) | (t0Ns_year == 2008) | (t0Ns_year == 2009))

t0Ns = t0Ns[iyNs]
t1Ns = t1Ns[iyNs]
lonmNs = lonmNs[iyNs]
ltmNs = ltmNs[iyNs]
hminsNs = hminsNs[iyNs]
hmaxsNs = hmaxsNs[iyNs]

s_allNs = s_allNs[iyNs]
N_allNs = N_allNs[iyNs]

#####
t0Nseq_datetime = pd.to_datetime(t0Nseq, unit='s')
t0Nseq_year = t0Nseq_datetime.year

iyNseq = np.where((t0Nseq_year == 2007) | (t0Nseq_year == 2008) | (t0Nseq_year == 2009))

t0Nseq = t0Nseq[iyNseq]
t1Nseq = t1Nseq[iyNseq]
lonmNseq = lonmNseq[iyNseq]
ltmNseq = ltmNseq[iyNseq]
hminsNseq = hminsNseq[iyNseq]
hmaxsNseq = hmaxsNseq[iyNseq]

s_allNseq = s_allNseq[iyNseq]
N_allNseq = N_allNseq[iyNseq]




####
t0S_datetime = pd.to_datetime(t0S, unit='s')
t0S_year = t0S_datetime.year

iyS = np.where((t0S_year == 2007) | (t0S_year == 2008) | (t0S_year == 2009))

t0S = t0S[iyS]
t1S = t1S[iyS]
lonmS = lonmS[iyS]
ltmS = ltmS[iyN]
hminsS = hminsS[iyS]
hmaxsS = hmaxsS[iyS]

s_allS = s_allS[iyS]
N_allS = N_allS[iyS]

####
t0Ss_datetime = pd.to_datetime(t0Ss, unit='s')
t0Ss_year = t0Ss_datetime.year

iySs = np.where((t0Ss_year == 2007) | (t0Ss_year == 2008) | (t0Ss_year == 2009))

t0Ss = t0Ss[iySs]
t1Ss = t1Ss[iySs]
lonmSs = lonmSs[iySs]
ltmSs = ltmSs[iyNs]
hminsSs = hminsSs[iySs]
hmaxsSs = hmaxsSs[iySs]

s_allSs = s_allSs[iySs]
N_allSs = N_allSs[iySs]


####
t0Sseq_datetime = pd.to_datetime(t0Sseq, unit='s')
t0Sseq_year = t0Sseq_datetime.year

iySseq = np.where((t0Sseq_year == 2007) | (t0Sseq_year == 2008) | (t0Sseq_year == 2009))

t0Sseq = t0Sseq[iySseq]
t1Sseq = t1Sseq[iySseq]
lonmSseq = lonmSseq[iySseq]
ltmSseq = ltmSseq[iyNseq]
hminsSseq = hminsSseq[iySseq]
hmaxsSseq = hmaxsSseq[iySseq]

s_allSseq = s_allSseq[iySseq]
N_allSseq = N_allSseq[iySseq]



#### Calculate mean altitudes in km
hmN = (hminsN + hmaxsN) / 2 / 1000
hmNs = (hminsNs + hmaxsNs) / 2 / 1000
hmNseq = (hminsNseq + hmaxsNseq) / 2 / 1000


hmS = (hminsS + hmaxsS) / 2 / 1000
hmSs = (hminsSs + hmaxsSs) / 2 / 1000
hmSseq = (hminsSseq + hmaxsSseq) / 2 / 1000


################################################################################
### GOCE
################################################################################
path_GOCE = 'C:/Poblet/IAP/JAPAN2023/research_topic/codes/read_outs_GOCE/'

t0_t1_lon_lt_hmin_hmax_GOCEN = np.load(path_GOCE+'t0_t1_lon_lt_hmin_hmax_u_lat-0-60N_10111213.npy', allow_pickle=True)

t0_GOCEN = t0_t1_lon_lt_hmin_hmax_GOCEN[0]
t1_GOCEN = t0_t1_lon_lt_hmin_hmax_GOCEN[1]
lonm_GOCEN = t0_t1_lon_lt_hmin_hmax_GOCEN[2]
ltm_GOCEN = t0_t1_lon_lt_hmin_hmax_GOCEN[3]
hmins_GOCEN = t0_t1_lon_lt_hmin_hmax_GOCEN[4]
hmaxs_GOCEN = t0_t1_lon_lt_hmin_hmax_GOCEN[5]

N_all_GOCEN = np.load(path_GOCE+'N_u_lat-0-60N_10111213.npy', allow_pickle=True)
s_all_GOCEN = np.load(path_GOCE+'s_u_lat-0-60N_10111213.npy', allow_pickle=True)
D_q_all_GOCEN = np.load(path_GOCE+'D_3_u_lat-0-60N_LT_10111213.npy', allow_pickle=True)
D_q_all_GOCEN_median = np.load(path_GOCE+'D_3_umedian_lat-0-60N_LT_10111213.npy', allow_pickle=True)


t0_t1_lon_lt_hmin_hmax_GOCENs = np.load(path_GOCE+'t0_t1_lon_lt_hmin_hmax_u_lat-30-60N_10111213.npy', allow_pickle=True)

t0_GOCENs = t0_t1_lon_lt_hmin_hmax_GOCENs[0]
t1_GOCENs = t0_t1_lon_lt_hmin_hmax_GOCENs[1]
lonm_GOCENs = t0_t1_lon_lt_hmin_hmax_GOCENs[2]
ltm_GOCENs = t0_t1_lon_lt_hmin_hmax_GOCENs[3]
hmins_GOCENs = t0_t1_lon_lt_hmin_hmax_GOCENs[4]
hmaxs_GOCENs = t0_t1_lon_lt_hmin_hmax_GOCENs[5]

N_all_GOCENs = np.load(path_GOCE+'N_u_lat-30-60N_10111213.npy', allow_pickle=True)
s_all_GOCENs = np.load(path_GOCE+'s_u_lat-30-60N_10111213.npy', allow_pickle=True)
D_q_all_GOCENs = np.load(path_GOCE+'D_3_u_lat-30-60N_LT_10111213.npy', allow_pickle=True)
D_q_all_GOCENs_median = np.load(path_GOCE+'D_3_umedian_lat-30-60N_LT_10111213.npy', allow_pickle=True)


t0_t1_lon_lt_hmin_hmax_GOCENseq = np.load(path_GOCE+'t0_t1_lon_lt_hmin_hmax_u_lat-0-30N_10111213.npy', allow_pickle=True)

t0_GOCENseq = t0_t1_lon_lt_hmin_hmax_GOCENseq[0]
t1_GOCENseq = t0_t1_lon_lt_hmin_hmax_GOCENseq[1]
lonm_GOCENseq = t0_t1_lon_lt_hmin_hmax_GOCENseq[2]
ltm_GOCENseq = t0_t1_lon_lt_hmin_hmax_GOCENseq[3]
hmins_GOCENseq = t0_t1_lon_lt_hmin_hmax_GOCENseq[4]
hmaxs_GOCENseq = t0_t1_lon_lt_hmin_hmax_GOCENseq[5]

N_all_GOCENseq = np.load(path_GOCE+'N_u_lat-0-30N_10111213.npy', allow_pickle=True)
s_all_GOCENseq = np.load(path_GOCE+'s_u_lat-0-30N_10111213.npy', allow_pickle=True)
D_q_all_GOCENseq = np.abs(np.load(path_GOCE+'D_3_u_lat-0-30N_LT_10111213_cyc.npy', allow_pickle=True))
D_q_all_GOCENseq_median = np.abs(np.load(path_GOCE+'D_3_umedian_lat-0-30N_LT_10111213_cyc.npy', allow_pickle=True))



t0_t1_lon_lt_hmin_hmax_GOCES = np.load(path_GOCE+'t0_t1_lon_lt_hmin_hmax_u_lat-0-60S_10111213.npy', allow_pickle=True)

t0_GOCES = t0_t1_lon_lt_hmin_hmax_GOCES[0]
t1_GOCES = t0_t1_lon_lt_hmin_hmax_GOCES[1]
lonm_GOCES = t0_t1_lon_lt_hmin_hmax_GOCES[2]
ltm_GOCES = t0_t1_lon_lt_hmin_hmax_GOCES[3]
hmins_GOCES = t0_t1_lon_lt_hmin_hmax_GOCES[4]
hmaxs_GOCES = t0_t1_lon_lt_hmin_hmax_GOCES[5]

N_all_GOCES = np.load(path_GOCE+'N_u_lat-0-60S_10111213.npy', allow_pickle=True)
s_all_GOCES = np.load(path_GOCE+'s_u_lat-0-60S_10111213.npy', allow_pickle=True)
D_q_all_GOCES = np.load(path_GOCE+'D_3_u_lat-0-60S_LT_10111213.npy', allow_pickle=True)
D_q_all_GOCES_median = np.load(path_GOCE+'D_3_umedian_lat-0-60S_LT_10111213.npy', allow_pickle=True)


t0_t1_lon_lt_hmin_hmax_GOCESs = np.load(path_GOCE+'t0_t1_lon_lt_hmin_hmax_u_lat-30-60S_10111213.npy', allow_pickle=True)

t0_GOCESs = t0_t1_lon_lt_hmin_hmax_GOCESs[0]
t1_GOCESs = t0_t1_lon_lt_hmin_hmax_GOCESs[1]
lonm_GOCESs = t0_t1_lon_lt_hmin_hmax_GOCESs[2]
ltm_GOCESs = t0_t1_lon_lt_hmin_hmax_GOCESs[3]
hmins_GOCESs = t0_t1_lon_lt_hmin_hmax_GOCESs[4]
hmaxs_GOCESs = t0_t1_lon_lt_hmin_hmax_GOCESs[5]

N_all_GOCESs = np.load(path_GOCE+'N_u_lat-30-60S_10111213.npy', allow_pickle=True)
s_all_GOCESs = np.load(path_GOCE+'s_u_lat-30-60S_10111213.npy', allow_pickle=True)
D_q_all_GOCESs = np.load(path_GOCE+'D_3_u_lat-30-60S_LT_10111213.npy', allow_pickle=True)
D_q_all_GOCESs_median = np.load(path_GOCE+'D_3_umedian_lat-30-60S_LT_10111213.npy', allow_pickle=True)


t0_t1_lon_lt_hmin_hmax_GOCESseq = np.load(path_GOCE+'t0_t1_lon_lt_hmin_hmax_u_lat-0-30S_10111213.npy', allow_pickle=True)

t0_GOCESseq = t0_t1_lon_lt_hmin_hmax_GOCESseq[0]
t1_GOCESseq = t0_t1_lon_lt_hmin_hmax_GOCESseq[1]
lonm_GOCESseq = t0_t1_lon_lt_hmin_hmax_GOCESseq[2]
ltm_GOCESseq = t0_t1_lon_lt_hmin_hmax_GOCESseq[3]
hmins_GOCESseq = t0_t1_lon_lt_hmin_hmax_GOCESseq[4]
hmaxs_GOCESseq = t0_t1_lon_lt_hmin_hmax_GOCESseq[5]

N_all_GOCESseq = np.load(path_GOCE+'N_u_lat-0-30S_10111213.npy', allow_pickle=True)
s_all_GOCESseq = np.load(path_GOCE+'s_u_lat-0-30S_10111213.npy', allow_pickle=True)
D_q_all_GOCESseq = np.load(path_GOCE+'D_3_u_lat-0-30S_LT_10111213_cyc.npy', allow_pickle=True)
D_q_all_GOCESseq_median = np.load(path_GOCE+'D_3_umedian_lat-0-30S_LT_10111213_cyc.npy', allow_pickle=True)



#### Calculate mean altitudes in km
hm_GOCEN = (hmins_GOCEN + hmaxs_GOCEN) / 2 / 1000
hm_GOCENs = (hmins_GOCENs + hmaxs_GOCENs) / 2 / 1000
hm_GOCENseq = (hmins_GOCENseq + hmaxs_GOCENseq) / 2 / 1000

hm_GOCES = (hmins_GOCES + hmaxs_GOCES) / 2 / 1000
hm_GOCESs = (hmins_GOCESs + hmaxs_GOCESs) / 2 / 1000
hm_GOCESseq = (hmins_GOCESseq + hmaxs_GOCESseq) / 2 / 1000

#################################################################################
### 1D plots of of D_3 s vs altitude- disgregating information
#################################################################################

bin_edges = np.logspace(np.log10(s_min), np.log10(s_max), num_bins + 1)

i_bin_edges = np.where( (bin_edges > 70) & (bin_edges <= 7000) )
bin_edges_plot = bin_edges[i_bin_edges]

### HERE WE SELECT NORTH OR SOUTH
Comp_N = [N_allN,N_allN,N_allS,N_allS,N_all_GOCEN,N_all_GOCEN,N_all_GOCES,N_all_GOCES]
Comp_Ns =[N_allNs,N_allNs,N_allSs,N_allSs,N_all_GOCENs,N_all_GOCENs,N_all_GOCESs,N_all_GOCESs]
Comp_Nseq =[N_allNseq,N_allNseq,N_allSseq,N_allSseq,N_all_GOCENseq,N_all_GOCENseq,N_all_GOCESseq,N_all_GOCESseq]

Comp_s = [s_allN,s_allN,s_allS,s_allS, s_all_GOCEN,s_all_GOCEN,s_all_GOCES,s_all_GOCES]
Comp_ss = [s_allNs,s_allNs,s_allSs,s_allSs,s_all_GOCENs,s_all_GOCENs,s_all_GOCESs,s_all_GOCESs]
Comp_sseq = [s_allNseq,s_allNseq,s_allSseq,s_allSseq,s_all_GOCENseq,s_all_GOCENseq,s_all_GOCESseq,s_all_GOCESseq]

Comp_D3 = [D_q_allN,D_q_allN,D_q_allS,D_q_allS,D_q_all_GOCEN,D_q_all_GOCEN,D_q_all_GOCES,D_q_all_GOCES]
Comp_D3s = [D_q_allNs,D_q_allNs,D_q_allSs,D_q_allSs,D_q_all_GOCENs,D_q_all_GOCENs,D_q_all_GOCESs,D_q_all_GOCESs]
Comp_D3seq = [D_q_allNseq,D_q_allNseq,D_q_allSseq,D_q_allSseq,D_q_all_GOCENseq,D_q_all_GOCENseq,D_q_all_GOCESseq,D_q_all_GOCESseq]

Comp_D3_median = [D_q_allN_median,D_q_allN_median,D_q_allS_median,D_q_allS_median,D_q_all_GOCEN_median,D_q_all_GOCEN_median,D_q_all_GOCES_median,D_q_all_GOCES_median]
Comp_D3s_median = [D_q_allNs_median,D_q_allNs_median,D_q_allSs_median,D_q_allSs_median,D_q_all_GOCENs_median,D_q_all_GOCENs_median,D_q_all_GOCESs_median,D_q_all_GOCESs_median]
Comp_D3seq_median = [D_q_allNseq_median,D_q_allNseq_median,D_q_allSseq_median,D_q_allSseq_median,D_q_all_GOCENseq_median,D_q_all_GOCENseq_median,D_q_all_GOCESseq_median,D_q_all_GOCESseq_median]

Comp_hm = [hmN,hmN,hmS,hmS,hm_GOCEN,hm_GOCEN,hm_GOCES,hm_GOCES]
Comp_hms = [hmNs,hmNs,hmSs,hmSs,hm_GOCENs,hm_GOCENs,hm_GOCESs,hm_GOCESs]
Comp_hmseq = [hmNseq,hmNseq,hmSseq,hmSseq,hm_GOCENseq,hm_GOCENseq,hm_GOCESseq,hm_GOCESseq]

panel = ['A1','B1','C1','D1','A2','B2','C2','D2',]
panel_inf_sup = ['sup', 'inf', 'sup', 'inf', 'sup', 'inf', 'sup', 'inf']
panel_title = ['(a) CHAMP - North','(c) CHAMP - North','(e) CHAMP - South','(g) CHAMP - South',
               '(b) GOCE - North','(d) GOCE - North','(f) GOCE - South','(h) GOCE - South']

hvalues_inf = [332, 310, 340, 318, 248, 228, 249, 236]
hvalues_sup = [354, 332, 362, 340, 268, 248, 277, 249]

fig, axd = plt.subplot_mosaic([['A1','A2'],['B1','B2'],['C1','C2'],['D1','D2']], constrained_layout=True)    
i = 0
D_3_all =[]
D_3s_all = []
D_3seq_all = []
D_3_all_median =[]
D_3s_all_median = []
D_3seq_all_median = []

D_3_all_GOCE = []
D_3s_all_GOCE = []
D_3seq_all_GOCE = []
D_3_all_GOCE_median = []
D_3s_all_GOCE_median = []
D_3seq_all_GOCE_median = []



for comp_s, comp_ss, comp_sseq, comp_N, comp_Ns, comp_Nseq, comp_D3, comp_D3s, comp_D3seq, comp_D3_median, comp_D3s_median, comp_D3seq_median, comp_hm, comp_hms, comp_hms, tpanel, lpanel, title_panel in zip(Comp_s,Comp_ss, Comp_sseq, Comp_N, Comp_Ns, Comp_Nseq, Comp_D3, Comp_D3s, Comp_D3seq, Comp_D3_median, Comp_D3s_median, Comp_D3seq_median, Comp_hm, Comp_hms, Comp_hmseq, panel_inf_sup, panel, panel_title):

    print(i+1)
    Sc = np.empty((len(bin_edges))) * np.nan
    Scs = np.empty((len(bin_edges))) * np.nan
    Scseq = np.empty((len(bin_edges))) * np.nan

    
    Nc = np.empty((len(bin_edges))) * np.nan
    Ncs = np.empty((len(bin_edges))) * np.nan
    Ncseq = np.empty((len(bin_edges))) * np.nan
    
    D_3c = np.empty((len(bin_edges))) * np.nan
    D_3cs = np.empty((len(bin_edges))) * np.nan
    D_3cseq = np.empty((len(bin_edges))) * np.nan
    
    D_3c_median = np.empty((len(bin_edges))) * np.nan
    D_3cs_median = np.empty((len(bin_edges))) * np.nan
    D_3cseq_median = np.empty((len(bin_edges))) * np.nan
    
    D_3std = np.empty((len(bin_edges))) * np.nan
    D_3stds = np.empty((len(bin_edges))) * np.nan
    D_3stdseq = np.empty((len(bin_edges))) * np.nan

    hminv = hvalues_inf[i]
    hmaxv = hvalues_sup[i]   

    # print(hminv, hmaxv)

    ihh = np.where((comp_hm >= hminv) & (comp_hm < hmaxv))

    Sloop = comp_s[ihh]
    Sloops = comp_ss[ihh]
    Sloopseq = comp_sseq[ihh]
    
    Nloop = comp_N[ihh]
    Nloops = comp_Ns[ihh]
    Nloopseq = comp_Nseq[ihh]
    
    D3loop = comp_D3[ihh]
    D3loops = comp_D3s[ihh]
    D3loopseq = comp_D3seq[ihh]

    D3loop_median = comp_D3_median[ihh]
    D3loops_median = comp_D3s_median[ihh]
    D3loopseq_median = comp_D3seq_median[ihh]

    
    # print(len(Sloop))
    # print(len(Nloop))
    # print(len(D3loop))
    
    for b in range(num_bins):
        S = []
        Ss = []
        Sseq = []
        Ns = []
        Nss = []
        Nsseq = []
        D3 = []
        D3s = []
        D3seq = []
        D3_median = []
        D3s_median = []
        D3seq_median = []

        for sloop, sloops, sloopseq, nloop, nloops, nloopseq, d3loop, d3loops, d3loopseq, d3loop_median, d3loops_median, d3loopseq_median in zip(Sloop, Sloops, Sloopseq, Nloop, Nloops, Nloopseq, D3loop, D3loops, D3loopseq, D3loop_median, D3loops_median, D3loopseq_median):
            # Ensure all arrays have the same size by trimming the longer ones
            ### full 0-60 sub-legs
            min_length = min(len(sloop), len(nloop), len(d3loop))
    
            sloop = sloop[:min_length]
            nloop = nloop[:min_length]
            d3loop = d3loop[:min_length]
            
            ###
            ifinal = np.where((sloop >= bin_edges[b]) & (sloop < bin_edges[b + 1]))

            if ifinal[0].size > 0:
                S.append(sloop[ifinal])
                Ns.append(nloop[ifinal])
                D3.append(d3loop[ifinal])
                D3_median.append(d3loop_median[ifinal])

            ### midlatitudes sub-legs
            min_lengths = min(len(sloops), len(nloops), len(d3loops))
    
            sloops = sloops[:min_lengths]
            nloops = nloops[:min_lengths]
            d3loops = d3loops[:min_lengths]
            d3loops_median = d3loops_median[:min_lengths]
        
            ###
            ifinals = np.where((sloops >= bin_edges[b]) & (sloops < bin_edges[b + 1]))

            if ifinals[0].size > 0:
                Ss.append(sloops[ifinals])
                Nss.append(nloops[ifinals])
                D3s.append(d3loops[ifinals])
                D3s_median.append(d3loops_median[ifinals])

            ### equatorial sub-legs
            min_lengthseq = min(len(sloopseq), len(nloopseq), len(d3loopseq))
    
            sloopseq = sloopseq[:min_lengthseq]
            nloopseq = nloopseq[:min_lengthseq]
            d3loopseq = d3loopseq[:min_lengthseq]
            d3loopseq_median = d3loopseq_median[:min_lengthseq]
        
            ###
            ifinalseq = np.where((sloopseq >= bin_edges[b]) & (sloopseq < bin_edges[b + 1]))

            if ifinalseq[0].size > 0:
                Sseq.append(sloopseq[ifinalseq])
                Nsseq.append(nloopseq[ifinalseq])
                D3seq.append(d3loopseq[ifinalseq])
                D3seq_median.append(d3loopseq_median[ifinalseq])


        if Ns:
            S_combined = np.concatenate(S)
            Ns_combined = np.concatenate(Ns)
            D3_combined = np.concatenate(D3)
            D3_combined_median = np.concatenate(D3_median)
            
            # Calculate the median and standard deviation
            median_S = np.nanmedian(S_combined)
            median_D3 = np.nanmedian(D3_combined)
            median_D3_median = np.nanmedian(D3_combined_median)
           
            std_D3 = np.nanstd(D3_combined)
            std_D3_median = np.nanstd(D3_combined_median)
            
            # Exclude values outside three sigma from the mean
            ## I keep this just to exclude some unrealistic values
            ## but it does not make any difference in practice

            valid_indices = np.where(np.abs(D3_combined - median_D3) < 3 * std_D3)
            valid_indices_median = np.where(np.abs(D3_combined_median - median_D3_median) < 3 * std_D3_median)
            
            # Apply the filter to S_combined, Ns_combined and D3_combined
            filtered_S = S_combined[valid_indices]
            filtered_Ns = Ns_combined[valid_indices]

            filtered_D3 = D3_combined[valid_indices]
            filtered_D3_median = D3_combined_median[valid_indices_median]
            
            Sc[b] = np.nanmedian(filtered_S)
            Nc[b] = np.nansum(filtered_Ns)
            D_3c[b] = np.nanmedian(filtered_D3)
            D_3c_median[b] = np.nanmedian(filtered_D3_median)

            D_3std[b] = np.nanstd(filtered_D3) / (Nc[b] - 1)
        

        if Nss:
            Ss_combined = np.concatenate(Ss)
            Nss_combined = np.concatenate(Nss)
            D3s_combined = np.concatenate(D3s)
            D3s_combined_median = np.concatenate(D3s_median)
            
            # Calculate the median and standard deviation
            median_Ss = np.nanmedian(Ss_combined)
            median_D3s = np.nanmedian(D3s_combined)
            median_D3s_median = np.nanmedian(D3s_combined_median)
            
            std_D3s = np.nanstd(D3s_combined)
            std_D3s_median = np.nanstd(D3s_combined_median)
            
            # Exclude values outside three sigma from the mean: 
            ## I keep this just to exclude some unrealistic values
            ## but it does not make any difference in practice
            valid_indicess = np.where(np.abs(D3s_combined - median_D3s) < 3 * std_D3s)
            valid_indicess_median = np.where(np.abs(D3s_combined_median - median_D3s_median) < 3 * std_D3s_median)
            
            # Apply the filter to Ns_combined and D3_combined
            filtered_Ss = Ss_combined[valid_indicess]
            filtered_Nss = Nss_combined[valid_indicess]
            filtered_D3s = D3s_combined[valid_indicess]
            filtered_D3s_median = D3s_combined_median[valid_indicess_median]
            
            Scs[b] = np.nanmedian(filtered_Ss)
            Ncs[b] = np.nansum(filtered_Nss)
            D_3cs[b] = np.nanmedian(filtered_D3s)
            D_3cs_median[b] = np.nanmedian(filtered_D3s_median)
            
            D_3stds[b] = np.nanstd(filtered_D3s) / (Ncs[b] - 1)

        if Nsseq:
            Sseq_combined = np.concatenate(Sseq)
            Nsseq_combined = np.concatenate(Nsseq)
            D3seq_combined = np.concatenate(D3seq)
            D3seq_combined_median = np.concatenate(D3seq_median)
            
            # Calculate the median and standard deviation
            median_Sseq = np.nanmedian(Sseq_combined)
            median_D3seq = np.nanmedian(D3seq_combined)
            median_D3seq_median = np.nanmedian(D3seq_combined_median)
            
            std_D3seq = np.nanstd(D3seq_combined)
            std_D3seq_median = np.nanstd(D3seq_combined_median)
            
            # Exclude values outside three sigma from the mean: 
            ## I keep this just to exclude some unrealistic values
            ## but it does not make any difference in practice
            valid_indicesseq = np.where(np.abs(D3seq_combined - median_D3seq) < 3 * std_D3seq)
            valid_indicesseq_median = np.where(np.abs(D3seq_combined_median - median_D3seq_median) < 3 * std_D3seq_median)
            
            # Apply the filter to Ns_combined and D3_combined
            filtered_Sseq = Sseq_combined[valid_indicesseq]
            filtered_Nsseq = Nsseq_combined[valid_indicesseq]
            filtered_D3seq = D3seq_combined[valid_indicesseq]
            filtered_D3seq_median = D3seq_combined_median[valid_indicesseq_median]
            
            Scseq[b] = np.nanmedian(filtered_Sseq)
            Ncseq[b] = np.nansum(filtered_Nsseq)
            D_3cseq[b] = np.nanmedian(filtered_D3seq)
            D_3cseq_median[b] = np.nanmedian(filtered_D3seq_median)
            
            D_3stdseq[b] = np.nanstd(filtered_D3seq) / (Ncseq[b] - 1)


    # Check that the short sublegs do not show artificial values outside s ranges
    iSvalid = np.where((Scs > 4000))
    D_3cs[iSvalid] =  np.nan
    D_3cs_median[iSvalid] =  np.nan

    # Panels discrimination            
    if lpanel in ['A1','B1','C1','D1']:
        D_3_all.append(D_3c)
        D_3s_all.append(D_3cs)
        D_3seq_all.append(D_3cseq)

        D_3_all_median.append(D_3c_median)
        D_3s_all_median.append(D_3cs_median)
        D_3seq_all_median.append(D_3cseq_median)

    else:
        D_3_all_GOCE.append(D_3c)
        D_3s_all_GOCE.append(D_3cs)
        D_3seq_all_GOCE.append(D_3cseq)

        D_3_all_GOCE_median.append(D_3c_median)
        D_3s_all_GOCE_median.append(D_3cs_median)
        D_3seq_all_GOCE_median.append(D_3cseq_median)
   

    # Positive - Negative discrimination
    D_3cN = np.empty((len(D_3c)))*np.nan
    D_3cP = np.empty((len(D_3c)))*np.nan
    iN = np.where(D_3c < 0)
    iP = np.where(D_3c >=0)
    D_3cN[iN] = D_3c[iN]
    D_3cP[iP] = D_3c[iP]
    
    D_3csN = np.empty((len(D_3cs)))*np.nan
    D_3csP = np.empty((len(D_3cs)))*np.nan
    isN = np.where(D_3cs < 0)
    isP = np.where(D_3cs >=0)
    D_3csN[isN] = D_3cs[isN] 
    D_3csP[isP] = D_3cs[isP]

    D_3cseqN = np.empty((len(D_3cseq)))*np.nan
    D_3cseqP = np.empty((len(D_3cseq)))*np.nan
    iseqN = np.where(D_3cseq < 0)
    iseqP = np.where(D_3cseq >=0)
    D_3cseqN[iseqN] = D_3cseq[iseqN] 
    D_3cseqP[iseqP] = D_3cseq[iseqP]

    D_3c_medianN = np.empty((len(D_3c_median)))*np.nan
    D_3c_medianP = np.empty((len(D_3c_median)))*np.nan
    i_medianN = np.where(D_3c_median < 0)
    i_medianP = np.where(D_3c_median >=0)
    D_3c_medianN[i_medianN] = D_3c_median[i_medianN] 
    D_3c_medianP[i_medianP] = D_3c_median[i_medianP]

    D_3cs_medianN = np.empty((len(D_3cs_median)))*np.nan
    D_3cs_medianP = np.empty((len(D_3cs_median)))*np.nan
    is_medianN = np.where(D_3cs_median < 0)
    is_medianP = np.where(D_3cs_median >=0)
    D_3cs_medianN[is_medianN] = D_3cs_median[is_medianN] 
    D_3cs_medianP[is_medianP] = D_3cs_median[is_medianP]    

    D_3cseq_medianN = np.empty((len(D_3cseq_median)))*np.nan
    D_3cseq_medianP = np.empty((len(D_3cseq_median)))*np.nan
    iseq_medianN = np.where(D_3cseq_median < 0)
    iseq_medianP = np.where(D_3cseq_median >=0)
    D_3cseq_medianN[iseq_medianN] = D_3cseq_median[iseq_medianN] 
    D_3cseq_medianP[iseq_medianP] = D_3cseq_median[iseq_medianP]    


    axd[lpanel].loglog(bin_edges, D_3cP, '.', color ='dodgerblue', label=r'$\langle \delta u^3\rangle$ 0-60°')
    axd[lpanel].loglog(bin_edges, D_3csP, 'x', color='blue', label=r'$\langle \delta u^3\rangle$ 30-60°')
    axd[lpanel].loglog(bin_edges, D_3cseqP, '<', markerfacecolor='none', color='cornflowerblue', label=r'$\langle \delta u^3\rangle$ 0-30°')

    axd[lpanel].loglog(bin_edges, D_3c_medianP, '>', markerfacecolor='none', color ='dodgerblue', label=r'$\langle \delta u^3\rangle_{med}$ 0-60°')
    axd[lpanel].loglog(bin_edges, D_3cs_medianP, '+', color='blue', label=r'$\langle \delta u^3\rangle_{med}$ 30-60°')
    axd[lpanel].loglog(bin_edges, D_3cseq_medianP, 'o', markerfacecolor='none', color='cornflowerblue', label=r'$\langle \delta u^3\rangle_{med}$ 0-30°')

    axd[lpanel].loglog(bin_edges, np.abs(D_3cN), '.', color ='crimson', label=r'$\langle \delta u^3\rangle$ 0-60°')
    axd[lpanel].loglog(bin_edges, np.abs(D_3csN), 'x', color='red', label=r'$\langle \delta u^3\rangle$ 30-60°')
    axd[lpanel].loglog(bin_edges, np.abs(D_3cseqN), '<', markerfacecolor='none', color='tomato', label=r'$\langle \delta u^3\rangle$ 0-30°')

    axd[lpanel].loglog(bin_edges, np.abs(D_3c_medianN), '>', markerfacecolor='none', color ='crimson', label=r'$\langle \delta u^3\rangle_{med}$ 0-60°')
    axd[lpanel].loglog(bin_edges, np.abs(D_3cs_medianN), '+', color='red', label=r'$\langle \delta u^3\rangle_{med}$ 30-60°')
    axd[lpanel].loglog(bin_edges, np.abs(D_3cseq_medianN), 'o', markerfacecolor='none', color='tomato', label=r'$\langle \delta u^3\rangle_{med}$ 0-30°')

    # axd[lpanel].loglog(bin_edges_plot, fix_coef_s2 *(bin_edges_plot)**2, '-k', linewidth=2)
    # axd[lpanel].loglog(bin_edges, fix_coef_s2 *0.0001* (bin_edges)**3, '-k', linewidth=2)

    axd[lpanel].set_xlim([s_min, s_max])
    axd[lpanel].set_ylim(ylimm_nc)
    axd[lpanel].grid()
    axd[lpanel].set_title(title_panel+' H: %02d-%02d km' % (hminv, hmaxv), fontsize=fs, loc='left')
    
    # axd[lpanel].set_ylabel('$<\delta u^3>$(m$^3$ s$^{-3})$', fontsize=fs)
    axd[lpanel].tick_params(labelsize=fs )
    # axd[lpanel].legend(loc='upper left', fontsize=fs - 1)

    # axd[lpanel].set_xlabel('$s$ (km)', fontsize=fs )

    if lpanel in ['A1','B1','C1','D1']:
        axd[lpanel].set_ylabel('(m$^3$ s$^{-3})$', fontsize=fs)
        axd[lpanel].loglog(bin_edges_plot, fix_coef_s2 *(bin_edges_plot)**2, '-k', linewidth=2, label='$\sim s^2$')
        axd[lpanel].loglog(bin_edges_plot, fix_coef_s2* 0.00001 *(bin_edges_plot)**3, '--k', linewidth=2, label='$\sim s^3$')

    else:
        axd[lpanel].loglog(bin_edges_plot, fix_coef_s2 * 0.2 *(bin_edges_plot)**2, '-k', linewidth=2, label='$\sim s^2$')
        axd[lpanel].loglog(bin_edges_plot, fix_coef_s2 * 0.00001 *(bin_edges_plot)**3, '--k', linewidth=2, label='$\sim s^3$')
        axd[lpanel].set_ylabel('')

    if lpanel not in ['D1', 'D2']:
        axd[lpanel].set_xlabel('')
    else:
        axd[lpanel].set_xlabel('$s$ (km)', fontsize=fs)
        
    handles, labels = axd[lpanel].get_legend_handles_labels()   
    i = i + 1 
fig.legend(handles, labels, loc='lower center', fontsize=fs, framealpha=0, ncol=5)
plt.show()

# # plt.savefig("C:/Poblet/IAP/JAPAN2023/research_topic/codes/D3_CHAMP_GOCE_severalH_50bins_mean_median_eqandmid.jpg", dpi=400)








######################################################################################
# A1: Lindborg and Cho (2001) - COMPARISON WITH AVERAGED CURVE
######################################################################################
dfp = pd.read_csv('C:/Poblet/IAP/JAPAN2023/research_topic/codes/D_3T_tropo_positive_lindborgcho2001.csv')
dfn = pd.read_csv('C:/Poblet/IAP/JAPAN2023/research_topic/codes/D_3T_tropo_negative_lindborgcho2001.csv')

dfpe = pd.read_csv('C:/Poblet/IAP/JAPAN2023/research_topic/codes/D_3T_strato_positive_lindborgcho2001.csv')
dfne = pd.read_csv('C:/Poblet/IAP/JAPAN2023/research_topic/codes/D_3T_strato_negative_lindborgcho2001.csv')

# Use the .to_numpy() method to convert the selected columns to a numpy array
# s = df[1].to_numpy()
DqTtropop = dfp.to_numpy()
DqTtropon = dfn.to_numpy()

DqTstratop = dfpe.to_numpy()
DqTstraton = dfne.to_numpy()

















### PLOT WITH AVERAGE CURVES 
#######################################################################################
# fig, axd = plt.subplot_mosaic([['A1','A2'],
                               # ['A3','A4']], constrained_layout=True)    
fig, axl = plt.subplot_mosaic([['A2','A3','A1']], constrained_layout=True)    

# Power scaling

axl['A1'].loglog(DqTtropop[:,0],DqTtropop[:,1],'.b',markersize=msize, label = r'$\langle \delta u_{T}^3 \rangle $, trop.')
axl['A1'].loglog(DqTtropon[:,0],DqTtropon[:,1],'.r',markersize=msize, label = r'$\langle \delta u_{T}^3 \rangle $, trop.')

axl['A1'].loglog(1.97907587845612,1.18723612000513,'xb',markersize=msize, label = r'$\langle \delta u_{T}^3 \rangle$, strat.')# this has only one value
axl['A1'].loglog(DqTstraton[:,0],DqTstraton[:,1],'xr',markersize=msize, label = r'$\langle \delta u_{T}^3 \rangle $, strat.')


# axd['A1'].loglog(DqTstraton[:,0],np.abs(fitd2l),'-k',linewidth=2)
axl['A1'].loglog(DqTstraton[:,0],0.1*(DqTstraton[:,0])**2,'-k',linewidth=2, label = '$\sim s^2$')


axl['A1'].set_title('(c) Cho & Lindborg (2001)', fontsize=fs-1, loc='left')
axl['A1'].set_xlabel('$s$ (km)', fontsize=fs)
axl['A1'].set_ylabel('m$^3$ s$^{-3}$', fontsize=fs)

axl['A1'].set_xlim([1, s_max])  # Se
axl['A1'].set_ylim([0.1, 1e6])
axl['A1'].grid()

axl['A1'].tick_params(labelsize=fs-1)

axl['A1'].legend(ncol=1, loc='upper left', fontsize=fs-2)  # Set fontsize to the desired value (e.g., fs)


#########################################################################################

D_3_all_mean = np.nanmedian(D_3_all, axis=0)
D_3s_all_mean = np.nanmedian(D_3s_all, axis=0)

D_3_all_GOCE_mean = np.nanmedian(D_3_all_GOCE, axis=0)
D_3s_all_GOCE_mean = np.nanmedian(D_3s_all_GOCE, axis=0)


axl['A2'].loglog(bin_edges, D_3_all_mean, '.', color ='blue', label='CHAMP 0-60°')
axl['A2'].loglog(bin_edges, D_3s_all_mean, 'x', color='purple', label='CHAMP 30-60°')

axl['A2'].loglog(bin_edges, D_3_all_GOCE_mean, '.', color ='darkgreen', label='GOCE 0-60°')
axl['A2'].loglog(bin_edges, D_3s_all_GOCE_mean, 'x', color='darkorange', label='GOCE 30-60°')

axl['A2'].loglog(bin_edges_plot, fix_coef_s2 * (bin_edges_plot)**2, '-k', linewidth=2, label = '$\sim s^2$')
# axl['A2'].loglog(bin_edges_plot, fix_coef_s2 * 0.00001 * (bin_edges_plot)**3, '--k', linewidth=2)

axl['A2'].set_title(r'(a) Averaged $\langle \delta u^3 \rangle$', fontsize=fs-1, loc='left')
axl['A2'].set_xlabel('$s$ (km)', fontsize=fs)
axl['A2'].set_ylabel('m$^3$ s$^{-3}$', fontsize=fs)

axl['A2'].set_xlim([1, s_max])  # Se
axl['A2'].set_ylim([0.1, 1e6])
axl['A2'].grid()

axl['A2'].tick_params(labelsize=fs-1)
axl['A2'].legend(ncol=1,  loc='upper left', fontsize=fs-2)  # Set fontsize to the desired value (e.g., fs)


#########################################################################################

D_3_all_mean_median = np.nanmedian(D_3_all_median, axis=0)
D_3s_all_mean_median = np.nanmedian(D_3s_all_median, axis=0)

D_3_all_GOCE_mean_median = np.nanmedian(D_3_all_GOCE_median, axis=0)
D_3s_all_GOCE_mean_median = np.nanmedian(D_3s_all_GOCE_median, axis=0)


# enstrophy flux estimation
Ens = [D_3_all_mean_median, D_3s_all_mean_median, D_3_all_GOCE_mean_median, D_3s_all_GOCE_mean_median] # The last two values of D_3_all_GOCE_mean_median are negative for some reason
Popt = []

for ens in Ens:
    # Find NaN indices
    inan = np.isnan(ens) | (ens < 0) # Remove any negative value (if there is one)
    
    # Remove NaN values from ens and corresponding bin_edges
    ens = ens[~inan]
    # print(ens)
    bin_edges_fit = bin_edges[~inan]  # Make sure bin_edges is defined in your code
    
    # Perform curve fitting only if there are non-NaN values
    if len(ens) > 0:
        # popt, pcov = curve_fit(enstrophy, bin_edges_fit, ens)
        popt_log10, pcov_log10 = curve_fit(enstrophy_linear, np.log10(bin_edges_fit*1000), np.log10(ens))
        popt = 10**(popt_log10)
        
        Popt.append(popt)
    else:
        # Handle the case where all values are NaN
        Popt.append(None)  # or any other suitable handling    
Popt = np.squeeze(Popt)
print(Popt)
print(np.mean(Popt))
print(np.std(Popt))


axl['A3'].loglog(bin_edges, D_3_all_mean_median, '.', color ='blue', label='CHAMP 0-60°')
axl['A3'].loglog(bin_edges, D_3s_all_mean_median, 'x', color='purple', label='CHAMP 30-60°')

axl['A3'].loglog(bin_edges, D_3_all_GOCE_mean_median, '.', color ='darkgreen', label='GOCE 0-60°')
axl['A3'].loglog(bin_edges, D_3s_all_GOCE_mean_median, 'x', color='darkorange', label='GOCE 30-60°')

# axl['A2'].loglog(bin_edges_plot, fix_coef_s2 * (bin_edges_plot)**2, '-k', linewidth=2)
axl['A3'].loglog(bin_edges_plot, fix_coef_s2 * 0.0001 * (bin_edges_plot)**3, '--k', linewidth=2, label = '$\sim s^3$')

# axl['A3'].loglog(bin_edges_plot, 1/8 *Popt[0] * (bin_edges_plot)**3, '-r', linewidth=2, label = '$\sim s^3$')
# axl['A3'].loglog(bin_edges_plot, 1/8 *Popt[1] * (bin_edges_plot)**3, '-r', linewidth=2, label = '$\sim s^3$')
# axl['A3'].loglog(bin_edges_plot, 1/8 *Popt[2] * (bin_edges_plot)**3, '-r', linewidth=2, label = '$\sim s^3$')
# axl['A3'].loglog(bin_edges_plot, 1/8 *Popt[3] * (bin_edges_plot)**3, '-r', linewidth=2, label = '$\sim s^3$')

axl['A3'].set_title(r'(b) Averaged $\langle \delta u^3 \rangle_{med}$', fontsize=fs-1, loc='left')
axl['A3'].set_xlabel('$s$ (km)', fontsize=fs)
axl['A3'].set_ylabel('m$^3$ s$^{-3}$', fontsize=fs)

axl['A3'].set_xlim([1, s_max])  # Se
axl['A3'].set_ylim([0.1, 1e6])
axl['A3'].grid()

axl['A3'].tick_params(labelsize=fs-1)
axl['A3'].legend(ncol=1,  loc='upper left', fontsize=fs-2)  # Set fontsize to the desired value (e.g., fs)


# # plt.savefig("C:/Poblet/IAP/JAPAN2023/research_topic/codes/D3_CHAMP_GOCE_total_50bins_mean_median.jpg", dpi=400)



 















