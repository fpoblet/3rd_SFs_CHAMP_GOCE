# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 16:54:40 2023

Combined plot of s vs altitude for CHAMP and GOCE. 

Number of pairs per s analysis 

@author: ma042
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.dates as mdates

def create_discrete_cmap(colormap, num_colors):
    # Create a ListedColormap with the specified number of colors
    cmap = plt.get_cmap(colormap, num_colors)
    discrete_cmap = mcolors.ListedColormap(cmap(np.arange(num_colors)))
    return discrete_cmap



def plot_pcolor_norm(ax, num_colors, vminl, vmaxl, fs, bin_edges, hvalues, tref1, tref2, D_qn_binned, D_qp_binned, titlet='0° - 60°'):
    cmap_p = discrete_oranges  # Use the 'Oranges' colormap
    # cmap_p = 'Oranges'  # Use the 'Oranges' colormap
    norm_p = mcolors.LogNorm(vmin=vminl, vmax=vmaxl)
    
    # Set x-axis to log scale
    # ax.set_xscale('log')
    
    imp = ax.pcolor(bin_edges, hvalues, D_qp_binned.T, shading='auto', cmap=cmap_p, norm=norm_p)

    
    cmap_n = discrete_greens  # Use the 'Greens' colormap
    # cmap_n = 'Greens'# Use the 'Greens' colormap
    norm_n = mcolors.LogNorm(vmin=vminl, vmax=vmaxl)
    imn = ax.pcolor(bin_edges, hvalues, D_qn_binned.T, shading='auto', cmap=cmap_n, norm=norm_n)
  
    # Create a colorbar for D_qp_binned inside the plot
    cbar_p = plt.colorbar(imp, ax=ax, shrink=0.5)
    cbar_p.ax.yaxis.set_tick_params(labelsize=fs-1)
    cbar_p.ax.set_ylabel('GOCE # of pairs', fontsize=fs-1)
    cbar_p.ax.yaxis.set_label_position('left')
    
    # Create a colorbar for D_qn_binned outside the plot
    cbar_n = plt.colorbar(imn, ax=ax, shrink=0.5, pad=0.1)
    cbar_n.ax.yaxis.set_tick_params(labelsize=fs-1)
    cbar_n.ax.set_ylabel('CHAMP # of pairs', fontsize=fs-1)
    cbar_n.ax.yaxis.set_label_position('left')
    
    ax.set_xlabel("$s$ (km)", fontsize=fs-1)
    ax.set_ylabel("Height (km)", fontsize=fs-1)
    
    ax.set_title(titlet, fontsize=fs-1, loc='center')
    ax.tick_params(labelsize=fs-1)


################################################################################
### PARAMETERS
################################################################################

#### Altitude ranges
hmin = 220
hmax = 380
step = 1
hvalues = np.arange(hmin, hmax, step)

#### s ranges
s_min = 70
s_max = 1e4
num_bins = 50

## FontSize for every label
fs = 16

ylimm_nc = [10, 1e6]
# ylimm_nc = [1, 1e2]


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

t0_t1_lon_lt_hmin_hmaxS = np.load(path+'t0_t1_lon_lt_hmin_hmax_u_lat-0-60S_07080910.npy', allow_pickle=True)

t0S = t0_t1_lon_lt_hmin_hmaxS[0]
t1S = t0_t1_lon_lt_hmin_hmaxS[1]
lonmS = t0_t1_lon_lt_hmin_hmaxS[2]
ltmS = t0_t1_lon_lt_hmin_hmaxS[3]
hminsS = t0_t1_lon_lt_hmin_hmaxS[4]
hmaxsS = t0_t1_lon_lt_hmin_hmaxS[5]

N_allS = np.load(path+'N_u_lat-0-60S_07080910.npy', allow_pickle=True)
s_allS = np.load(path+'s_u_lat-0-60S_07080910.npy', allow_pickle=True)


### USE only 2007-2009 ############################################################
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

#### Calculate mean altitudes in km
hmN = (hminsN + hmaxsN) / 2 / 1000
hmS = (hminsS + hmaxsS) / 2 / 1000

plt.figure()
plt.plot((hmaxsS - hminsS)/1000,'or')
plt.plot((hmaxsN - hminsN)/1000,'ob')


#### Time for reference in plot
t_refN = np.empty((len(hvalues)))*np.nan
t_refS = np.empty((len(hvalues)))*np.nan
i = 0
for h in hvalues:
    itN = np.where((hmN >= h) & (hmN < h + step))
    itS = np.where((hmS >= h) & (hmS < h + step))
    if len(itN[0]) > 0:  # Check if it is not an empty array
        try:
            t_refN[i] = np.nanmean((t0N[itN] + t1N[itN]) / 2)
        except Exception as e:
            print(f"An error occurred: {e}")
    if len(itS[0]) > 0:  # Check if it is not an empty array
        try:
            t_refS[i] = np.nanmean((t0S[itS] + t1S[itS]) / 2)
        except Exception as e:
            print(f"An error occurred: {e}")
    
    i += 1

# # Now, explicitly set NaN for values of h within the range 318 to 325
# maskn = (hvalues >= 318) & (hvalues <= 325)
# t_refN[maskn] = np.nan
tref_datetimeN = pd.to_datetime(t_refN, unit='s')
# masks = (hvalues >= 328) & (hvalues <= 335)
# t_refS[masks] = np.nan
tref_datetimeS = pd.to_datetime(t_refS, unit='s')




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


t0_t1_lon_lt_hmin_hmax_GOCES = np.load(path_GOCE+'t0_t1_lon_lt_hmin_hmax_u_lat-0-60S_10111213.npy', allow_pickle=True)

t0_GOCES = t0_t1_lon_lt_hmin_hmax_GOCES[0]
t1_GOCES = t0_t1_lon_lt_hmin_hmax_GOCES[1]
lonm_GOCES = t0_t1_lon_lt_hmin_hmax_GOCES[2]
ltm_GOCES = t0_t1_lon_lt_hmin_hmax_GOCES[3]
hmins_GOCES = t0_t1_lon_lt_hmin_hmax_GOCES[4]
hmaxs_GOCES = t0_t1_lon_lt_hmin_hmax_GOCES[5]

N_all_GOCES = np.load(path_GOCE+'N_u_lat-0-60S_10111213.npy', allow_pickle=True)
s_all_GOCES = np.load(path_GOCE+'s_u_lat-0-60S_10111213.npy', allow_pickle=True)

#### Calculate mean altitudes in km
hm_GOCEN = (hmins_GOCEN + hmaxs_GOCEN) / 2 / 1000
hm_GOCES = (hmins_GOCES + hmaxs_GOCES) / 2 / 1000

plt.figure()
plt.plot((hmaxs_GOCEN - hmins_GOCEN)/1000,'or')
plt.plot((hmaxs_GOCES - hmins_GOCES)/1000,'ob')



#### Time for reference in plot
t_ref_GOCEN = np.empty((len(hvalues)))*np.nan
t_ref_GOCES = np.empty((len(hvalues)))*np.nan
i = 0
for h in hvalues:
    itN = np.where((hm_GOCEN >= h) & (hm_GOCEN < h + step))
    itS = np.where((hm_GOCES >= h) & (hm_GOCES < h + step))
    if len(itN[0]) > 0:  # Check if it is not an empty array
        try:
            t_ref_GOCEN[i] = np.nanmean((t0_GOCEN[itN] + t1_GOCEN[itN]) / 2)
        except Exception as e:
            print(f"An error occurred: {e}")
    if len(itS[0]) > 0:  # Check if it is not an empty array
        try:
            t_ref_GOCES[i] = np.nanmean((t0_GOCES[itS] + t1_GOCES[itS]) / 2)
        except Exception as e:
            print(f"An error occurred: {e}")
    i += 1

tref_datetime_GOCEN = pd.to_datetime(t_ref_GOCEN, unit='s')
tref_datetime_GOCES = pd.to_datetime(t_ref_GOCES, unit='s')

#################################################################################
### 2D histogram of counts s vs altitude
#################################################################################
### Here we choose north or south
s_all = s_allN
N_all = N_allN

s_all_GOCE = s_all_GOCEN
N_all_GOCE = N_all_GOCEN

hm = hmN
hm_GOCE = hm_GOCEN

tref_datetime = tref_datetimeN
tref_datetime_GOCE = tref_datetime_GOCEN


# bin_edges = np.logspace(np.log10(s_min), np.log10(s_max), num_bins + 1)
# bin_edges = np.arange(s_min, s_max, num_bins+1)
bin_edges = np.round(s_all[0])
delta_bin = round(np.diff(bin_edges)[0]/2)

bin_edgesG = np.round(s_all_GOCE[0])
delta_binG = round(np.diff(bin_edgesG)[0]/2)

Nc = np.empty((len(bin_edges), len(hvalues))) * np.nan
Nc_GOCE = np.empty((len(bin_edgesG), len(hvalues))) * np.nan
for ih in range(0, len(hvalues)):
    ihh = np.where((hm >= hvalues[ih]) & (hm < hvalues[ih] + 1))
    Sloop = s_all[ihh]
    Nloop = N_all[ihh]

    ihhG = np.where((hm_GOCE >= hvalues[ih]) & (hm_GOCE < hvalues[ih] + 1))
    SloopG = s_all_GOCE[ihhG]
    NloopG = N_all_GOCE[ihhG]

    # for b in range(num_bins):
    for b in range(len(bin_edges)):
        Ns = []
        for sloop, nloop in zip(Sloop, Nloop):
            ifinal = np.where(np.abs(sloop-bin_edges[b]) <= delta_bin)
            Ns.append(nloop[ifinal])

        if Ns:
            Ns_combined = np.concatenate(Ns)
            Nc[b, ih] = np.nansum(Ns_combined)
           
    for b in range(len(bin_edgesG)):
        NsG = []
        for sloopG, nloopG in zip(SloopG, NloopG):
            ifinalG = np.where(np.abs(sloopG-bin_edgesG[b]) <= delta_binG)
            NsG.append(nloopG[ifinalG])

        if NsG:
            NsG_combined = np.concatenate(NsG)
            Nc_GOCE[b, ih] = np.nansum(NsG_combined)



vminl = ylimm_nc[0]
vmaxl = ylimm_nc[1]
fig, axes = plt.subplots(1, 1)

# Ensure axes is a list if it's a single subplot
if not isinstance(axes, (list, np.ndarray)):
    axes = [axes]


num_colors = 20  # Specify the number of colors you want
discrete_oranges = create_discrete_cmap('Oranges', num_colors)
discrete_greens = create_discrete_cmap('Greens', num_colors)

plot_pcolor_norm(axes[0], 
                 num_colors, 
                 vminl, 
                 vmaxl, 
                 fs, 
                 bin_edges, 
                 hvalues, 
                 tref_datetime,
                 tref_datetime_GOCE,
                 Nc, 
                 Nc_GOCE, 
                 titlet='Sampling 0-60°N')

plt.show()

# # plt.savefig("C:/Poblet/IAP/JAPAN2023/research_topic/codes/sampling_N_CHAMP_GOCE.jpg", dpi=400)

plt.figure()
plt.plot(bin_edgesG,Nc_GOCE, '-o')

plt.figure()
plt.plot(bin_edges,Nc, '-o')




#################################################################################
### t0 vs Altitude 
#################################################################################
t0N_datetime = pd.to_datetime(t0N, unit='s')
t0_GOCEN_datetime = pd.to_datetime(t0_GOCEN, unit='s')
t0S_datetime = pd.to_datetime(t0S, unit='s')
t0_GOCES_datetime = pd.to_datetime(t0_GOCES, unit='s')

fs = 14

# Assuming you have your data already defined, e.g., t0_datetime, hmins, hmaxs, t0_GOCE_datetime, hmins_GOCE, hmaxs_GOCE

# Set the window size for the rolling average
window_size = 5

ms = 3
# Smooth the data using rolling mean
smoothed_dataN = pd.Series((hminsN + hmaxsN) / 2 / 1000).rolling(window=window_size, center=True).mean()
smoothed_data_GOCEN = pd.Series((hmins_GOCEN + hmaxs_GOCEN) / 2 / 1000).rolling(window=window_size, center=True).mean()

# Create the plot
plt.figure()
plt.plot(t0N_datetime, (hminsN + hmaxsN) / 2 / 1000, 'o',markersize=ms, color='darkgreen', label='CHAMP 0-60°N')
# plt.plot(t0_datetime, smoothed_data, '-', color='green', label='Smoothed Data')
plt.plot(t0_GOCEN_datetime, (hmins_GOCEN + hmaxs_GOCEN) / 2 / 1000, 'o',markersize=ms, color='darkorange', label='GOCE 0-60°N')
# plt.plot(t0_GOCE_datetime, smoothed_data_GOCE, '-', color='orange', label='Smoothed Data (GOCE)')

plt.plot(t0S_datetime, (hminsS + hmaxsS) / 2 / 1000, 'o',markersize=ms, color='darkblue', label='CHAMP 0-60°S')
# plt.plot(t0_datetime, smoothed_data, '-', color='green', label='Smoothed Data')
plt.plot(t0_GOCES_datetime, (hmins_GOCES + hmaxs_GOCES) / 2 / 1000, 'o',markersize=ms, color='darkred', label='GOCE 0-60°S')
# plt.plot(t0_GOCE_datetime, smoothed_data_GOCE, '-', color='orange', label='Smoothed Data (GOCE)')


plt.xlabel('$t_0$ (year)', fontsize=fs + 2)
plt.ylabel('mean sub-leg altitude (km)', fontsize=fs)
plt.title('(b) Altitude evolution of sub-legs', fontsize=fs, loc='left')  # Changed from ylabel to title
plt.legend(fontsize=fs-2)
plt.grid()

# Set the tick label font size for both x and y axes
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)


# plt.savefig("C:/Poblet/IAP/JAPAN2023/research_topic/codes/altitude_t0_CHAMP_GOCE.jpg", dpi=400)