# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 18:08:36 2023

This code calculate high-order structure functions using CHAMP observations.



@author: ma042
"""

import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np
import pandas as pd

def extract_winds(ipath, year, month, days, daysn, uu, dip = 'True'):     
    T = []
    H = []
    Lon = []
    Lat = []
    LT = []
    Wind_component = []
   
    for y in year:
        for m in month:
            dd = 0
            for d in range(1, daysn[dd] + 1): 
            # for d in days:
                try:
                    filename = ipath + '%04d%02d/wind_%04d%02d%02d.txt' % (y, m, y, m, d)
                    # filename = ipath + '/wind_%04d%02d%02d.txt' % (y, m, d)
                    # print(filename)
                    # Your code to process the file goes here
                    
                    data = []
                    # Open the file
                    with open(filename, 'r') as file:
                        for line in file:
                            # Skip comment lines starting with #
                            if line.startswith('#'):
                                continue
                
                            # Split the line into columns using whitespace as the delimiter
                            columns = line.split()
                
                            # Append the columns as a list to the data list
                            data.append(columns)
                    
                    # extract info                    
                    date = [row[0] for row in data]
                    time = [row[1] for row in data]
                    timestamps = []
                    
                    ### this should be changed so timestamp is a GPS time timestamp and not a UT.
                    ### The differences are small, but just to be complete.
                    for date_str, time_str in zip(date, time):
                        # Combine the date and time strings into a single datetime string
                        datetime_str = date_str + ' ' + time_str
                    
                        # Parse the datetime string into a datetime object
                        dt = datetime.strptime(datetime_str, "%Y-%m-%d %H:%M:%S.%f")
                    
                        # Convert the datetime object to a numerical timestamp (e.g., Unix timestamp)
                        timestamp = dt.timestamp()
                    
                        # Append the timestamp to the list
                        timestamps.append(timestamp)                    
                    
                    h =   np.array([float(row[3]) for row in data])
                    if dip:
                        # print('dip')
                        lon = np.array([float(row[7]) for row in data]) # Geo lon = 4, Geomag lon = 7 
                        lat = np.array([float(row[8]) for row in data]) # Geo lat = 5, Geomag lat = 8
                    else:
                        lon = np.array([float(row[4]) for row in data]) # Geo lon = 4, Geomag lon = 7 
                        lat = np.array([float(row[5]) for row in data]) # Geo lat = 5, Geomag lat = 8

                    lt = np.array([float(row[6]) for row in data])
                    u =   np.array([float(row[uu]) for row in data])
                    v =   np.array([float(row[uu+1]) for row in data])
                    
                    # u = u**2+v**2
                    
                    # stack info
                    
                    T =   np.concatenate((T, np.array(timestamps)))
                    H =   np.concatenate((H, h))
                    Lon = np.concatenate((Lon, lon))
                    Lat = np.concatenate((Lat, lat))
                    LT = np.concatenate((LT, lt))
                    Wind_component = np.concatenate((Wind_component, u))  
                    
                except FileNotFoundError:
                    print(f"File not found for {y}{m}{d}. Skipping...")
                    continue
                
                dd = dd + 1

    return T, H, Lon, Lat, LT, Wind_component


def great_circle(lon0, lat0, lon1, lat1, height=90.0):
    """ measure distance to lat0,lon0 """
    # spherical earth approximation!!!
    return( (6371.0+height) * (
        np.arccos(np.sin(np.pi*lat0/180.0) * np.sin(np.pi*lat1/180.0) + np.cos(np.pi*lat0/180.0) * np.cos(np.pi*lat1/180.0) * np.cos(np.pi*(lon0 - lon1)/180.0))
    ))



def calculate_structure_function_q(wind_speed_data, q, lat, lon, h):
    """ Estimate s values and D_u for user-given q values"""

    n = len(wind_speed_data)
    max_s = n - 1
    
    D_q = []  # Initialize as a Python list
    D_q_med = []  # Initialize as a Python list
    D_q_std = np.zeros(max_s)
    Ns = np.zeros(max_s)
    s = np.zeros(max_s)
    s_std = np.zeros(max_s)
    
    for ri in range(1, max_s + 1):               
        # s values
        latsa = lat[ri:]
        latsb = lat[:n - ri]

        lonsa = lon[ri:]
        lonsb = lon[:n - ri]
        
        # Determining the distances s        
        nc = len(latsa)
        sii = np.zeros(nc)
        for ci in range(0, nc):
            sii[ci] = great_circle(lonsa[ci], latsa[ci], lonsb[ci], latsb[ci], height=np.nanmean((h[:n - ri] + h[ri:]) / 2 / 1000))
            
        # s and sampling
        s[ri - 1] = np.nanmean(sii)
        s_std[ri - 1] = np.nanstd(sii)
        Ns[ri - 1] = len(sii)
        
        # Dq values
        min_length = min(len(wind_speed_data[:n - ri]), len(wind_speed_data[ri:]))

        # delta_v = wind_speed_data[:min_length] - wind_speed_data[ri:min_length + ri]
        delta_v = -(wind_speed_data[:min_length] - wind_speed_data[ri:min_length + ri])

        D_q.append(np.nanmean(delta_v ** q))
        D_q_med.append(np.nanmedian(delta_v ** q))

        D_q_std[ri - 1] = np.nanstd(delta_v ** q)

    # print(len(s), len(D_q), len(Ns))
    return s, np.array(D_q), np.array(D_q_med), Ns




def cube(x):
    if x >= 0:
        return x**(1/3)
    elif x < 0:
        return -(abs(x)**(1/3))



# EXTRACT WINDS FOR PARTICULAR TIME PERIOD
## path
ipath = 'C:/Poblet/IAP/JAPAN2023/research_topic/Data_CHAMP/'
# ipath = 'C:/Poblet/IAP/JAPAN2023/research_topic/05092023/'

# ## dates
year = [2007,2008,2009,2010]
months = [1,2,3,4,5,6,7,8,9,10,11,12]
days = [25,26,27,28] ## days of the month
dayssn = [31,28,31,30,31,30,31,31,30,31,30,31] ## how many days in each month
uu = 10 # wind component in the files: 10 crosstack winds, 12 zonal wind from HWM07 model

# year = [2007]
# months = [6]
# days = [30] ## days of the month
# dayssn = [30] ## how many days in each month
# uu = 10 # wind component in the files
          
## call the routine - check th edata
ts, h, lon, lat, lt, u = extract_winds(ipath, year, months, days, dayssn, uu, dip = 'True')

# plt.figure()
# plt.plot(pd.to_datetime(ts, unit='s'),h/1000,'.')

# plt.grid()

# ## plots to check the data
# plt.figure()
# plt.plot(lon,lat,'.')
# plt.grid()

# fig, ax = plt.subplots()
# hhh = ax.hist2d(lon, lat, bins=500)
# fig.colorbar(hhh[3], ax=ax)


## Filter data to separate individual paths at certain latitude ranges
### EQUATORIAL data
ts_eq = np.empty((len(lat)))*np.nan
lat_eq = np.empty((len(lat)))*np.nan
lt_eq = np.empty((len(lat)))*np.nan
lon_eq = np.empty((len(lat)))*np.nan
h_eq = np.empty((len(lat)))*np.nan
u_eq = np.empty((len(lat)))*np.nan

# ieq = np.where(((lat <= 60) & (lat >30)) | ((lat >= -60) & (lat < -30)))
ieq = np.where((lat < 0) & (lat >= -30))

ts_eq[ieq] = ts[ieq]
lat_eq[ieq] = lat[ieq]
lt_eq[ieq] = lt[ieq]
lon_eq[ieq] = lon[ieq]
h_eq[ieq] = h[ieq]
u_eq[ieq] = u[ieq]

# plt.figure()
# plt.plot(pd.to_datetime(ts_eq, unit='s'),h_eq/1000,'.')
# plt.grid()


# Find indices of "nan" values
nan_indices = np.where(np.isnan(lat_eq))[0] ## only look nan on this variable because they are located in the same position in the other variables

# Initialize the start index
start_idx = 0

ts_eq_subarrays = []
lon_eq_subarrays = []
lat_eq_subarrays = []
lt_eq_subarrays = []
h_eq_subarrays = []
u_eq_subarrays = []
len_eq = []
for nan_idx in nan_indices:
    subarray = lat_eq[start_idx:nan_idx]
  
    # Check if the subarray is not empty
    # print(len(subarray))
    if (len(subarray) >= 42) & (len(subarray) < 50): # cutoff for every subleg size (to improve statistics) 89-95: for 0-60 degrees. 42-50: for 30-60 degrees.
        # print(len(subarray))
        lat_eq_subarrays.append(subarray)
        lt_eq_subarrays.append(lt_eq[start_idx:nan_idx])
        ts_eq_subarrays.append(ts_eq[start_idx:nan_idx])
        lon_eq_subarrays.append(lon_eq[start_idx:nan_idx])
        h_eq_subarrays.append(h_eq[start_idx:nan_idx])
        u_eq_subarrays.append(u_eq[start_idx:nan_idx])
        len_eq.append(len(subarray))
    
    # Update the start index for the next subarray
    start_idx = nan_idx + 1

# # mean altitude difference between consecutive points
# print(np.mean(np.diff(h_eq_subarrays[1]/1000)))
# print(np.std(np.diff(h_eq_subarrays[1]/1000)))



# ESTIMATE THE STRUCTURE FUNCTIONS 
q = 3  # Replace with the desired q value

## ESTIMATE D_q FOR EVERY SATELLITE SUBLEG AND ESTIMATE IDENTIFICATION PARAMETERS
s = []
Ns = []
D_q = []
D_q_med = []
ts_desde_subleg =  []
ts_hasta_subleg =  []
lon_mean_subleg =  []
lt_mean_subleg =  []
h_min_subleg = []
h_max_subleg = []
i = 0
for lati in lat_eq_subarrays:
    # # # Anticyclonic zonal wind decomposition: North
    # # if lati[0] > lati[-1]:
    # #     u_eq_subarrayslt =  -u_eq_subarrays[i]
    # # else:
    # #     u_eq_subarrayslt =  u_eq_subarrays[i]
        
    # # Anticyclonic zonal wind decomposition: South
    # if lati[0] > lati[-1]:
    #     u_eq_subarrayslt =  u_eq_subarrays[i]
    # else:
    #     u_eq_subarrayslt = -u_eq_subarrays[i]

    # # Cyclonic zonal wind decomposition: North
    # if lati[0] > lati[-1]:
    #     u_eq_subarrayslt =  u_eq_subarrays[i]
    # else:
    #     u_eq_subarrayslt =  -u_eq_subarrays[i]
        
    # Cyclonic zonal wind decomposition: South
    if lati[0] > lati[-1]:
        u_eq_subarrayslt =  -u_eq_subarrays[i]
    else:
        u_eq_subarrayslt = u_eq_subarrays[i]

    sl, D_ql, D_ql_med, ns = calculate_structure_function_q(u_eq_subarrayslt,
                                                  q,
                                                  # max_s,
                                                  lati,
                                                  lon_eq_subarrays[i],
                                                  h_eq_subarrays[i])
    
    ts_desde_subleg.append(np.min(ts_eq_subarrays[i]))
    ts_hasta_subleg.append(np.max(ts_eq_subarrays[i]))
    lon_mean_subleg.append(np.nanmean(lon_eq_subarrays[i]))
    lt_mean_subleg.append(np.nanmean(lt_eq_subarrays[i]))
    h_min_subleg.append(np.nanmin(h_eq_subarrays[i]))
    h_max_subleg.append(np.nanmax(h_eq_subarrays[i]))
    
    
    ## print esto
    dt_desde = datetime.utcfromtimestamp(np.min(ts_eq_subarrays[i]))
    dt_hasta = datetime.utcfromtimestamp(np.max(ts_eq_subarrays[i]))

    print('Subleg time from: '+dt_desde.strftime('%Y-%m-%d') +' '+dt_desde.strftime('%H:%M:%S.%f')[:-3]+' GPS time to '+dt_hasta.strftime('%Y-%m-%d')+' '+dt_hasta.strftime('%H:%M:%S.%f')[:-3]+' GPS time')
    
    s.append(sl)
    Ns.append(ns)
    D_q.append(D_ql)
    D_q_med.append(D_ql_med)
    
    i = i + 1

# SAVE data into files

t0_t1_lon_lt_hmin_hmax = np.vstack((np.array(ts_desde_subleg,dtype=object),
                         np.array(ts_hasta_subleg,dtype=object), 
                         np.array(lon_mean_subleg,dtype=object),
                         np.array(lt_mean_subleg,dtype=object),
                         np.array(h_min_subleg,dtype=object),
                         np.array(h_max_subleg,dtype=object)))

s = np.array(s,dtype=object)
Ns = np.array(Ns,dtype=object)
D_q = np.array(D_q,dtype=object)
D_q_med = np.array(D_q_med,dtype=object)

np.save('t0_t1_lon_lt_hmin_hmax_u_diplat-0-30S_07080910.npy', t0_t1_lon_lt_hmin_hmax)
np.save('s_u_diplat-0-30S_07080910.npy', s)
np.save('N_u_diplat-0-30S_07080910.npy', Ns)
np.save('D_3_u_diplat-0-30S_LT_07080910_cyc.npy', D_q)
np.save('D_3_umedian_diplat-0-30S_LT_07080910_cyc.npy', D_q_med)
