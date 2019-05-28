
# coding: utf-8

# In[2]:


import sys
sys.path.insert(0,'/home/noah/Desktop/large_analysis/ca_analysis/')
from matplotlib.pyplot import*
import fileread as fr
import predefined as pd
import wave_ops as wo
import time
import numpy as np
from scipy.optimize import curve_fit
import os
from scipy.constants import *

if len(sys.argv) < 2:
    print( 'Arguments: Run number (required) \n'+' '*len('Arguments: ')+'Coincidence Window (*optional* defaults to 100 timebins)')
    exit()
elif len(sys.argv) ==2:
    run= int(sys.argv[1])
    twindow=100
else:
    run= int(sys.argv[1])
    twindow=int(sys.argv[2])
print('Summing multipixel events for Run {0:0.0f} with a coincidence window of {1:0.0f} timebins'.format(run,twindow))


calibration=np.load('/home/noah/Desktop/large_analysis/ca_analysis/simulation_comparison/calibration.npy')
calibration=calibration.view(np.recarray)
def calibrate(energy_type,board,channel): 
    bdch=int(board*8+channel) 
    if bdch ==6: 
        m,b=1/calibration.slope[3],calibration.offset[3]
    elif bdch==11: 
        m,b=1/calibration.slope[0],calibration.offset[0]
    elif bdch==12: 
        m,b=1/calibration.slope[1],calibration.offset[1]
    elif bdch==35: 
        m,b=1/calibration.slope[2],calibration.offset[2]
    else: 
        m,b=0,0   
    return (energy_type-b)*m
vec_calibrate=np.vectorize(calibrate)

path='/home/noah/Desktop/large_analysis/ca_analysis/cur_data/'


data=fr.gen_output(path+'Run_'+str(run)+'-all.dat')[0]
data=data[pd.good_timestamps(data)]
data.sort(order='timestamp')
beg=time.time()
energy_type='energy'
multi=data.copy()
multi[energy_type]=vec_calibrate(multi[energy_type],multi['board'],multi['channel'])
multi=multi[multi['t0']>600]
multi.sort(order='timestamp')
multi=multi[pd.lor(multi['pilediff']<twindow,multi['pileup']<2)]
multi=multi[pd.doubles(multi,energy_type)]
print('Cut double events in {:0.0f} s'.format(time.time()-beg))

cut=multi.copy()

beg=time.time()
i,j=0,0
while i<len(multi)-1 and j<len(multi)-1:
    j=i+1
    backscattering=multi['timestamp'][j]-multi['timestamp'][i] < twindow
    energy=multi[energy_type][i]
    while backscattering and j<len(multi)-1:
        energy+=multi[energy_type][j]
        multi[energy_type][j]=-10
        j+=1
        backscattering=multi['timestamp'][j]-multi['timestamp'][i] < twindow
    multi[energy_type][i]=energy
    i=j
print('Summed multipixel events in {:0.0f} s'.format(time.time()-beg))
np.save(file=path+'multi_'+str(run)+'_window-'+str(twindow),arr=multi)       ####This is the array of interest
np.save(file=path+'pre-multi_cut_'+str(run)+'_window-'+str(twindow),arr=cut)


# In[4]:




