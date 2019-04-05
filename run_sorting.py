import sys
lib_path='/nics/d/home/nwbirge/large_analysis/ca_analysis/'
sys.path.insert(0,lib_path)
import fileread as fr
import predefined as pd
import numpy as np
from scipy.optimize import curve_fit

def gauss(x,a,mu,sigma):
    return a*np.exp(-(x-mu)**2./(2.*sigma**2.))

source_list=[31,32,48,49,50,51,52,53,82,83,113,114,115,131,132,133,134,135,136,137]
runlog=np.zeros(200,dtype=[('mu','3f'),('muerror','3f'),('sigma','3f'),('sigmaerror','3f'),])

for run in range(80,154):
    try:
        if run in source_list:
            continue
        loc='/lustre/haven/gamma/neutrons/ca45_data/2017/disk0/'
        if run>63:
            loc='/lustre/haven/gamma/neutrons/ca45_data/2017/disk1/'
        data=fr.trig(loc+'Run_'+str(run)+'.trig')[0]
        data=data[data['energy']>6000]
        i=0
        for pixel in [11,12,35]:
            bd,ch=int(pixel/8),int(pixel%8)
            h,b=np.histogram(pd.single_pixel(data,bd,ch)['energy'],bins=305,range=[0,20000])
            b=pd.cbins(b)
            mx=np.amax(h)
            mu=b[np.argmax(h)]
            window=1000
            beg,end=mu-window,mu+window
            fithist=h[pd.land(b>beg,b<end)]
            fitbins=b[pd.land(b>beg,b<end)]
            pars,vrs=curve_fit(gauss,fitbins,fithist,p0=[mx,mu,200],bounds=[0,np.inf])
            vrs=np.sqrt(np.diag(vrs))
            runlog[run]['mu'][i]=pars[1]
            runlog[run]['muerror'][i]=vrs[1]
            runlog[run]['sigma'][i]=pars[2]
            runlog[run]['sigmaerror'][i]=vrs[2]
            i+=1
    except FileNotFoundError:
        continue

np.save('/nics/d/home/nwbirge/job_activity/run_log',runlog)
