import numpy as np
import fileread as fr
import predefined as pd
from scipy.optimize import curve_fit
#DONT FORGET TO IMPORT CONDA ON THE ACF IF RUNNING THERE

def gauss(x,*pars):
    a,mu,sigma = pars
    return a*np.exp(-(x-mu)**2./(2.*sigma**2.))

def lingauss(x,*pars):
    a,mu,sigma,m,b=pars
    return a*np.exp(-(x-mu)**2./(2.*sigma**2.))+m*x+b


ebins,erange=1000,[0,5000]
bd,ch=4,3
sigs=np.zeros((361,8))
psigs = np.zeros((361,8))
count=0
lst = np.concatenate((np.linspace(10,100,10,dtype=int),np.linspace(200,1000,9,dtype=int)))
for i in lst:
    for j in lst:
        sigs[count,0:2]=i,j
        name = str(i)+'-'+str(j)+'-all.dat'
        fname='./Run_131-'
        data=fr.gen_output(fname+name)[0]
        data=pd.precuts(data)
        hist,bins= np.histogram(pd.pixcut(data,'energy',bd,ch),ebins,erange)
        bins=pd.cbins(bins)
        start=len(bins)-np.sum(bins>1500)
        amax= np.argmax(hist[start:])
        window=60
        beg,end = bins[start+amax]-window,bins[start+amax]+window
        fitbins=bins[pd.land(bins>beg,bins<end)]
        fithist = hist[pd.land(bins>beg,bins<end)]
        pars=[hist[start+amax],bins[start+amax],20,0,0]
        pars=curve_fit(lingauss,fitbins,fithist,p0=pars,ftol=0.001)[0]
        chisq= np.sqrt(np.sum((lingauss(fitbins,*pars)-fithist)**2.))
        sigs[count,2:7],sigs[count,7]=pars,chisq/(2*window+1-len(pars))


        data=fr.gen_output(fname+name)[0]
        data=data[data['energy']>3000]
        hist,bins= np.histogram(pd.pixcut(data,'energy',bd,ch),ebins,erange)
        bins=pd.cbins(bins)
        start=len(bins)-np.sum(bins>1500)
        amax= np.argmax(hist[start:])
        window=60
        beg,end = bins[start+amax]-window,bins[start+amax]+window
        fitbins=bins[pd.land(bins>beg,bins<end)]
        fithist = hist[pd.land(bins>beg,bins<end)]
        pars=[hist[start+amax],bins[start+amax],20,0,0]
        pars=curve_fit(lingauss,fitbins,fithist,p0=pars,ftol=0.001)[0]
        chisq= np.sqrt(np.sum((lingauss(fitbins,*pars)-fithist)**2.))
        psigs[count,2:7],psigs[count,7]=pars,chisq/(2*window+1-len(pars))



	print i,j,sigs[count,4],psigs[count,4]
        count+=1

np.save('shaping_time_scan-363pk-lingauss',sigs)
np.save('shaping_time_scan-pulser-lingauss',psigs)
