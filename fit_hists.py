import numpy as np
import fileread as fr
import predefined as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib import cm
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import sys

def gauss(x,*pars):
    a,mu,sigma = pars
    return a*np.exp(-(x-mu)**2./(2.*sigma**2.))

def lingauss(x,*pars):
    a,mu,sigma,m,b = pars
    return a*np.exp(-(x-mu)**2./(2.*sigma**2.))+m*x+b

def scanfit(x,*pars):
    a,b,c = pars
    return a*1.*x+b*1./x+c

fil='./shaping_time_scan-363pk-lingauss.npy'
col=4

sigs=np.load(fil)
x,y=np.concatenate((np.linspace(10,100,10),np.linspace(200,1000,9))),np.concatenate((np.linspace(10,100,10),np.linspace(200,1000,9)))
x,y=np.meshgrid(x,y)
z=np.zeros_like(x)
for i in range(len(sigs)):
    if col ==4:   
        z[i/19,i%19]= abs(sigs[i,col])
    elif col ==3:
        z[i/19,i%19]= abs((sigs[i,col])*0.157)
    else:
        z[i/19,i%19]= sigs[i,col]
    #print i/5,i%5,sigs[i,4]
z=np.transpose(z)

pts =np.concatenate((np.linspace(10,100,10,dtype=int),np.linspace(200,1000,9,dtype=int)))

with PdfPages('rise,top=200,100.pdf') as pdf:
    counts=0
    big,small=1000,3000
    for i in pts:
        for j in pts:
            if i>00 or j>00:
                arry=np.load(str(i)+'-'+str(j)+'-363-histbins.npy')
                hist=arry[:,0]
                bins=arry[:,1]
                amp,mx = sigs[counts,2:4]
                start=len(bins)-np.sum(bins>1500)
                amax= np.argmax(hist[start:])
                window=12
                beg,end = bins[start+amax]-window*5*0.5,bins[start+amax]+window*5*0.5
                fitbins=bins[pd.land(bins>beg,bins<end)]
                s,e=fitbins[0],fitbins[len(fitbins)-1]
                rng=np.linspace(s,e,200)
                fithist = hist[pd.land(bins>beg,bins<end)]
                plt.figure(figsize=(30,20))
                plt.errorbar(bins,hist,yerr=np.sqrt(hist),fmt='o',label='rise='+str(i)+' top='+str(j))
                plt.plot(rng,lingauss(rng,*sigs[counts,2:7]),'r',label=r'$\chi^2/DoF$= %0.1f $\sigma$= %0.1f '%(sigs[counts,7],sigs[counts,4]))
                plt.plot(bins[start+amax],hist[start+amax],'y*',label='Max hist val.',markersize=25)
                plt.legend(fontsize=25)
                #plt.yscale('log')
                plt.tick_params(labelsize=25)
                plt.xlabel('Energy (ADC)',fontsize=25)
                plt.ylabel('Counts',fontsize=25)
                plt.xlim((2000,2500))
                pdf.savefig()
                plt.close()
            counts+=1

print max(sigs[:,3]),min(sigs[sigs[:,3]>100][:,3])
'''
with PdfPages('top=100.pdf') as pdf:
    counts=0
    for i in pts:
        for j in pts:
            if i==200 or j==100:
                window=60
                arry=np.load(str(i)+'-'+str(100)+'-363-histbins.npy')
                hist=arry[:,0]
                bins=arry[:,1]
                mxbin=bins[len(bins[bins<1500])+np.argmax(hist[bins>1500])]
                mxhist=hist[len(bins[bins<1500])+np.argmax(hist[bins>1500])]
                amp,mx = sigs[counts,2:4]
                plotbins=np.linspace(mx-60,mx+60,300)
                plt.figure(figsize=(30,20))
                plt.plot(bins,hist,ls='steps',label='rise='+str(i)+' top='+str(100))
                plt.plot(plotbins,lingauss(plotbins,*sigs[counts,2:5]),'r',label=r'$\chi^2/DoF$= %0.1f'%(sigs[counts,7]))
                plt.plot(mxbin,mxhist,'y*',label='Max hist val.',markersize=25)
                plt.legend(fontsize=25)
                plt.yscale('log')
                plt.tick_params(labelsize=25)
                plt.xlim((2000,2500))
                pdf.savefig()
                plt.close()
                counts+=1'''
