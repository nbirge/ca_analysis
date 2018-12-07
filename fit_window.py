import numpy as np
import predefined as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import erf
from matplotlib.backends.backend_pdf import PdfPages
import sys


def gauss(x,*pars):
    a,mu,sigma = pars
    return a*np.exp(-(x-mu)**2./(2.*sigma**2.))

def lingauss(x,*pars):
    a,mu,sigma,m,b = pars
    return a*np.exp(-(x-mu)**2./(2.*sigma**2.))+m*x+b

def erfgauss(x,*pars):
    a,mu,sigma,b,c=pars
    return a*np.exp(-(x-mu)**2./(2.*sigma**2.))+b*erf(x-mu)+c

def offgauss(x,*pars):
    a,mu,sigma,b=pars
    return a*np.exp(-(x-mu)**2./(2.*sigma**2.))+b

fittype=int(sys.argv[1])

fil='200-100-363-histbins.npy'
x=np.load(fil)

hist,bins= x[:,0],x[:,1]
start=len(bins)-np.sum(bins>1500)
amax= np.argmax(hist[start:])
window=12  #in bins

if fittype==0:
    with PdfPages(fil[:-12]+'-gauss.pdf') as pdf:
        mini,minj,minchi=100.,100.,100.
        for i in np.linspace(0.2,1.,9):
            for j in np.linspace(0.2,1.,5):
                beg,end = bins[start+amax]-window*5*i,bins[start+amax]+window*5*j
                fitbins=bins[pd.land(bins>beg,bins<end)]
                fithist = hist[pd.land(bins>beg,bins<end)]
                fithist[fithist==0]=1
                weights=np.sqrt(fithist)
                pars=[hist[start+amax],bins[start+amax],30]
                pars=curve_fit(gauss,fitbins,fithist,p0=pars,sigma=weights,maxfev=120000)[0]
                chisq= np.sum(((gauss(fitbins,*pars)-fithist)/weights)**2.)
                chisq=chisq/(len(fitbins)-len(pars))
                if chisq<minchi and chisq>0.8:
                    mini,minj,minchi=i,j,chisq
                plt.figure(figsize=(30,20))
                plt.errorbar(bins+0.0,hist,yerr=np.sqrt(hist),fmt='o',label='Beg.='+str(i)+' End='+str(j))
                plt.plot(np.linspace(min(fitbins),max(fitbins),200),gauss(np.linspace(min(fitbins),max(fitbins),200),*pars),'r',label=r'$\chi^2/DoF$= %0.1f'%(chisq)+'\n'+r'$\sigma =$'+str(pars[2]))
                #plt.plot(np.arange(1000,5000),gauss(np.arange(1000,5000),0,pars[1],1,*pars[3:]))
                plt.plot(pars[1],hist[start+amax],'go')
                plt.legend(fontsize=25)
                #plt.yscale('log')
                plt.tick_params(labelsize=25)
                plt.xlim((2100,2300))
                pdf.savefig()
                plt.close()
        print('gauss ',mini,minj,minchi)

if fittype==1:
    with PdfPages(fil[:-12]+'-offgauss.pdf') as pdf:
        mini,minj,minchi=100.,100.,100.
        for i in np.linspace(0.2,1.,9):
            for j in np.linspace(0.2,1.,5):
                beg,end = bins[start+amax]-window*5*i,bins[start+amax]+window*5*j
                fitbins=bins[pd.land(bins>beg,bins<end)]
                fithist = hist[pd.land(bins>beg,bins<end)]
                fithist[fithist==0]=1
                weights=np.sqrt(fithist)
                pars=[hist[start+amax],bins[start+amax],40,1]
                pars=curve_fit(offgauss,fitbins,fithist,p0=pars,sigma=weights,maxfev=120000)[0]
                chisq= np.sum(((offgauss(fitbins,*pars)-fithist)/weights)**2.)
                chisq=chisq/(len(fitbins)-len(pars))
                if chisq<minchi and chisq>0.8:
                    mini,minj,minchi=i,j,chisq
                plt.figure(figsize=(30,20))
                plt.errorbar(bins+0.0,hist,yerr=np.sqrt(hist),fmt='o',label='Beg.='+str(i)+' End='+str(j))
                plt.plot(np.linspace(min(fitbins),max(fitbins),200),offgauss(np.linspace(min(fitbins),max(fitbins),200),*pars),'r',label=r'$\chi^2/DoF$= %0.1f'%(chisq)+'\n'+r'$\sigma =$'+str(pars[2]))
                #plt.plot(np.arange(1000,5000),offgauss(np.arange(1000,5000),0,pars[1],1,*pars[3:]))
                plt.plot(pars[1],hist[start+amax],'go',markersize=24)
                plt.legend(fontsize=25)
                #plt.yscale('log')
                plt.tick_params(labelsize=25)
                plt.xlim((2100,2300))
                pdf.savefig()
                plt.close()
        print('offgauss ',mini,minj,minchi)


if fittype==2:
    with PdfPages(fil[:-12]+'-lingauss.pdf') as pdf:
        mini,minj,minchi=100.,100.,100.
        for i in np.linspace(0.1,1.,10):
            for j in np.linspace(0.1,1.,10):
                print(i,j)
                beg,end = bins[start+amax]-window*5*i,bins[start+amax]+window*5*j
                fitbins=bins[pd.land(bins>beg,bins<end)]
                fithist = hist[pd.land(bins>beg,bins<end)]
                fithist[fithist==0]=1
                weights=np.sqrt(fithist)
                pars=[hist[start+amax],bins[start+amax],30,1,1]
                pars=curve_fit(lingauss,fitbins,fithist,p0=pars,sigma=weights,bounds=([0,2170,0,-np.inf,-np.inf],[np.inf,2200,np.inf,np.inf,np.inf]),maxfev=120000)[0]
                chisq= np.sum(((lingauss(fitbins,*pars)-fithist)/weights)**2.)
                chisq=chisq/(len(fitbins)-len(pars))
                if chisq<minchi and chisq>0.8:
                    mini,minj,minchi=i,j,chisq
                plt.figure(figsize=(30,20))
                plt.errorbar(bins+0.0,hist,yerr=np.sqrt(hist),fmt='o',label='Beg.='+str(i)+' End='+str(j))
                plt.plot(np.linspace(min(fitbins),max(fitbins),200),lingauss(np.linspace(min(fitbins),max(fitbins),200),*pars),'r',label=r'$\chi^2/DoF$= %0.1f'%(chisq)+'\n'+r'$\sigma =$'+str(pars[2]))
                #plt.plot(np.arange(1000,5000),lingauss(np.arange(1000,5000),0,pars[1],1,*pars[3:]))
                plt.plot(pars[1],hist[start+amax],'go',markersize=24)
                plt.legend(fontsize=25)
                #plt.yscale('log')
                plt.tick_params(labelsize=25)
                plt.xlim((2100,2300))
                pdf.savefig()
                plt.close()
        print('lingauss ',mini,minj,minchi)

if fittype==3:
    with PdfPages(fil[:-12]+'-erfgauss.pdf') as pdf:
        mini,minj,minchi=100.,100.,100.
        for i in np.linspace(0.2,1.,9):
            for j in np.linspace(0.2,1.,5):
                beg,end = bins[start+amax]-window*5*i,bins[start+amax]+window*5*j
                fitbins=bins[pd.land(bins>beg,bins<end)]
                fithist = hist[pd.land(bins>beg,bins<end)]
                fithist[fithist==0]=1
                weights=np.sqrt(fithist)
                pars=[hist[start+amax],bins[start+amax],30,1,1]
                pars=curve_fit(erfgauss,fitbins,fithist,p0=pars,sigma=weights,maxfev=120000)[0]
                chisq= np.sum(((erfgauss(fitbins,*pars)-fithist)/weights)**2.)
                chisq=chisq/(len(fitbins)-len(pars))
                if chisq<minchi and chisq>0.8:
                    mini,minj,minchi=i,j,chisq
                plt.figure(figsize=(30,20))
                plt.errorbar(bins+0.0,hist,yerr=np.sqrt(hist),fmt='o',label='Beg.='+str(i)+' End='+str(j))
                plt.plot(np.linspace(min(fitbins),max(fitbins),200),erfgauss(np.linspace(min(fitbins),max(fitbins),200),*pars),'r',label=r'$\chi^2/DoF$= %0.1f'%(chisq)+'\n'+r'$\sigma =$'+str(pars[2]))
                #plt.plot(np.arange(1000,5000),erfgauss(np.arange(1000,5000),0,pars[1],1,*pars[3:]))
                plt.plot(pars[1],hist[start+amax],'go',markersize=24)
                plt.legend(fontsize=25)
#                #plt.yscale('log')
                plt.tick_params(labelsize=25)
                plt.xlim((2100,2300))
                pdf.savefig()
                plt.close()
        print('erfgauss ',mini,minj,minchi)
