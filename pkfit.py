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

savehist=1


ebins,erange=4000,[0,20000]
for pixel in range(48):
    bd,ch=pixel/8,pixel%8
    print 'BOARD = '+str(bd)+' CHANNEL = '+str(ch)
    sigs=np.zeros((361,8))
    psigs = np.zeros((361,8))
    count=0
    lst = np.concatenate((np.linspace(10,100,10,dtype=int),np.linspace(200,1000,9,dtype=int)))
    path='./run48/'
    for i in lst:
        for j in lst:
            if i > 00 or j > 00:
                pars=[0,0,0]
                sigs[count,0:2]=i,j
                name = str(i)+'-'+str(j)+'-all.dat'
                fname='Run_48-'
                data=fr.gen_output(path+fname+name)[0]
                #data=pd.precuts(data) #uncomment this line for nonpulser data
                hist,bins= np.histogram(pd.pixcut(data,'energy',bd,ch),ebins,erange)
                bins=pd.cbins(bins)
                start=len(bins)-np.sum(bins>1000)
                amax= np.argmax(hist[start:])
                window=12
                beg,end = bins[start+amax]-window*5*0.5,bins[start+amax]+window*5*0.5
                fitbins=bins[pd.land(bins>beg,bins<end)]
                fithist = hist[pd.land(bins>beg,bins<end)]
                weights=np.sqrt(fithist)
                weights[weights==0]=1
                pars=[hist[start+amax],bins[start+amax],30,1,1]
                pars=curve_fit(lingauss,fitbins,fithist,p0=pars,sigma=weights,bounds=([0,bins[start+amax]-30,0,-np.inf,-np.inf],[np.inf,bins[start+amax]+30,np.inf,np.inf,np.inf]),max_nfev=120000)[0]
                chisq= np.sum(((lingauss(fitbins,*pars)-fithist)/weights)**2.)
                sigs[count,2:7],sigs[count,7]=pars,chisq/(len(fitbins)-len(pars))

                if savehist==1 and i ==400 and j==70:
                    x=np.zeros((len(hist),2))
                    x[:,0]=hist
                    x[:,1]=bins
                    np.save(path+str(i)+'-'+str(j)+'-48-histbins-'+str(pixel),x)
                if count%20 ==0:
                    print i,j,bd,ch,'\n',sigs[count,2:8]
                count+=1
    print sigs
    np.save('./run48/shaping_time_scan-48-'+str(pixel),sigs)

'''

            data=fr.gen_output(fname+name)[0]
            data=data[data['energy']>2800]
            hist,bins= np.histogram(pd.pixcut(data,'energy',bd,ch),ebins,erange)
            bins=pd.cbins(bins)
            start=len(bins)-np.sum(bins>1000)
            amax= np.argmax(hist[start:])
            beg,end = bins[start+amax]-window,bins[start+amax]+window
            fitbins=bins[pd.land(bins>beg,bins<end)]
            fithist = hist[pd.land(bins>beg,bins<end)]
            weights=np.sqrt(fithist)
            weights[weights==0]=1
            pars=[hist[start+amax],bins[start+amax],30,1,1]
            pars=curve_fit(lingauss,fitbins,fithist,p0=pars,sigma=weights,maxfev=120000)[0]
            chisq= np.sum(((lingauss(fitbins,*pars)-fithist)/weights)**2.)
            psigs[count,2:7],psigs[count,7]=pars,chisq/(len(fitbins)-len(pars))

            if savehist==1:
                x=np.zeros((len(hist),2))
                x[:,0]=hist
                x[:,1]=bins
                np.save(str(i)+'-'+str(j)+'-pulser-histbins',x)

            print i,j,sigs[count,2:8],count
        count+=1

np.save('shaping_time_scan-363pk-lingauss',sigs)
np.save('shaping_time_scan-pulser-lingauss',psigs)'''
