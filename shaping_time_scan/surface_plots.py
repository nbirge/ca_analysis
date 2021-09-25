
# coding: utf-8

# In[1]:


import numpy as np
import fileread as fr
import predefined as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib import cm
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
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


# ebins,erange=1000,[0,5000]
# bd,ch=4,3
# sigs=np.zeros((100,6))
# count=0
# for i in range(100,1100,100):
#     for j in range(100,1100,100):
# #i,j = 200,200
#         sigs[count,0:2]=i,j
#         name = str(i)+'-'+str(j)+'-all.dat'
#         fname='./Run_131-'
#         data=fr.gen_output(fname+name)[0]
#         data=pd.precuts(data)
#         hist,bins= np.histogram(pd.pixcut(data,'energy',bd,ch),ebins,erange)
#         bins=pd.cbins(bins)
#         start=len(bins)-np.sum(bins>1500)
#         amax= np.argmax(hist[start:])
#         window=60
#         beg,end = bins[start+amax]-window,bins[start+amax]+window
#         fitbins=bins[pd.land(bins>beg,bins<end)]
#         fithist = hist[pd.land(bins>beg,bins<end)]
#         pars=[hist[start+amax],bins[start+amax],20]
#         pars=curve_fit(gauss,fitbins,fithrist,p0=pars,ftol=0.001)[0]
#         chisq= np.sqrt(np.sum((gauss(fitbins,*pars)-fithist)**2.))
#         sigs[count,2:5],sigs[count,5]=pars,chisq/(2*window+1-len(pars))
#         count+=1
# 
# np.save('shaping_time_scan',sigs)

# In[6]:
fsize=20
with PdfPages('/home/noah/Desktop/large_analysis/ca_analysis/shaping_time_scan/testing/sn-shape_scan.pdf') as pdf:
    for pixel in [35]:#]range(48):
        fil='shaping_time_scan-131-'+str(pixel)+'.npy'
        col=4
        if col == 5 or col == 7:
            zlabel=r'$\chi ^{2}$/DoF'
        else:
            zlabel = r'$\sigma$ (ADC)'
        fig = plt.figure(figsize=(30,20))
        ax = fig.gca()
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
        mins=np.zeros((19,2),dtype=float)

        plt.figure(figsize=(30,20),edgecolor='r')
        plt.imshow((z[:,1:]),interpolation='none',origin='lower',cmap=cm.seismic,vmin=0,vmax=30)
#        plt.imshow((z[:,1:]),interpolation='none',origin='lower',cmap=cm.coolwarm)
        plt.xticks(np.arange(18),np.concatenate((np.linspace(20,100,9,dtype=int),np.linspace(200,1000,9,dtype=int))),)
        plt.yticks(np.arange(19),np.concatenate((np.linspace(10,100,10,dtype=int),np.linspace(200,1000,9,dtype=int))))
        cbar=plt.colorbar()
        #surf = plt.contourf(x,y,z,np.linspace(10,50,25),cmap=cm.coolwarm)
        #cbar=plt.colorbar()
        cbar.ax.tick_params(labelsize='xx-large')
        cbar.set_label(zlabel,fontsize=fsize,rotation=270, labelpad = 20)
        plt.xlabel('Rise (shaping) time (4 ns timebins)', fontsize=fsize)
        plt.ylabel('Top (4 ns timebins)',fontsize=fsize)
        #x.set_zlim(np.min(z),np.max(z))
        #ax.set_title(,fontsize=fsize)
        #ax.set_yscale('log')
        #ax.set_xscale('log')
        plt.tick_params(labelsize=fsize)
        plt.title('Pixel '+pd.pixel(pixel/8,pixel%8),fontsize=fsize)
        if col ==4:
            for i in range(19):
                plt.plot(np.argmin(z[i,1:]),i,'k*',markersize=20,label='Minimum Shaping Time for fixed top')
                if i ==0:
                    plt.legend(loc=1,fontsize=fsize)
            #plt.legend(fontsize=fsize)
        else:
            for i in range(19):
                plt.plot(np.argmin(zz[i,1:]),i,'k*',markersize=15,label='Minimum Shaping Time for fixed top')
            plt.legend(fontsize=fsize)
        #ax.zaxis.set_major_locator(LinearLocator(5))
        #ax.zaxis.set_major_formatter(FormatStrFormatter('%.01f'))
        #fig.colorbar(surf, label=zlabel)
        #ax.tick_params(labelsize=fsize) 
        print fil[:-3]
        #ax.legend(fontsize=fsize)
        if col == 4:
            nm='-sigma'
        elif col ==5 or col ==7:
            nm='-chisq'
        pdf.savefig(bbox_inches='tight')
        plt.close('all')

'''
plt.figure(figsize=(30,20))
for top in [30,70,200]:
    x=fr.gen_output('./testing/Run_131-200-'+str(top)+'-all.dat')[0]
    x=pd.precuts(x)
    x=pd.single_pixel(x,board=4,channel=3)
    plt.hist(x['energy'],ebins,erange,histtype='step',label='Rise= %0.0f top= %0.0f' %(200,top))
plt.yscale('log')
plt.xlabel('Energy (ADC)',fontsize=fsize)
plt.xlim((100,2450))
plt.legend(loc='lower left',fontsize=fsize)
plt.title(r'Ca$^{45} +$ Sn$^{113}$ Spectrum',fontsize=fsize)
plt.tick_params(labelsize=fsize)
plt.savefig('./testing/spectra',bbox_inches='tight',format='pdf')
plt.close()


pars=np.zeros((3,5))

pars[0,:]=5571.9994187055,2178.7152448736,12.6063814266,-8.1357325576,18141.1958720878
pars[1,:]=5218.4396403656,2183.1811473664,13.2081262183,-7.9260969042,17735.7998353576
pars[2,:]=4555.0355465584,2193.0288492362,14.8697686843,-7.729296634,17475.0722229956

plt.figure(figsize=(30,20))
fmt=['C0','C1','C2']
tops = [30,70,200]
for i in range(3):
    x=fr.gen_output('./testing/Run_131-200-'+str(tops[i])+'-all.dat')[0]
    x=pd.precuts(x)
    x=pd.single_pixel(x,board=4,channel=3)
    hist,bins=np.histogram(x['energy'],ebins,erange)
    bins=pd.cbins(bins)
    t=np.linspace(2140,2250,501)
    plt.plot(bins,hist,fmt[i],ls='steps',label='Rise= %0.0f top= %0.0f' %(200,tops[i]))
    plt.plot(t,gauss(t,*pars[i,:]),fmt[i]+'--',label='Fit: Rise= %0.0f top= %0.0f' %(200,tops[i]))
plt.yscale('log')
plt.xlabel('Energy (ADC)',fontsize=fsize)
plt.legend(loc='lower left',fontsize=fsize)
plt.title(r'Ca$^{45} +$ Sn$^{113}$ Spectrum',fontsize=fsize)
plt.tick_params(labelsize=fsize)
plt.savefig('spectra')
plt.xlim((2100,2300))
plt.savefig('./testing/spectra+fits')

# In[43]:


best=np.argmin(sigs[:,4])
bd,ch=4,3
fname='./testing/Run_131-'
for i in range(100,600,100):
    name=str(i)+'-'+str(i)+'-all.dat'
    data=fr.gen_output(fname+name)[0]
    data=pd.precuts(data)
    hist,bins= np.histogram(pd.pixcut(data,'energy',bd,ch),ebins,erange)
    bins=pd.cbins(bins)
    start=len(bins)-np.sum(bins>1500)
    amax= np.argmax(hist[start:])
    window=60
    beg,end = bins[start+amax]-window,bins[start+amax]+window
    fitbins=bins[pd.land(bins>beg,bins<end)]
    t=np.linspace(fitbins[0],fitbins[-1],4*window)
    fithist = hist[pd.land(bins>beg,bins<end)]
    pars=[hist[start+amax],bins[start+amax],20]
    pars=curve_fit(gauss,fitbins,fithist,p0=pars,ftol=0.001)[0]
    plt.figure(figsize=(20,20))
    plt.plot(fitbins,fithist,ls='steps',label=name)
    plt.plot(t,gauss(t,*pars),'r',label=r'$\mu$='+str(pars[1])+r',$\sigma$='+str(pars[2]))
    plt.legend()
    #plt.plot(fitbins,pars[-2]*fitbins+pars[-1],'g-')
    plt.savefig('./testing/plots/Run_131-'+name[:-4]+'.pdf')
print pars
#plt.plot(t,gauss(t,10,2000,1))
#plt.xlim(1990,2010)


# In[42]:


j=200
plt.figure(figsize=(20,20))
for i in range(100,600,100):
    name = str(i)+'-'+str(i)+'-all.dat'
    fname='./testing/Run_131-'
    data=fr.gen_output(fname+name)[0]
    data=pd.precuts(data)
    hist,bins= np.histogram(pd.pixcut(data,'energy',bd,ch),ebins,erange)
    bins=pd.cbins(bins)
    plt.plot(bins,hist,ls='steps',label=name)
plt.xlim((2000,2600))
plt.yscale('log')
plt.legend()
plt.show()
'''
