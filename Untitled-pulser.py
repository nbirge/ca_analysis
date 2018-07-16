
# coding: utf-8

# In[2]:


import numpy as np
import fileread as fr
import predefined as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import sys


def gauss(x,*pars):
    a,mu,sigma = pars
    return a*np.exp(-(x-mu)**2./(2.*sigma**2.))


# In[3]:

fil=sys.argv[1]
sigs=np.load(fil)
#print sigs.shape
fig=plt.figure(figsize=(20,20))
ax=fig.add_subplot(111,projection='3d')
x,y=np.linspace(100,1000,10),np.linspace(100,1000,10)
x,y=np.meshgrid(x,y)
z=np.zeros_like(x)
for i in range(len(sigs)):
    z[i/10,i%10]=np.log(sigs[i,4])
    #print i/5,i%5,sigs[i,4]
ax.plot_surface(x,y,z,cmap='Spectral')
plt.show()


# In[6]:


'''ebins,erange=1000,[0,5000]
plt.figure(figsize=(20,20))
bd,ch=4,3
sigs=np.zeros((25,6))
count=0
for i in range(100,600,100):
    for j in range(100,600,100):
#i,j = 200,200
        sigs[count,0:2]=i,j
        name = str(i)+'-'+str(j)+'-all.dat'
        fname='./testing/Run_131-'
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
        pars=[hist[start+amax],bins[start+amax],20]
        pars=curve_fit(gauss,fitbins,fithist,p0=pars,ftol=0.001)[0]
        chisq= np.sqrt(np.sum((gauss(fitbins,*pars)-fithist)**2.))
        sigs[count,2:5],sigs[count,5]=pars,chisq/(2*window+1-len(pars))
        count+=1
fig=plt.figure(figsize=(20,20))
ax=fig.add_subplot(111,projection='3d')
x,y=np.linspace(100,500,5),np.linspace(100,500,5)
x,y=np.meshgrid(x,y)
z=np.zeros_like(x)
for i in range(len(sigs)):
    z[i/5,i%5]=sigs[i,4]
ax.plot_surface(x,y,z,cmap='Spectral')
plt.show()
np.savetxt('sigmas.txt',sigs,delimiter=',')


# In[7]:


bd,ch=4,3
fname='./testing/Run_131-'
name='100-500-all.dat'
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
pars=[hist[start+amax],bins[start+amax],20]
pars=curve_fit(gauss,fitbins,fithist,p0=pars,ftol=0.001)[0]
plt.plot(bins,hist,ls='steps')
plt.yscale('log')
plt.xlim((3000,3500))
plt.show()
print np.sum(hist)
#plt.figure(figsize=(20,20))
#plt.plot(fitbins,fithist,ls='steps')
#plt.plot(fitbins,gauss(fitbins,*pars),'r')
#plt.show()
#plt.plot(t,gauss(t,10,2000,1))
#plt.xlim(1990,2010)


# In[11]:


j=200
plt.figure(figsize=(20,20))
for i in range(100,600,300):
    name = str(i)+'-'+str(j)+'-all.dat'
    fname='./testing/Run_131-'
    data=fr.gen_output(fname+name)[0]
    data=pd.precuts(data)
    hist,bins= np.histogram(pd.pixcut(data,'energy',bd,ch),ebins,erange)
    bins=pd.cbins(bins)
    plt.plot(bins,hist,ls='steps',label=name)
plt.xlim((2100,2600))
plt.yscale('log')
plt.legend()
plt.show()
'''
