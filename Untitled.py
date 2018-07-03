
# coding: utf-8

# In[1]:


import numpy as np
import fileread as fr
import predefined as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib import cm
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import sys

def gauss(x,*pars):
    a,mu,sigma = pars
    return a*np.exp(-(x-mu)**2./(2.*sigma**2.))

def gauss2(x,*pars):
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
fsize=25
rsize=35
if len(sys.argv[:]) < 3:
	print 'Input: file+location column'
	sys.exit()
else:
	fil,col = sys.argv[1:3]
col=int(col)
if col == 5 or col == 7:
	zlabel=r'$\chi ^{2}$/DoF'
else:
	zlabel = r'$\sigma$ (ADC)'
sigs=np.load(fil)
#print sigs.shape
fig = plt.figure(figsize=(30,20))
ax = fig.gca()
x,y=np.concatenate((np.linspace(10,100,10),np.linspace(200,1000,9))),np.concatenate((np.linspace(10,100,10),np.linspace(200,1000,9)))
x,y=np.meshgrid(x,y)
z=np.zeros_like(x)
#np.savetxt('./testing/sigs.txt',sigs, delimiter=',')

for i in range(len(sigs)):
	if col ==4:   
		z[i/19,i%19]= abs(sigs[i,col])
	elif col ==3:
		z[i/19,i%19]= abs((sigs[i,col])*0.157)
	else:
		z[i/19,i%19]= sigs[i,col]
    #print i/5,i%5,sigs[i,4]
z=np.transpose(z)
mins=np.zeros((16,2),dtype=float)
plt.figure(figsize=(40,30))
for i in np.linspace(2,len(z[:,1])-4,len(z[:,1])-3-2,dtype=int):
	if i%3 == 0:
		fmt='-'
	elif i%3 == 1:
		fmt= '--'
	else:
		fmt=':'
	pars=curve_fit(scanfit,x[i,:],z[i,:],p0=[1.,1.,1.],ftol=0.0001)[0]
	plt.plot(x[i,:],z[i,:],'ko')
	chisq=np.sum((z[i,:]-scanfit(x[i,:],*pars))**2.)
	t=np.linspace(min(x[i,:]),max(x[i,:]),100)
	mins[i,:]=(pars[1]/pars[0])**0.5,y[i,i]
	plt.plot(t,scanfit(t,*pars),fmt,label=r'top= %0.0f | $\Sigma r^2=$ %0.1f | min = %0.1f' %(y[i,i],chisq,mins[i,0]),markersize=10)
plt.legend(fontsize=rsize)
plt.tick_params(labelsize=rsize)
plt.ylabel('363 keV Fitted Peak Width (ADC)',fontsize=rsize)
plt.xlabel('Rise (shaping) time (4ns timebins)',fontsize=rsize)
plt.savefig('./testing/residuals')
plt.close()

surf = plt.contourf(x,y,z,18,cmap=cm.coolwarm)
cbar=plt.colorbar()
cbar.ax.tick_params(labelsize='x-large')
cbar.set_label(zlabel,fontsize=fsize,rotation=270, labelpad = 50)
ax.set_xlabel('Rise (shaping) time (4 ns timebins)', fontsize=fsize)
ax.set_ylabel('Top (4 ns timebins)',fontsize=fsize)
#x.set_zlim(np.min(z),np.max(z))
#ax.set_title(,fontsize=fsize)
ax.set_yscale('log')
ax.set_xscale('log')
ax.tick_params(labelsize=fsize)
plt.plot(mins[:,0],mins[:,1],'y*',markersize=15,label='Minimum Shaping Time for fixed top')
#ax.zaxis.set_major_locator(LinearLocator(5))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.01f'))
#fig.colorbar(surf, label=zlabel)
ax.tick_params(labelsize=fsize) 
print fil[:-3]
ax.legend(fontsize=fsize)
plt.savefig('./testing/Sn363contour')
#np.savetxt('./testing/z.txt',z,delimiter = ',')

plt.close()

ebins,erange=1000,[0,5000]
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
plt.savefig('./testing/spectra')
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
	plt.plot(t,gauss2(t,*pars[i,:]),fmt[i]+'--',label='Fit: Rise= %0.0f top= %0.0f' %(200,tops[i]))
plt.yscale('log')
plt.xlabel('Energy (ADC)',fontsize=fsize)
plt.legend(loc='lower left',fontsize=fsize)
plt.title(r'Ca$^{45} +$ Sn$^{113}$ Spectrum',fontsize=fsize)
plt.tick_params(labelsize=fsize)
plt.savefig('spectra')
plt.xlim((2100,2300))
plt.savefig('./testing/spectra+fits')

# In[43]:

'''
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
