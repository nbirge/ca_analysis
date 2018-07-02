import matplotlib.pyplot as plt
import numpy as np
import fileread as fr
import predefined as pd
bd,ch=4,3
ebins,erange=1000,[0,5000]
rbins,rrange=100,[0,100]
x=fr.gen_output('./testing/Run_131-200-70-all.dat')[0]
x=pd.precuts(x)
x=pd.single_pixel(x,board=bd,channel=ch)
fsize=25

plt.figure(figsize=(20,20))
for i in np.linspace(0,2000,5,dtype=int):
	beg,end=i+0.,i+500.
	hist,bins=np.histogram(x['risetime'][pd.land(x['energy']>beg,x['energy']<end)],rbins,rrange)
	bins=pd.cbins(bins)
	plt.plot(bins,hist,ls='steps',label='%0.0d < E (ADC) < %0.0d' %(beg,end) )
plt.legend(fontsize=fsize)
plt.title(pd.pixel(bd,ch),fontsize=fsize)
plt.tick_params(labelsize=fsize)
plt.xlabel('Risetime (4ns timebins)',fontsize=fsize)
plt.xlim((0,100))
plt.yscale('log')
plt.savefig('./testing/rises'+pd.pixel(bd,ch))
	
