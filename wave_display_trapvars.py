import numpy as np
import fileread as fr
import matplotlib.pyplot as plt
import wave_ops as wo
import predefined as pd
from scipy import signal
import sys

if len(sys.argv[:]) < 6:
	print 'Input: run part path numwaves startrow'
	sys.exit()
else:
	run,part,path,numwaves,startrow = sys.argv[1:6]
if not path.endswith('/'):
	path+='/'
numwaves=int(numwaves)
startrow=int(startrow)
length=3500

#for i in range(4,len(sys.argv[4:])):		#PARSE ARGS LIKE IN ca_analysis.py

name= 'Run_'+str(run)+'_'+str(part)+'.bin'

data = fr.raw(path+name,length=length,numwaves=numwaves,row=startrow)
bd,ch=4,3
data=pd.single_pixel(data,bd,ch)
wo.baseline_restore(data,pretrigger=600)

tbins=np.arange(length)

#Making trap to plot on top of Waveform
rise,top,fall=300,200,1073
trapar= np.zeros(len(tbins))
wo.trap(trapar,rise,top,fall)
#

for i in range(len(data)):
	if i%2 == 0:	
		i+=1
	for rise in range(100,1100,300):
		wo.trap(trapar,rise,rise,fall)
		convolution=signal.fftconvolve(data[i]['wave'],trapar)[0:length]/(rise*int(fall))
		if np.max(data[i]['wave'])>3000:
			plt.plot(tbins,data['wave'][i,0:length],'b-',label='Raw Wave')
			plt.plot(tbins,trapar*1.,'r-',label='Convolved shape')
			convolution=signal.fftconvolve(data[i]['wave'],trapar)[0:length]/(rise*int(fall))
			plt.plot(tbins,convolution[0:len(tbins)],'g-',label='Trapezoid')
			plt.title('Rise, top = '+str(rise))
			plt.legend()
			plt.show()
