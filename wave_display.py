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
rise,top,fall=100,250,1100
trapar= np.zeros(len(tbins))
wo.trap(trapar,rise,top,fall)
#

ecut=3000

for i in range(len(data)):
	convolution=signal.fftconvolve(data[i]['wave'],trapar)[0:length]/(rise*int(fall))
	if np.max(data[i]['wave']) > ecut:
		plt.plot(tbins,data['wave'][i,0:length],'b-',label='Raw Wave')
		plt.plot(tbins,trapar*.1,'r-',label='Convolved shape (scaled by 1/10)')
		convolution=signal.fftconvolve(data[i]['wave'],trapar)[0:length]/(rise*int(fall))
		plt.plot(tbins,convolution[0:len(tbins)],'g-',label='Trapezoid')
		plt.legend()
    Title = str('Raw_Wave'+'_'+path[39:43]+'_'+path[44:49]+'_'+str(run)+'_'+str(part)+'_'+str(int(startrow+i))+'_'+str(startrow))
    plt.plot(tbins,data['wave'][i,0:length],'b-',label=Title)
    plt.xlabel('Time /(4 ns)')
    plt.ylabel('Arbitrary Units')
    plt.savefig('plots/wave_display/'+str(Title)+'.png')
    plt.show()
