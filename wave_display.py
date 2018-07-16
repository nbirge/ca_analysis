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

#for i in range(4,len(sys.argv[4:])):        #PARSE ARGS LIKE IN ca_analysis.py

name= 'Run_'+str(run)+'_'+str(part)+'.bin'

data = fr.raw(path+name,length=length,numwaves=numwaves,row=startrow)
bd,ch=1,3
data=pd.single_pixel(data,bd,ch)
wo.baseline_restore(data,pretrigger=600)

tbins=np.arange(length)

#Making trap to plot on top of Waveform
rise,top,fall=300,100,250
trapar= np.zeros(len(tbins))
wo.trap(trapar,rise,top,fall)
#

emin,emax=0000.,100.

for i in range(len(data)):
    convolution=signal.fftconvolve(data[i]['wave'],trapar)[0:length]/(rise*int(fall))
    if np.max(data[i]['wave']) > emin and np.max(data[i]['wave'])<emax:
        plt.figure(figsize=(7,5))
        plt.plot(tbins,data['wave'][i,0:length],'b-',label='Raw Waveform')
#        plt.plot(tbins,trapar*.1,'r-',label='Convolved shape (scaled by 1/10)')
#        convolution=signal.fftconvolve(data[i]['wave'],trapar)[0:length]/(rise*int(fall))
#        plt.plot(tbins,convolution[0:len(tbins)],'g-',label='Trapezoid')
        plt.ylabel('ADC Bins')
        plt.xlabel('Timebins (4 ns/bin)')
        plt.legend()
        plt.show()