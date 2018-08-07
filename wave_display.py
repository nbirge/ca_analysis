import numpy as np
import fileread as fr
import matplotlib.pyplot as plt
import wave_ops as wo
import predefined as pd
from scipy import signal
import sys

def idealpulse(t,*pars):
    t0,amp,decay=pars
    return np.heaviside(t-t0,1.)*amp*np.exp(-(t-t0)/decay)

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
bd,ch=0,6
data=pd.single_pixel(data,bd,ch)
wo.baseline_restore(data,pretrigger=600)

tbins=np.linspace(0,length-1,length)

#Making trap to plot on top of Waveform
rise,top,fall=70,400,250
trapar= np.zeros(len(tbins))
wo.trap(trapar,rise,top,fall)
#

emin,emax=100.,10000.

for i in range(len(data)):
    convolution=signal.fftconvolve(data[i]['wave'],trapar)[0:length]/(rise*int(fall))
    if np.max(data[i]['wave']) > emin and np.max(data[i]['wave'])<emax:
        plt.figure(figsize=(7,5))
        plt.plot(tbins,data['wave'][i,0:length],'b-',label='Raw Waveform')
        mx=np.argmax(data['wave'][i,0:length])
        amp=np.mean(data['wave'][i,mx-10:mx],dtype=float)
        print data['wave'][i,1020:1030]
        plt.plot(tbins,idealpulse(tbins,*[1010,amp,fall]),'r-',label='Ideal electronic response') 
#        plt.plot(tbins,trapar*.1,'r-',label='Convolved shape (scaled by 1/10)')
        convolution=signal.fftconvolve(data[i]['wave'],trapar)[0:length]/(rise*int(fall))
        plt.plot(tbins,convolution[0:len(tbins)],'g-',label='Trapezoid')
        plt.ylabel('ADC Bins')
        plt.xlabel('Timebins (4 ns/bin)')
        plt.legend()
        plt.show()
