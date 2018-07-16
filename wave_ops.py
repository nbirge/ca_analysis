import numpy as np
from scipy.optimize import curve_fit
from scipy import signal

#Average fall times for run 131c pixel by pixel
#means=[1000, 1031.3367, 1086.8575, 1217.0291, 1041.5563, 1000, 1230.2096, 1188.8999, 1000, 1263.1642, 1233.1743, 1056.3289, 1213.4717, 1112.0769, 1049.4534, 1219.0482, 1000, 1000, 1077.4932, 1157.1627, 1000, 1163.2235, 1000, 1000, 1000, 1027.103, 1111.1212, 1033.5468, 1109.469, 1022.693, 1929.7336, 1000, 1000, 1124.478, 1073.1306, 1040.2197, 1100.4457, 1045.0566, 1135.8975, 1073.1854, 1000, 1000, 1087.187, 1133.1069, 1005.3494, 1000, 1000, 1000] # UNCOMMENT THIS FOR PRODUCTION DATA

means=np.ones(48)*250.

def baseline_restore(wave,pretrigger):
    '''Bitwise & operation and range restore (-=16384) with a baseline restoration \n Use: >> baseline_restore(wavetorestorebaseline,pretrig) '''
    numwaves = len(wave['wave'])
    wave['wave']=wave['wave']&16383
    wave['wave'][wave['wave']>8192]-=16384
    wave['wave']=np.subtract(wave['wave'][0:numwaves],np.mean(wave['wave'][0:numwaves,0:pretrigger],axis=1).reshape((numwaves,1)))


def trap(arr,rise,top,fall):
    '''A trapezoid filter convolving function is stored in arr with parmeters rise,top,fall \n Use >> trap(arr,rise,top,fall) '''
    r,t,d = int(rise),int(top),int(fall)                                                                                                                                                                 
    length=len(arr)     
    x=np.arange(length)               
    arr[0:r]= fall+x[0:r]         
    arr[r+t:r+r+t]=r-d-x[0:r]
    arr[r:r+t]=r                          
    arr[r+r+t:length] = 0

def maxes(waves,startpoint,wavelength,maxamps,maxlocs):
    '''Maxamps and maxlocs should have the shape (1,numwaves). Each entry will contain the max(amplitude/time) for a given wave in waves'''
    numwaves = len(waves)
    maxamps[0:numwaves] =  np.amax(waves[0:len(waves),startpoint:wavelength],axis=1)
    maxlocs[0:numwaves] = np.argmax(waves[0:len(waves),startpoint:wavelength],axis=1)+startpoint


def rises(waves,maxamps,maxlocs,rises):
    t=np.arange(len(waves[0]))
    rises[0:len(rises)]=-1.
    for i in range(len(waves)):
        t1=np.logical_and(waves[i][int(maxlocs[i]-100):int(maxlocs[i])]<0.1*maxamps[i],waves[i][int(maxlocs[i]+1-100):int(maxlocs[i]+1)]>0.1* maxamps[i])
        t2=np.logical_and(waves[i][int(maxlocs[i]-100):int(maxlocs[i])]<0.9*maxamps[i],waves[i][int(maxlocs[i]+1-100):int(maxlocs[i]+1)]>0.9* maxamps[i])
        if np.sum(t1)>0 and np.sum(t2)>0:
            rises[i]=t[int(maxlocs[i]-100):int(maxlocs[i])][t2][0]-t[int(maxlocs[i]-100):int(maxlocs[i])][t1][-1]

def trap_energy(traps,length,output):
    '''Will return an array of trap energies'''
    numwaves=len(traps)
    t=np.arange(length)
    maxamps=np.amax(traps[0:numwaves,0:length],axis=1).reshape((numwaves,1))
    t1=np.logical_and(traps[0:numwaves,0:length-1]-0.8*maxamps<0,traps[0:numwaves,1:length]-0.8*maxamps>0)    #Here should be a flag to tag multiple crappy t1's
    t2=np.logical_and(traps[0:numwaves,0:length-1]-0.8*maxamps>0,traps[0:numwaves,1:length]-0.8*maxamps<0)
    #Here could be that check (not sure how to do it without a for loop...
    #Fix shit this
    for i in range(numwaves):
        if np.sum(t1[i])>=1 and np.sum(t2[i])>=1:
            output[i]=traps[i][(t[0:length-1][t2[i]][0]-t[0:length-1][t1[i]][0])/2+t[0:length-1][t1[i]][0]]
        else:
            output[i]=-1.

def tail_fit(data,output):
    '''Data should be smoothed beforehand. AND data is column of waveforms (not structured array)'''
    length= len(data[0])
    t= np.arange(length)
    fitpars=np.zeros(3)
    for i in range(len(data)):
        maxbin=np.argmax(data[i])
        if maxbin > length-800:
            maxbin=1800
        tail = lambda t,a,b: a*np.exp(b*t[maxbin+200:length])
        fitpars = [data[i][maxbin],-1./1000]
        fitpars = curve_fit(tail,t,data[i][maxbin+200:length],p0=fitpars,ftol=0.1)[0]
        output[i]=1./fitpars[1]

#OLD taifit code
#    length = len(data[0])
#    t = np.arange(length)
#    line = lambda t,a,b: a+b*t[1200:2000]
#    guess = [0.,1000.]
#    for i in range(len(data)):
#        if np.any(data[i][1200:2000]<0.):
#            output[i] = np.array([-0,0],dtype='f')[0]
#        else:
    #        output[i]=1./curve_fit(line,t,np.log(data[i][1200:2000]),p0=(np.log(np.amax(data[i])),1000.))[0][1]



def pileup(data,workarr,thresh):
    '''Data should be smoothed beforehand'''
    length=len(data[0])
    for i in range(len(data)):
        workarr[i]=np.sum(data[i][signal.argrelmax(data[i,0:length],order=20)[0]]>thresh)

def pixel_traps(workarr,rise,top):
    '''48xLENGTH array of traps for each individual pixel | taken via mean of distribution of falltimes for run 137'''
    length=len(workarr[0])
    for i in range(len(means)):
        means[i]=int(means[i])
    for i in range(len(means)):
        trap(workarr[i],rise,top,means[i])


def apply_trap(rise,data,trap,output):
    numwaves=len(data)
    length=len(data[0]['wave'])
    for i in range(len(data)):
        bdch=data[i]['board']*8+data[i]['channel']
        output[i][0:length]=(signal.fftconvolve(data[i]['wave'], trap[bdch], mode='full')[0:length])/(rise*means[bdch])


def fitted_trap(data,rise,top,fall,output):
    #fall should be an array of length numwaves
    numwaves=len(data)
    length=len(data[0]['wave'])
    trp = np.zeros(length)
    for i in range(numwaves):
        trap(trp,rise,top,int(fall[i]))
        output[i][0:length]=signal.fftconvolve(data[i]['wave'],trp)[0:length]/(rise*int(fall[i]))

def find_t0(data,output):
    length = len(data)
    output[0:length]=np.argmax(data['wave'],axis=1)

    























