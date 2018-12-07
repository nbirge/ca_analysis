import numpy as np
import fileread as fr
import predefined as pd
import wave_ops as wo
import sys,os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

run = sys.argv[1]
lastpart = sys.argv[2]
timestampfile = sys.argv[3]
inpath = sys.argv[4]
portion = sys.argv[5]
length = 3500 

if not inpath.endswith('/'):
    inpath+='/'


dt=[('result', 'u1'), ('evID', '<i4'), ('board', '<i4'), ('channel', '<i4'), ('timestamp', '<u8'), ('requesttime', '<u8'), ('risetime', '<i4'), ('energy', '<f4'), ('falltime', '<f4'), ('t0', '<i4')]

timestamps = np.fromfile(timestampfile,dtype=dt,count = os.stat(timestampfile).st_size/1+12+16+16)

print(len(timestamps))
chunk = 10000
temp= fr.raw(inpath+'Run_'+run+'_0.bin',length=length,row=0,numwaves=1)
temp=np.zeros(len(timestamps),dtype=temp.dtype)
count = 0
for i in range(int(lastpart)+1):
    name = 'Run_'+run+'_'+str(i)+'.bin'
    print('Checking '+name+' for timestamps')
    numwaves= (os.stat(inpath+name).st_size-8)/(1+12+16+4+2*length)
    for i in range(numwaves/chunk):
        rem = 0
        if i == numwaves/chunk:
            rem = numwaves%chunk
        tot = chunk+rem
        data = fr.raw(inpath+name,length=length,numwaves=tot,row = i*chunk)
        for j in range(len(timestamps)):
            x= data[pd.land(pd.land(data['timestamp']==timestamps[j]['timestamp'],data['board']==timestamps[j]['board']),data['channel']==timestamps[j]['channel'])]
            count+=len(x)
            if len(x) >0:
#                print x['board'],x['channel'],x['timestamp']
                temp[count-1]=x
print(temp['board'])
wo.baseline_restore(temp,600)
tbins=np.arange(3500)
#for i in range(len(temp)):
with PdfPages('fallwaves.pdf') as pdf:
    for j in range(int(portion),len(timestamps)):
    #        print temp[temp['timestamp']==timestamps[j]['timestamp']]['wave'][0]
        if timestamps[j]['board']>2:
            plt.plot(tbins,temp[temp['timestamp']==timestamps[j]['timestamp']]['wave'][0], label='falltime = '+str(timestamps[j]['falltime'])+'\n '+pd.pixel(timestamps[j]['board'],timestamps[j]['channel']))
            plt.legend(fontsize=20)
            plt.tight_layout()
            pdf.savefig()
            plt.close()
        
