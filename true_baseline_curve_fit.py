from mpi4py import MPI
import numpy as np
import matplotlib.pyplot as plt
import sys,os
sys.path.insert(0, '/nics/d/home/nwbirge/large_analysis/ca_analysis/')
import fileread as fr
from scipy.signal import fftconvolve
import predefined as pd
import wave_ops as wo
from scipy import signal
from scipy.optimize import curve_fit
import time

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def background(t,*pars):
    a,b,offset,omega=pars
    return a*np.sin(omega*t)+b*np.cos(omega*t)+offset*np.ones(len(t))

means=np.array([1000, 1031.3367, 1086.8575, 1217.0291, 1041.5563, 1000, 1230.2096, 1188.8999,                1000, 1263.1642, 1233.1743, 1056.3289, 1213.4717, 1112.0769, 1049.4534, 1219.0482,                1000, 1000, 1077.4932, 1157.1627, 1000, 1163.2235, 1000, 1000,                1000, 1027.103, 1111.1212, 1033.5468, 1109.469, 1022.693, 1929.7336, 1000,                1000, 1124.478, 1073.1306, 1040.2197, 1100.4457, 1045.0566, 1135.8975, 1073.1854,                1000, 1000, 1087.187, 1133.1069, 1005.3494, 1000, 1000, 1000])


if len(sys.argv)==1:
    lowerbound,upperbound=8.E-4,8.E-3
elif len(sys.argv)==2:
    lowerbound,upperbound=float(sys.argv[1]),8.E-3
else:
    lowerbound,upperbound=float(sys.argv[1]),float(sys.argv[2])

runs=['120']
for r in runs:
    if int(r)>99:
        start=-13
    else:
        start=-12
    if int(r)>64:
        path='/lustre/haven/gamma/neutrons/ca45_data/2017/disk1/'
#        path='/home/noah/Desktop/large_analysis/ca_analysis/'
    else:
        path='/lustre/haven/gamma/neutrons/ca45_data/2017/disk0/'

    run=path+'Run_'+r+'_0.bin'
    name=run[start:-6]

    if rank>0:
        print('Starting on '+name)
        chunk=int((os.stat(run).st_size-8)/7033)
        chunk=int(chunk/(size-1))
        row=(rank-1)*chunk

        numwaves=chunk
        if rank == size-1:
            numwaves=chunk+int((os.stat(run).st_size-8)/7033)%(size-1)

        data=fr.raw(run,length=3500,row=row,numwaves=numwaves)
        numwaves=len(data)

        wo.baseline_restore(data,600)
        length=3500

        t=np.arange(length); out=np.zeros(length); trap=np.zeros(length)
        fall=0
        maxouts=0
        omega=2*np.pi/length
        pars=[]
        beg=time.time()
        sys.stdout.flush()
        for i in range(numwaves):
            if i%int(numwaves/100)==0 and rank==1:
                string='Done with {:0.2f}%'.format(100*i/numwaves)
                sys.stdout.write(string)
                sys.stdout.write('\b'*len(string))
                sys.stdout.flush()
            try:
                bd,ch=data['board'][i],data['channel'][i]
                fall=means[data['board'][i]*8+data['channel'][i]]
                wo.trap(arr=trap,rise=100,top=70,fall=fall)
                out=signal.convolve(trap,data['wave'][i])[0:length]/(100*fall)
                maxout=np.amax(out)
                if maxout<50:
                    #p,v=curve_fit(background,t,data['wave'][i],p0=[maxout,maxout,0,omega],bounds=([-100,-100,-300,8E-4],[100,100,300,8E-3]))       #UNCOMMENT THIS TO DO A REGULAR FIT (NOT ESPECIALLY BOUNDED)
                    p,v=curve_fit(background,t,data['wave'][i],p0=[maxout,maxout,0,omega],bounds=([-100,-100,-300,lowerbound],[100,100,300,upperbound]))
                    chisq=np.sum(np.power(background(t,*p)-data['wave'][i],2.))/(length-p.shape[0])
                    pars.append(np.concatenate((p,np.sqrt(np.diag(v)),[chisq,bd,ch,i])));

            except RuntimeError:
                continue

        pars=np.array(pars)
        end=time.time()
        print(name,rank,end-beg)
        np.save('/nics/d/home/nwbirge/baselinefits/True-Fit-pars_'+name+'-'+str(rank)+'.npy',pars)
        comm.send(rank,dest=0)

    if rank==0:
        for check in np.arange(1,size,1):
            print(check == comm.recv(source=check),check)
            sys.stdout.flush()
        data=np.load('/nics/d/home/nwbirge/baselinefits/True-Fit-pars_'+name+'-1.npy')
        os.system('rm /nics/d/home/nwbirge/baselinefits/True-Fit-pars_'+name+'-1.npy')
        for part in np.arange(2,size,1):
            data=np.concatenate((data,np.load('/nics/d/home/nwbirge/baselinefits/True-Fit-pars_'+name+'-'+str(part)+'.npy')))
            os.system('rm /nics/d/home/nwbirge/baselinefits/True-Fit-pars_'+name+'-'+str(part)+'.npy')
        np.save('/nics/d/home/nwbirge/baselinefits/True-Fit-pars_'+name+'.npy',data)
        print('Finished '+name+' woohoo')
