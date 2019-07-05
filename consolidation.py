import fileread as fr
import numpy as np
import sys,os

def land(x,y):
    return np.logical_and(x,y)

path,run = str(sys.argv[1]),str(sys.argv[2])
#path='./'
beg= int(sys.argv[1])
if len(sys.argv[:])>2:
    end = int(sys.argv[2])+1
else:
    end = beg+1
runs = np.arange(beg,end,1))

for i in runs:
    path='./'+str(i)+'/'
    print('combining parts of run: '+str(i))
    run=str(i)
    name='Run_'+run
    x=[x for x in os.listdir(path) if x.startswith('Run_'+str(run)) and x.endswith('-comb.bin') and x!='Run_'+str(run)+'_0-comb.bin']
    x=sorted(x)
    print(x)
    header= np.zeros(1,dtype=[('theader','Q'),('formats','10i')])

    data,header['theader'],header['formats'] = fr.gen_output(path+name+'_0-comb.bin')
    for fyle in x:
        print(fyle)
        if os.stat(path+fyle).st_size>1000:
            data=np.concatenate((fr.gen_output(path+fyle)[0],data))
#
#       CORRUPTION REJECTION
    print('Removing event corruption for board:')
    data=np.sort(data,order='timestamp')
    cut=float(len(data))
    y=[]
    size=np.zeros(6,dtype=int)
    for j in range(6):
        print(j)
        bd=j
        x=data[data['board']==bd]
        t=np.zeros(len(x),dtype=bool)
        counter=0
        truth1=np.zeros(len(x),dtype=bool)
        truth2=np.zeros(len(x),dtype=bool)
        for i in range(len(x)):
            if i<99:
                beg = i
            xx=x[i-beg:i+99]    #Maybe just assign xx=x[blah]['requesttime']
            truth1=land(xx['requesttime']<x['timestamp'][i],x['timestamp'][i]<xx['requesttime']+60000/4) # I think that this is a 60 us = 60,000/ us/4ns timebins
            truth2 = land(xx['requesttime']<x['timestamp'][i]+3500,x['timestamp'][i]<xx['requesttime']+60000/4)
            t[i] = not np.any(land(truth1,truth2))
        y.append(x[t])
        size[j]=len(y[j])
    data=np.zeros(np.sum(size),dtype=data.dtype)
    for i in range(len(size)):
        if i ==0:
            data[0:size[0]]=y[0]
        elif i>0:
            data[int(np.sum(size[0:i])):np.sum(size[0:i])+size[i]]=y[i]
    data=np.sort(data,order='timestamp')
    print('Removed %0.2f %% of %0.0f' %((1.-len(data)/cut)*100.,cut))

    output=path+'Run_'+run+'-all.dat'
    print('\n'+output)
    with open(path+'Run_'+run+'-all.dat','wb') as f:
        header.tofile(f)
        data.tofile(f)









