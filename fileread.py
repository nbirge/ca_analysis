import numpy as np
import os


def trapnfit(path):
    '''This will load the output (generally -all.fin files) of the trap and fit processing and return a structured numpy array'''
#dline = struct.pack('<BiiiQQii%sf' %len(set1),result,eventid,board,channel,timestamp,writetimestamp,part,count,*set1) len(set1) should be 10
#set1 = [risetime,Amp of Trap,fitted fall of trap, chi^2/dof,tau1,tau2,v0,t0,chi^2/dof,variance of pretrigger,]
    size = os.stat(path).st_size/(1+3*4+2*8+2*4+4*11)
    data=[]
    with open(path) as f:
        data= np.core.records.fromfile(f,formats='B,i,i,i,Q,Q,i,i,f,f,f,f,f,f,f,f,f,f,f',shape=size, names='result,evID,board,channel,timestamp,writetime,part,count,rise,trap,tail,tchi,midbin,tau1,tau2,v0,t0,echi,var',byteorder='<')
#        energy = data['v0']*((data['tau1']/data['tau2'])**(-1.*data['tau2']/(data['tau1']-data['tau2']))-(data['tau1']/data['tau2'])**(-1.*data['tau1']/(data['tau1']-data['tau2'])))
#    return data,energy
    return data

def raw(path,length=3500,numwaves=10000,row=1000):
    '''This will read numwaves number of waveforms of length length into a structured numpy array *** Starts with row 0!'''
    fsize=os.stat(path).st_size
    totrows = (fsize-8)/(33+length*2)
    data=[]
    pos = 8+row*(33+length*2)
    with open(path,'rb') as f:
        f.seek(pos)
        data= np.core.records.fromfile(f, formats= 'B,i4,i4,i4,Q,Q,i4,(3500)i2',shape=numwaves,\
            names='result,evID,board,channel,timestamp,requesttime,length,wave',byteorder='<')
    return data

def analysis_output(fname):
    with open(fname) as f:
        data=[]
        formats='B,i,i,i,Q,Q,i,f,f'
        names='result,evID,board,channel,timestamp,requesttime,risetime,energy,falltime'
        numwaves = os.stat(fname).st_size/(41)
        data=np.core.records.fromfile(f,formats=formats,shape=numwaves,names=names,byteorder='<')
    return data

def gen_output(fname):
    with open(fname,'rb') as f:
        data=[]
        formatting = 'B,i,i,i,Q,Q,i,f,f'
        names='result,evID,board,channel,timestamp,requesttime,risetime,energy,pretrigrms'
        file_timestamp=np.fromfile(f,dtype='Q',count=1)[0]
        f.seek(8)
        formats = np.fromfile(f,dtype='10i',count=1)[0]
        data_byte_size=(1+3*4+2*8+4+4+4)
        if formats[0]==1:            #For sets with fall times calculated
            data_byte_size+=4
            formatting+=',f'
            names+=',falltime'
        if formats[1] >= 1:            #For data sets where pileup is determined
            data_byte_size+=4+4+4
            formatting+=',i,i,f'
            names+=',pileup,pilediff,pileamp'
        if formats[2] == 1:
            data_byte_size+=4
            formatting+=',f'
            names+=',fitenergy'
        if formats[3] == 1:
            data_byte_size+=4
            formatting+=',i'
            names+=',t0'
        if formats[4] >0:
            data_byte_size+=int(4*formats[4])
            formatting+=',4f,4f,f,f'
            names+=',osc_amps,osc_errors,s2,osc_energy'
        if formats[5] >0:
            data_byte_size+=int(3*4)
            formatting+=',i,i,i'
            names+=',pfilter0,pfilter1,pfilter2'
        numwaves=int((os.stat(fname).st_size-(8+4*10))/data_byte_size)
        f.seek(8+4*10)
        data=np.core.records.fromfile(f,formats=formatting,shape=numwaves,names=names,byteorder='<')
        return data,file_timestamp,formats


def simulation(fname):
    dtype={'names':('entry','detector','pixel','timestamp','energy'),'formats':('i4','U1','i4','f8','f8')}
    data= np.loadtxt(fname,dtype=dtype,delimiter=' ')
    data['energy']=data['energy']/1.e3
    return data

def initial_simulation(fname):
    dtype={'names':('entry','particle','t0','energy','x','y','z','vx','vy','vz'),'formats':('i4','i4','f8','f8','f8','f8','f8','f8','f8','f8')}
    data= np.loadtxt(fname,dtype=dtype,delimiter=' ')
    data['energy']=data['energy']/1.e3
    return data

def trig(fname):
    with open(fname,'rb') as f:
        formatting='i4,i1,i8,i2'
        names='board,channel,timestamp,energy'
        timestamp=np.fromfile(f,dtype='Q',count=1)[0]
        data_byte_size=4+1+8+2
        numwaves=int((os.stat(fname).st_size-8)/data_byte_size)
        f.seek(8)
        data=np.core.records.fromfile(f,formats=formatting,shape=numwaves,names=names,byteorder='<')
    return data,timestamp

def par(fname):
    dtype={'names':('threshold','falltime','top','rise'),'formats':('>i4','>i8','>i2','>i2')}
    return np.fromfile(fname,dtype=dtype)

'''def file_consolidation(path,runnumber)
    x=filter(lambda x: x.startswith('Run_'+str(runnumber)) and x.endswith('-comb.bin') and x!='Run_'+str(runnumber)+'_0-comb.bin',os.listdir(path))
    name='Run_'+str(runnumber)
    data=[]
    formatting = 'B,i,i,i,Q,Q,i,f'
    names='result,evID,board,channel,timestamp,requesttime,risetime,energy'
    file_timestamp=np.fromfile(path+name,dtype='Q',count=1)[0]
    f.seek(8)
    formats = np.fromfile(f,dtype='10B',count=1)[0]
    data_byte_size=(1+3*4+2*8+4+4)
    if formats[0]==1:            #For sets with fall times calculated
        data_byte_size+=4
        formatting+=',f'
        names+=',falltime'
    if formats[1] == 1:            #For data sets where pileup is determined
        data_byte_size+=4
        formatting+=',i'
        names+=',pileup'
    numwaves=(os.stat(path+name).st_size-(8+1*10))/data_byte_size
    f.seek(8+1*10)
    data=np.core.records.fromfile(f,formats=formatting,shape=numwaves,names=names,byteorder='<')

def temp_gen_output(fname):
    with open(fname,'rb') as f:
        data=[]
        formatting = 'B,i,i,i,Q,Q,i,f'
        names='result,evID,board,channel,timestamp,requesttime,risetime,energy'
        file_timestamp=np.fromfile(f,dtype='Q',count=1)[0]
        f.seek(8)
        formats = np.fromfile(f,dtype='10i',count=1)[0]
        data_byte_size=(1+3*4+2*8+4+4)
        if formats[0]==1:            #For sets with fall times calculated
            data_byte_size+=4
            formatting+=',f'
            names+=',falltime'
        if formats[1] >= 1:            #For data sets where pileup is determined
            data_byte_size+=4+4+4
            formatting+=',i,i,f'
            names+=',pileup,pilediff,pileamp'
        if formats[2] == 1:
            data_byte_size+=4
            formatting+=',f'
            names+=',fitenergy'
        if formats[3] == 1:
            data_byte_size+=4
            formatting+=',i'
            names+=',t0'
        numwaves=(os.stat(fname).st_size-(8+4*10))/data_byte_size
        f.seek(8+4*10)
        data=np.core.records.fromfile(f,formats=formatting,shape=numwaves,names=names,byteorder='<')
        return data,file_timestamp,formats
'''










