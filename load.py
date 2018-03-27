import numpy as np
import os


def help():
	print 'To use load, just call the load function. Ex: load.load(/path/to/run_84.fin)'
	print 'This program returns the original data array and a corresponding energy array'
	print 'The data are arranged in columns as follows: \n Result	EvID	Board	Channel	Timestamp	tau1 (falltime)	tau2 (risetime)	V0	t0	Residual	Energy'
	return

def load0(path):
#For older format
	size = os.stat(path).st_size/41
	data =[]
	with open(path) as f:
#		f.seek(8)
		data = np.core.records.fromfile(f,formats='B,i,i,i,Q,f,f,f,f,f',shape=size,names='Result,EvID,Board,Channel,Timestamp,tau1,tau2,v0,t0,Residual',byteorder='<')
		energy = data['v0']*((data['tau1']/data['tau2'])**(-1.*data['tau2']/(data['tau1']-data['tau2']))-(data['tau1']/data['tau2'])**(-1.*data['tau1']/(data['tau1']-data['tau2'])))
	return data,energy

def load1(path):
#May 2017 format
        size = os.stat(path).st_size/(33+4*4)
        data =[]
        with open(path) as f:
#               f.seek(8)
                data = np.core.records.fromfile(f,formats='B,i,i,i,Q,Q,f,f,f,f,f',shape=size,names='Result,EvID,Board,Channel,Timestamp,Writetime,tau1,tau2,v0,t0,Residual',byteorder='<')
                energy = data['v0']*((data['tau1']/data['tau2'])**(-1.*data['tau2']/(data['tau1']-data['tau2']))-(data['tau1']/data['tau2'])**(-1.*data['tau1']/(data['tau1']-data['tau2'])))
        return data,energy

def load2(path):
#May 2017 format - includes entry number
        size = os.stat(path).st_size/(33+4+4+4*4) # includes part and count now
        data =[]
        with open(path) as f:
#               f.seek(8)
                data = np.core.records.fromfile(f,formats='B,i,i,i,Q,Q,i,i,f,f,f,f,f',shape=size,names='result,evID,board,channel,timestamp,writetime,part,count,tau1,tau2,v0,t0,Residual',byteorder='<')
                energy = data['v0']*((data['tau1']/data['tau2'])**(-1.*data['tau2']/(data['tau1']-data['tau2']))-(data['tau1']/data['tau2'])**(-1.*data['tau1']/(data['tau1']-data['tau2'])))
        return data,energy

def loadtrap(path):
#	dline = struct.pack('<BiiiQQiiiffff',result,eventid,board,channel,timestamp,writetimestamp,part,count,risetime,maxguess,var,sol.x[1],wav[s])
	size= os.stat(path).st_size/(33-4+4*7)
	data=[]
	with open(path) as f:
		data = np.core.records.fromfile(f,formats='B,i,i,i,Q,Q,i,i,i,f,f,f,f',shape=size,names='result,evID,board,channel,timestamp,writetime,part,count, risetime,max,variance,falltime,trap',byteorder='<')
	return data

def loadtrapnfit(path):
#dline = struct.pack('<BiiiQQii%sf' %len(set1),result,eventid,board,channel,timestamp,writetimestamp,part,count,*set1) len(set1) should be 10
#set1 = [risetime,Amp of Trap,fitted fall of trap, chi^2/dof,tau1,tau2,v0,t0,chi^2/dof,variance of pretrigger,]
	size = os.stat(path).st_size/(1+3*4+2*8+2*4+4*10)
	data=[]
	with open(path) as f:
		data= np.core.records.fromfile(f,formats='B,i,i,i,Q,Q,i,i,f,f,f,f,f,f,f,f,f,f',shape=size, names='result,evID,board,channel,timestamp,writetime,part,count,rise,trap,tail,tchi,tau1,tau2,v0,t0,echi,var',byteorder='<')
		energy = data['v0']*((data['tau1']/data['tau2'])**(-1.*data['tau2']/(data['tau1']-data['tau2']))-(data['tau1']/data['tau2'])**(-1.*data['tau1']/(data['tau1']-data['tau2'])))
	return data,energy

def loadtrapnfit1(path):
#dline = struct.pack('<BiiiQQii%sf' %len(set1),result,eventid,board,channel,timestamp,writetimestamp,part,count,*set1) len(set1) should be 10
#set1 = [risetime,Amp of Trap,fitted fall of trap, chi^2/dof,tau1,tau2,v0,t0,chi^2/dof,variance of pretrigger,]
	size = os.stat(path).st_size/(1+3*4+2*8+2*4+4*11)
	data=[]
	with open(path) as f:
		data= np.core.records.fromfile(f,formats='B,i,i,i,Q,Q,i,i,f,f,f,f,f,f,f,f,f,f,f',shape=size, names='result,evID,board,channel,timestamp,writetime,part,count,rise,trap,tail,tchi,midbin,tau1,tau2,v0,t0,echi,var',byteorder='<')
#		energy = data['v0']*((data['tau1']/data['tau2'])**(-1.*data['tau2']/(data['tau1']-data['tau2']))-(data['tau1']/data['tau2'])**(-1.*data['tau1']/(data['tau1']-data['tau2'])))
#	return data,energy
	return data

