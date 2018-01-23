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
#		energy = data['v0']*((data['tau1']/data['tau2'])**(-1.*data['tau2']/(data['tau1']-data['tau2']))-(data['tau1']/data['tau2'])**(-1.*data['tau1']/(data['tau1']-data['tau2'])))
#	return data,energy
	return data

def raw(path,length,numwaves,row):
	'''This will read numwaves number of waveforms of length length into a structured numpy array *** Starts with row 0!'''
	fsize=os.stat(path).st_size
	totrows = (fsize-8)/(33+length*2)
	data=[]
	pos = 8+row*(33+length*2)
	with open(path) as f:
		f.seek(pos)
		data= np.core.records.fromfile(f,formats='B,i4,i4,i4,Q,Q,i4,(3500)i2',shape=numwaves,names='result,evID,board,channel,timestamp,requesttime,length,wave',byteorder='<')
	return data

def analysis_output(fname):
	with open(fname) as f:
		data=[]
		formats='B,i,i,i,Q,Q,i,f,f'
		names='result,evID,board,channel,timestamp,requesttime,risetime,energy,falltime'
		numwaves = os.stat(fname).st_size/(41)
		data=np.core.records.fromfile(f,formats=formats,shape=numwaves,names=names,byteorder='<')
	return data
