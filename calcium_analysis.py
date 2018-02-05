#GIT IS NOT COOPERATING!!!

from mpi4py import MPI
import numpy as np
import wave_ops as wo
import fileread as fr
from scipy import signal
import sys,os,time


comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


run,part,inpath,outpath = int(sys.argv[1]),int(sys.argv[2]),str(sys.argv[3]),str(sys.argv[4])
if inpath[-1] != '/':
	inpath+='/'
if outpath[-1] != '/':
	outpath+='/'
fname='Run_'+str(run)+'_'+str(part)+'.bin'
fsize = os.stat(inpath+fname).st_size

fitting=1
pileup=1

length = -1.
if rank == 0:
	with open(inpath+fname) as f:
		theader = np.core.records.fromfile(f,formats='Q',shape=1,byteorder='<')[0][0]
		f.seek(8+1+3*4+2*8)
		length=np.core.records.fromfile(f,formats='i4',shape=1,byteorder='<')[0][0]
		for i in range(1,size):
			comm.send(length,dest=i)
	if (fsize-8)%(33+2*length) != 0:
		print 'File format incorrect or file corrupted.', fsize
else:
	length = comm.recv(source=0)

####
numrows = (fsize-8)/(33+2*length)
datachunk = numrows/(size-1)

if rank ==1 :
	datachunk = numrows/(size-1)-1000
	row = 1000
	if datachunk<0:
		datachunk=numrows/(size-1)		#THIS isn't right
		row = 0
elif rank>1 and rank< size-1:
	row = (rank-1)*datachunk
	datachunk= numrows/(size-1)
elif rank == size-1:
	row = (rank-1)*datachunk
	datachunk=numrows/(size-1)+numrows%(size-1)
	if datachunk < 0:
		datachunk = 0


piece = int(10000) # 120,000 waveforms in memory < 1GB, there will be datachunk/piece iterations per core

begin=time.time()


dtype=[('result','B'),('evID','i'),('board','i'),('channel','i'),('timestamp','Q'),('requesttime','Q'),('risetime','i'),('energy','f')]
fformat=np.zeros(10,'B')

if fitting ==1:
	dtype.append(('falltime','f'))
	fformat[0]=1

if pileup ==1:
	dtype.append(('pileup','i'))
	fformat[1]=1
	
writebuffer=np.zeros(piece+datachunk%piece,dtype=dtype)
if rank>0:
	trap = np.zeros((48,length))
	rise,top,fall=300,200,1100		#not going to use fall
	wo.pixel_trap(workarr=trap,rise=rise,top=top)
	b,a = signal.bessel(5,0.05,btype='low',analog=False)
	maxamps,maxlocs,risetimes=np.zeros(piece+datachunk%piece),np.zeros(piece+datachunk%piece),np.zeros(piece+datachunk%piece)
	traps=np.zeros((piece+datachunk%piece,length))
	if pileup == 1:
		liltrap=np.zeros((48,length))
		fast_rise,fast_top=10,4
		wo.pixel_trap(arr=liltrap,rise=fast_rise,top=fast_top)
	with open(outpath+fname[:-4]+'-'+str(rank)+'.part','w') as f:
		for i in range(datachunk/piece):

			if i == datachunk/piece-1:
				rem = datachunk%piece
			else:
				rem = 0
			print i,row+i*piece+rem,piece,rank
			try:
				writebuffer[0:piece+rem]=0
				numwaves = piece+rem
				data = fr.raw(path=inpath+fname,length=length,numwaves=piece+rem,row=row+i*piece)
				writebuffer[0:piece+rem]['result'], writebuffer[0:piece+rem]['evID'], writebuffer[0:piece+rem]['board'], writebuffer[0:piece+rem]['channel'], writebuffer[0:piece+rem]['timestamp'], writebuffer[0:piece+rem]['requesttime'] = data['result'], data['evID'], data['board'], data['channel'], data['timestamp'], data['requesttime']
				wo.baseline_restore(data,600)		#restores baseline and performs necessary preformatting of the data (data & 16383...)
				smooth_wave= signal.filtfilt(b,a,data['wave'])
				wo.maxes(waves=smooth_wave,startpoint=500,wavelength=length,maxamps=maxamps[0:piece+rem],maxlocs=maxamps[0:piece+rem])
				wo.rises(smooth_wave,maxamps[0:piece+rem],maxamps[0:piece+rem],risetimes[0:piece+rem])
				writebuffer[0:piece+rem]['risetime']=risetimes[0:piece+rem]

				if fitting ==1:
					wo.tail_fit(data=smooth_wave,output=maxamps[0:piece+rem])
					writebuffer[0:piece+rem]['falltime']=-1.*maxamps[0:piece+rem]
				if pileup ==1:
#					traps= np.apply_along_axis(lambda m: signal.fftconvolve(m, liltrap, mode='full'), axis=1, arr=data['wave'])/(fast_rise*fall)		#gotta now smooth the waves and then look for peaks #FUUUUUUUUUCK NOT A GOOD WAY TO DO THIS FOR A SPECIFIC PIXEL!!!!
					wo.apply_trap(workarr=data,trap=liltrap,output=traps)
					wo.pileup(data=traps,workarr=maxamps,thresh=125)
					writebuffer[0:piece+rem]['pileup']=maxamps[0:piece+rem]

#				traps= np.apply_along_axis(lambda m: signal.fftconvolve(m, trap, mode='full'), axis=1, arr=data['wave'])/(rise*fall)	#FUUUUUUUUUCK NOT A GOOD WAY TO DO THIS FOR A SPECIFIC PIXEL!!!!
				wo.apply_trap(data=data,trap=trap,output=traps)
				wo.trap_energy(traps,length=length,output=maxamps[0:piece+rem])	#need these maxamps for fitting later...
				writebuffer[0:piece+rem]['energy']=maxamps[0:piece+rem]
			
				writebuffer[0:piece+rem].tofile(f)
			except ValueError:
				print 'Fuckup occurred here:'
				print rank,i,row+i*piece+rem,piece+rem

	end=time.time()
	print 'Rank ',rank,' finished in ',end-begin,' seconds'
	comm.send(rank,dest=0)


if rank == 0:
	check = np.arange(1,size,1)
	header= np.zeros(1,dtype=[('theader','Q'),('formats','10B')])
	header['theader'][0]=theader
	header['formats'][0:10]=fformat[0:10]
	header.tofile(outpath+'Run_'+str(run)+'_'+str(part)+'_0.part')
	for i in np.arange(1,size,1):
		print check[i-1]==comm.recv(source=i),i
	os.system('cat '+outpath+'Run_'+str(run)+'_'+str(part)+'_0.part '+outpath+'Run_'+str(run)+'_'+str(part)+'-*.part > '+outpath+'Run_'+str(run)+'_'+str(part)+'-comb.bin')
	os.system('rm '+outpath+'Run_'+str(run)+'_'+str(part)+'*.part')

	#File consolidation should go here!















