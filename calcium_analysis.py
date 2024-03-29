#GIT IS NOT COOPERATING!!!
import mpi4py
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

fitting,pileup,trapNfit,findt0,osc_removal,backscatter,simulation = 0,0,0,0,0,0,0
for i in range(len(sys.argv[:])):
    if sys.argv[i] == '-f':
        fitting=1
    if sys.argv[i] == '-p':
        pileup=1
        pile_thresh=int(sys.argv[i+1])
    if sys.argv[i] == '-tNf':
        trapNfit=1
        fitting=1
    if sys.argv[i] == '-t0':
        findt0=1
    if sys.argv[i] == '-osc_removal':
        osc_removal=1
    if sys.argv[i] == '-backscatter':
        backscatter=1
    if sys.argv[i] == '-simulation':
        simulation=1

length = int(3500)
if rank == 0:
    with open(inpath+fname,'rb') as f:
        theader = np.core.records.fromfile(f,formats='Q',shape=1,byteorder='<')[0][0]
        f.seek(8+1+3*4+2*8)
        length=np.core.records.fromfile(f,formats='i4',shape=1,byteorder='<')[0][0]
        for i in range(1,size):
            comm.send(length,dest=i)
    if (fsize-8)%(33+2*length) != 0:
        print('File format incorrect or file corrupted.', fsize)
else:
    length = comm.recv(source=0)

####
numrows = (fsize-8)/(33+2*length)
datachunk = int(numrows/(size-1))

if rank< size-1:
    row = (rank-1)*datachunk
    datachunk= numrows/(size-1)
elif rank == size-1:
    row = (rank-1)*datachunk
    datachunk=numrows/(size-1)+numrows%(size-1)

datachunk=int(datachunk)

piece = int(100000) # 120,000 waveforms in memory < 1GB, there will be datachunk/piece iterations per core
begin=time.time()


dtype=[('result','B'),('evID','i'),('board','i'),('channel','i'),('timestamp','Q'),('requesttime','Q'),('risetime','i'),('energy','f'),('pretrigrms','f')]
fformat=np.zeros(10,'i')

if fitting ==1:
    dtype.append(('falltime','f'))
    fformat[0]=1

if pileup ==1:
    liltrap=np.zeros((48,length))
    fast_rise,fast_top=10,4
    wo.pixel_traps(workarr=liltrap,rise=fast_rise,top=fast_top)
    for i in [('pileup','i'),('pilediff','i'),('pileamp','f')]:
        dtype.append(i)
    fformat[1]=pile_thresh


if trapNfit == 1:
    dtype.append(('fitenergy','f'))
    fformat[2]=1

if findt0==1:
    dtype.append(('t0','i'))
    fformat[3]=1

if osc_removal==1:
    fitpars=np.zeros((piece+datachunk%piece,9),dtype=float)
    dtype.append(('osc_amps','4f'))
    dtype.append(('osc_errors','4f'))
    dtype.append(('s2','f'))
    dtype.append(('osc_energy','f'))
    fformat[4]=10           # Equal to the number of sine/cosine functions currently using

if backscatter ==1:
    filters=np.zeros((piece+datachunk%piece,3),dtype=int)
    for i in [('pfilter0','i'),('pfilter1','i'),('pfilter2','i')]:
        dtype.append(i)
    fformat[5]=1


writebuffer=np.zeros(piece+datachunk%piece,dtype=dtype)
if rank>0:
    trap = np.zeros((48,length))
    rise,top,fall=300,100,1100        #not going to use fall
    wo.pixel_traps(workarr=trap,rise=rise,top=top)
#    b,a = signal.bessel(5,0.05,btype='low',analog=False)
    maxamps,maxlocs,risetimes=np.zeros(piece+datachunk%piece),np.zeros(piece+datachunk%piece),\
np.zeros(piece+datachunk%piece)
    traps=np.zeros((piece+datachunk%piece,length))


    with open(outpath+fname[:-4]+'-'+str(rank)+'.part','w') as f:
        for i in range(int(datachunk/piece)):
            if i == int(datachunk/piece)-1:
                rem = datachunk%piece
            else:
                rem = 0
            print(i,row+i*piece+rem,piece,rank)
            sys.stdout.flush()
            try:
                writebuffer[0:piece+rem]=0
                numwaves = piece+rem
                data = fr.raw(path=inpath+fname,length=length,numwaves=piece+rem,row=int(row+i*piece))
                writebuffer[0:piece+rem]['result'], writebuffer[0:piece+rem]['evID'], writebuffer[0:piece+rem]['board'], writebuffer[0:piece+rem]['channel'], writebuffer[0:piece+rem]['timestamp'], writebuffer[0:piece+rem]['requesttime'] = data['result'], data['evID'], data['board'], data['channel'], data['timestamp'], data['requesttime']
                if simulation == 1:
                    wo.sim_baseline_restore(data,600)
                else:
                    wo.baseline_restore(data,600)        #restores baseline and performs necessary preformatting of the data (data & 16383...)
                wo.pretrigger_rms(wave=data['wave'],workarr=maxamps[0:piece+rem],pretrig_timebin=600)
                writebuffer[0:piece+rem]['pretrigrms']=maxamps[0:piece+rem].copy()
#                smooth_wave= signal.filtfilt(b,a,data['wave'])
                wo.maxes(waves=data['wave'],startpoint=500,wavelength=length,maxamps=maxamps[0:piece+rem],maxlocs=maxamps[0:piece+rem])
                wo.rises(data['wave'],maxamps[0:piece+rem],maxamps[0:piece+rem],risetimes[0:piece+rem])
                writebuffer[0:piece+rem]['risetime']=risetimes[0:piece+rem].copy()

                if fitting ==1:
                    wo.tail_fit(data=data['wave'],output=maxamps[0:piece+rem])
                    writebuffer[0:piece+rem]['falltime']=maxamps[0:piece+rem].copy()
                if trapNfit ==1 and fitting ==1:
                    wo.fitted_trap(data=data,rise=rise,top=top,fall=maxamps[0:piece+rem],output=traps[0:piece+rem])
                    wo.trap_energy(traps=traps[0:piece+rem],length=length,output=maxamps[0:piece+rem])
                    writebuffer[0:piece+rem]['fitenergy'] = maxamps[0:piece+rem].copy()
                if pileup ==1:
                    wo.apply_trap(rise=fast_rise,data=data,trap=liltrap,output=traps)
                    wo.pileup(data=traps[0:piece+rem],thresh=pile_thresh,amplitudes=maxamps[0:piece+rem],tdiff=risetimes[0:piece+rem],numpeaks=maxlocs[0:piece+rem])
                    writebuffer[0:piece+rem]['pileup']=maxlocs[0:piece+rem].copy()
                    writebuffer[0:piece+rem]['pilediff']=risetimes[0:piece+rem].copy()
                    writebuffer[0:piece+rem]['pileamp']=maxamps[0:piece+rem].copy()

#               UNADJUSTED TRAP OUTPUT
                wo.apply_trap(rise=rise,data=data,trap=trap,output=traps)
                wo.trap_energy(traps[0:piece+rem],length=length,output=maxamps[0:piece+rem])    #need these maxamps for fitting later...
                writebuffer[0:piece+rem]['energy']=maxamps[0:piece+rem].copy()

                if osc_removal==1:
                    wo.osc_removal(data,fitpars)
                    writebuffer[0:piece+rem]['osc_amps']=fitpars[0:piece+rem,0:4]
                    writebuffer[0:piece+rem]['osc_errors']=fitpars[0:piece+rem,4:8]
                    writebuffer[0:piece+rem]['s2']=fitpars[0:piece+rem,8]
                    wo.apply_trap(rise=rise,data=data,trap=trap,output=traps)
                    wo.trap_energy(traps[0:piece+rem],length=length,output=maxamps[0:piece+rem])    #need these maxamps for fitting later...
                    writebuffer[0:piece+rem]['osc_energy']=maxamps[0:piece+rem].copy()

                if backscatter ==1:
                    wo.backscatter(data=data,output=filters)
                    writebuffer[0:piece+rem]['pfilter0']=filters[0:piece+rem,0]
                    writebuffer[0:piece+rem]['pfilter1']=filters[0:piece+rem,1]
                    writebuffer[0:piece+rem]['pfilter2']=filters[0:piece+rem,2]

                if findt0==1:
                    wo.find_t0(traps=data['wave'],output=maxamps)
                    writebuffer[0:piece+rem]['t0']=maxamps[0:piece+rem].copy()
                writebuffer[0:piece+rem].tofile(f)
            except ZeroDivisionError:
                print('F.up occurred here:')   #lol
                print(rank,i,row+i*piece+rem,piece+rem)
                sys.stdout.flush()

    end=time.time()
    print('Rank ',rank,' finished in ',end-begin,' seconds')
    comm.send(rank,dest=0)


if rank == 0:
    check = np.arange(1,size,1)
    header= np.zeros(1,dtype=[('theader','Q'),('formats','10i')])
    header['theader'][0]=theader
    header['formats'][0:10]=fformat[0:10]
    with open(outpath+'Run_'+str(run)+'_'+str(part)+'_0.part','wb') as f:
        header.tofile(f)
        f.close()
        print('Created '+'Run_'+str(run)+'_'+str(part)+'_0.part')
        sys.stdout.flush()
    name=''
    for i in np.arange(1,size,1):
        print(check[i-1]==comm.recv(source=i),i)
        sys.stdout.flush()
        name+=outpath+'Run_'+str(run)+'_'+str(part)+'-'+str(i)+'.part '
    os.system('cat '+outpath+'Run_'+str(run)+'_'+str(part)+'_0.part '+name+' > '+outpath+'Run_'+str(run)+'_'+str(part)+'-comb.bin')
    os.system('rm '+outpath+'Run_'+str(run)+'_'+str(part)+'_0.part '+name)
    print('Successfully finished')

    #File consolidation should go here!















