import fileread as fr
import numpy as np
import sys,os

#path,run = str(sys.argv[1]),str(sys.argv[2])
path='./'
beg= int(sys.argv[1])
if len(sys.argv[:])>2:
    end = int(sys.argv[2])+1
else:
    end = beg+1
runs = range(beg,end,1)

for i in runs:
	path='./'+str(i)+'/'
	print 'combining parts of run: '+str(i)
	run=str(i)
	name='Run_'+run
	x=filter(lambda x: x.startswith('Run_'+str(run)) and x.endswith('-comb.bin') and x!='Run_'+str(run)+'_0-comb.bin',os.listdir(path))
	x=sorted(x)
	print x
	header= np.zeros(1,dtype=[('theader','Q'),('formats','10i')])

	data,header['theader'],header['formats'] = fr.gen_output(path+name+'_0-comb.bin')
	for fyle in x:
		print fyle
		if os.stat(path+fyle).st_size>1000:
			data=np.concatenate((fr.gen_output(path+fyle)[0],data))
#
	print data
	output=path+'Run_'+run+'-all.dat'
	print '\n'+output
	with open(path+'Run_'+run+'-all.dat','wb') as f:
		header.tofile(f)
		data.tofile(f)
