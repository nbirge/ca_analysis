import fileread as fr
import numpy as np
import sys

path,run = str(sys.argv[1]),str(sys.argv[2])

name='Run_'+run
x=filter(lambda x: x.startswith('Run_'+str(run)) and x.endswith('-comb.bin') and x!='Run_'+str(run)+'_0-comb.bin',os.listdir(path))
x=sorted(x)

header= np.zeros(1,dtype=[('theader','Q'),('formats','10i')])

data,header['theader'],header['formats'] = fr.gen_output(path+name+'_0-comb.bin')
for fyle in x:
	data=np.concatenate((fr.gen_output(path+fyle)[0],data))


with open(path+'Run'+run+'-all.dat','wb') as f:
	header.tofile(f)
	data.tofile(f)







