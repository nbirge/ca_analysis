import numpy as np
import fileread as fr
import matplotlib.pyplot as plt
import wave_ops as wo
import sys

run,part,path,numwaves,startrow = sys.argv[1:6]
if not path.endswith('/'):
	path+='/'
numwaves=int(numwaves)
startrow=int(startrow)
length=3500

#for i in range(4,len(sys.argv[4:])):		#PARSE ARGS LIKE IN ca_analysis.py

name= 'Run_'+str(run)+'_'+str(part)+'.bin'

data = fr.raw(path+name,length=length,numwaves=numwaves,row=startrow)
wo.baseline_restore(data,pretrigger=600)

tbins=np.arange(length)
for i in range(len(data)):
	plt.plot(tbins,data['wave'][i,0:length],'b-',label='Raw Wave')
	plt.show()
