import sys
sys.path.insert(0,'/home/noah/Desktop/large_analysis/ca_analysis/')
from matplotlib.pyplot import*
import fileread as fr
import predefined as pd

dataloc='/home/noah/Desktop/large_analysis/ca_analysis/cur_data/'

runs=[120]

for run in runs:
    data=fr.gen_output(dataloc+'Run_'+str(run)+'-all.dat')
