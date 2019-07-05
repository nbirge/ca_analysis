import os
import sys

if len(sys.argv)<2:
    print('Enter run number')
    exit()

threshes=[50,40,30,20,10,0]

folder='/home/noah/Desktop/large_analysis/ca_analysis/cur_data/'
run=int(sys.argv[1])
coinc_window=100

for thresh in threshes:
    os.system('python multipixel_sum.py {:d} -t {:d} -c {:d}'.format(run,thresh,coinc_window))

