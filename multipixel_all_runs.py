import os

path='/home/noah/Desktop/large_analysis/ca_analysis/cur_data'

runs=[run for run in os.listdir(path) if run.endswith('all.dat')]

for run in runs:
    r=run[4:len(run)-8]
    if r == 0 :
        continue
    os.system('python multipixel_sum.py '+r)
