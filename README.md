# ca_analysis
Ca analysis toolkit (2017 data)


load.py is deprecated and its utility is encompassed in fileread.py

sample data available here: https://drive.google.com/file/d/1mmF6zKcPnwSyYhD-OYlQbGIy4Y1afj-P/view


Full data is located on the ACF
/lustre/haven/gamma/neutrons/ca45\_data
Python version is 2\*

Scripts:
wave\_display.py:
Arguments:
run,part,path,numwaves,startrow 
path = path to .bin files, where files are named like Run\_60\_1.bin
run = for this example 60 is run number
part = for this example 1 is part number
numwaves = integer from 1 to ??
startrow = integer from 1 to ??

Output:
Creates a wave plot, or multiple wave plots on screen, saves all plots to the plots/ directory

Example:
python wave\_display.py 60 1 /lustre/haven/gamma/neutrons/ca45\_data/2016/disk0/Calcium45/ 1 40

calcium_analysis.py
Arguments:
run,part,inpath,outpath
inpath = path to .bin files, where files are named like Run\_60\_1.bin
run = for this example 60 is run number
part = for this example 1 is part number
outpath = ./ most times ??
-f enable fitting
-t0 enable T_0 determination ??
-t
-Tn??

Output:
Several .part files ??
?? ??

Example: 
mpiexec -n 4 python calcium_analysis.py 131 0 /lustre/haven/gamma/neutrons/ca45_data/2017/disk1/ ./ -f -t0
