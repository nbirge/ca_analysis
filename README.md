# ca_analysis
Ca analysis toolkit (2017 data)


load.py is deprecated and its utility is encompassed in fileread.py

sample data available here: https://drive.google.com/file/d/1mmF6zKcPnwSyYhD-OYlQbGIy4Y1afj-P/view


Full data is located on the ACF
/lustre/haven/gamma/neutrons/ca45_data
Python version is 2\*

Scripts:
wave_display.py:
Arguments:
run,part,path,numwaves,startrow 
path = path to .bin files, where files are named like Run_60_1.bin
run = for this example 60 is run number
part = for this example 1 is part number
numwaves = integer from 1 to ??
startrow = integer from 1 to ??

Output:
Creates a wave plot, or multiple wave plots on screen, saves all plots to the plots/wave_display/ directory

Example:
python wave_display.py 60 1 /lustre/haven/gamma/neutrons/ca45_data/2016/disk0/Calcium45/ 1 40

calcium_analysis.py
Arguments:
run,part,inpath,outpath
inpath = path to .bin files, where files are named like Run_60_1.bin
run = for this example 60 is run number
part = for this example 1 is part number
outpath = ./ most times ??
-f enable fitting
-t0 enable T_0 determination ??
-t
-Tn??

Output:
Several .part files which are consolidated upon successful completion to a .bin file for that run. The number of part files is determined by
the number of cores which is passed to the python multiplexer


Example: 
mpiexec -n 4 python calcium_analysis.py 131 0 /lustre/haven/gamma/neutrons/ca45_data/2017/disk1/ ./ -f -t0

consolidation.py
Arguments
begin,end,path 
begin = run list begin, Run must have a file 0 'Run_130_0-comb.bin' for some reason, seems to be able to handle missing intermediate files
end = the last run to use
path = an absolute path for now, will test relative paths later ??

Output:
.dat file, In this example output is something like Run_130-all.dat

Example:
python consolidation.py 130 130 /lustre/haven/user/griley4/CalciumExperiment/ca_analysis/
