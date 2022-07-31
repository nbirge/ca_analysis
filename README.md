# ca_analysis
Ca analysis toolkit (2017 data)<br>

This contains the digital signal processing (DSP) toolkit I used for processing the raw Calcium-45 waveforms (typically .py files) as well as many jupyter notebooks used to analyze the processed waveform data


load.py is deprecated and its utility is encompassed in fileread.py <br>

sample data available here: https://drive.google.com/file/d/1mmF6zKcPnwSyYhD-OYlQbGIy4Y1afj-P/view <br>


Full data is located on the UTK's ACF <br>
/lustre/haven/gamma/neutrons/ca45_data <br>
Python version is 3.6.7\*

<b>Main analysis scripts:</b>

<u>wave_display.py:</u> <br>
Arguments:<br>
run,part,path,numwaves,startrow <br>
path = path to .bin files, where files are named like Run_60_1.bin<br>
run = for this example 60 is run number<br>
part = for this example 1 is part number<br>
numwaves = integer number of waves to show<br>
startrow = integer row or waveform number in file on which to start processing<br>

Output:
Creates a wave plot, or multiple wave plots on screen, saves all plots to the plots/wave_display/ directory <br>

Example: <br>
python wave_display.py 60 1 /lustre/haven/gamma/neutrons/ca45_data/2016/disk0/Calcium45/ 1 40 <br>

<u>calcium_analysis.py</u> <br>
Arguments: <br>
run,part,inpath,outpath <br>
inpath = path to .bin files, where files are named like Run_60_1.bin <br>
run = for this example 60 is run number <br>
part = for this example 1 is part number <br>
outpath = ./  <br>
-f enable fitting <br>
-t0 enable T_0 determination for waveforms <br> <br>

Output:
Several .part files which are consolidated upon successful completion to a .bin file for that run. The number of part files is determined by
the number of cores which is passed to the python multiplexer <br> <br>


Example:  <br>
mpiexec -n 4 python calcium_analysis.py 131 0 /lustre/haven/gamma/neutrons/ca45_data/2017/disk1/ ./ -f -t0

<u>consolidation.py </u> <br>
Arguments <br>
begin,end,path <br> 
begin = run list begin, Run must have a file 0 'Run_130_0-comb.bin' for some reason, seems to be able to handle missing intermediate files <br>
end = the last run to use <br>
path = an absolute path <br>

Output:
.dat file, In this example output is something like Run_130-all.dat <br>

Example: <br>
python consolidation.py 130 130 /lustre/haven/user/griley4/CalciumExperiment/ca_analysis/
