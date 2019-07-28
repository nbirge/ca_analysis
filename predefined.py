import numpy as np
from scipy.constants import *
from scipy.optimize import curve_fit
import os
cwd=os.getcwd()
#import load

pix_map={ 
39:(0,1),
40:(0,2),
41:(0,3),
50:(0,4),
51:(0,5),
52:(0,6),
53:(0,7),

62:(1,1),
63:(1,2),
64:(1,3),
65:(1,4),
66:(1,5),
75:(1,6),
76:(1,7),

77:(2,1),
78:(2,2),
87:(2,3),
88:(2,4),
89:(2,5),

}
#---------------------------------------------------------------------------------------
#               Taking simulation 'PixelDetector' --> Board,Channel
#---------------------------------------------------------------------------------------
def pixel_map(pixel):
    return pix_map[pixel]

v_pixel_map=np.vectorize(pixel_map)

def detector_board(detector):
    if detector == 'W':
        return 0
    else:
        return 3
v_detector_board=np.vectorize(detector_board)

def vec_pixel_map(sim):
    return v_pixel_map(sim)

def pixel_to_bdch(sim):
    bdch=v_pixel_map(sim['pixel'])
    bdch=np.array(bdch,dtype=int)
    bdch[0]+=v_detector_board(sim['detector'])
    return bdch
#---------------------------------------------------------------------------------------
try:
    calibration=np.load('/home/noah/Desktop/large_analysis/ca_analysis/simulation_comparison/calibration.npy')
except FileNotFoundError:
    calibration=np.load('/nics/d/home/nwbirge/large_analysis/ca_analysis/simulation_comparison/calibration.npy')

calibration=calibration.view(np.recarray)
def calibrate(energy_type,board,channel): 
    '''Use vec_calibrate(energy,board,channel) to return an array of calibrated energies which
           depend on board + ch'''
    bdch=int(board*8+channel) 
    if bdch ==6: 
        m,b=1/calibration.slope[3],calibration.offset[3]
    elif bdch==11: 
        m,b=1/calibration.slope[0],calibration.offset[0]
    elif bdch==12: 
        m,b=1/calibration.slope[1],calibration.offset[1]
    elif bdch==35: 
        m,b=1/calibration.slope[2],calibration.offset[2]
    elif bdch==34:
        m,b=1/calibration.slope[4],calibration.offset[4]
    else: 
        m,b=0,0   
    return (energy_type-b)*m
vec_calibrate=np.vectorize(calibrate)

def land(x,y):
    return np.logical_and(x,y)

def lor(x,y):
    return np.logical_or(x,y)

def cbins(x):
    return 0.5*(x[:-1]+x[1:])

def gauss(t,mu=0,sigma=1):
    '''Returns a normalize gaussian distribution centered at mu and with sigma width'''
    return (2*3.14159*sigma**2.)**(-0.5)*np.exp(-(t-t[mu])**2./(2.*sigma**2.))

def offdblgauss(x,*pars):
    a,mu,sigma,b,mu1,sigma1,offset=pars
    return a*np.exp(-1.*(x-mu)**2./(2.*sigma**2.))+b*np.exp(-1.*(x-mu1)**2./(2.*sigma1**2.))+offset

def dblgauss(x,*pars):
    a,mu,sigma,b,mu1,sigma1=pars
    return a*np.exp(-1.*(x-mu)**2./(2.*sigma**2.))+b*np.exp(-1.*(x-mu1)**2./(2.*sigma1**2.))

def single_pixel(x,board,channel):
    '''Returns a subset of x that only contains event with board == board and channel == channel'''
    return x[land(x['board']==board,x['channel']==channel)]

def pixel(board,channel):
    names=['Psr','39W','40W','41W','50W','51W','52W','53W','Psr','62W','63W','64W','65W','66W','75W','76W','Psr','77W','78W','87W','88W','89W','N/A','N/A','Psr','39E','40E','41E','50E','51E','52E','53E','Psr','62E','63E','64E','65E','66E','75E','76E','Psr','77E','78E','87E','88E','89E','N/A','N/A']
    bdch = board*8+channel
    return names[bdch]

def pixel_bdch(bdch): 
    names=['Psr','39W','40W','41W','50W','51W','52W','53W','Psr','62W','63W','64W','65W','66W','75W','76W','Psr','77W','78W','87W','88W','89W','N/A','N/A','Psr','39E','40E','41E','50E','51E','52E','53E','Psr','62E','63E','64E','65E','66E','75E','76E','Psr','77E','78E','87E','88E','89E','N/A','N/A']
    return names[bdch]

def sim_single_pixel(x,board,channel):
    pix= pixel(board,channel)
    return x[land(x['detector']==pix[-1],x['pixel']==int(pix[:-1]))]

def precuts(x,twindow=250,pilewindow=100,energy=100,remove_double=False):
    '''Returns portion of data with > twindow between separate neighboring events,
        <pilewindow for detected pileups, and have an energy > energy'''
    if remove_double==True:
        x=x[doubles(x)]
    x=x[x['energy']>energy]
    x.sort(order='timestamp')
    t0,t1,t2=x['timestamp'][0:-2],x['timestamp'][1:-1],x['timestamp'][2:]
    trutharray=land(t2-t1>twindow,t1-t0>twindow)
    x=x[1:-1][trutharray]
    x=x[lor(x['pileup']<2,np.abs(x['pilediff'])<pilewindow)]
    x=x[x['t0']>100]
    return x

def doubles(data,etype='energy'): 
    '''Returns a logical array with False -> 2nd+ recorded event'''
    timewindow=3500
    Ediff=1
    length=len(data)
    trutharray=np.ones(length,dtype=bool)
    multi=0
    i=0
    bdch=0
    while i < length:
        bdch=data['board'][i]*8+data['channel'][i]
        multi=len(data[i:i+100][(data[i:i+100]['board']*8+data[i:i+100]['channel'] == bdch) \
                                *(np.abs(data[i:i+100][etype]-data[i][etype])<Ediff) \
                                *(np.abs(data[i:i+100]['timestamp']-data[i:i+100]['timestamp'])<timewindow)])
        if multi>1 :
            trutharray[i]=False
        i+=1
    return trutharray


#(data[i:i+100]['board']*8+data[i:i+100]['channel'] == bdch)*(np.abs(data[i:i+100][etype]-data[i][etype])<Ediff)*(np.abs(data[i:i+100]['timestamp']-data[i:i+100]['timestamp'])<timewindow)



def precuts_multipixel_twindow(x,twindow=250,pilewindow=100,energy=100):
    x=x[x['energy']>100]
    x.sort(order='timestamp')
    x=x[doubles(x)]
    t0,t1,t2=x['timestamp'][0:-2],x['timestamp'][1:-1],x['timestamp'][2:]
    trutharray=land(t2-t1>250,t1-t0>250)
    x=x[1:-1][trutharray]
    x=x[lor(x['pileup']<2,np.abs(x['pilediff'])<pilewindow)]
    x=x[x['t0']>100]
    return x

def good_timestamps(data,time_cut=2*3600/4e-9):
    trutharray=data['timestamp']<time_cut
    return trutharray

def only_cal_pixels(data):
    trutharray=np.zeros(len(data))
    for i in [11,12,35]:
        bd,ch=int(i/8),int(i%8)
        trutharray=lor(trutharray,\
                    land(data['board']==bd,data['channel']==ch))
    return trutharray


def standard_cuts(x,etype):
    return 0

def pixcut(x,attr,board,channel):
    '''Returns x['attr'][logical_and(x[board]==board,x[channel]==channel)
        Ex >>> subarray = predefined.pixcut(data,'trap',0,6)'''
    return x[attr][land(x['board']==board,x['channel']==channel)]


def sim_spixel_cut(x):
    '''Returns x[g] where all events of x[g] are isolated to a single pixel. I.e. backscattering
        between on a single pixel is allowed. Any other is'''
    length=len(x)
    g=np.ones(len(x),dtype=bool)
    i=0
    while i<length-1:
        if x['entry'][i]!=x['entry'][i+1]:  # if next entry is different keep event
            i+=1
            continue
        elif x['detector'][i]!=x['detector'][i+1] or x['pixel'][i]!=x['pixel'][i+1] or not g[i]:
            g[i]=False
            g[i+1]=False
            i+=1
            continue
        else:
            i+=1
    return x[g]

def sim_comb_single_pixel(x):
    '''Returns an array where all energies corresponding to the the same event number are summed.
            This expects that x has already been reduced to single pixel data.'''
    comb=np.zeros_like(x)
    comb['energy']+=-10
    i=0
    j=0
    while i < len(x)-1 and j<len(x)-1:
        j=i+1
        backscattering=x[i]['entry']==x[j]['entry']
        energy= x['energy'][i]
        while backscattering and j<len(x)-1:
            energy+=x['energy'][j]
            j+=1
            backscattering=x[i]['entry']==x[j]['entry']
        comb[i]=x[i]
        comb[i]['energy']=energy
        i=j
    return comb

def sim_single_event(x):
    '''Returns a trutharray of len(x), where True entries correspond to entries with
        only one event. I.e. backscattering is cut.'''
    trutharray=np.zeros(len(x),dtype=bool)
    i,j=0,1
    while i<len(x)-1 and j<len(x)-1:
        backscattering=x['entry'][i]==x['entry'][j]
        while backscattering:
            j+=1
            backscattering =x['entry'][i]==x['entry'][j]
        if j-i == 1:
            trutharray[i:j]=True
        else:
            trutharray[i:j]=False
        i=j
        j=i+1
    return trutharray


def sim_thresh(x,thresh=15):
    length=len(x)
    trutharray=ones(length,dtype=bool)
    thresh=15.
    i=0
    while i < length-2:
        j=i+1
        backscattering=x['entry'][i]==x['entry'][j]
        energy=x['energy'][i]
        while backscattering and j<length-2:
            if x['detector'][j]== x['detector'][i] and x['pixel'][j]== x['pixel'][i]:
                energy+=x['energy'][j]
            j+=1
            backscattering=x['entry'][i]==x['entry'][j]
        if energy < thresh:
            trutharray[i]=False
            j=i+1
            backscattering=x['entry'][i]==x['entry'][j]
            while backscattering and j <length-2:
                if x['energy'][j]< thresh:
                    trutharray[j]=False
                j+=1
                backscattering=x['entry'][i]==x['entry'][j]
        i=j
    return trutharray

def board(pixel):
    if int(pixel)<54:
        return 0
    elif int(pixel)<77:
        return 1
    else:
        return 2

def channel(pixel):
    ch=int(pixel)
    if ch==39 or ch==62 or ch==77:
        return 1
    elif ch==40 or ch==63 or ch==78:
        return 2
    elif ch==41 or ch==64 or ch==87:
        return 3
    elif ch==50 or ch==65 or ch==88:
        return 4
    elif ch==51 or ch==66 or ch==89:
        return 5
    elif ch==52 or ch==75:
        return 6
    elif ch==53 or ch==76:
        return 7
    else:
        return 0

def east_west(detector):
    if detector=='W':
        return 0
    else:
        return 3


vectorboard=np.vectorize(board,otypes=[int])
vectorchannel=np.vectorize(channel,otypes=[int])
vectordetector=np.vectorize(east_west,otypes=[int])

def sim_restructure(simdata):
#    vectorboard=np.vectorize(board,otypes=[int])
#    vectorchannel=np.vectorize(channel,otypes=[int])
#    vectordetector=np.vectorize(east_west,otypes=[int])
    dtype=[('entry', '<i4'), ('board', '<i4'), ('channel', '<i4'), ('timestamp', '<f4'), ('energy', '<f4')]
    output=np.zeros(len(simdata),dtype=dtype)
    output['board']=vectorboard(simdata['pixel'])+vectordetector(simdata['detector'])
    output['channel']=vectorchannel(simdata['pixel'])
    output['entry']=simdata['entry']
    output['timestamp']=simdata['timestamp']
    output['energy']=simdata['energy']
    return output

def Fierz_fit(histogram,bins,beg,end,normbin,sigma=1):
    ''' histogram is generally a ratio of data to simulation w/b=0 (should be a normalized ratio),
        bins are the energybins for histogram, beg and end are the beginning 
        and ending to the fit range, and normbin is the bin of energybins to 
        which you're normalizing'''
    trutharray=land(bins>beg,bins<end)
    fitbins=bins[trutharray].astype(float)
    normE=bins[normbin]
    normarray=fitbins==bins[normbin]
    fithist=histogram[trutharray].astype(float)
    if isinstance(sigma,int):
        sigma=np.ones(len(fithist),dtype=float)
    sigma=sigma[trutharray]
    m_e_kev=m_e*c**2./(kilo*eV)
    shape=lambda x,b: (1+b*m_e_kev/(m_e_kev+x))/(1+b*m_e_kev/(m_e_kev+normE))
    pars,vrs=curve_fit(shape,fitbins,fithist,p0=-0.01,sigma=sigma)
    vrs=np.sqrt(np.diag(vrs))
    return pars[0],vrs[0]

def Fierz_arb_norm(X,a,b):
    '''X=tuple of simulation, bins & trutharray'''
    histo,bins,trutharray = X
    trutharray=trutharray.astype(bool)
    redHisto=histo[trutharray].astype(float)
    redBins=bins[trutharray].astype(float)
    return a*redHisto*(1+b*m_e*c**2./(m_e*c**2.+redBins*kilo*eV))

def Fierz_arb_norm_fit(simulation,data,bins,beg=100,end=200,sigma=1):
    '''Fits simulation to data. I.e. Data= Fierz_arb_norm((simulation,bins,truth),a,b) | 
        To plot the fit, Simulation is scaled by fit pars to data '''
    trutharray=land(bins>beg,bins<end)
    guess=[simulation[trutharray][0]/data[trutharray][0],0]
    bounds=[(0,-10),(np.inf,10)]
    X=(simulation,bins,trutharray)
    if isinstance(sigma,int):
        sigma=np.ones(np.sum(len(bins)),dtype=float)
    pars,vrs=curve_fit(Fierz_arb_norm,X,data[trutharray],sigma=sigma[trutharray],p0=guess)#,bounds=bounds)
    vrs=np.sqrt(np.diag(vrs))
    return pars,vrs


#def combsets(runs):
#    '''Returns the entered runs as a single numpy array. 
#    Ex: >>> data = predefined.combsets([131,132,133])'''
#    dat=[]
#    for run in runs:
#        print run
#        dat.append(precuts(load.loadtrapnfit1('./new/Run_'+str(run)+str('-all.fin'))))
#    dat=np.concatenate(dat[:])
#    return dat




