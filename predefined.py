import numpy as np
#import load

def land(x,y):
    return np.logical_and(x,y)

def lor(x,y):
    return np.logical_or(x,y)

def cbins(x):
    return 0.5*(x[:-1]+x[1:])

def gauss(x,*pars):
    a,mu,sigma=pars
    return a*np.exp(-(x-mu)**2./(2.*sigma**2.))

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

def sim_single_pixel(x,board,channel):
    pix= pixel(board,channel)
    return x[land(x['detector']==pix[-1],x['pixel']==int(pix[:-1]))]

def precuts(x):
    x=x[x['energy']>100]
    x.sort(order='timestamp')
    t0,t1,t2=x['timestamp'][0:-2],x['timestamp'][1:-1],x['timestamp'][2:]
    trutharray=land(t2-t1>250,t1-t0>250)
    x=x[1:-1][trutharray]
    x=x[lor(x['pileup']<2,x['pilediff']<60)]
    x=x[x['t0']>100]
    return x

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
    comb['energy']+=-5000
    i=0
    while i < len(x)-2:
        j=i+1
        backscattering=x[i]['entry']==x[j]['entry']
        energy= x['energy'][i]
        while backscattering:
            energy+=x['energy'][j]
            j+=1
            if j<len(x)-2:
                backscattering=x[i]['entry']==x[j]['entry']
        comb[i]=x[i]
        comb[i]['energy']=energy
        i=j
    return comb

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


def sim_restructure(simdata):
    vectorboard=np.vectorize(board,otypes=[int])
    vectorchannel=np.vectorize(channel,otypes=[int])
    vectordetector=np.vectorize(east_west,otypes=[int])
    dtype=[('entry', '<i4'), ('board', '<i4'), ('channel', '<i4'), ('timestamp', '<f4'), ('energy', '<f4')]
    output=np.zeros(len(simdata),dtype=dtype)
    output['board']=vectorboard(simdata['pixel'])+vectordetector(simdata['detector'])
    output['channel']=vectorchannel(simdata['pixel'])
    output['entry']=simdata['entry']
    output['timestamp']=simdata['timestamp']
    output['energy']=simdata['energy']
    return output



#def combsets(runs):
#    '''Returns the entered runs as a single numpy array. 
#    Ex: >>> data = predefined.combsets([131,132,133])'''
#    dat=[]
#    for run in runs:
#        print run
#        dat.append(precuts(load.loadtrapnfit1('./new/Run_'+str(run)+str('-all.fin'))))
#    dat=np.concatenate(dat[:])
#    return dat

def calibrate(bins,pixel):
    '''Returns calibrated bins for the given pixel (or bdch).
        Ex: >>> histbins = predefined.calibrate(histbins,'52W')
                or
            >>> histbins = predefined.calibrate(histbins,6)'''
    if pixel=='52W'or pixel == 6:
        return bins/6.2844 + 19.1/6.2844


