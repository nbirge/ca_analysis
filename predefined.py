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

def precuts(x):
    x=x[x['energy']>100]
    x.sort(order='timestamp')
    t0,t1,t2=x['timestamp'][0:-2],x['timestamp'][1:-1],x['timestamp'][2:]
    trutharray=land(t2-t1>250,t1-t0>250)
    x=x[1:-1][trutharray]
    for field in x.dtype.names:
        if field == 'falltime':
            x=x[land(x['falltime']>200,x['falltime']<2000)]
    return x

def pixcut(x,attr,board,channel):
    '''Returns x['attr'][logical_and(x[board]==board,x[channel]==channel)
        Ex >>> subarray = predefined.pixcut(data,'trap',0,6)'''
    return x[attr][land(x['board']==board,x['channel']==channel)]

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

def pixel(board,channel):
    names=['Psr','39W','40W','41W','50W','51W','52W','53W','Psr','62W','63W','64W','65W','66W','75W','76W','Psr','77W','78W','87W','88W','89W','N/A','N/A','Psr','39E','40E','41E','50E','51E','52E','53E','Psr','62E','63E','64E','65E','66E','75E','76E','Psr','77E','78E','87E','88E','89E','N/A','N/A']
    bdch = board*8+channel
    return names[bdch]

