import scipy as sc
output=sc.test('all',raise_warnings='release')
import pylab as gr
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider, Button, RadioButtons
import scipy.io as sio
import sympy as sy
from analysis1 import GetLocalExtrema
gr.ion()


x, p= sy.symbols('x p')
a, pa, r, xr= sy.symbols('a pa r xr')
h= sy.symbols('h')

dpAP = a*(pa-p) + h*(1-p)
dxAP = r*x*(xr - x) + p*(1-x)
dp = a*(p - pa)
dx = r*x*(xr - x)

# Stimulus
pars={'timeStart':0.0, 'timeMax': 500.0, 'timeStep':1e-3,
      'tauX':75.0,  'tauP':75.0,
      'asynX':1.0, 'asynP':0.3,
      'kP':0.0, 'kX':1.0,
      'jumpP':0.2,
      'stimHz':20.0, 'stimStart':10.0}

#sampTimes=sc.arange(0,pars['timeMax'],pars['timeStep'])
#strHz=r'$\omega_{stim}$=%.0f Hz'%(pars['stimHz'])
#preSpikes=sc.arange(pars['stimStart'], pars['timeMax'], 1/pars['stimHz'])

def alphaCurve(x, A=1.0, tau=1.0, downAccel=1.0):
    aa= sc.int32(x>0)
    xovertau = x/tau
    return A* xovertau * sc.exp( downAccel*(1 - xovertau))

def trainAlpha(samplingTimes, pulseTimes, tauTrain=1.0,downAccel=1.0):
    """
    train= trainGauss(samplingTimes, pulseTimes,bellSpread=1.0)
    """
    nPts= len(samplingTimes)
    train = sc.zeros(nPts)
    alpha=alphaCurve(samplingTimes,A=1.0,tau=tauTrain,downAccel=downAccel)
    for n in range(len(pulseTimes)):
        nn=gr.find(samplingTimes<pulseTimes[n]).max()
        train[nn:] = train[nn:] + alpha[0: nPts-nn]
    return train


def pJump(p,h):
    return h + (1-h) * p

def xJump(x,p):
    #return x - (1-x)*p
    return x*(1+p) -p

def sy_pSolution(t, p0, pAsyn, a):
    return pAsyn- (pAsyn-p0)*sy.exp(-t*a)

def sy_xSolution(t, x0, xAsyn, a):
    denom = x0 + (xAsyn-x0) * sy.exp(-t*xAsyn*a )
    return xAsyn * x0 / denom

def xJump(x,p):
    #return x - (1-x)*p
    return x*(1+p) -p


def bellCurve(x, A=1.0, mu=0.0, sigma=1.0):
    aa= ((x-mu)**2)
    bb= 2* (sigma**2)
    cc= sc.exp(-aa/bb)
    return A * cc

def trainGauss(samplingTimes, pulseTimes,bellSpread=1.0):
    """
    train= trainGauss(samplingTimes, pulseTimes,bellSpread=1.0)
    """
    #nPts= len(samplingTimes)
    train = 0 #sc.zeros(nPts)
    amp = 1/(bellSpread*sc.sqrt(2*sc.pi))
    for n in range(len(pulseTimes)):
        train=train + bellCurve(samplingTimes,A=amp, mu=pulseTimes[n],sigma=bellSpread)
    return train

def plasticity(y, y0, sigma=0.5, eta=1.0):
    e1= sc.exp( (y-y0)*eta)
    ee = e1**(sigma)
    return 2/(ee + ee/e1)


#
def stimTrains(pars, tauTrain=0.1):
    #print 'Obtaining presynaptic spike times and stimulus train'
    pars['sampleTimes']=sc.arange(0,pars['timeMax'],pars['timeStep'])
    pars['preAPs']=sc.arange(pars['stimStart'], pars['timeMax'], 1000.0/pars['stimHz'])
    alphaTrainP=trainAlpha(samplingTimes=pars['sampleTimes'], pulseTimes=pars['preAPs'][1:],tauTrain=1.0,downAccel=1.0)
    alphaTrainX=trainAlpha(samplingTimes=pars['sampleTimes'], pulseTimes=pars['preAPs'],tauTrain=1.0,downAccel=1.0)
    pars['alphaTrainP']= alphaTrainP/alphaTrainP.max()
    pars['alphaTrainX']= alphaTrainX/alphaTrainX.max()
    pars['spikesCa']=lambda t: sc.interp(t,pars['sampleTimes'],pars['alphaTrainP'])
    pars['spikesX']=lambda t: sc.interp(t,pars['sampleTimes'],pars['alphaTrainX'])
    return pars
#
def stFD(U,t,pa):
    x,p,dxp=U
    spikeX = pa['spikesX'](t)
    spikeP = pa['spikesCa'](t)
    dx = (x**pa['kX']) * ( pa['asynX']-x)/pa['tauX'] - spikeX * x * p
    dp = (p**pa['kP']) * ( pa['asynP']-p)/pa['tauP'] + spikeP * pa['jumpP'] * (1-p)
    dxp= p*dx + x*dp
    #if  (spike):
    #    print t, (1-p), pa['jumpP'], (1-p)* pa['jumpP']
    return dx, dp, dxp

def simulateFD(pa):
    # Done setting up the stimulus
    U0=[pa['asynX'], pa['asynP'], pa['asynX']*pa['asynP']]
    tCrit=pa['preAPs']
    orbit= sc.integrate.odeint(stFD,U0,pa['sampleTimes'],args=(pa,), tcrit=tCrit).transpose()
    return orbit


def graphSimulation(ax, pa, orbit):
    xOrbit=orbit[0]; pOrbit=orbit[1]; xpOrbit=orbit[2]
    nStim= len(pa['preAPs'])
    peakInds= sc.zeros(nStim,'int32')
    for nn in range(0,nStim):
        peakInds[nn] = 1+gr.find(pa['sampleTimes']< pa['preAPs'][nn]+tauTr).max()
    firstPeakInd=peakInds[0]
    relAmp= xpOrbit[peakInds[0]]
    xpNorm=xpOrbit[peakInds]/relAmp
    #
    gr.ioff()
    ax.plot(pa['preAPs'], 1.1*sc.ones(len(pa['preAPs'])), 'k|', ms=4, lw=2, label=r'$t_1,...,t_n$')
    ax.plot(pa['sampleTimes'], xOrbit, 'k:', lw=1, alpha=0.5, label=r'$(t,x)$')
    ax.plot(pa['sampleTimes'], pOrbit, 'k--', lw=1, alpha=0.5, label=r'$(t,p)$')
    ax.plot(pa['sampleTimes'], xpOrbit, 'k', lw=1, alpha=0.5,  label=r'$(t,xp)$')
    ax.plot([pa['sampleTimes'][firstPeakInd],], [xpNorm[0],], 'ro', ms=3, alpha=0.99)
    ax.plot([0,pa['timeMax']], sc.ones(2), 'r', lw=1, alpha=1)#, label=r'$t_1,...,t_n$')
    ax.plot( pa['sampleTimes'][peakInds], xpNorm, 'b', lw=2, alpha=0.8)
    ax.plot( pa['sampleTimes'][peakInds], xpNorm, 'wo', ms=5, alpha=1,  label=r'$(t, xp(t))/xp_1$')
    ax.set_ylim(-0.1,sc.maximum(1.5,1.3*xpNorm.max()))
    ax.set_xlim(0,pa['timeMax'])
    ax.set_xlabel('time (ms)')
    ax.legend(loc='upper center',ncol=5,fontsize=10)
    gr.ion(); gr.draw()
    return xpNorm

# Check plot sim
if 1:
    pars['tauX']=150.0
    pars['tauP']=75.0
    pars['asynX']=0.9
    pars['asynP']=0.3
    pars['jumpP']=0.9
    pars['timeMax']=500.
    tauTr=0.1
    pars= stimTrains(pars, tauTrain=tauTr)
    orbit=simulateFD(pars)
    r=1; c=1
    f0=gr.figure(figsize=(7,5))
    ax1=f0.add_subplot(r,c,1)
    xpNorm=graphSimulation(ax=ax1, pa=pars, orbit=orbit)

# Roll  over taus
def rollOverTaus(pars, orbit, tauMin=25.0, tauMax=500.0, tauStep=25.0):
    #tauTr=0.1
    #pars= stimTrains(pars, tauTrain=tauTr)
    tauXs= sc.arange(tauMin, tauMax, tauStep)
    tauPs= sc.arange(tauMin, tauMax, tauStep)
    tauXs[0]=1.0; tauPs[0]=1.0
    mTau= len(tauXs); nTau= len(tauPs)
    figs=list()
    nStim= len(pars['preAPs'])
    peakInds= sc.zeros(nStim,'int32')
    peakInds[0]=1+gr.find(pars['sampleTimes']< pars['preAPs'][0]+tauTr).max()
    for nn in range(1,nStim):
        peakInds[nn] = 1+gr.find(pars['sampleTimes']< pars['preAPs'][nn]+tauTr).max()

    for m in range(mTau):
        pars['tauX']=tauXs[m]
        f=gr.figure(figsize=(23,11))
        gr.ioff()
        r=5; c=sc.ceil(sc.float32(nTau)/r)
        ax1=list()
        figStr1=r'$\tau_x$=%.0f, $\tau_p \in $ [%.0f,%.0f], $h$=%.1f, $x_{\infty}$=%.1f, $p_{\infty}$=%.1f'%(tauXs[m], tauPs[0], tauPs[-1], pars['jumpP'], pars['asynX'], pars['asynP'])
        f.suptitle(figStr1, fontsize=13)
        for n in range(nTau):
            ax1.append(f.add_subplot(r,c,n+1))
            figStr2=r'$(\tau_x,\tau_p)$'
            figStr3=r'(%.0f, %.0f)'%(tauXs[m], tauPs[n])
            figName='taux=%.0f_taup=%.0f_h=%.1f_xa=%.1f_pa=%.1f'%(tauXs[m], tauPs[n], pars['jumpP'], pars['asynX'], pars['asynP'])
            print figName
            pars['tauP']=tauPs[n]
            orbit=simulateFD(pars)
            xpNorm=graphSimulation(ax=ax1[n], pa=pars, orbit=orbit)
        f.savefig('jj/'+figName+'.png')
        gr.ion(); gr.draw()
        gr.pause(2)
        gr.close('all')
    return

# Roll taus and jumps
if 1:
    hArray=sc.arange(0.1,1,0.2)
    pars['asynX']=0.9
    pars['asynP']=0.3
    pars['jumpP']=0.9
    pars['timeMax']=1000.0
    tauTr=0.1
    pars= stimTrains(pars, tauTrain=tauTr)
    for thisH in hArray:
        pars['jumpP']=thisH
        orbit=simulateFD(pars)
        rollOverTaus(pars, orbit, tauMin=0.0, tauMax=100.0, tauStep=10.0)
