# -*- coding: utf-8 -*-
import scipy as sc
output=sc.test('all',raise_warnings='release')
import pylab as gr
gr.ion()
# Importar los datos en archivos para matlab
import scipy.io as sio

def bellCurve(x, A=1.0, mu=0.0, sigma=1.0):
    aa= ((x-mu)**2)
    bb= 2* (sigma**2)
    cc= sc.exp(-aa/bb)
    return A * cc

def trainGauss(samplingTimes, pulseTimes,bellSpread=0.5):
    """
    train= trainGauss(samplingTimes, pulseTimes,bellSpread=1.0)
    """
    #nPts= len(samplingTimes)
    train = 0 #sc.zeros(nPts)
    for n in range(len(pulseTimes)):
        amp = 1/(bellSpread*sc.sqrt(2*sc.pi))
        train=train + bellCurve(samplingTimes,A=amp, mu=pulseTimes[n]+3*bellSpread,sigma=bellSpread)
    return train

def plasticity(h, h0, sigma=0.5, eta=1.0):
    e1= sc.exp( -(h-h0)*eta)
    ee = e1**(sigma)
    return 2/(ee + ee/e1)

def stf(U,t,pars):
    tau_x, tau_p, xinf, pinf, q,  spikes =pars
    x,p=U
    #q = plasticity(h=20.0, h0=40.0)
    spike = trainGauss(samplingTimes=t,pulseTimes=spikes)
    #print spike
    dx = x*(xinf-x)/tau_x - spike * p * x 
    #dp = (pinf-p)/tau_p + spike * q * (1-p)
    dp = (pinf-p)/tau_p + spike * q *(1-p) #dice marco que te vayas para arriba y no estés chingando
    if spike:
        print p,1-p,q, spike * q * (1-p)
    return dx, dp

#
stimHz=20.0
strCa=r'$[Ca]_{pre}$'
strHz=r'$\omega_{stim}$=%.0f Hz'%(stimHz)
firstSpike=10.0; period=1000.0/stimHz
timeStart=0.0; timeMax = 500.0; timeStep=1e-3
sampTimes=sc.arange(timeStart,timeMax,timeStep)
preSpikes=sc.arange(firstSpike, timeMax, period)
bellTrain= trainGauss(samplingTimes=sampTimes, pulseTimes=preSpikes,bellSpread=0.1)
#gr.figure(); gr.plot(sampTimes, bellTrain)
tauFac=500.0; pinf=1; xinf=1.0; q=0.2; pinf2=pinf*10; 
U0=[1,0.2]
tCrit=preSpikes

def tauRecSweep(tauRecs=5.0*sc.arange(1.0,2.00,0.2), tauFaci=500.0, pinf=1.0, xinf=1.0, q=0.2,  stimHz=20.0, firstSpike=10.0, timeStart=0.0, timeMax = 500.0, timeStep=1e-3, U0=[1,0.2],embedLaTeX=0,saveFig=0):
    period=1000.0/stimHz
    sampTimes=sc.arange(timeStart,timeMax,timeStep)
    preSpikes=sc.arange(firstSpike, timeMax, period)
    # bellTrain= trainGauss(samplingTimes=sampTimes, pulseTimes=preSpikes,bellSpread=1.0)
    orbits=list(); figs=list()
    strCa=r'$[Ca]_{pre}$'
    strFac=r'$\tau_{fac}$=%2.0f msec'%(tauFaci)
    strHz2='omegaStim%2.0fHz'%(stimHz)
    strFac2='tauFac%2.0fmsec'%(tauFaci)
    strPinf=r'$\pinf2$=%2.0f'%(pinf)
    strPinf2='pinf2%2.0f'%(pinf)
    strPinf3 ='pinf{0}'.format(pinf)
    print 'Performing simulations'
    #print tauRecs
    nSimulations=len(tauRecs)
    strRecFirst='tauRecs[0]%2.0fmsec'%(tauRecs[0])
    strRecLast='tauRecs[-1]%2.0fmsec'%(tauRecs[-1])
    nPulses=len(preSpikes)
    for m in range(nSimulations):
        print '%d of %d'%(m+1,nSimulations)
        tau_rec=tauRecs[m]
        strRec=r'$\tau_{rec}$=%2.0f msec'%(tau_rec)
        strRec2='tauRec%2.0fmsec'%(tau_rec)
        pars= [tau_rec, tauFaci, xinf, pinf, q, preSpikes]
        orbit= sc.integrate.odeint(stf,U0,sampTimes,args=(pars,), tcrit=tCrit).transpose()
        orbits.append(orbit)
        xOrbit=orbit[0]
        pOrbit=orbit[1]
        xp=xOrbit*pOrbit
        startInd=1+gr.find(sampTimes<firstSpike).max()
        periodInd=1+gr.find(sampTimes<period).max()
        pulseInds= startInd+ sc.arange(0,nPulses*periodInd, periodInd)
        xpNorm = xp[pulseInds]/xp[startInd]
        #
        thisFig=gr.figure(figsize=(11,5))
        figs.append(thisFig)
        r=2; c=1
        ax1=thisFig.add_subplot(r,c,1)
        ax2=thisFig.add_subplot(r,c,2)
        gr.ioff()
        ax2.plot(sampTimes[pulseInds], xpNorm,'k',alpha=1.0, label='')
        ax1.plot(sampTimes, 0.2*bellTrain - 0.1,'k',alpha=1.0, label=strCa+', '+strHz)
        ax1.plot(sampTimes, orbit[0],'b.', ms=1.0, alpha=0.8, label='Occupancy ($x$),' + ' ' +strRec )
        ax1.plot(sampTimes, orbit[1],'r.', ms=1.0, alpha=0.8, label='P(release) ($p$),' + ' ' +strFac )
        ax1.plot(sampTimes, orbit[0]*orbit[1],'k', lw=2.0, alpha=1.0, label='max release $(xp)$')
        ax1.plot(sampTimes, pinf*sc.ones(len(sampTimes)),'k:', lw=1.0, label='initial max release $(x_0 p_0)$')
        ax1.set_ylim(-0.11,1.1)
        ax1.set_xlabel('time (msec)')
        ax1.legend(loc='upper right',fontsize=10.0)
        gr.ion(); gr.draw()
        figDir='C:\Users\Janet Barroso\Dropbox\Modelo\BarrosoHerrera\Figuras\\'
        #baseName='q0dot%2d'%(100*q)+strHz2+'_'+'_'+strFac2+strRec2
        baseName=strFac2+strRec2+strPinf3
        # sio.savemat('sampTimes'+'q0dot%2d'%(100*q)+strHz2+'_'+'_'+'strFac2,{'sampTimes', sampTimes})
        figName=figDir+baseName+'.png'
        matDir='C:\Users\Janet Barroso\Dropbox\Modelo\BarrosoHerrera\Matrices\\'
        matlabName=matDir+baseName+'.mat'
        # Opciones: sio.loadmat, sio.savemat, sio.whosmat
        # Abrir archivo para escritura en un objeto f
        f = open('matlabName', 'w')
        #Hacer carpeta para las matrices
        #sio.savemat(matlabName, {'orbit':orbit})
        f.close()

        if embedLaTeX:
            #print r'\newpage'
            print r'\begin{figure}[h]'
            print r'\includegraphics[width=\textwidth]{%s}'%figName
            print r'\end{figure}'
        if saveFig:
            gr.savefig(figName)
            gr.close('all')
    timeName='sampTimes'+'q0dot%2d'%(100*q)+strHz2+'_'+'_'+strFac2+strRecFirst+'-'+strRecLast
    #sio.savemat(timeName,{'sampTimes',sampTimes})
    return orbits, figs

#print "Ran from LaTeX"


if 0:
    orbits, fig100= tauRecSweep(q=0.2, tauRecs=500.0*sc.arange(1.0,1.001,0.2), tauFaci=5.0, pinf=0.2, stimHz=20.0, firstSpike=10.0, timeStart=0.0, timeMax = 500.0, timeStep=1e-3, U0=[1,0.2], embedLaTeX=0, saveFig=0)

# Subió y bajó (poquito)
if 0:
    orbits, fig100= tauRecSweep(q=1.29, tauRecs=75.0*sc.arange(1.0,1.001,0.2), tauFaci=100.0, pinf=0.6, stimHz=20.0, firstSpike=10.0, timeStart=0.0, timeMax = 500.0, timeStep=1e-3, U0=[1,0.2], embedLaTeX=0, saveFig=0)

#Deprime
#orbits, fig100= tauRecSweep(q=0.5, tauRecs=1000.0*sc.arange(1.0,1.001,0.2), tauFaci=10.0, pinf=0.2, stimHz=20.0, firstSpike=10.0, timeStart=0.0, timeMax = 500.0, timeStep=1e-3, U0=[1,0.2], embedLaTeX=0, saveFig=0)
#Biphasic
#orbits, fig100= tauRecSweep(q=0.5, tauRecs=100.0*sc.arange(1.0,1.001,0.2), tauFaci=50.0, pinf=0.6, stimHz=20.0, firstSpike=10.0, timeStart=0.0, timeMax = 500.0, timeStep=1e-3, U0=[1,0.2], embedLaTeX=0, saveFig=0)
#La que sigue facilita
orbits, fig100= tauRecSweep(q=0.5, tauRecs=10.0*sc.arange(1.0,1.001,0.2), tauFaci=100.0, pinf=0.3, stimHz=20.0, firstSpike=10.0, timeStart=0.0, timeMax = 500.0, timeStep=1e-3, U0=[1,0.2], embedLaTeX=0, saveFig=0)

if 0:
    qs= sc.arange(0.2,1.0,0.3)
    orbits=list()

    for m in range(len(qs)):
        orbits.append(list())
        orbits[m], figs= tauRecSweep(tauRecs=200.0*sc.arange(1.0,1.001,0.2), tauFaci=10.0, pinf=0.4, q=qs[m], stimHz=20.0, firstSpike=10.0, timeStart=0.0, timeMax = 500.0, timeStep=1e-3, U0=[1,0.2], embedLaTeX=0, saveFig=1)



