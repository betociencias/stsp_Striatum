import scipy as sc
import matplotlib.pylab as gr

cPre={'c_PreInfty_muM':0.0001,'tau_c_Pre':50.0,'c_conv':0.0001,'c_halfAct_r_muM':0.3}
rPre={'alpha_r':2.0,'c_coop_r':4.0}
qPre={'q_Infty':1.0,'tau_q':25.0}
p=dict(); p.update(cPre);p.update(rPre);p.update(qPre)
print(p)

def RK2_Autonomous(f,pars, eParNames=[], eParList=[]):
    """Second-order Runge-Kutta method to solve x' = f(x,t) with U(t[0]) = U0.
    NOTES:
        This version is based on the algorithm presented in "Numerical
        Analysis", 6th Edition, by Burden and Faires, Brooks-Cole, 1997.
    """
    U=sc.zeros((pars['nSteps'], sc.prod(sc.shape(pars['ic']))),"float64")
    U[0]=pars['ic']
    for i in range(pars['nSteps']-1):
        for k in range(len(eParNames)):
            pars[eParNames[k]]= eParList[k][i]
        #t = stepSize * i
        k1 = pars['stepSize'] * f( U[i], pars) / 2.0
        U[i+1] = U[i] + pars['stepSize'] * f( U[i] + k1, pars)
    return U.transpose()

def DiracComb(p,spikeTrain):
    pulseInds = sc.int32(spikeTrain[spikeTrain<p['tMax']]/p['stepSize']);
    comb = sc.zeros(p['nSteps']);  comb[pulseInds]=p['pulseAmp']
    return comb

def gammaTrain(meanISI, shape=2.0, nSpikes=100):
    isiSample= sc.zeros(nSpikes)
    isiSample=sc.random.gamma(shape, scale=meanISI/shape, size=nSpikes)
    spikeTrain= isiSample.cumsum()
    return spikeTrain

def Hill(u,c,n=2):
    uu = u**n
    return uu/(uu+ c**n)

def presynCalcium(c,p):
    dc = (p['c_PreInfty_muM'] - c)/p['tau_c_Pre'] + p['c_PreIn']/p['stepSize']
    return sc.array([dc])

def presynSTSP(U,p):
    c,r,q= U
    rSS = Hill(c,p['c_halfAct_r_muM'],p['c_coop_r'])
    tau_r = rSS/p['alpha_r']
    dc = (p['c_PreInfty_muM'] - c)/p['tau_c_Pre'] + p['c_PreIn']/p['stepSize']
    dr = (rSS-r)/tau_r
    dq = (p['q_Infty']-q)/p['tau_q'] - r*q/p['stepSize']
    return sc.array([dc,dr,dq])

def presynSTSP_2D(U,p):
    c,q= U
    rSS = Hill(c,p['c_halfAct_r_muM'],p['c_coop_r'])
    dc = (p['c_PreInfty_muM'] - c)/p['tau_c_Pre'] + p['c_PreIn']/p['stepSize']
    dq = (p['q_Infty']-q)/p['tau_q'] - rSS*q/p['stepSize']
    return sc.array([dc,dq])

# Function to plot the 2D dynamics
def solve2DDynamics(p,spikeTrain,dispFig=1):
    p['ic'] = sc.array([0.0001,p['q_Infty']]);
    pulses=DiracComb(p,spikeTrain)
    c2,q2=RK2_Autonomous(f=presynSTSP_2D,pars=p, eParNames=['c_PreIn'], eParList=[pulses])
    rSS2 = Hill(c2,p['c_halfAct_r_muM'],p['c_coop_r'])
    nt=q2*rSS2
    if dispFig==1:
        f1=gr.figure(figsize=(15,7)); rows=3; cols=1; ax=list()
        ax=[f1.add_subplot(rows,cols,n+1) for n in range(rows*cols)]
        ax[0].plot(p['sampTimes'],c2,'orange'+'.',label='$c(t)$')
        ax[0].plot(p['sampTimes'],rSS2,'b.',lw=3,alpha=0.4,label='$r_{\infty}(c)$')
        ax[1].plot(p['sampTimes'],q2,'k.',label='$q(t)$')
        ax[2].plot(p['sampTimes'],nt,'k.',label='$r_{\infty}(c)q(t)$')
        ax[0].plot(spikeTrain, sc.maximum(c2.max(),rSS2.max())*sc.ones(len(spikeTrain)),'r|', ms=10)
        ax[1].plot(spikeTrain, sc.ones(len(spikeTrain)),'r|', ms=10)
        ax[2].plot(spikeTrain, -0.01*sc.ones(len(spikeTrain)),'r|', ms=10)
        [ax[n].legend() for n in range(rows*cols)]
    return c2,q2,rSS2
