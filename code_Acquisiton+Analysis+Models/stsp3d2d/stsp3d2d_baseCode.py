# MAHV. 20180909
# Short-term plasticity model

import scipy as sc
import matplotlib.pylab as gr

def RK2_Autonomous(f,pars, eParNames=[], eParList=[]):
    """Second-order Runge-Kutta method to solve x' = f(x,t) with U(t[0]) = U0.
    NOTES:
        This version is based on the algorithm presented in "Numerical
        Analysis", 6th Edition, by Burden and Faires, Brooks-Cole, 1997.
    """
    U=sc.zeros((pars['nSteps'], sc.prod(sc.shape(pars['ic']))),"float64")
    U[0]=pars['ic']
    for i in range(pars['nSteps']-1):
        for k in sc.arange(len(eParNames)):
            pars[eParNames[k]]= eParList[k][i]
        #t = stepSize * i
        k1 = pars['stepSize'] * f( U[i], pars) / 2.0
        U[i+1] = U[i] + pars['stepSize'] * f( U[i] + k1, pars)
    return U.transpose()

def linearDynamics(x, pp):
    return pp['forcing'] + (pp['ss']-x)/pp['tau']

def logisticDynamics(x, pp):
    return pp['forcing'] + (x**pp['le'])*(pp['ss']-x)/pp['tau']

def stsp3D(U, pp):
    c,p,x=U
    dc = pp['Flux_Ca'] + (pp['ss_Ca']-c)/pp['tau_Ca']
    #dp = pp['alpha_p']* c * (1-p) - pp['beta_p']*p
    dp = pp['alpha_p']* ( c * (1-p) - pp['beta/alpha_p'] * p)
    dx = (x**pp['le_x'])*(pp['ss_x']-x)/pp['tau_x'] - x*p
    return sc.array([dc,dp,dx])

def ss_p(c,pp):
    return c/(c+ pp['beta/alpha_p'])

def stsp2D(U, pp):
    c,x=U
    dc = pp['Flux_Ca'] + (pp['ss_Ca']-c)/pp['tau_Ca']
    dx = (x**pp['le_x'])*(pp['ss_x']-x)/pp['tau_x'] - x*ss_p(c,pp)
    return sc.array([dc,dx])


#
def stsp3d_profile(pa):
    lO= RK2_Autonomous(f= stsp3D, pars=pa, eParNames=['Flux_Ca'],eParList=[ff])
    c=lO[0]; p=lO[1]; x=lO[2]
    #
    f3d = gr.figure(figsize=(11,7)); gr.ioff()
    rows = 3; cols = 2; ax=list()
    for n in range(rows*cols):
        ax.append(f3d.add_subplot(rows,cols,n+1))
    #
    ax[0].plot(timeSamples, c,'.',ms=1.0, color='darkgreen', label=r'$p$')
    ax[2].plot(timeSamples, p,'b.',ms=1.0, label=r'$p$')
    ax[2].plot(timeSamples, p*x,'k.',ms=1.0, label=r'$px$')
    ax[2].plot(timeSamples, ss_p(c,pa),'b', alpha=0.3, lw=3.0, label=r'$p_{\infty}$')
    ax[4].plot(timeSamples, x,'.',ms=1, color='orange',label=r'$x$')
    ax[1].plot(c, p,'b.',ms=1.0, label=r'$(c,p)$')
    ax[1].plot(c, p*x,'k.',ms=1.0, label=r'$(c,px)$')
    ax[1].plot(c, x,'.',ms=1,color='orange',label=r'$(c,x)$')
    ax[3].plot(c, p*x,'k.',ms=1.0, label=r'$(c,px)$')
    ax[3].plot(c, p,'b.',ms=1.0, label=r'$(c,p)$')
    ax[5].plot(c, x,'.',ms=1.0, color='orange', label=r'$(c,x)$')
    #
    ax[0].set_ylim(0,x.max()); ax[1].set_ylim(0,x.max())
    for n in range(rows*cols):
        ax[n].legend()
    gr.ion(); gr.draw()
    return c,p,x,f3d

tau_p = lambda c,pp: 1/(c*pp['alpha_p'] + pp['beta_p'])
U0 = sc.array([0.02,0.01,0.9])
pa={'Flux_Ca':0.0,'tau_Ca':20.0,'ss_Ca':0.01,'alpha_p':1.0,'beta_p':10.0,'ss_x':1.0, 'tau_x':50.0, 'le_x':1.0}
pa['beta/alpha_p']= pa['beta_p']/pa['alpha_p']
pinfty = U0[0] /( U0[0] + pa['beta/alpha_p'])
print('p_infty = %g'%pinfty)
timeStep=1/40.0; timeMax= 550.
nSteps = sc.int32(timeMax/timeStep)
timeSamples= sc.arange(0,timeMax,timeStep)
stimInterval = 50.0
stimTimes = stimInterval * sc.arange(1,11)
stimIdx = sc.int32( stimTimes/ timeStep)
print(stimIdx)
ff = sc.zeros(len(timeSamples))
ff[stimIdx]= 200.0 * timeStep
pa['ic']= U0; pa['stepSize']=timeStep; pa['timeMax']= timeMax; pa['nSteps']= nSteps
c,p,x,f3d= stsp3d_profile(pa)



# Runs on 2D
U0 = sc.array([0.05,0.5])
pa['ic']=U0
lO= RK2_Autonomous(f= stsp2D, pars=pa, eParNames=['Flux_Ca'],eParList=[ff])
c=lO[0]; x=lO[1]
p = ss_p(c,pa)
#
f0 = gr.figure(figsize=(11,7)); gr.ioff()
rows = 2; cols = 2; ax=list()
for n in range(rows*cols):
    ax.append(f0.add_subplot(rows,cols,n+1))

ax[0].plot(timeSamples, c,'.',ms=1.0, color='darkgreen', label=r'$p$')
ax[0].plot(timeSamples, p,'b.',ms=1.0, label=r'$p$')
ax[0].plot(timeSamples, p*x,'k.',ms=1.0, label=r'$px$')
ax[0].plot(timeSamples, x,'.',ms=1, color='orange',label=r'$x$')
ax[1].plot(lO[0], p,'b.',ms=1.0, label=r'$(c,p)$')
ax[1].plot(lO[0], p*x,'k.',ms=1.0, label=r'$(c,px)$')
ax[1].plot(lO[0], x,'.',ms=1,color='orange',label=r'$(c,x)$')
ax[2].plot(timeSamples, p*x,'k.',ms=1.0, label=r'$(c,px)$')
ax[2].plot(timeSamples, p,'b.',ms=1.0, label=r'$p$')
ax[2].plot(timeSamples, ss_p(c,pa),'b', alpha=0.3, lw=3.0, label=r'$p$')
ax[3].plot(lO[0], p*x,'k.',ms=1.0, label=r'$(c,px)$')
ax[3].plot(lO[0], p,'b.',ms=1.0, label=r'$(c,p)$')

ax[0].set_ylim(0,x.max()); ax[1].set_ylim(0,x.max())
for n in range(rows*cols):
    ax[n].legend()

gr.ion(); gr.draw()
