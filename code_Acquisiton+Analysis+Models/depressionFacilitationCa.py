# coding: utf-8

# # Ca-dependent short-term synaptic plasticity
#
# Presynaptic release
# \begin{eqnarray}
# \partial_t c &=& \frac{c_m - c}{\tau_c}  - \phi(t)
# \\
# \partial_t p &=& \frac{p_{\infty}(c) - p}{\tau_p}
# \\
# \partial_t q &=& q \left( \frac{q_{\infty} - q}{\tau_q} \right) (1- p)
# \end{eqnarray}
#
# Stimulation functions
# \begin{eqnarray}
# \phi(t) &=& \sum_{i} \alpha(t-t_i; \tau)
# \\
# \alpha(s) &=& H(s)\frac{s}{\tau} \exp\left( 1-\frac{s}{\tau}\right),
# \end{eqnarray}
# with $H(s)= 1$ if $s\geq 0$ and $0$ if $s<0$.
#

# In[11]:

import scipy as sc
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import matplotlib.pylab as gr
#get_ipython().magic('matplotlib inline')


# Stimulation functions

def alphaFunction(x, A=1.0, tau=1.0, downAccel=1.0):
    """
    alphaFunction creates an alpha function with amplitude A, time constant tau,
    and downward acceleration downAccel.
    Example:
    alphaFunction(x, A=1.0, tau=1.0, downAccel=1.0)
    """
    aa= sc.int32(x>0)
    xovertau = x/tau
    return A* aa* xovertau * sc.exp( downAccel*(1 - xovertau))

def trainAlpha(samplingTimes, pulseTimes, stimAmp=1.0, tauTrain=1.0,downAccel=1.0):
    """
    Example:
    a= trainAlpha(samplingTimes, pulseTimes, tauTrain=1.0,downAccel=1.0)
    """
    nPts= len(samplingTimes)
    train = sc.zeros(nPts)
    alpha=alphaFunction(samplingTimes,A=stimAmp,tau=tauTrain,downAccel=downAccel)
    for n in range(len(pulseTimes)):
        #nn=gr.find(samplingTimes<pulseTimes[n]).max()
        nn=sc.where(samplingTimes<pulseTimes[n])[0].max()
        train[nn:] = train[nn:] + alpha[0: nPts-nn]
    return train

def F(x):
    ee= sc.exp(x)
    return ee/(1+ee)


# Parameters
p={"g_pR": 100.0, "q_tau": 30.0, "c_tau":50.0, "tau_pulse":0.2,
   "q_ss":0.9, "c_ss": 1e-3, "stim_Amp": 0.005,
   "c_pR":0.25, "tau_pR": 1.0}
#
p["timeMax"]=600.0; p["timeStep"]=1e-2
sampTimes=sc.arange(0,p["timeMax"],p["timeStep"])
p["sampTimes"] = sampTimes
#
p["stimStart"]=60.0; p["stimEnd"]=500.0; p["stimInterval"]=50.0;
p["stimTimes"] = sc.arange(p["stimStart"], p["stimEnd"], p["stimInterval"])

def createTrainFunction(p):
    p["stimTrain"] = trainAlpha(samplingTimes=p["sampTimes"],
        pulseTimes=p["stimTimes"],
        stimAmp=p["stim_Amp"],
        tauTrain=p["tau_pulse"])
    return lambda t: sc.interp(t, p["sampTimes"], p["stimTrain"])

# Graph for the train of alpha pulses
p["phi"]=createTrainFunction(p)
f1=gr.figure(figsize=(15,6))
r=1; c=2; gr.ioff()
ax1=f1.add_subplot(r,c,1)
ax2=f1.add_subplot(r,c,2)
ca=sc.arange(0,1,0.001)
ax1.plot(ca, F((ca-p["c_pR"])*p["g_pR"]))
ax2.plot(sampTimes, p["phi"](sampTimes))
ax1.set_xlabel("$[Ca]_i$ ($\mu$ M)")
ax1.set_ylabel("$p_{\infty}([Ca]_i)$ ($\mu$ M)")
ax2.set_xlabel("time (msecs)")
ax2.set_ylabel("$\mu$ M / msecs")
gr.ion(); gr.draw()
gr.show()

# System dynamics.
# Phi represents the increments in Ca from the activation of Ca-currents
def preSynCa(u,t,par):
    dc = (par["c_ss"]-u)/par["c_tau"] + par["phi"](t)
    #print(t,u,dc)
    return dc

def preSynCaSTSP2(U,t,par):
    c,q=U
    dc = (par["c_ss"]-c)/par["c_tau"] + par["phi"](t)
    dq =  (q * (par["q_ss"] - q) /par["q_tau"]) - (par["q_ss"]-q) * F((c-p["c_pR"])*p["g_pR"])
    return sc.array([dc, dq])

def RK4(f,p):
    return lambda y, t, dt: (
            lambda dy1: (
            lambda dy2: (
            lambda dy3: (
            lambda dy4: (dy1 + 2*dy2 + 2*dy3 + dy4)/6
            )( dt * f( y + dy3, t + dt  , p   ) )
        )( dt * f( y + dy2/2, t + dt/2, p ) )
        )( dt * f( y + dy1/2, t + dt/2, p ) )
        )( dt * f( y, t, p ) )

def solveRK4(p):
    U=sc.zeros((p["nSteps"],p["nDim"]),"float32")
    dU = RK4(p["rhs"],p)
    U[0]=p["ic"]
    for n in sc.arange(1,p["nSteps"]-1):
        U[n] = U[n-1] + dU( U[n-1], p["sampTimes"][n-1], p["timeStep"])
        if (p["sampTimes"][n]%1000)<0.001:
            print("t=%g, y=%g, w=%g"%(p["sampTimes"][n], U[n][0],U[n][1]))
    return U.transpose()

# use tain of alpha pulses for the
p["c_ss"]=1e-3 #in muM
p["c_tau"]=50.0; p["tau_pulse"]=0.2 #in muSec, gives pulses of about 1 msec
p["stim_Amp"]=0.2
p["q_ss"]= 0.9; p["q_tau"]=100.0
p["c_pR"]=0.175; 
p["g_pR"]=70.0
p["pRss"]= lambda c: F((c-p["c_pR"])*p["g_pR"])
p["timeMax"]=600.0; p["timeStep"] = 1e-3
sampTimes=sc.arange(0,p["timeMax"],p["timeStep"])
p["sampTimes"] = sampTimes
p["sampTimes"]= sc.arange(0,p["timeMax"],p["timeStep"])
p["nSteps"]=len(p["sampTimes"])
p["ic"]=[0.01,0.9]
p["nDim"]=len(p["u0"])
p["stimTimes"] = sc.arange(p["stimStart"], p["stimEnd"], p["stimInterval"])
p["phi"] = createTrainFunction(p)
p["rhs"] = preSynCaSTSP2
cOrbit,qOrbit=solveRK4(p)
#CaOrbit1= odeint(preSynCa,0.01,p["sampTimes"][:-1],(p,),tcrit=p["stimTimes"]).transpose()[0]
#tCrit= sc.hstack((p["stimTimes"],p["stimTimes"]-(1e-5)))
#tCrit.sort()
#cOrbit,qOrbit= odeint(p["rhs"],p["ic"], p["sampTimes"], (p,), tcrit=tCrit,rtol=1e-9,atol=1e-3).transpose()
pOrbit = F((cOrbit-p["c_pR"])*p["g_pR"])
ca=sc.arange(0,1,0.001)
#
f2= gr.figure(figsize=(11,9)); gr.ioff()
axCa=f2.add_subplot(2,2,1)
axPR=f2.add_subplot(2,2,2)
axqP=f2.add_subplot(2,2,3)
ax4=f2.add_subplot(2,2,4)
axCa.plot(p["sampTimes"],cOrbit,label=r"$(t,[Ca]_i)$")
axPR.plot(pOrbit,cOrbit, label=r"$(t,P_{rel}([Ca]_i))$")
axPR.plot( F((ca-p["c_pR"])*p["g_pR"]),ca)
axqP.plot(p["sampTimes"],pOrbit,label=r"$(t,P_{rel}([Ca]_i))$")
axqP.plot(p["sampTimes"],qOrbit,label=r"$(t,q)$")
axqP.plot(p["sampTimes"],pOrbit*qOrbit,label=r"$(t,q P_{rel})$")
axCa.set_ylim(-0.01,1.01)
axqP.set_ylim(-0.01,1.01)
axCa.legend(); axqP.legend()
gr.ion(); gr.draw();  gr.show()


# In[19]:
p["tau_pR"]= 1.0; p["tau_q"]=30.0; p["tau_c"]=20.0; p["tau_syn"]=10.0,
p["g_pR"]= 1.0/100.0;
p["q_ss"]=0.9; p["c_ss"]= 1e-3; p["stim_Amp"]= 0.001,
p["c_p"]=0.5
p["ic"]= sc.array([0.001, 0.1, 0.6 ])
p["nDim"]=3; p["nSteps"]=len(p["sampTimes"])
p["timeMax"]=200.0; p["timeStep"]=1e-3
p["sampTimes"] = sc.arange(0,p["timeMax"],p["timeStep"])
p["rhs"]=preSynCaSTSP2
#
orbit2=solveRK4(p)
ca1=orbit1[0];ca2=orbit2[0]
pR1=orbit1[1];pR2=orbit2[1]
q1=orbit1[2];q2=orbit2[2]
pRq1= pR1*q1;pRq2= pR2*q2
#
print(len(p["sampTimes"]),ca1.shape, pR1.shape, q1.shape, pRq1.shape)


# In[10]:

f2=gr.figure(figsize=(15,11))
r=3; c=2; gr.ioff()
ax=list()
labs=[r"$c$", r"$p$", r"$q$"]
for n in sc.arange(r*c):
    ax.append(f2.add_subplot(r,c,n+1))

for n in sc.arange(r):
    #ax[c*n].plot(p["sampTimes"],orbit[n],label=labs[n])
    ax[c*n].plot(orbit1[n],label=labs[n])
    ax[c*n].plot(orbit2[n],label=labs[n])

ax[2].plot(p["sampTimes"], pRq,label=r"$pq$")
ax[4].plot(p["sampTimes"], pRq,label=r"$pq$")

ax[1].plot(pR, ca, label=r"$(p,c)$")
ax[3].plot(ca, pR, 'b', label=r"$(c,p)$")
ax[3].plot(ca, pRq, 'g', label=r"$(c,qp)$")
ax[5].plot(ca, q, 'b', label=r"$(c,q)$")
ax[5].plot(ca, pRq, 'g', label=r"$(c,qp)$")


for n in sc.arange(r):
    ax[c*n +1].set_xlim(0,1)

for n in sc.arange(r*c):
    ax[n].set_ylim(0,1)
    ax[n].legend()

gr.ion(); gr.draw();


# In[ ]:
