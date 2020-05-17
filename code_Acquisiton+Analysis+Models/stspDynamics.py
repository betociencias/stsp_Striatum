from lowDExcitable import * 
gr.ion()

p=dict()
# Auxiliary functions
ssP= lambda y,pp: nExpSigmoid(pp['gainP']*(y-pp['yHalfP']))
rateP= lambda y,pp: pp['rateP'] * nExpSum(pp['gainP']*(y-pp['yHalfP']),s=pp['biasP'])
#p["ssP"]= ssP; p["rateP"]= rateP


def STSP1(U,t,par):
    dU = (U ** par["expP"]) * (ssP(y,par) - U ) * rateP(y,par)
    return dU

def STSP2D(U,t,par,y):
    p,q=U
    uu= y(t)
    dp = (p ** par["expP"]) * (ssP(uu,par) - p) * rateP(uu,par)
    #dq = (q ** par["expQ"]) * (par["ssQ"] - q * (1+p) ) * par["rateQ"] 
    dq = (q ** par["expQ"]) * (par["ssQ"] - q) * par["rateQ"] - p * q
    return dp,dq

# Analysis of static trayectories for the system
def staticSTSP(p,timeMax=800.0, vC = sc.arange(-70,11,10.0), vT=26.7):
    nSims=len(vC)
    pOrbits= list(); qOrbits= list()
    sampTimes=sc.arange(0,timeMax,p["timeStep"])
    vss0 = sc.arange(vC.min(),vC.max(),0.1)
    yss0= vss0/vT; yC= vC/vT
    for n in sc.arange(nSims):
        y=lambda t: yC[n]
        orbit= sc.integrate.odeint(STSP2D, [0.15, 0.9], sampTimes, args=(p,y)).transpose()
        pOrbits.append(orbit[0])
        qOrbits.append(orbit[1])

    fig=gr.figure(figsize=(15,11)); gr.ioff()
    ax=list(); rows=1; cols=2
    for n in sc.arange(rows*cols):
        ax.append(fig.add_subplot(rows,cols,n+1))
    for n in sc.arange(nSims):
        ax[0].plot(sampTimes, pOrbits[n],label=r"p(t,%g)"%(vss1[n]))
        ax[0].plot(sampTimes, qOrbits[n],label=r"q(t,%g)"%(vss1[n]))
        ax[0].plot(sampTimes, pOrbits[n]*qOrbits[n],label=r"$pq$")
        ax[1].plot(pOrbits[n], qOrbits[n],label=r"(p,q)")
        ax[1].plot(pOrbits[n], pOrbits[n]*qOrbits[n], 'r', label=r"(p,qp)")

    ax[1].plot(vss0, ssP(yss0,p),label=r"$F_p(v)$")
    ax[1].plot(vss0, 1/rateP(yss0,p),label=r"$\tau_p(v)$")
    ax[1].plot(vC, ssP(yC,p),'ko')
    ax[1].plot(vC, 1/rateP(yC,p),'ko')

    for n in sc.arange(rows*cols):
        ax[n].legend(ncol=4,fontsize=10,loc="lowerleft")
    #ax[0].set_ylim(-0.1,1.1)
    gr.ion(); gr.draw()
    return pOrbits, qOrbits


#staticSTSP()

def dynamicSTSP(p, ax, ic=[0.2,0.9], timeMax=600.0, vT=26.7, pulses=100+sc.arange(0,500,50.0),stri=""):
    sampTimes=sc.arange(0,timeMax,p["timeStep"])
    yy= trainAlpha(samplingTimes=sampTimes, pulseTimes=pulses, tauTrain=1.0,downAccel=1.0) 
    y = lambda t: 3*sc.interp(t, sampTimes,yy) - 60.0/vT
    orbits= sc.integrate.odeint(STSP2D, [0.25, 0.8], sampTimes, args=(p,y)).transpose()
    pOrbit=orbits[0]; qOrbit=orbits[1]
    #
    ax[0].plot(sampTimes, (yy-yy.min())/(yy.max()-yy.min()),'k',alpha=0.3)
    ax[0].plot(sampTimes, pOrbit,'r',label=r"$p(t)$")
    ax[0].plot(sampTimes, qOrbit,'b',label=r"$q(t)$")
    ax[0].plot(sampTimes, pOrbit*qOrbit,'k', label=r"$pq$")
    ax[1].plot(pOrbit, qOrbit,'b', label=r"$(p,q)$")
    ax[1].plot(pOrbit, pOrbit*qOrbit, 'k', label=r"$(p,qp)$")
    ax[0].legend(ncol=1,fontsize=10,loc="upper right")
    ax[1].legend(ncol=1,fontsize=10,loc="upper right")
    ax[0].text(10,0.2,r"$\omega$=%g Hz"%(p["stimFreq"]))
    ax[0].text(10,0.4,stri)
    ax[0].set_ylim(0,1); ax[1].set_ylim(0,1)
    ax[1].set_xlim(0,1)
    return orbits


# -----------------------------------------------------------------
# Facilitation to depression based on the steepness of the threshold for the probability of release
# -----------------------------------------------------------------
p={"gainP":3.0,"vHalfP":-10.0, "biasP":1.0, "rateP":1/20.0, "expP":1.0,
    "ssQ": 0.9, "rateQ": 1/1.0, "expQ":1.0,
    "stimFreq": 20.0, "timeStep":0.01}
p["yHalfP"]= p["vHalfP"]/26.7
orbits=list(); ax=list()
gs= sc.arange(1,5.01,1)
nSims= len(gs)
fig=gr.figure(figsize=(17,13)); gr.ioff()
rows=nSims; cols=2
for n in sc.arange(nSims*cols):
    ax.append(fig.add_subplot(rows,cols,n+1))
print(len(ax))
for n in sc.arange(nSims):
    p["gainP"]= gs[n]
    str0=r"$\tau_P=$%gms, $\tau_Q=$%gms, $v_P=$%gmV, $g_P=$%g, $b_P=$%g"%(1/p["rateP"],1/p["rateQ"],p["vHalfP"],p["gainP"],p["biasP"])
    a=2*n; ll= ax[a:a+2]
    oo=dynamicSTSP(p,ax[a:a+2], ic=[0.2,0.9],timeMax=600.0, vT=26.7, pulses=100+sc.arange(0,400,1000/p["stimFreq"]),stri=str0)
    orbits.append(oo)
gr.ion(); gr.draw()
fig.subplots_adjust(bottom=0.05, left=0.05, top=0.98, right=0.98, hspace=0.1, wspace=0.1)
fig.savefig("stspFigures/stsp_%gHz_tauP%gms_tauQ%gms_vP%gmV_gP%g.png"%(p["stimFreq"], 1/p["rateP"], 1/p["rateQ"],p["vHalfP"],p["gainP"]))


# -----------------------------------------------------------------
# Facilitation to depression based on the rate of recovery of the readily releasable pool
# Soft spots for switching
# -----------------------------------------------------------------
p={"gainP":3.0,"vHalfP":-10.0, "biasP":1.0, "rateP":1/20.0, "expP":0.0,
    "ssQ": 0.99, "rateQ": 1/1.0, "expQ":1.0,
    "stimFreq": 20.0, "timeStep":0.01}
p["yHalfP"]= p["vHalfP"]/26.7
orbits=list(); ax=list()
rQs= 1/sc.arange(0.25,3.01,0.5)
nSims= len(rQs)
fig=gr.figure(figsize=(17,13)); gr.ioff()
rows=nSims; cols=2;
print(len(ax))
p["vHalfP"]= -20.0
p["gainP"]= 2.0
p["rateP"]= 1/40.0
for n in sc.arange(nSims):
    ax.append(fig.add_subplot(rows,cols,2*n+1))
    ax.append(fig.add_subplot(rows,cols,2*n+2))
    p["rateQ"]= rQs[n]
    str0=r"$\tau_P=$%gms, $\tau_Q=$%gms, $v_P=$%gmV, $g_P=$%g, $b_P=$%g"%(1/p["rateP"],1/p["rateQ"],p["vHalfP"],p["gainP"],p["biasP"])
    a=2*n; ll= ax[a:a+2]
    oo=dynamicSTSP(p,ax[a:a+2], ic=[0.1,1.0],timeMax=800.0, vT=26.7, pulses=100+sc.arange(0,400,1000/p["stimFreq"]),stri=str0)
    orbits.append(oo)
gr.ion(); gr.draw()
fig.subplots_adjust(bottom=0.05, left=0.05, top=0.98, right=0.98, hspace=0.1, wspace=0.1)
fig.savefig("stspFigures/stsp_%gHz_tauP%gms_tauQ%gms_vP%gmV_gP%g.png"%(p["stimFreq"], 1/p["rateP"], 1/p["rateQ"],p["vHalfP"],p["gainP"]))

# -----------------------------------------------------------------
# Facilitation to depression based on the rate of decay of the probability of release
# -----------------------------------------------------------------
p={"gainP":3.0,"vHalfP":-10.0, "biasP":1.0, "rateP":1/20.0, "expP":0.0,
    "ssQ": 1.0, "rateQ": 1/1.0, "expQ":1.0,
    "stimFreq": 20.0, "timeStep":0.005}
p["yHalfP"]= p["vHalfP"]/26.7
orbits=list(); ax=list()
rPs= 1/sc.arange(10,40.01,5.0)
nSims= len(rPs)
fig=gr.figure(figsize=(17,13)); gr.ioff()
rows=nSims; cols=2
for n in sc.arange(nSims*cols):
    ax.append(fig.add_subplot(rows,cols,n+1))
print(len(ax))
p["vHalfP"]= 0.0
p["gainP"]= 2.0
p["rateQ"]= 1/0.5
for n in sc.arange(nSims):
    p["rateP"]= rPs[n]
    str0=r"$\tau_P=$%gms, $\tau_Q=$%gms, $v_P=$%gmV, $g_P=$%g, $b_P=$%g"%(1/p["rateP"],1/p["rateQ"],p["vHalfP"],p["gainP"],p["biasP"])
    a=2*n; ll= ax[a:a+2]
    oo=dynamicSTSP(p,ax[a:a+2],ic=[0.2,0.9], timeMax=600.0, vT=26.7, pulses=100+sc.arange(0,400,1000/p["stimFreq"]),stri=str0)
    orbits.append(oo)
gr.ion(); gr.draw()
fig.subplots_adjust(bottom=0.1, left=0.05, top=0.98, right=0.98, hspace=0.1, wspace=0.1)
fig.savefig("stspFigures/stsp_%gHz_tauP%gms_tauQ%gms_vP%gmV_gP%g.png"%(p["stimFreq"], 1/p["rateP"], 1/p["rateQ"],p["vHalfP"],p["gainP"]))

# -----------------------------------------------------------------
# Facilitation to depression based on the rate of recovery of the readily releasable pool
# -----------------------------------------------------------------
p={"gainP":3.0,"vHalfP":-10.0, "biasP":1.0, "rateP":1/20.0, "expP":0.0,
    "ssQ": 1.0, "rateQ": 1/1.0, "expQ":1.0,
    "stimFreq": 20.0, "timeStep":0.01}
p["yHalfP"]= p["vHalfP"]/26.7
orbits=list(); ax=list()
vHPs= sc.arange(-30,10.1,10.0)
nSims= len(vHPs)
fig=gr.figure(figsize=(17,13)); gr.ioff()
rows=nSims; cols=2
for n in sc.arange(nSims*cols):
    ax.append(fig.add_subplot(rows,cols,n+1))
print(len(ax))

p["rateP"]= 1/20.0
p["rateQ"]= 1/2.0
for n in sc.arange(nSims):
    p["vHalfP"]= vHPs[n]
    str0=r"$\tau_P=$%gms, $\tau_Q=$%gms, $v_P=$%gmV, $g_P=$%g, $b_P=$%g"%(1/p["rateP"],1/p["rateQ"],p["vHalfP"],p["gainP"],p["biasP"])
    a=2*n; ll= ax[a:a+2]
    oo=dynamicSTSP(p,ax[a:a+2],ic=[0.2,0.9], timeMax=600.0, vT=26.7, pulses=100+sc.arange(0,400,1000/p["stimFreq"]),stri=str0)
    orbits.append(oo)
gr.ion(); gr.draw()
fig.subplots_adjust(bottom=0.1, left=0.05, top=0.98, right=0.98, hspace=0.1, wspace=0.1)
fig.savefig("stspFigures/stsp_%gHz_tauP%gms_tauQ%gms_vP%gmV_gP%g.png"%(p["stimFreq"], 1/p["rateP"], 1/p["rateQ"],p["vHalfP"],p["gainP"]))


# -----------------------------------------------------------------
p={"gainP":3.0,"vHalfP":-10.0, "biasP":1.0, "rateP":1/20.0, "expP":0.0,
    "ssQ": 1.0, "rateQ": 1/1.0, "expQ":1.0,
    "stimFreq": 20.0, "timeStep":0.005}
p["yHalfP"]= p["vHalfP"]/26.7
orbits=list(); ax=list()
b= sc.arange(0.5,1.01,0.1)
nSims= len(b)
fig=gr.figure(figsize=(17,13)); gr.ioff()
rows=nSims; cols=2
for n in sc.arange(nSims*cols):
    ax.append(fig.add_subplot(rows,cols,n+1))
print(len(ax))
p["vHalfP"]= -20.0
p["gainP"]= 2.0
p["rateQ"]= 1/4.0
p["rateP"]= 1/60.0
for n in sc.arange(nSims):
    p["biasP"]= b[n]
    str0=r"$\tau_P=$%gms, $\tau_Q=$%gms, $v_P=$%gmV, $g_P=$%g, $b_P=$%g"%(1/p["rateP"],1/p["rateQ"],p["vHalfP"],p["gainP"],p["biasP"])
    a=2*n; ll= ax[a:a+2]
    oo=dynamicSTSP(p,ax[a:a+2],ic=[0.2,0.9], timeMax=600.0, vT=26.7, pulses=100+sc.arange(0,400,1000/p["stimFreq"]),stri=str0)
    orbits.append(oo)
gr.ion(); gr.draw()
fig.subplots_adjust(bottom=0.1, left=0.05, top=0.98, right=0.98, hspace=0.1, wspace=0.1)
fig.savefig("stspFigures/stsp_%gHz_tauP%gms_tauQ%gms_vP%gmV_gP%g.png"%(p["stimFreq"], 1/p["rateP"], 1/p["rateQ"],p["vHalfP"],p["gainP"]))




