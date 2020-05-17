from stdste1 import *
import channelsTransporters as tp

#
CmMicroFaradsPerCm2=0.75
#CmMicroFaradsPerCm2=1.0
CmNanoFaradsPerCm2=CmMicroFaradsPerCm2*1e3
CmNanoFaradsPerMu2=CmNanoFaradsPerCm2/(10**8)
CmPicoFaradsPerCm2=CmMicroFaradsPerCm2*1e6
CmPicoFaradsPerMu2=CmPicoFaradsPerCm2/(10**8)
# Compartment dimensions
r=5.0; L= 10.0; 
areaCylinder=2*r*sc.pi*L; 
areaSphere=4*sc.pi*(r**2);
CmCylinder = areaCylinder* CmPicoFaradsPerMu2
CmSphere = areaSphere* CmPicoFaradsPerMu2
print("radius=%f um, length=%f um"%(r,L))
print("Cm Cylinder=%f pF, Cm Sphere=%f pF"%(CmCylinder, CmSphere))

# Create a random stimulus from a renewal process. Done.
# Translate the random stimuli into a continuous trace that mimmics the calcium transients
# Use the "calcium transients" to generate a trace of neurotransmitter release 
# Examine the voltage transients resulting from the neurotransmitter-dependent activations of fast excitatory receptors (AMPA/kainate)

pars={'timeMin':0.0, 'timeMax': 1000.0, 'timeStep':1e-3,
      'tauX':75.0,  'tauP':25.0, 'tauE':0.775, 'tauI':1.4, 'tauAlpha':1.0,
      'asynX':1.0, 'asynP':0.1, 
      'kP':1.0, 'kX':1.0, 'kX':0.0, 
      'jumpP':0.1, 
      'stimHz':20.0, 'stimStart':10.0,
      'vE':0.0,'vI':-70.0,'vNa':70.0, 'vK':-90.0,'vATP':-450.0,
      'preFR': 10.0,
      'Cm': 2.356,
      'aBarNaK': 50.0, 'aBarE': 1.0, 'aBarI': 0.0, 
      'tempC': 37.0,
    }

pars['vT'] = tp.koverq * (273.15 + pars['tempC'])
pars['v2T'] = 2*pars['vT']
# 
pars['vNaK']= pars['vATP'] + 3*pars['vNa'] - 2*pars['vK'] 
sampleTimes=sc.arange(pars['timeMin'],pars['timeMax'],pars['timeStep'])
spikesPre=RenewalProcess(tmin=pars['timeMin'], tmax=pars['timeMax'], fr=pars['preFR'], refract=2.0)
nSpikesPre=len(spikesPre)
trainX=trainAlpha(samplingTimes=sampleTimes, pulseTimes=spikesPre, tauTrain=pars['tauAlpha'],downAccel=1.0)
trainV=trainAlpha(samplingTimes=sampleTimes, pulseTimes=spikesPre+3*pars['tauAlpha'], tauTrain=pars['tauAlpha'],downAccel=1.0)
trainP= sc.zeros(len(sampleTimes))
indSecondSpike= gr.find(sampleTimes<spikesPre[1]).max()
trainP[indSecondSpike:]= trainX[indSecondSpike:]
pars['spikesV']=lambda t: sc.interp(t, sampleTimes, trainV)
pars['spikesCa']=lambda t: sc.interp(t, sampleTimes, trainX)
pars['spikesP']=lambda t: sc.interp(t, sampleTimes, trainP)

tCrit=spikesPre
U0=[pars['asynX'], pars['asynP'], 0, -65.0 ]
orbit = sc.integrate.odeint(stFDPrePostE,U0,sampleTimes,args=(pars,), tcrit=tCrit).transpose()
orbit[3]=orbit[3]-pars['vNaK']
labels=[r'$x$', r'$p$', r'$q_E$', r'$v$']
colors=['b','g','k:','k']
fig=gr.figure(figsize=(11,3))
r=1; c=1;
axSyn=fig.add_subplot(r,c,1)
gr.ioff()
axSyn.plot(spikesPre, sc.ones(nSpikesPre), 'k|', ms=20, label=r'$v_{pre}$')
axSyn.plot(sampleTimes,orbit[0]*orbit[1], 'r', label=r'$xp$')
for n in range(len(orbit)): 
        axSyn.plot(sampleTimes,orbit[n], colors[n], label=labels[n])

axSyn.set_ylim(-0.1,1.1)
axSyn.legend()
gr.ion(); gr.draw()

