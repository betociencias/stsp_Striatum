from stdste1 import *
plt.ion()

# Stimulus
pars={'timeStart':0.0, 'timeMax': 500.0, 'timeStep':1e-3,
      'tauX':75.0,  'tauP':75.0, 
      'asynX':1.0, 'asynP':0.2, 
      'kP':0.0, 'kX':1.0, 
      'jumpP':0.2, 
      'stimHz':20.0, 'stimStart':10.0}

strHz=r'$\omega_{stim}$=%.0f Hz'%(pars['stimHz'])
pars= stimTrain(pars)
p0=pars

# Initial simulation
orbit=simulateFD(pa=p0)
xOrbit=orbit[0]
pOrbit=orbit[1]
xpOrbit=orbit[2]
#
r=1; c=1
f0, ax1=plt.subplots()
#plt.ioff()
plt.subplots_adjust(bottom=0.35)
pLine,=plt.plot(p0['sampleTimes'],pOrbit)
xLine,=plt.plot(p0['sampleTimes'],xOrbit)
pxLine,=plt.plot(p0['sampleTimes'],xpOrbit)
#
axcolor = 'lightgoldenrodyellow'
#
axTauX = plt.axes([0.1, 0.05, 0.65, 0.03], axisbg=axcolor)
axTauP = plt.axes([0.1, 0.1, 0.65, 0.03], axisbg=axcolor)
axAsynX = plt.axes([0.1, 0.15, 0.65, 0.03], axisbg=axcolor)
axAsynP = plt.axes([0.1, 0.2, 0.65, 0.03], axisbg=axcolor)
axJumpP = plt.axes([0.1, 0.25, 0.65, 0.03], axisbg=axcolor)
# 
tauXSlider = Slider(axTauX, r'$\tau_x$', 1.0, 500.0, valinit=p0['tauX'])
tauPSlider = Slider(axTauP, r'$\tau_p$', 1.0, 500.0, valinit=p0['tauP'])
asynXSlider = Slider(axAsynX, r'$x_{\infty}$', 0.0, 1.0, valinit=p0['asynX'])
asynPSlider = Slider(axAsynP, r'$p_{\infty}$', 0.0, 1.0, valinit=p0['asynP'])
jumpPSlider = Slider(axJumpP, r'$h_p$', 0.0, 1.0, valinit=p0['jumpP'])

def updateSliders(val):
    pars['tauX'] = tauXSlider.val
    pars['tauP'] = tauPSlider.val
    pars['asynX'] = asynXSlider.val
    pars['asynP'] = asynPSlider.val
    pars['jumpP'] = asynPSlider.val
    #pars['stimF'] = stimFreqSlider.val
    orbit=simulateFD(pa=pars)
    xLine.set_ydata(orbit[0])
    pLine.set_ydata(orbit[1])
    xpLine.set_ydata(orbit[2])
    f0.canvas.draw_idle()

tauXSlider.on_changed(updateSliders)
tauPSlider.on_changed(updateSliders)
asynXSlider.on_changed(updateSliders)
asynPSlider.on_changed(updateSliders)
jumpPSlider.on_changed(updateSliders)
#stimFSlider.on_changed(update)


#plt.ion()
#plt.draw()
plt.show()
