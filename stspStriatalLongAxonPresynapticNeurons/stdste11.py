from stdste1 import *

orbits, figs= tauRecSweep(tauRecs=75.0*sc.arange(0.2,3.001,0.4),
tauFaci=75.0, p0=0.2, q=0.2, stimHz=20.0, firstSpike=10.0,
timeStart=0.0, timeMax = 1000.0, timeStep=1e-3, U0=[1,0.2],
embedLaTeX=1, saveFig=1)


