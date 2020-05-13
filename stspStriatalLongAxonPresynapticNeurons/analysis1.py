import sympy as sy
import scipy as sc

t,p,x= sy.symbols('t p x')
pinf, taup, h = sy.symbols('pinf taup h')
xinf, taux= sy.symbols('xinf taux')

def pSolution(t, t0, p0):
    return pinf/ (1 - (1 - pinf/p0)* sy.exp( pinf*(t0-t)/taup))

def xSolution(t, t0, x0):
    return xinf/ (1 - (1 - xinf/x0)* sy.exp( xinf*(t0-t)/taux))

def pPulse(p,h):
    return p *(1-h) + h

def xPulse(p,x):
    return x * (1 - p*x)


def createPulseSymbols(nPulses=3):
    pulsos=list()
    for n in range(0,nPulses):
        str1= 't%d'%n
        print(str1)
        pulsos.append(sy.symbols(str1))
    return pulsos

def pRecursivePulses(pulseTimes, transientSolution, pulseJump, pStart=pinf):
    pTransient=pStart
    print(pTransient)
    for n in range(0,len(pulseTimes)):
        p0 = pulseJump(pTransient,h).subs('t',pulseTimes[n])
        pTransient=transientSolution(t0=pulseTimes[n],t=t, p0=p0)
        print(p0)
        print(pTransient)
        print('\n')
    return pTransient

pulsos=createPulseSymbols(nPulses=4)
if 1: 
    pSol1= pRecursivePulses(pulseTimes=pulsos, transientSolution=pSolution, pulseJump=pPulse, pStart=pinf)


def recursivePulses2(pulseTimes, pTransientSolution, xTransientSolution, pPulseJump, xPulseJump, pStart=pinf, xStart=xinf):
    pTransient=pStart; xTransient=xStart
    for n in range(0,len(pulseTimes)):
        print n
        pAtPulse = pPulseJump(pTransient,h).subs('t',pulseTimes[n])
        xAtPulse = xPulseJump(pTransient,xTransient).subs('t',pulseTimes[n])
        pTransient=pTransientSolution(t0=pulseTimes[n],t=t, p0=pAtPulse)
        xTransient=xTransientSolution(t0=pulseTimes[n],t=t, x0=xAtPulse)
    return pTransient,xTransient


pars={'pinf':0.3, 'xinf':0.8, 'taup':10.0, 'taux':50.0, 'h':0.4}
pulsos=createPulseSymbols(nPulses=2)
p,x =recursivePulses2(pulseTimes=pulsos, pTransientSolution=pSolution, xTransientSolution=xSolution, pPulseJump=pPulse, xPulseJump=xPulse, pStart=pinf, xStart=xinf)

xx= x.subs(pars)
pp= p.subs(pars)
pulseTimes=sc.arange(10.0, 200.0, 50.0)
for m in range(len(pulsos)):
    str1= 't%d'%m
    pp=pp.subs(str1,pulseTimes[m])
    xx=xx.subs(str1,pulseTimes[m])



 
