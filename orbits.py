#RV orbit modeling code based on RVLIN by Jason Wright (IDL) and orbits.f by Alex Wolszczan (FORTRAN77)

#from __future__ import print_function #interface Python 2/3, need this?
import numpy as np
import read_rdb_harpsn as rr
from pwkit import lsqmdl
from utils import * #just fan for now

#Given a target name, does everything!
def orbits_test(targname='HD209458',norbits=1,nterms=0,jitter=0.0,modelstart=0,modelrange=0,modelstep=0.1,nboot=1000):
    
    #read RV data 
    #jdb, rv, srv, labels = rr.process_all(targname)
 
    #demo!
    
    sfile = open('/home/sgettel/Dropbox/cfasgettel/py.lib/sgcode/rvorbits/209458.vel.txt')
    jdb, rv, srv = np.loadtxt(sfile,unpack=True,usecols=(0,1,2))
   

    npts = jdb.size
    tstart = jdb[0]

    

    #read Mstar
    Mstar = 1.0
    #Mstar = 1.69

    #process offsets - very not implemented 

    #read/process parameter guesses from somewhere
    #format [(T0, per, K1, ecc, omega, omegadot)*nplanets,voffset]
    #guess = np.array([53460.9375,2093.060,18.242,0.40179,206.896,0.000,-1.704])
    #guess = np.array([54292.4256,746.139,91.406,0.37450,108.235,0.0,83.342])
        
    guesspars = np.array([3.524733, 2452836.1, 0.0, 336.5415, 85.49157, -1.49, 0.0])#HD209458
    #p = orbel[0+i*7]
    #tp = orbel[1+i*7]
    #ecc = orbel[2+i*7]
    #om = orbel[3+i*7] *np.pi/180. #degrees to radians
    #k = orbel[4+i*7]
    #gamma = orbel[5+i*7]
    #dvdt = orbel[6+i*7]
    #curv = 0 


    m = rvfit_lsqmdl(guesspars, jdb, rv, srv, jitter=jitter)
    
    
    #calculate mpsini

    #display some output - customize
    m.print_soln()
    
    
    #make plots
    tmod = np.linspace(np.min(jdb),np.max(jdb),1000)
    model_init = rv_drive(guesspars,tmod)
    model_final = rv_drive(m.params,tmod)
     
    #scramble residuals & call bootstrapping
    bootpar, sigpar = bootstrap_rvs(m.params, jdb, rv, srv,nboot=nboot,jitter=jitter)

    #call mass estimate

    return m, bootpar #jdb, rv, nsrv

#def mass_estimate


#this should set limits and call lsqmdl, should be callable by bootstrapper...
def rvfit_lsqmdl(orbel,jdb,rv,srv,jitter=0, param_names=0):

    if jitter > 0.0: #this should happen in rvfit_lsqmdl
        nsrv = np.sqrt(srv**2 + jitter**2)
    else:
        nsrv = srv

    if param_names == 0: #do something more clever here later
        param_names = ['Per', 'Tp', 'ecc', 'om', 'K1', 'gamma', 'dvdt'] 

    m = lsqmdl.Model(None, rv, 1./nsrv) #m is a class...
    m.set_func(rv_drive,param_names, args=[jdb])

    #make some reasonable limits
    m.lm_prob.p_limit(0, lower=0.1, upper=(np.max(jdb)-np.min(jdb))) #per no longer than range of survey
    m.lm_prob.p_limit(1, lower=orbel[1]-orbel[0]/2., upper=orbel[1]+orbel[0]/2.) #T0 within one period guess
    m.lm_prob.p_limit(2, lower=0.0, upper=0.99) #ecc must be physical
    m.lm_prob.p_limit(3, lower=0.0, upper=360.0)
    m.lm_prob.p_limit(4, lower=0.0, upper=1.0e5) #K must be physical  
    m.lm_prob.p_value(6, 0.0, fixed=True) #fix dvdt

    m.solve(orbel)
    #m.print_soln()
  
    return m

def rv_drive(orbel, t):
    #From rv_drive.pro in RVLIN by JTW

    rv = np.zeros_like(t)
    nplanets = orbel.size/7 

    phase = np.zeros((rv.size,nplanets))
    #print orbel
    for i in range(nplanets):
        p = orbel[0+i*7]
        tp = orbel[1+i*7]
        ecc = orbel[2+i*7]
        om = orbel[3+i*7] *np.pi/180. #degrees to radians
        k = orbel[4+i*7]
        gamma = orbel[5+i*7]
        dvdt = orbel[6+i*7]
        curv = 0
        if i == 0 and nplanets*7 == orbel.size-1: # if 1 too many elements,
            curv = orbel[-1]                      # last is curvature

        #Error checking
        if p < 0 or ecc < 0 or ecc >= 1 or k < 0:
            print 'Bad inputs to rv_drive'
            
            if p < 0:
                p = 1e-2
            if ecc < 0:
                ecc = 0
            if ecc >= 1:
                ecc = 0.99
            if k < 0:
                k = 1e-2
        
        #calculate the approximate eccentric anomaly, E1, from the mean anomaly, M
        M = 2.0*np.pi*( ((t-tp)/p) - np.floor((t-tp)/p) ) #phase in radians
        
        E1 = kepler(M,ecc) #returns a matrix, why?
        E1 = E1[0,:]
        #calculate true anomaly
        n1 = 1.0 + ecc
        n2 = 1.0 - ecc
        theta = 2.0*np.arctan(np.sqrt(n1/n2)*np.tan(E1/2.0))

        phase0 = theta + om - np.pi/2.0
        under = np.squeeze(np.where(phase0 < -np.pi)) 
        phase0[under] += np.pi
        over = np.squeeze(np.where(phase0 >= np.pi))
        phase0[over] -= np.pi
        
        phase[:,i] = phase0

        #calculate radial velocity
        epoch = t[0]
        rv = rv + k*(np.cos(theta + om) + ecc*np.cos(om)) + gamma + dvdt*(t - epoch) + curv*(t - epoch)**2 #default epoch = t[0]

    return rv


#Iteratively solve for E (anomoly) given M (mean anomoly) and e (eccentricity)
def kepler(inM,inecc):
    #Port of kepler.pro by JTW, from Murray & Dermott p. 35, from Danby (1988)
    #SJG Aug 2014

    #make floats into arrays
    if len(np.array(inM).shape) == 0:
        marr = np.array([inM])
    else:
        marr = inM
    if len(np.array(inecc).shape) == 0:
        eccarr = np.array([inecc])
    else:
        eccarr = inecc
    
    nm = marr.size
    nec = eccarr.size
    if nec == 1 and nm > 1:   #Do I need this?
        eccarr = fan(eccarr,nm)
    if nec > 1 and nm == 1:
        marr = fan(marr,nec)

    conv = 1e-12 #convergence criteria
    k = 0.85     #scale factor?
    
    mphase = marr/(2.0*np.pi)
    #print mphase.size, mphase.shape

    #apply the relevant part of restrict.pro? Missing some error checking
    under = np.squeeze(np.where(mphase < 0))
    mphase[under] +=1
    over = np.squeeze(np.where(mphase >= 1))
    mphase[over] -= 1

    ssm = np.sign(np.sin(marr))

    #make first guess at E
    Earr = marr + ssm*k*eccarr
    fiarr = (Earr - eccarr*np.sin(Earr) - marr) #E - ecc*sin(E) - M, converge to 0

    convd = np.squeeze(np.where(np.abs(fiarr) > conv)) #which indices converged?
    nd = convd.size
    #print nd
    count = 0
    while nd > 0:          #while unconverged elements exist...
        count += 1
        M = marr[convd]    #just keep the unconverged elements
        ecc = eccarr[convd]
        E = Earr[convd]
        
        fi = fiarr[convd]        #fi = E - e*sin(E)-M,  should go to 0
        fip = 1 - ecc*np.cos(E)  #d/dE(fi), fi^(\prime)
        fipp = ecc*np.sin(E)     #d/dE(d/dE(fi)), fi^(\prime\prime)
        fippp = 1 - fip          #d/dE(d/dE(d/dE(fi))),  fi^(\prime\prime\prime)

        d1 = -fi/fip                          #first order correction to E
        d2 = -fi/(fip + d1*fipp/2.0)          #second order correction to E
        d3 = -fi/(fip + d2*fipp/2.0 + d2*d2*fippp/6.0) #third order correction
        E += d3
        Earr[convd] = E
        fiarr = (Earr - eccarr*np.sin(Earr) - marr)
        convd = np.squeeze(np.where(np.abs(fiarr) > conv))
        nd = convd.size
        if count > 100:
            print 'WARNING!  Keplers equation not solved!!!'
            nd = 0

    return Earr

def bootstrap_rvs(bestpar, jdb, rv, srv,nboot=1000,jitter=0):
    #based on the general method of bootpar.pro by SXW
    #Wow, this is slow! Make faster

    bestfit = rv_drive(bestpar, jdb)
    print rv.size, bestfit.size, jdb.size, bestpar.size
    resid = rv - bestfit

    nboot = np.max([nboot,resid.size*np.log10(resid.size)**2]).astype(int) #why this limit?
    print 'nboot = ',nboot

    bootpar = np.zeros((nboot,bestpar.size))
    
    for i in range(nboot):
        if i%100 == 0:
            print i/nboot*100,'%'
         
        scramble = np.random.randint(jdb.size,size=jdb.size)#random indices with replacement
        tmprv = resid[scramble] + bestfit            #scramble residuals
        tmperr = srv[scramble]                      #error goes with individual residual

        mi = rvfit_lsqmdl(bestpar, jdb, tmprv, tmperr, jitter=jitter)
        bootpar[i,:] = mi.params

    
    sigpar = np.stddev(bootpar,axis=0)

    return bootpar, sigpar
