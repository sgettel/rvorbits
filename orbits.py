#RV orbit modeling code based on RVLIN by Jason Wright (IDL) and orbits.f by Alex Wolszczan (FORTRAN77)

#from __future__ import print_function #interface Python 2/3, need this?
import emcee
import socket
import triangle
import numpy as np
import read_rdb_harpsn as rr
import matplotlib.pyplot as plt
import utils as ut
from pwkit import lsqmdl

#
# TO DO:
# - print output to file
# - histogram plots
# - fix TT for eccentric orbits
# - multiple planets
# - read orbit params from somewhere
# - add offset terms
# - MCMC


#Given a target name, do everything!
def orbits_test(targname='K00273',jitter=0.0,nboot=1000,epoch=2.455e6,circ=0,maxrv=1e6,minrv=-1e6,maxsrv=5, webdat='no', nwalkers=200, pfix=1):

    transit = np.zeros(1)
    
    circ = ut.arrayify(circ)
    pfix = ut.arrayify(pfix)
    print pfix

    if socket.gethostname() == 'sara-gettels-macbook-2.local':
    	home = '/Users/Sara/'
    else:
    	home = '/home/sgettel/'

    #read RV data
    if targname == 'HD209458':
        #demo!
        sfile = open(home+'Dropbox/cfasgettel/py.lib/sgcode/rvorbits/209458.vel.txt')
        jdb, rv, srv = np.loadtxt(sfile,unpack=True,usecols=(0,1,2))

    else:
        print targname
        jdb, rv, srv, labels = rr.process_all(targname,maxsrv=maxsrv,maxrv=maxrv,minrv=minrv)
        jdb += 2.4e6 #because process_all gives truncated JDBs
     
    #adjust values to be sensible
    #if jdb[0] > epoch:
    #    print 'truncating dates'
    tnorm = jdb - epoch #truncate
    rvnorm = rv - np.median(rv)
    print np.min(rvnorm),np.max(rvnorm),np.min(rv),np.max(rv)
    
  
    #read Mstar
    
    #Mstar = 1.69

    #process offsets - very not implemented 

    #read/process parameter guesses from somewhere
    #p = orbel[0+i*7]
    #tp = orbel[1+i*7]
    #ecc = orbel[2+i*7]
    #om = orbel[3+i*7] *np.pi/180. #degrees to radians
    #k = orbel[4+i*7]
    #gamma = orbel[5+i*7]
    #dvdt = orbel[6+i*7]
    #curv = orbel[-1], optional
    

    if targname == 'HD209458':
        guesspars = np.array([3.524733, 2452826.628514, 0.0, 336.5415, 85.49157+30, -1.49-50, 0.01])#HD209458
        transit = np.array([2452826.628514]) 
        mstar = 1.0
    
    if targname == 'K00273':
        guesspars = np.array([10.57377, 2455008.066, 0.0, 90.0, 1.7358979, -3398.0498, 1.1889011, 0.0]) #K00273
        transit = np.array([2455008.066,0.0]) 
        mstar = 1.07
        #2 planets...
        #guesspars = np.array([10.57377, 2455008.066, 0.0, 90.0, 1.7358979, -3398.0498, 0.0, 1400.0, 2455008.066, 0.0, 90.0, 100.0, 0.0, 0.0])
    
    if targname == 'K00069':
        guesspars = np.array([4.72673978, 2454944.29227, 0.0, 90.0, 1.733, -91.08, 0.0329])
        transit = np.array([2454944.29227])
        mstar = 0.887

    norbits = guesspars.size/7
    ip = np.arange(norbits)
    guesspars[1+ip*7] -= epoch

    for i in range(transit.size):  #do something more clever here...
        if not transit[i] == 0.0:
            transit[i] -= epoch 
            print 'set transit time: ', transit
 
    if not guesspars[6] == 0: #only set for 1st planet
        trend = 1
    else:
        trend = 0

    npars = guesspars.size
    if npars % 7 == 1: #1 too many params
        curv = 1
        print 'parabolic term allowed'
    else:
        curv = 0
    
    
    m = rvfit_lsqmdl(guesspars, tnorm, rvnorm, srv, jitter=jitter,trend=trend,circ=circ, curv=curv,tt=transit,epoch=epoch,pfix=pfix)
    m0 = np.copy(m)

    
    #display initial fit - want to show fixed params too
    m.print_soln()  #read this in detail
    
    #mass estimate
    mpsini = mass_estimate(m, mstar, norbits=norbits)
    print 'mp*sin(i):         ',str(mpsini)
    
    #make plots
    plot_rv(targname,tnorm,rvnorm,srv,guesspars,m,nmod=200,home=home)
    
    #call bootstrapping
    if nboot > 0:
        bootpar, meanpar, sigpar = bootstrap_rvs(m.params, tnorm, rvnorm, srv,nboot=nboot,jitter=jitter,circ=circ,trend=trend,curv=curv,tt=transit,pfix=pfix)
   
        mpsini, mparr_all = mass_estimate(m,mstar,norbits=norbits,bootpar=bootpar)
       
        print_errs(meanpar, sigpar, mpsini, mparr_all, norbits=norbits,curv=curv)

        #return m0, bootpar,sigpar, mparr_all #jdb, rv, nsrv

    #call MCMC    
    if nwalkers > 0:

        m, flt, samples = setup_emcee(targname, m, tnorm, rvnorm, srv, circ=circ, trend=trend, curv=curv, tt=transit, jitter=jitter, nwalkers=nwalkers)
    
        return m, flt, samples# bestpars, varpars, flt, pnames # 
    else:
        return

def plot_rv(targname,jdb,rv,srv,guesspars,m,nmod=1000,norbits=1,home='/home/sgettel/'):
   
    tmod = np.linspace(np.min(jdb),np.max(jdb),nmod)

    #guesspars2 = np.copy(guesspars)
    model_init = rv_drive(guesspars,tmod) # something wrong if trend
    model_final = rv_drive(m.params,tmod)
    
    #unphased data
    plt.figure(1)
    plt.errorbar(jdb,rv,yerr=srv,fmt='bo')
    plt.plot(tmod,model_final,'r-')
    #plt.plot(tmod,model_init,'g-')
    plt.savefig(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+'_autoplot.png')
    plt.close(1)

    #phase at first period - only taking 1 period right now - LEFT OFF HERE
    
    pars = np.copy(m.params)  #for model
    pars[6] = 0 #remove trend
    if pars.size % 7 == 1: #remove curve if it is used
        pars[-1] = 0
    pars3 = np.copy(m.params) #for obs
    pars3[4] = 0.0 #trend only
    pars3[5] = 0.0
    phase = (jdb - pars[1])/pars[0] % 1.0
    
    rvt = rv_drive(pars3,jdb)

    
    #print rvt[0:10]

    plt.figure(2)
    plt.errorbar(phase, rv-rvt, yerr=srv,fmt='bo')
    plt.plot((tmod - pars[1])/pars[0] % 1.0, rv_drive(pars, tmod),'r.')
    #plt.plot((tmod - guess))
    plt.savefig(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+'_phase_autoplot.png')
    plt.close(2)
    
#def plot_pars(targname,bootpars,mparr_all,norbits=1):
#    print targname
#    plot histograms!

def print_errs(meanpar,sigpar, mpsini, mparr_all,norbits=1,curv=0):
    
    for i in range(norbits):
        print '*****Planet ',str(i+1),' errors:*****'
        print 'Per: ', str(meanpar[i*7]),'+/-',str(sigpar[i*7])
        print 'Tp: ', str(meanpar[i*7+1]),'+/-',str(sigpar[i*7+1])
        print 'ec: ', str(meanpar[i*7+2]),'+/-',str(sigpar[i*7+2])
        print 'om: ', str(meanpar[i*7+3]),'+/-',str(sigpar[i*7+3])
        print 'K: ', str(meanpar[i*7+4]),'+/-',str(sigpar[i*7+4])
        print 'gamma: ', str(meanpar[i*7+5]),'+/-',str(sigpar[i*7+5])
        print 'dvdt: ', str(meanpar[i*7+6]),'+/-',str(sigpar[i*7+6])
        if curv == 1:
            print 'curv: ',str(meanpar[i*7+7]),'+/-',str(sigpar[i*7+7]) #BAD, only true for 1 planet
        print 'mp*sin(i): ',str(np.mean(mparr_all[i,:])),'+/-',str(np.std(mparr_all[i,:]))
        print 'mass error:', str(np.std(mparr_all[i,:])/mpsini*100),'%'
    return

def mass_estimate(m,mstar,norbits=1,bootpar=-1):
   
    #some constants, maybe use astropy here
    msun = 1.9891e30
    mearth = 5.97219e24
    G = 6.673e-11 #m^3 kg^-1 s^-2
    etoj = 317.83
    
    ip = np.array(range(norbits))
    
    pers = m.params[ip*7]
    eccs = m.params[ip*7+2]
    amps = m.params[ip*7+4]
    
    #mass estimate
    fm = (1 - eccs*eccs)**(1.5)*amps**3*(pers*86400.0)/(2*np.pi*G) #kg
    mpsini = ((mstar*msun)**2*fm)**(1./3.)/mearth

    #calculate error on mass if bootpar array is input
    if len(np.array(bootpar).shape) > 0:

        #print bootpar.shape, m.params.shape
        mparr_all = np.zeros((norbits,bootpar.shape[0]))

        for i in ip:
            fmarr = (1 - bootpar[:,i*7+2]**2)**(1.5)*bootpar[:,i*7+4]**3*(bootpar[:,i*7]*86400.0)/(2.0*np.pi*G)
            mparr = ((mstar*msun)**2*fmarr)**(1./3.)/mearth
            mparr_all[i,:] = mparr

        return mpsini, mparr_all

    else:
        return mpsini
    
#this should set limits and call lsqmdl, should be callable by bootstrapper...
def rvfit_lsqmdl(orbel,jdb,rv,srv,jitter=0, param_names=0,trend=0,circ=0,curv=0, tt=np.zeros(1),epoch=2.455e6,pfix=1):
    
    if jitter > 0.0: #this should happen in rvfit_lsqmdl
        nsrv = np.sqrt(srv**2 + jitter**2)
    else:
        nsrv = srv

    norbits = orbel.size/7
    ip = np.arange(norbits)

    if param_names == 0: #do something more clever here later
        param_names = ['Per', 'Tp', 'ecc', 'om', 'K1', 'gamma', 'dvdt']*norbits 
        if curv == 1:
            param_names.extend(['curv']) 
    
    m = lsqmdl.Model(None, rv, 1./nsrv) #m is a class
    m.set_func(rv_drive,param_names, args=[jdb] )

    #print pfix, ' = pfix'
    #make some reasonable limits - don't need the loop?
    for i in range(norbits):

        #apply normal limits first
        #fix/limit period
        if pfix[i] == 1:
            m.lm_prob.p_value(0+i*7, orbel[0+i*7], fixed=True)
        else:
            m.lm_prob.p_limit(0+i*7, lower=0.1, upper=(np.max(jdb)-np.min(jdb))*10) #per no longer than 10x range of survey 
        
        #limit Tp
        m.lm_prob.p_limit(1+i*7, lower=orbel[1+i*7]-orbel[0+i*7]/4., upper=orbel[1+i*7]+orbel[0+i*7]/4.) #T0 within one period guess

        #fix/limit ecc
        if circ[i] == 1:
            m.lm_prob.p_value(2+i*7, 0.0, fixed=True) 
        else:
            m.lm_prob.p_limit(2+i*7, lower=0.0, upper=0.99) #ecc must be physical 
        #limit omega
        m.lm_prob.p_limit(3+i*7, lower=0.0, upper=360.0)

        #limit K
        m.lm_prob.p_limit(4+i*7, lower=0.0, upper=1.0e5) #K must be physical
  
        #fix gamma except first
        if i > 0:
            m.lm_prob.p_value(5+i*7, 0.0, fixed=True)
        
        #optionally fix trend
        if i == 0:              #use a different logic statement...
            if trend == 0:
                m.lm_prob.p_value(6+i*7, 0.0, fixed=True) 
            #else set some reasonable limits?
        else:
            m.lm_prob.p_value(6+i*7, 0.0, fixed=True)
      
        #now include known transit effects
        if not tt[i] == 0:
            if circ[i] == 1:  #by convention tt=tp & omega=90
                m.lm_prob.p_value(1+i*7, tt[i], fixed=True)
                m.lm_prob.p_value(3+i*7, 90.0, fixed=True)
            #else:
                #tie omega to ecc

    #if orbel.size % 7 == 1:
            #set some reasonable limits?

    m.solve(orbel)
    #m.print_soln()
    return m

def rv_drive(orbel, t):
    #From rv_drive.pro in RVLIN by JTW

    if orbel.size > 20:
        print 'Warning: large number of params - are you sure you put the args in the right order?'

    rv = np.zeros_like(t)
    nplanets = orbel.size/7 

    phase = np.zeros((rv.size,nplanets))
    
    for i in range(nplanets):  #get rid of the for loops...
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
        
        E1 = kepler(M,ecc) #returns a matrix, because fan
        #E1 = E1[0,:]       #1 planet at a time
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
        epoch = 0.0 #Keep this
        
        rv = rv + k*(np.cos(theta + om) + ecc*np.cos(om)) + gamma + dvdt*(t - epoch) + curv*(t - epoch)**2 

    return rv



#Iteratively solve for E (anomoly) given M (mean anomoly) and e (eccentricity)
#@profile
def kepler(inM,inecc):
    #Port of kepler.pro by JTW, from Murray & Dermott p. 35, from Danby (1988)
    #SJG Aug 2014

    #make floats into arrays
    marr = ut.arrayify(inM)
    ecc = ut.arrayify(inecc)
    
#    if len(np.array(inM).shape) == 0:
#        marr = np.array([inM])
#    else:
#        marr = np.array(inM)
#    if len(np.array(inecc).shape) == 0:
#        ecc = np.array([inecc])
#    else:
#        ecc = np.array(inecc)
    
    nm = marr.size    #nm = nrv, ~100
    nec = ecc.size #nec = 1
    #if nec == 1 and nm > 1:   #Get rid of this...
    #    eccarr = fan(eccarr,nm)
    #if nec > 1 and nm == 1:
    #    marr = fan(marr,nec)

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
    Earr = marr + ssm*k*ecc
    fiarr = (Earr - ecc*np.sin(Earr) - marr) #E - ecc*sin(E) - M, converge to 0

    convd = np.squeeze(np.where(np.abs(fiarr) > conv)) #which indices converged?
    nd = convd.size
    #print nd
    count = 0
    while nd > 0:          #while unconverged elements exist...
        count += 1
        M = marr[convd]    #just keep the unconverged elements
        #ecc = eccarr[convd]
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
        fiarr = (Earr - ecc*np.sin(Earr) - marr)
        convd = np.squeeze(np.where(np.abs(fiarr) > conv))
        nd = convd.size
        if count > 100:
            print 'WARNING!  Keplers equation not solved!!!'
            nd = 0

    return Earr

def bootstrap_rvs(bestpar, jdb, rv, srv,nboot=1000,jitter=0,circ=0,trend=0,curv=0,tt=np.zeros(1),pfix=1):
    #based on the general method of bootpar.pro by SXW
    #Wow, this is slow! Make faster

    bestfit = rv_drive(bestpar, jdb)
    #print bestpar.size, bestpar.shape
    resid = rv - bestfit

    #nboot = np.max([nboot,resid.size*np.log10(resid.size)**2]).astype(int) #why this limit?
    print 'nboot = ',nboot

    bootpar = np.zeros((nboot,bestpar.size))
    
    for i in range(nboot): #get rid of the for loop...
        if i%100 == 0:
            print (i+0.0)/nboot*100,'%'
        #print i 
        scramble = np.random.randint(jdb.size,size=jdb.size)#random indices with replacement
        tmprv = resid[scramble] + bestfit            #scramble residuals
        tmperr = srv[scramble]                      #error goes with individual residual

        mi = rvfit_lsqmdl(bestpar, jdb, tmprv, tmperr, jitter=jitter,circ=circ,curv=curv,trend=trend,tt=tt,pfix=pfix)
        bootpar[i,:] = mi.params

    meanpar = np.mean(bootpar,axis=0)
    sigpar = np.std(bootpar,axis=0)
    
    return bootpar, meanpar, sigpar

#to run as ...orbits.py
if __name__ == '__main__':
    orbits_test(webdat='yes',nboot=10)


#begin MCMC setup - following emcee line-fitting demo
def lnprior(theta, fullpars, flt, pnames):
    
    #some pretty basic limits for flat priors
    low = np.array([0.1, 0.0, 0.0, 0.0, 0, -1e6, -1e3, -1])
    high = np.array([10*365.25, 3e6, 0.99, 360.0, 1e4, 1e6, 1e3, 1])


    #figure out which params are varied and the associated limits
    pfloat = pnames[flt.nonzero()]
    lfloat = low[flt.nonzero()]
    hfloat = high[flt.nonzero()]


    if (theta > lfloat).all() and (theta < hfloat).all():
        return 0.0
    else:
        return -np.inf
    

def lnprob(theta, jdb, rv, srv, fullpars, flt, pnames):
    lp = lnprior(theta, fullpars, flt, pnames)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, jdb, rv, srv, fullpars, flt)

def lnlike(theta, jdb, rv, srv, fullpars, flt):
    
    #TODO: combine theta with fixed pars!
    newpars = np.copy(fullpars)
    newpars[flt.nonzero()] = theta
    
    model = rv_drive(newpars, jdb)
    
    return -0.5*np.sum((rv - model)**2/srv**2) #chisq

def setup_emcee(targname, m, jdb, rv, srv, nwalkers=200, circ=0, trend=0, curv=0, tt=np.zeros(1),jitter=0): 

    bestpars = m.params

    if jitter > 0.0: #this should happen in rvfit_lsqmdl
        nsrv = np.sqrt(srv**2 + jitter**2)
    else:
        nsrv = srv

    #p = orbel[0+i*7]
    #tp = orbel[1+i*7]
    #ecc = orbel[2+i*7]
    #om = orbel[3+i*7] *np.pi/180. #degrees to radians
    #k = orbel[4+i*7]
    #gamma = orbel[5+i*7]
    #dvdt = orbel[6+i*7]
    #curv = 0 

    #How many planets?
    npars = bestpars.size
    
    norbits = npars/7

    ip = np.arange(norbits) #index for planets
    if norbits > 1:
        ip2 = np.arange(norbits-1)+1 #index for planets after 1st one

    #separate the params being varied from the full list
    flt = np.ones(npars) #all params float now, turn off individually

    #force fix period for now
    flt[0+ip*7] = 0
    
    if circ > 0:   
        flt[2+ip*7] = 0  #fix ecc - this fixes all planets, BAD
        if not tt[0] == 0:
            flt[1+ip*7] = 0 #if tt known, fix tp - also all planets
            flt[3+ip*7] = 0 #also fix omega - also all planets

    #want some testing that the output of LM is sensible!
    #need to consider eccentric transiting planets...

    if norbits > 1:
        flt[5+ip2*7] = 0 #always force fix all gammas except 1st
        flt[6+ip2*7] = 0 #same for dvdt
    
    #if curv = 1, flt[-1] = 1, but this already happens
    print flt
    varpars = bestpars[flt.nonzero()]
    ndim = varpars.size
    pnames = np.copy(m.pnames)
    print 'MCMC params: ',pnames[flt.nonzero()] 
    
    samples = run_emcee(targname, bestpars, varpars, flt, jdb, rv, srv, pnames, ndim, nwalkers=nwalkers)

    return m, flt, samples
   # return bestpars, varpars, flt, pnames


def run_emcee(targname, bestpars, varpars, flt, jdb, rv, srv, pnames, ndim, nwalkers=200, nsteps=1000):
    
    #Initialize walkers in tiny Gaussian ball around MLE results
    #number of params comes from varpars
    #pos = [varpars + 1e-3*np.random.randn(ndim) for i in range(nwalkers)] #??
    pos = [varpars + 1e-2*varpars*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(jdb,rv,srv,bestpars,flt, pnames))

    #Run MCMC
    sampler.run_mcmc(pos, nsteps)

    #It takes a number of iterations to spread walkers throughout param space
    #This is 'burning in'
    samples = sampler.chain[:, 200:, :].reshape((-1, ndim))

    fig = triangle.corner(samples, labels=pnames[flt.nonzero()], truths=bestpars[flt.nonzero()]) #makes nice plots...
    fig.savefig('/home/sgettel/Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+'_triangle.png')
    plt.close(fig)

    return sampler.chain
