#RV orbit modeling code based on RVLIN by Jason Wright (IDL) and orbits.f by Alex Wolszczan (FORTRAN77)

#
# Branch master
#


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
# - jitter term
# - plot trend if there is one
# - histogram plots
# - read orbit params from somewhere
# - add offset terms
# - test MCMC


#Given a target name, do everything! Eventually.
def orbits_test(targname='K00273',jitter=0.0,nboot=1000,epoch=2.455e6,circ=0,maxrv=1e6,minrv=-1e6,maxsrv=5, webdat='no', nwalkers=200, pfix=1,norbits=1,npoly=0):

    if npoly > 4:
        print 'Must have <= 4th order polynomial'


    transit = np.zeros(1)
    
    circ = ut.arrayify(circ)
    pfix = ut.arrayify(pfix)
    

    if socket.gethostname() == 'sara-gettels-macbook-2.local':
    	home = '/Users/Sara/'
    else:
    	home = '/home/sgettel/'

    #read RV data
    if targname == 'HD209458':
        #demo!
        sfile = open(home+'Dropbox/cfasgettel/py.lib/sgcode/rvorbits/209458.vel.txt')
        jdb, rv, srv = np.loadtxt(sfile,unpack=True,usecols=(0,1,2))
  
        #fake two telescopes
        telvec = np.zeros_like(jdb)
        telvec[0:10] = 1
        rv[0:10] += 30.0 + np.random.randn(10) #noiser data with offset

    else:
        print targname
        jdb, rv, srv, labels, fwhm, contrast, bis_span, rhk, sig_rhk = rr.process_all(targname,maxsrv=maxsrv,maxrv=maxrv,minrv=minrv)
        jdb += 2.4e6 #because process_all gives truncated JDBs
     
    #adjust values to be sensible
    #if jdb[0] > epoch:
    #    print 'truncating dates'
    tnorm = jdb - epoch #truncate
    rvnorm = rv - np.median(rv)
    #print np.min(rvnorm),np.max(rvnorm),np.min(rv),np.max(rv)
    print np.min(tnorm),np.max(tnorm)

    if jitter > 0.0: 
        nsrv = np.sqrt(srv**2 + jitter**2)
        print 'Adding ',str(jitter),' m/s fixed jitter'
    else:
        nsrv = srv



#    #process offsets  
#    ntel = np.unique(telvec).size
#    if ntel > 1:
#        print 'Including ',ntel-1,' offset terms'
#        offset = np.zeros(ntel-1)



    #read/process parameter guesses from somewhere
    #p = orbel[0+i*6]
    #tp = orbel[1+i*6]
    #ecc = orbel[2+i*6]
    #om = orbel[3+i*6] *np.pi/180. #degrees to radians
    #k = orbel[4+i*6]
    #gamma = orbel[5+i*6]

    #dvdt = orbel[0+norbits*6] - optional polynomial fit
    #quad = orbel[1+norbits*6]
    #cubic = orbel[2+norbits*6]
    #quart = orbel[3+norbits*6]

    if targname == 'HD209458':
        guesspars = np.array([3.524733, 2452826.628514, 0.0, 336.5415, 85.49157+30, -1.49-50, 0.01])#HD209458
        transit = np.array([2452826.628514]) 
        mstar = 1.0
    
    if targname == 'K00273':
        guesspars = np.array([10.573769, 2455008.06601, 0.0, 90.0, 1.7358979, -3398.0498, 1.1889011, 0.010])#, 0.0, 0.0]) #K00273
        transit = np.array([2455008.06601,0.0]) 
        mstar = 1.07
        #2 planets...
        #guesspars = np.array([10.573769, 2455008.06601, 0.0, 90.0, 1.7358979, -3398.0498,  530.0, 2455008.066, 0.0, 90.0, 100.0,0.0])

    
    if targname == 'K00069':
        guesspars = np.array([4.72673978, 2454944.29227, 0.0, 90.0, 1.733, -91.08, 0.0329])
        transit = np.array([2454944.29227])
        mstar = 0.887

    
    ip = np.arange(norbits)
    guesspars[1+ip*6] -= epoch

    for i in range(transit.size):  #do something more clever here...
        if not transit[i] == 0.0:
            transit[i] -= epoch 
            #print 'set transit time: ', transit
 


    
    m = rvfit_lsqmdl(guesspars, tnorm, rvnorm, nsrv, jitter=jitter,circ=circ, npoly=npoly,tt=transit,epoch=epoch,pfix=pfix,norbits=norbits)


    
    #display initial fit - want to show fixed params too
    m.print_soln()  
    
    
    #mass estimate
    mpsini = mass_estimate(m, mstar, norbits=norbits)
    print 'mp*sin(i):         ',str(mpsini)
    
    #make plots
    plot_rv(targname,tnorm,rvnorm,nsrv,guesspars,m,nmod=200,home=home,norbits=norbits,npoly=npoly)
    
    #call bootstrapping
    if nboot > 0:
        bootpars, meanpar, sigpar = bootstrap_rvs(m.params, tnorm, rvnorm, nsrv,nboot=nboot,jitter=jitter,circ=circ,npoly=npoly,tt=transit,pfix=pfix,norbits=norbits)
   
        mpsini, mparr_all = mass_estimate(m, mstar, norbits=norbits, bootpar=bootpars)
       
        print_boot_errs(meanpar, sigpar, mpsini, mparr_all, norbits=norbits,npoly=npoly)

        #return m0, bootpar,sigpar, mparr_all #jdb, rv, nsrv
    else:
        bootpars = -1
        mparr_all = -1

    #call MCMC    
    if nwalkers > 0:
        bestpars, pnames, flt, samples, mcpars = setup_emcee(targname, m, tnorm, rvnorm, nsrv, circ=circ, npoly=npoly, tt=transit, jitter=jitter, nwalkers=nwalkers, pfix=pfix)
#        m, flt, chain, samples, mcpars = setup_emcee(targname, m, tnorm, rvnorm, nsrv, circ=circ, npoly=npoly, tt=transit, jitter=jitter, nwalkers=nwalkers, pfix=pfix)
        mpsini, mparr_mc = mass_estimate(m, mstar, norbits=norbits, mcpar=mcpars)
        #print output from mass_estimate for mc
        print_mc_errs(mcpars, mpsini, mparr_mc,norbits=norbits,npoly=npoly)

        #make a nice triangle plot
        
        print pnames
        f = np.squeeze(flt.nonzero())
        fig = triangle.corner(samples, labels=pnames[f], truths=bestpars[f]) 
        fig.savefig('/home/sgettel/Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+'_triangle.png')
        plt.close(fig)
    else:
        mcpars = -1
        mparr_mc = -1

        #return m, flt, chain, samples, mcpars# bestpars, varpars, flt, pnames # 
    write_full_soln(m, targname, mpsini, bootpars=bootpars, mparr_all = mparr_all, mcpars=mcpars, mparr_mc=mparr_mc, norbits=norbits, npoly=npoly)


    return m

def plot_rv(targname,jdb,rv,srv,guesspars,m,nmod=1000,home='/home/sgettel/', norbits=1,npoly=0):
   
    tmod = np.linspace(np.min(jdb),np.max(jdb),nmod)

    
    model_init = rv_drive(guesspars,tmod,norbits,npoly) 
    model_final = rv_drive(m.params,tmod,norbits,npoly)
    
    if npoly > 0:
        parst =  np.copy(m.params)
        parst[4] = 0.0
        poly = rv_drive(parst,tmod,norbits,npoly)

    #unphased data
    plt.figure(1)
    plt.errorbar(jdb,rv,yerr=srv,fmt='bo')
    plt.plot(tmod,model_final,'r-')
    plt.plot(tmod,poly,'g-')
    #plt.plot(tmod,model_init,'g-')
    plt.savefig(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+'_autoplot.png')
    plt.close(1)


    #phase at 1st period
    pars = m.params[0:6]  #for model 
    pars[6:] = 0 #select first planet only
        
    parst = np.copy(m.params)
    parst[4] = 0.0 #other planets only
    parst[5] = 0.0
    phase = (jdb - pars[1])/pars[0] % 1.0
    rvt = rv_drive(parst,jdb,norbits,npoly)   
    plt.figure(2)
    plt.errorbar(phase, rv-rvt, yerr=srv,fmt='bo')
    plt.plot((tmod - pars[1])/pars[0] % 1.0, rv_drive(pars, tmod,1,0),'r.')
    #plt.plot((tmod - guess))
    plt.savefig(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+'_phase_autoplot.png')
    plt.close(2)

#def plot_pars(targname,bootpars,mparr_all,norbits=1):
#    print targname
#    plot histograms!

def write_full_soln(m,targname,mpsini, bootpars=-1, mparr_all=-1, mcpars=-1, mparr_mc=-1,norbits=1,npoly=0):
    
    poly_names = ['dvdt:  ','quad:  ', 'cubic: ','quart: ']
    
    f = open('/home/sgettel/Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+'_orbit.dat','w')
    for i in range(norbits):
        f.write('                                               \n') 
        f.write('*****Planet '+str(i+1)+' Solution:***** \n')
        f.write('Per: '+str(m.params[0+i*6])+'\n')
        f.write('Tp: '+str(m.params[1+i*6])+'\n')
        f.write('ecc: '+str(m.params[2+i*6])+'\n')
        f.write('om: '+str(m.params[3+i*6])+'\n')
        f.write('K1: '+str(m.params[4+i*6])+'\n')
        f.write('gamma: '+str(m.params[5+i*6])+'\n')

        f.write('mp*sin(i): '+str(mpsini[i])+'\n')
    
    for i in range(npoly):
        f.write(str(poly_names[i])+str(m.params[i+norbits*6]) +'\n')
    f.write(str('r chi sq: ')+str(m.rchisq)+'\n')

    if len(np.array(bootpars).shape) > 0: #print bootstrap errs
        for i in range(norbits):
            f.write('                                               \n')
            f.write('*****Planet '+str(i+1)+' Bootstrap Errors:***** \n')
            f.write('Per: '+ str(np.mean(bootpars[:,0+i*6]))+'+/-'+str(np.std(bootpars[:,0+i*6]))+'\n')
            f.write('Tp: '+ str(np.mean(bootpars[:,1+i*6]))+'+/-'+str(np.std(bootpars[:,1+i*6]))+'\n')
            f.write('ecc: '+ str(np.mean(bootpars[:,2+i*6]))+'+/-'+str(np.std(bootpars[:,2+i*6]))+'\n')
            f.write('om: '+ str(np.mean(bootpars[:,3+i*6]))+'+/-'+str(np.std(bootpars[:,3+i*6]))+'\n')
            f.write('K1: '+ str(np.mean(bootpars[:,4+i*6]))+'+/-'+str(np.std(bootpars[:,4+i*6]))+'\n')
            f.write('gamma: '+ str(np.mean(bootpars[:,5+i*6]))+'+/-'+str(np.std(bootpars[:,5+i*6]))+'\n')

            f.write('mp*sin(i): '+str(np.mean(mparr_all[i,:]))+'+/-'+str(np.std(mparr_all[i,:]))+'\n')
            f.write('mass error:'+ str(np.std(mparr_all[i,:])/mpsini[i]*100)+'%'+'\n')

        for i in range(npoly):
            f.write(str(poly_names[i])+str(np.mean(bootpars[:,i+norbits*6]))+'+/-'+str(np.std(bootpars[:,i+norbits*6])) +'\n') 

    if len(np.array(mcpars).shape) > 0: #print bootstrap errs
        for i in range(norbits):
            f.write('                                               \n')
            f.write('*****Planet '+str(i+1)+' MCMC Errors:***** \n')    
            f.write('Per: '+ str(np.mean(mcpars[:,0+i*7]))+'+/-'+str(np.std(mcpars[:,0+i*7]))+'\n')
            f.write('Tp: '+ str(np.mean(mcpars[:,1+i*7]))+'+/-'+str(np.std(mcpars[:,1+i*7]))+'\n')
            f.write('ecc: '+ str(np.mean(mcpars[:,2+i*7]))+'+/-'+str(np.std(mcpars[:,2+i*7]))+'\n')
            f.write('om: '+ str(np.mean(mcpars[:,3+i*7]))+'+/-'+str(np.std(mcpars[:,3+i*7]))+'\n')
            f.write('K1: '+ str(np.mean(mcpars[:,4+i*7]))+'+/-'+str(np.std(mcpars[:,4+i*7]))+'\n')
            f.write('gamma: '+ str(np.mean(mcpars[:,5+i*7]))+'+/-'+str(np.std(mcpars[:,5+i*7]))+'\n')
            
            f.write('mp*sin(i): '+str(np.mean(mparr_mc[i,:]))+'+/-'+str(np.std(mparr_mc[i,:]))+'\n')
            f.write('mass error:'+ str(np.std(mparr_mc[i,:])/mpsini[i]*100)+'%'+'\n') 

        for i in range(npoly):
            f.write(str(poly_names[i])+str(np.mean(mcpars[:,i+norbits*6]))+'+/-'+str(np.std(mcpars[:,i+norbits*6])) +'\n')
    f.close()

def print_boot_errs(meanpar,sigpar, mpsini, mparr_all,norbits=1,npoly=0):
 #   print mparr_all.shape
    poly_names = ['dvdt:  ','quad:  ', 'cubic: ','quart: ']
    for i in range(norbits):
        print '                                               '
        print '*****Planet ',str(i+1),' Bootstrap Errors:*****'
        print 'Per: ', str(meanpar[i*6]),'+/-',str(sigpar[i*6])
        print 'Tp: ', str(meanpar[i*6+1]),'+/-',str(sigpar[i*6+1])
        print 'ec: ', str(meanpar[i*6+2]),'+/-',str(sigpar[i*6+2])
        print 'om: ', str(meanpar[i*6+3]),'+/-',str(sigpar[i*6+3])
        print 'K: ', str(meanpar[i*6+4]),'+/-',str(sigpar[i*6+4])
        print 'gamma: ', str(meanpar[i*6+5]),'+/-',str(sigpar[i*6+5])

        print 'mp*sin(i): ',str(np.mean(mparr_all[i,:])),'+/-',str(np.std(mparr_all[i,:]))
        print 'mass error:', str(np.std(mparr_all[i,:])/mpsini[i]*100),'%'

    for i in range(npoly):
        print str(poly_names[i]), str(meanpar[i+norbits*6]), '+/-', str(sigpar[i+norbits*6])
    return

def print_mc_errs(mcpars, mpsini, mparr_all,norbits=1,npoly=0):
    poly_names = ['dvdt:  ','quad:  ', 'cubic: ','quart: ']
    for i in range(norbits):
        print '                                               '
        print '*****Planet ',str(i+1),' MCMC Errors:*****'  
        print 'Per: ', str(np.mean(mcpars[:,0+i*6])),'+/-',str(np.std(mcpars[:,0+i*6]))
        print 'Tp: ', str(np.mean(mcpars[:,1+i*6])),'+/-',str(np.std(mcpars[:,1+i*6]))
        print 'Ecc: ', str(np.mean(mcpars[:,2+i*6])),'+/-',str(np.std(mcpars[:,2+i*6]))
        print 'Om: ', str(np.mean(mcpars[:,3+i*6])),'+/-',str(np.std(mcpars[:,3+i*6]))
        print 'K: ', str(np.mean(mcpars[:,4+i*6])),'+/-',str(np.std(mcpars[:,4+i*6]))
        print 'gamma: ', str(np.mean(mcpars[:,5+i*6])),'+/-',str(np.std(mcpars[:,5+i*6]))

        print 'mp*sin(i): ',str(np.mean(mparr_all[i,:])),'+/-',str(np.std(mparr_all[i,:]))
        print 'mass error:', str(np.std(mparr_all[i,:])/mpsini[i]*100),'%'
    
    for i in range(npoly):
        print str(poly_names[i]), str(np.mean(mcpars[:,i+norbits*6])),'+/-',str(np.std(mcpars[:,i+norbits*6]))
    print 'jitter: ', str(np.mean(mcpars[:,-1])),'+/-',str(np.std(mcpars[:,-1]))
    return

def mass_estimate(m,mstar,norbits=1,bootpar=-1,mcpar=-1):
   
    #some constants, maybe use astropy here
    msun = 1.9891e30
    mearth = 5.97219e24
    G = 6.673e-11 #m^3 kg^-1 s^-2
    etoj = 317.83
    
    ip = np.array(range(norbits))
    
    pers = m.params[ip*6]
    eccs = m.params[ip*6+2]
    amps = m.params[ip*6+4]
    
    #mass estimate
    fm = (1 - eccs*eccs)**(1.5)*amps**3*(pers*86400.0)/(2*np.pi*G) #kg
    mpsini = ((mstar*msun)**2*fm)**(1./3.)/mearth

    #calculate error on mass if bootpar array is input
    if len(np.array(bootpar).shape) > 0:

        #print bootpar.shape, m.params.shape
        mparr_all = np.zeros((norbits,bootpar.shape[0]))

        for i in ip:
            fmarr = (1 - bootpar[:,i*6+2]**2)**(1.5)*bootpar[:,i*6+4]**3*(bootpar[:,i*6]*86400.0)/(2.0*np.pi*G)
            mparr = ((mstar*msun)**2*fmarr)**(1./3.)/mearth
            mparr_all[i,:] = mparr

        return mpsini, mparr_all

   # mass estimate for mcmc
    if len(np.array(mcpar).shape) > 0:
        mparr_mc = np.zeros((norbits,mcpar.shape[0]))

        for i in ip:
            fmarr = (1 - mcpar[:,2+i*6]**2)**(1.5)*mcpar[:,i*6+4]**3*(mcpar[:,i*6]*86400.0)/(2.0*np.pi*G)
            mparr = ((mstar*msun)**2*fmarr)**(1./3.)/mearth
            mparr_mc[i,:] = mparr
        return mpsini, mparr_mc
    else:
        return mpsini

   
#this should set limits and call lsqmdl, should be callable by bootstrapper...

def rvfit_lsqmdl(orbel,jdb,rv,srv,jitter=0, param_names=0,npoly=0,circ=0, tt=np.zeros(1),epoch=2.455e6,pfix=1,norbits=1):

#def rvfit_lsqmdl(orbel,jdb,rv,srv,jitter=0, param_names=0,npoly=0,circ=0, tt=np.zeros(1),epoch=2.455e6,pfix=1,norbits=1,noffsets=0,telvec=-1):

    
    

    ip = np.arange(norbits)

    if param_names == 0: #do something more clever here later
        param_names = ['Per', 'Tp', 'ecc', 'om', 'K1', 'gamma']*norbits
        poly_names = ['dvdt','quad','cubic','quart']
        param_names.extend(poly_names[:npoly]) 

    m = lsqmdl.Model(None, rv, 1./srv) #m is a class
    #m.set_func(rv_drive,param_names, args=[jdb] )
    m.set_func(rv_drive,param_names, args=(jdb,norbits,npoly) )

    #print pfix, ' = pfix'
    #make some reasonable limits - don't need the loop?
    for i in range(norbits):

        #apply normal limits first
        #fix/limit period
        if pfix[i] == 1:
            m.lm_prob.p_value(0+i*6, orbel[0+i*6], fixed=True)
        else:
            m.lm_prob.p_limit(0+i*6, lower=0.1, upper=(np.max(jdb)-np.min(jdb))*20) #per no longer than 10x range of survey 
        
        #limit Tp
        m.lm_prob.p_limit(1+i*6, lower=orbel[1+i*6]-orbel[0+i*6]/2., upper=orbel[1+i*6]+orbel[0+i*6]/2.) #T0 within one period guess

        #fix/limit ecc
        if circ[i] == 1:
            m.lm_prob.p_value(2+i*6, 0.0, fixed=True) 
        else:
            m.lm_prob.p_limit(2+i*6, lower=0.0, upper=0.99) #ecc must be physical 
        #limit omega
        m.lm_prob.p_limit(3+i*6, lower=0.0, upper=360.0)

        #limit K
        m.lm_prob.p_limit(4+i*6, lower=0.0, upper=1.0e5) #K must be physical
  
        #fix gamma except first
        if i > 0:
            m.lm_prob.p_value(5+i*6, 0.0, fixed=True)

        #now include known transit effects
        if not tt[i] == 0:
            if circ[i] == 1:  #by convention tt=tp & omega=90
                m.lm_prob.p_value(1+i*6, tt[i], fixed=True)
                m.lm_prob.p_value(3+i*6, 90.0, fixed=True)
            else:
                tiefunc = tie_omega_function(tt, i) #how does this know what orbel is?
                m.lm_prob.p_tie(3+i*6, tiefunc)

    #limit polynomial terms
    for i in range(npoly):
        
        m.lm_prob.p_limit(i+norbits*6, lower=-1e6, upper=1e6) #dvdt and higher


    m.solve(orbel)
    #m.print_soln()
    return m

#OO? magic from Peter...
def tie_omega_function(tt, i):

    def calculate_omega(orbel):
        
        p = orbel[0+i*6]
        tp = orbel[1+i*6]
        ecc = orbel[2+i*6]

        theta_tt = calc_true_anomaly(p, tp, ecc, tt[i]) #in radians
        omega = (np.pi/2.0 - theta_tt)*180.0/np.pi
         
        if omega < 0:     #do something more clever
            omega += 360.0
        if omega > 360.0:
            omega -= 360.0
        
        return omega

    return calculate_omega

def rv_drive(orbel, t, norbits, npoly):
    #From rv_drive.pro in RVLIN by JTW

    if orbel.size > 20:
        print 'Warning: large number of params - are you sure you put the args in the right order?'

    rv = np.zeros_like(t)
   
    phase = np.zeros((rv.size,norbits))
    
    for i in range(norbits):  #get rid of the for loops...
        p = orbel[0+i*6]
        tp = orbel[1+i*6]
        ecc = orbel[2+i*6]
        om = orbel[3+i*6] *np.pi/180. #degrees to radians
        k = orbel[4+i*6]
        gamma = orbel[5+i*6]


#        dvdt = orbel[6+i*7]
#        curv = 0
#        if i == 0 and nplanets*7 == orbel.size-1: # if 1 too many elements,
#            curv = orbel[-1]                      # last is curvature

        #Error checking
        if p < 0 or ecc < 0 or ecc >= 1 or k < 0:
            print 'Bad inputs to rv_drive'
            print p, ecc, k
            if p < 0:
                p = 1e-2
            if ecc < 0:
                ecc = 0
            if ecc >= 1:
                ecc = 0.99
            if k < 0:
                k = 1e-2
     
        theta = calc_true_anomaly(p, tp, ecc, t)


        phase0 = theta + om - np.pi/2.0
        under = np.squeeze(np.where(phase0 < -np.pi)) 
        phase0[under] += np.pi
        over = np.squeeze(np.where(phase0 >= np.pi))
        phase0[over] -= np.pi
        
        phase[:,i] = phase0

        #calculate radial velocity
        epoch = 0.0 #Yes, epoch corrected elsewhere
            
        rv = rv + k*(np.cos(theta + om) + ecc*np.cos(om)) + gamma #+ dvdt*(t - epoch) + curv*(t - epoch)**2

    #now add polynomial
    for i in range(npoly):
        rv = rv + orbel[i+norbits*6]*(t - epoch)**(i+1)

    return rv

def calc_true_anomaly(p, tp, ecc, t):
    #first calculate the approximate eccentric anomaly, E1, from the mean anomaly, M
    M = 2.0*np.pi*( ((t-tp)/p) - np.floor((t-tp)/p) ) #phase in radians
    E1 = kepler(M,ecc) 

    #calculate true anomaly
    n1 = 1.0 + ecc
    n2 = 1.0 - ecc
    theta = 2.0*np.arctan(np.sqrt(n1/n2)*np.tan(E1/2.0))

    return theta 

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

def bootstrap_rvs(bestpar, jdb, rv, srv,nboot=1000,jitter=0,circ=0,npoly=0,tt=np.zeros(1),pfix=1,norbits=1):
    #based on the general method of bootpar.pro by SXW
    
    bestfit = rv_drive(bestpar, jdb,norbits,npoly)
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

        mi = rvfit_lsqmdl(bestpar, jdb, tmprv, tmperr, jitter=jitter,circ=circ,norbits=norbits,npoly=npoly,tt=tt,pfix=pfix)
        bootpar[i,:] = mi.params

    meanpar = np.mean(bootpar,axis=0)
    sigpar = np.std(bootpar,axis=0)
    
    return bootpar, meanpar, sigpar

#to run as ...orbits.py
if __name__ == '__main__':
    orbits_test(webdat='yes',nboot=10)


#begin MCMC setup - following emcee line-fitting demo
def lnprior(theta, fullpars, flt, pnames, plo, phi):
    
    #some pretty basic limits for flat priors - I do this in setup_emcee
    #low = np.array([0.1, 0.0, 0.0, 0.0, 0, -1e6, -1e3, -1])
    #high = np.array([10*365.25, 3e6, 0.99, 360.0, 1e4, 1e6, 1e3, 1])
    #print theta

    #figure out which params are varied and the associated limits
    pfloat = pnames[flt.nonzero()]
    lfloat = plo[flt.nonzero()]
    hfloat = phi[flt.nonzero()]


    if (theta >= lfloat).all() and (theta < hfloat).all():
        return 0.0
    else:
        return -np.inf
    

def lnprob(theta, jdb, rv, srv, fullpars, flt, pnames, plo, phi, norbit, npoly):
    lp = lnprior(theta, fullpars, flt, pnames, plo, phi)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, jdb, rv, srv, fullpars, flt, norbit, npoly)

def lnlike(theta, jdb, rv, srv, fullpars, flt, norbit, npoly):
    
    newpars = np.copy(fullpars)
    newpars[flt.nonzero()] = theta
    #print newpars
    model = rv_drive(newpars, jdb, norbit, npoly)
    nsrv = np.sqrt(srv**2 + newpars[-1]**2) #add floating jitter term

    return -0.5*np.sum((rv - model)**2/nsrv**2) #chisq

def setup_emcee(targname, m, jdb, rv, srv_in, nwalkers=200, circ=0, npoly=0, norbits=1, tt=np.zeros(1),jitter=0, pfix=1,nburn=200): 

    bestpars = np.copy(m.params)
    bestpars = np.append(bestpars,0.001) #add placeholder for jitter
    pnames = np.copy(m.pnames)
    pnames = np.append(pnames,'jitter')

    #p = orbel[0+i*6]
    #tp = orbel[1+i*6]
    #ecc = orbel[2+i*6]
    #om = orbel[3+i*6] *np.pi/180. #degrees to radians
    #k = orbel[4+i*6]
    #gamma = orbel[5+i*6]
    #dvdt = orbel[0+norbits*6] - optional polynomial fit
    #quad = orbel[1+norbits*6]
    #cubic = orbel[2+norbits*6]
    #quart = orbel[3+norbits*6]
    #jitter = orbel[-1]

    npars = bestpars.size

    if jitter > 0:
        print 'removing ',str(jitter),' m/s fixed jitter'
        srv = np.sqrt(srv_in**2 - jitter**2)
    else:
        srv = srv_in

#vectorizing loop
#    ip = np.arange(norbits) #index for planets
#    if norbits > 1:
#        ip2 = np.arange(norbits-1)+1 #index for planets after 1st one

    #these will hold reasonable limits on priors
    plo = np.zeros(npars)
    phi = np.zeros(npars)

    #separate the params being varied from the full list
    flt = np.ones(npars) #all params float now, turn off individually  

    #probably doesn't need to be a loop...
    for i in range(norbits):
        
        #fix normal orbit params first...
        
        #fix/limit period:
        if pfix[i] == 1:
            flt[0+i*6] = 0
        plo[0+i*6] = 0.1
        phi[0+i*6] = (np.max(jdb)-np.min(jdb))*20

        #limit Tp
        plo[1+i*6] = -3e6
        phi[1+i*6] = 3e6

        #fix/limit ecc
        if circ[i] == 1:
            flt[2+i*6] = 0
        plo[2+i*6] = 0.0
        phi[2+i*6] = 0.99

        #limit omega
        plo[3+i*6] = 0.0
        phi[3+i*6] = 360.0

        #limit K
        plo[4+i*6] = 0.0
        phi[4+i*6] = 1.0e5

        #fix gamma except first
        if i > 0:
            flt[5+i*6] = 0
        plo[5+i*6] = -1e8
        phi[5+i*6] = 1e8

#        #optionally fix 1st trend, fix all after
#        if i == 0:              #use a different logic statement...
#            if trend == 0:
#                flt[6+i*7] = 0 
#        else:
#            flt[6+i*7] = 0 
#        plo[6+i*7] = -1e3
#        phi[6+i*7] = 1e3

        #now consider known transits...
        if not tt[i] == 0:
            if circ[i] == 1:
                flt[1+i*7] = 0
                bestpars[1+i*7] = tt[i] #this is not necessary?
                flt[3+i*7] = 0
                bestpars[3+i*7] = 90.0
            #else:
                #tie omega to ecc

    #limit polynomial terms
    for i in range(npoly):
        plo[i+norbits*6] = -1e6
        phi[i+norbits*6] = 1e6

    #limit jitter
    plo[-1] = 0.0
    phi[-1] = 0.01

    #want some testing that the output of LM is sensible!
       
    
    f = np.squeeze(flt.nonzero())
    varpars = bestpars[f]
    ndim = varpars.size
    
    print 'MCMC params: ',pnames[f] 
   

    chain = run_emcee(targname, bestpars, varpars, flt, plo, phi, jdb, rv, srv, pnames, ndim, nwalkers=nwalkers, norbits=norbits, npoly=npoly)
  
    #It takes a number of iterations to spread walkers throughout param space
    #This is 'burning in'
    samples = chain[:, nburn:, :].reshape((-1, ndim))
    
    #combine samples with best-fit values
    mcpars = ut.fan(bestpars,samples.shape[0])
    
    for i in range(ndim):
        mcpars[:,f[i]] = samples[:,i]

    #return m, flt, chain, samples, mcpars
    return bestpars, pnames, flt, samples, mcpars


def run_emcee(targname, bestpars, varpars, flt, plo, phi, jdb, rv, srv, pnames, ndim, nwalkers=200, nsteps=1000, norbits=1, npoly=0):
    
    #Initialize walkers in tiny Gaussian ball around MLE results
    #number of params comes from varpars
    #pos = [varpars + 1e-3*np.random.randn(ndim) for i in range(nwalkers)] #??
    pos = [varpars + 1e-3*varpars*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(jdb,rv,srv,bestpars,flt, pnames, plo, phi, norbits, npoly))

    #Run MCMC
    sampler.run_mcmc(pos, nsteps)
    print("Mean acceptance fraction: {0:.3f}"
                .format(np.mean(sampler.acceptance_fraction)))

    return sampler.chain


