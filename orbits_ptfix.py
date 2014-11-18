#RV orbit modeling code based on RVLIN by Jason Wright (IDL) and orbits.f by Alex Wolszczan (FORTRAN77)

#
# Branch bic
#


import emcee
import socket
import triangle
import numpy as np
import read_rdb_harpsn as rr
import matplotlib.pyplot as plt
import utils as ut
from pwkit import lsqmdl
from scipy.stats import norm

def orbits_test(targname='K00273',jitter=0.0,epoch=2.455e6,circ=0,maxrv=1e6,minrv=-1e6,maxsrv=5, webdat='no', nwalkers=200, pfix=1,norbits=1,npoly=0,keck='no',outer_loop='no',nsteps=1000,ttfloat='no'):


    if npoly > 4:
        print 'Must have <= 4th order polynomial'


    transit = np.zeros(1)
    telvec = np.zeros(1)
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

        s = np.argsort(jdb)
        jdb = jdb[s]
        rv = rv[s]
        srv = srv[s]

        #fake two telescopes
        telvec = np.zeros_like(jdb)
        telvec[80:] = 1
        rv[80:] += 80.0 #+ np.random.randn(telvec.size-80) #noiser data with offset
        rvnorm = rv - np.median(rv)

    elif targname == 'K00069':
        sfile = open(home+'Dropbox/cfasgettel/research/papers_collab/kep93_harpsn_data_4test.txt')
        jdb, rv, srv = np.loadtxt(sfile,unpack=True,usecols=(0,1,2),skiprows=1)
        rvnorm = rv - np.median(rv)

        jdb0 = np.copy(jdb)
        rv0 = np.copy(rv)
        srv0 = np.copy(srv)

        if keck == 'yes':
            sfile = open(home+'Dropbox/cfasgettel/research/papers_collab/kep93_hires_data_4test.txt')
            kjdb, krv, ksrv = np.loadtxt(sfile,unpack=True,usecols=(0,1,2),skiprows=1)
            krvnorm = krv - np.median(krv)
            ktel = np.ones_like(kjdb)
                
            jdb = np.append(jdb,kjdb)
            rvnorm = rv - np.median(rv)
            rvnorm = np.append(rvnorm,krvnorm)
            srv = np.append(srv,ksrv)
            telvec = np.append(telvec,ktel)
        
    else:


        print targname
        jdb, rv, srv, labels, fwhm, contrast, bis_span, rhk, sig_rhk = rr.process_all(targname,maxsrv=maxsrv,maxrv=maxrv,minrv=minrv)
        jdb += 2.4e6 #because process_all gives truncated JDBs
        telvec = np.zeros_like(jdb)
        print jdb.size,' HARPS-N obs'

        jdb0 = np.copy(jdb)
        rv0 = np.copy(rv)
        srv0 = np.copy(srv)

        if keck == 'yes':
            
            

            sfile = open(home+'Dropbox/cfasgettel/research/keck/'+targname+'.dat')
            kjdb, krv, ksrv = np.loadtxt(sfile,unpack=True,usecols=(2,3,4))
            kjdb = kjdb + 2.45e6 
            krvnorm = krv - np.median(krv)
            ktel = np.ones_like(kjdb)
                
            jdb = np.append(jdb,kjdb)
            rvnorm = rv - np.median(rv)
            rvnorm = np.append(rvnorm,krvnorm)
            srv = np.append(srv,ksrv)
            telvec = np.append(telvec,ktel)
            #print kjdb
        else:
            rvnorm = rv - np.median(rv)

    #adjust values to be sensible
    #if jdb[0] > epoch:
    #    print 'truncating dates'
    tnorm = jdb - epoch #truncate
  
    if jitter > 0.0: 
        nsrv = np.sqrt(srv**2 + jitter**2)
        print 'Adding ',str(jitter),' m/s fixed jitter'
    else:
        nsrv = srv



    #process offsets  
    ntel = np.unique(telvec).size
    if ntel > 1:
        print 'Including ',ntel-1,' offset terms'
        
        offset = np.zeros(ntel-1)
    


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

    #offset = orbel[norbits*6 + npoly + (0-3)] #up to 4 offset terms
        
    #ttfloat = orbel[norbits*6 + npoly + ntel-1] #up to norbits terms

    #jitter - final term, MCMC only
    rpl = 0.0

    if targname == 'HD209458':
        guesspars = np.array([3.524733, 2452826.628514, 0.0, 336.5415, 85.49157+10, -1.49-10])#HD209458
        transit = np.array([2452826.628514]) 
        mstar = 1.0
    
    if targname == 'K00273':

        #guesspars = np.array([10.573769, 2455008.06601, 0.0, 90.0, 1.7358979, -3398.0498, 1.1889011, 0.010, 0.0])#, 0.0]) #K00273
        guesspars = np.array([10.573737, 2455008.06787, 0.0, 90.0, 1.7358979, -3398.0498, 1.1889011, 0.010, 0.0])#, 0.0])
        #transit = np.array([2455008.06601,0.0])
        transit = np.array([2455008.06778,0.0]) 
        mstar = 1.07
        rpl = 1.82 #Me
        rple = 0.36

        #2 planets...       
        #guesspars = np.array([10.573769, 2455008.06601, 0.0, 90.0, 2.2, -16.0,  500.0, 2455041.9, 0.23, 40.0, 137.0,0.0])
        #guesspars = np.array([10.573738, 2455008.06778, 0.0, 90.0, 2.2, -16.0,  500.0, 2455041.9, 0.23, 40.0, 137.0,0.0])

    
    if targname == 'K00069':
        guesspars = np.array([4.72673978, 2454944.29227, 0.0, 90.0, 1.733, -91.08, 0.0329])
        transit = np.array([2454944.29227])
        mstar = 0.887

    
    ip = np.arange(norbits)
    guesspars[1+ip*6] -= epoch

    for i in range(transit.size):  
        if not transit[i] == 0.0:
            transit[i] -= epoch 
            #print 'set transit time: ', transit

    #append offsets if needed
    if ntel > 1:
        guesspars = np.append(guesspars,offset)
 
    if outer_loop == 'yes':
        #brute force search of outer period
        nboot = 0
        nwalkers = 0
        pguess = np.linspace(450,650,num=201)
        outchisq = np.zeros_like(pguess)
        outp = np.zeros_like(pguess)

        for i in range(pguess.size):
            guesspars[6] = pguess[i]
            if (pguess[i] == 487.0) or (pguess[i] == 520.0):
                continue
            print pguess[i]
            pfix = [1,1]
            m, flt = rvfit_lsqmdl(guesspars, tnorm, rvnorm, nsrv, jitter=jitter,circ=circ, npoly=npoly,tt=transit,epoch=epoch,pfix=pfix,norbits=norbits)
            m.print_soln()  
            mpsini = mass_estimate(m, mstar, norbits=norbits)
            print 'mp*sin(i):         ',str(mpsini)
            outp[i] = m.params[6]
            outchisq[i] = m.rchisq
        return m, pguess, outp, outchisq

    else:
        m, flt = rvfit_lsqmdl(guesspars, tnorm, rvnorm, nsrv, jitter=jitter,circ=circ, npoly=npoly,tt=transit,epoch=epoch,pfix=pfix,norbits=norbits,telvec=telvec)

    
        #display initial fit - want to show fixed params too
        m.print_soln()  
    
        #calc bayesian information criterion
        bic0 = calc_bic(m,flt,srv)
        print 'BIC:          ',str(bic0)

        #mass estimate
        mpsini = mass_estimate(m, mstar, norbits=norbits)
        print 'mp*sin(i):         ',str(mpsini)

        #density estimate
        dpl = density_estimate(mpsini,rpl)
        print 'density:           ',str(dpl),'g/cc'

        #correct offset before plotting
        tels = np.unique(telvec)
        #rvp = np.zeros_like(rvnorm)
        rvp = rvnorm
        par0 = np.copy(m.params)
        for i in range(ntel-1):
            a = np.squeeze(np.where(telvec == tels[i+1]))
            print 'offset: ', m.params[i+norbits*6+npoly]
            rvp[a] -= m.params[i+norbits*6+npoly]
            m.params[i+norbits*6+npoly] = 0


        #make plots
        plot_rv(targname,tnorm,rvnorm,nsrv,guesspars,m,nmod=200,home=home,norbits=norbits,npoly=npoly,telvec=telvec)
        m.params = par0


        #call MCMC    
        if nwalkers > 0:
            mcbest, bestpars, pnames, flt, samples, mcpars, chain = setup_emcee(targname, m, tnorm, rvnorm, nsrv, circ=circ, npoly=npoly, tt=transit, jitter=jitter, nwalkers=nwalkers, pfix=pfix, telvec=telvec, norbits=norbits, nsteps=nsteps)
            plot_rv(targname,tnorm,rvnorm,nsrv,mcbest,m,nmod=200,home=home,norbits=norbits,npoly=npoly,telvec=telvec,mc=1)
            mcmod = rv_drive(mcbest,tnorm,norbits,npoly,telvec)
            mcres = mcmod - rvnorm
            print m.resids
            print mcres
            mpsini, mparr_mc = mass_estimate(m, mstar, norbits=norbits, mcpar=mcpars)
            #dpl = density_estimate(mpsini,rpl, mcpar=mcpars, rple=rple)

            #print output from mass_estimate for mc
            print_mc_errs(mcpars, mpsini, mparr_mc,norbits=norbits,npoly=npoly)
            bic = calc_bic(m, flt, srv, mcresid=mcres)

            #print 'MC BIC:          ',str(bic)
            #make a nice triangle plot
            print pnames
            f = np.squeeze(flt.nonzero())
            fig = triangle.corner(samples, labels=pnames[f], truths=bestpars[f]) 
            fig.savefig('/home/sgettel/Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+'_triangle.png')
            plt.close(fig)
        else:
            mcpars = -1
            mparr_mc = -1


        
        write_full_soln(m, targname, mpsini, bic0, mcpars=mcpars, mparr_mc=mparr_mc, norbits=norbits, npoly=npoly, telvec=telvec)

    return m#, jdb0, rv0, srv0, fwhm, contrast, bis_span, rhk, sig_rhk, chain

def calc_bic(m, flt, srv, mcresid=-1):
    
    #calc BIC for LM fit, use original errors
    chisq = np.sum(m.resids**2/srv**2)
    
    k = np.sum(flt)
    n = srv.size
    print chisq, k, n
    bic = chisq + k*np.log(n) 

    if len(np.array(mcresid).shape) > 0:
      chisq = np.sum(mcresid**2/srv**2)
      print chisq
      bic = chisq + k*np.log(n)

    return bic
    

def plot_rv(targname,jdb,rv,srv,guesspars,m,nmod=1000,home='/home/sgettel/', norbits=1,npoly=0,telvec=-1,mc=0):

    #save uncorrected RVs, if needed
    if len(np.array(telvec).shape) > 0:
        ntel = np.unique(telvec).size

   
    tmod = np.linspace(np.min(jdb),np.max(jdb),nmod)

    
    model_final = rv_drive(m.params,tmod,norbits,npoly,telvec)
    
    if npoly > 0:
        parst =  np.copy(m.params)
        parst[4] = 0.0
        poly = rv_drive(parst,tmod,norbits,npoly,telvec)

    #unphased data
    plt.figure(1)
    plt.errorbar(jdb,rv,yerr=srv,fmt='bo')
    plt.plot(tmod,model_final,'r-')
    if npoly > 0:
        plt.plot(tmod,poly,'g-')
    
    if mc > 0:
       plt.savefig(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+'_autoplot_mc.png') 
    else:
        plt.savefig(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+'_autoplot.png')
    plt.close(1)

    #phase at 1st period
    pars = m.params[0:6]  #for model 
    pars[6:] = 0 #select first planet only
        
    parst = np.copy(m.params)
    parst[4] = 0.0 #other planets only
    parst[5] = 0.0
    phase = (jdb - pars[1])/pars[0] % 1.0

    rvt = rv_drive(parst,jdb,norbits,npoly,telvec)   
    telvec = np.zeros_like(tmod)

    plt.figure(2)
    plt.errorbar(phase, rv-rvt, yerr=srv,fmt='bo')
    plt.plot((tmod - pars[1])/pars[0] % 1.0, rv_drive(pars, tmod,1,0,telvec),'r.')
    #plt.plot((tmod - guess))
    plt.savefig(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+'_phase_autoplot.png')
    plt.close(2)


def write_full_soln(m,targname,mpsini, bic, mparr_all=-1, mcpars=-1, mparr_mc=-1,norbits=1,npoly=0,telvec=-1):
    
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
    
    for i in range(npoly):
        f.write(str(poly_names[i])+str(m.params[i+norbits*6]) +'\n')
    f.write('r chi sq: '+str(m.rchisq)+'\n')
    f.write('BIC: '+str(bic)+'\n')
    
    for i in range(norbits):
        f.write('mp*sin(i): '+str(mpsini[i])+'\n')
    
    if len(np.array(telvec).shape) > 0:
        ntel = np.unique(telvec).size

    
    if len(np.array(mcpars).shape) > 0: #print MCMC errs
        for i in range(norbits):
            f.write('                                               \n')
            f.write('*****Planet '+str(i+1)+' MCMC Errors:***** \n')    
            f.write('Per: '+ str(np.mean(mcpars[:,0+i*6]))+'+/-'+str(np.std(mcpars[:,0+i*6]))+'\n')
            f.write('Tp: '+ str(np.mean(mcpars[:,1+i*6]))+'+/-'+str(np.std(mcpars[:,1+i*6]))+'\n')
            f.write('ecc: '+ str(np.mean(mcpars[:,2+i*6]))+'+/-'+str(np.std(mcpars[:,2+i*6]))+'\n')
            f.write('om: '+ str(np.mean(mcpars[:,3+i*6]))+'+/-'+str(np.std(mcpars[:,3+i*6]))+'\n')
            f.write('K1: '+ str(np.mean(mcpars[:,4+i*6]))+'+/-'+str(np.std(mcpars[:,4+i*6]))+'\n')
            f.write('gamma: '+ str(np.mean(mcpars[:,5+i*6]))+'+/-'+str(np.std(mcpars[:,5+i*6]))+'\n')
            
            

        for i in range(npoly):
            f.write(str(poly_names[i])+str(np.mean(mcpars[:,i+norbits*6]))+'+/-'+str(np.std(mcpars[:,i+norbits*6])) +'\n')
    
        for i in range(norbits):
            f.write('mp*sin(i): '+str(np.mean(mparr_mc[i,:]))+'+/-'+str(np.std(mparr_mc[i,:]))+'\n')
            f.write('mass error:'+ str(np.std(mparr_mc[i,:])/mpsini[i]*100)+'%'+'\n') 
    f.close() 
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

def density_estimate(mpsini,rpl,mcpar=-1,rple=-1):
    
    rpl = ut.arrayify(rpl)
    nplanets = rpl.size

    re2cm = 1./6.378e8
    vol = 4./3.*np.pi*(rpl/re2cm)**3 #cm**3
    dpl = mpsini[0:nplanets]*5.973e27/vol

    return dpl


   

def rvfit_lsqmdl(orbel,jdb,rv,srv,jitter=0, param_names=0,npoly=0,circ=0, tt=np.zeros(1),epoch=2.455e6,pfix=1,norbits=1,telvec=-1):

    ip = np.arange(norbits)

    if len(np.array(telvec).shape) > 0:
        ntel = np.unique(telvec).size
    
    if param_names == 0: 
        param_names = ['Per', 'Tp', 'ecc', 'om', 'K1', 'gamma']*norbits
        poly_names = ['dvdt','quad','cubic','quart']
        param_names.extend(poly_names[:npoly]) 
        if ntel > 1:
            off_names = ['offset']*(ntel-1)
            param_names.extend(off_names)

    flt = np.ones_like(orbel) #Free param? Turn off individually

    m = lsqmdl.Model(None, rv, 1./srv) 
    m.set_func(rv_drive,param_names, args=(jdb,norbits,npoly,telvec) )

    #print pfix, ' = pfix'
    #make some reasonable limits - don't need the loop?
    for i in range(norbits):

        #apply normal limits first
        #fix/limit period
        if pfix[i] == 1:
            m.lm_prob.p_value(0+i*6, orbel[0+i*6], fixed=True)
            flt[0+i*6] = 0
        else:
            m.lm_prob.p_limit(0+i*6, lower=0.1, upper=(np.max(jdb)-np.min(jdb))*20) #per no longer than 10x range of survey 
        
        #limit Tp
        m.lm_prob.p_limit(1+i*6, lower=orbel[1+i*6]-orbel[0+i*6]/2., upper=orbel[1+i*6]+orbel[0+i*6]/2.) #T0 within one period guess

        #fix/limit ecc
        if circ[i] == 1:
            m.lm_prob.p_value(2+i*6, 0.0, fixed=True) 
            flt[2+i*6] = 0
        else:
            m.lm_prob.p_limit(2+i*6, lower=0.0, upper=0.99) #ecc must be physical 
        #fix/limit omega
        if circ[i] == 1:
            m.lm_prob.p_value(3+i*6, 0.0, fixed=True) #convention
            flt[3+i*6] = 0
        else:
            m.lm_prob.p_limit(3+i*6, lower=0.0, upper=360.0)

        #limit K
        m.lm_prob.p_limit(4+i*6, lower=0.0, upper=1.0e5) #K must be physical
  
        #fix gamma except first
        if i > 0:
            m.lm_prob.p_value(5+i*6, 0.0, fixed=True)
            flt[5+i*6] = 0

        #now include known transit effects
        if not tt[i] == 0:
            if circ[i] == 1:  #by convention tt=tp & omega=90
                m.lm_prob.p_value(1+i*6, tt[i], fixed=True)
                m.lm_prob.p_value(3+i*6, 90.0, fixed=True)
                flt[1+i*6] = 0
                flt[3+i*6] = 0
            else:
                tiefunc = tie_omega_function(tt, i) #how does this know what orbel is?
                m.lm_prob.p_tie(3+i*6, tiefunc)
                flt[3+i*6] = 0

    #limit polynomial terms
    for i in range(npoly):
        
        m.lm_prob.p_limit(i+norbits*6, lower=-1e6, upper=1e6) #dvdt and higher


    #limit offset terms
    for i in range(ntel-1):
        m.lm_prob.p_limit(i + norbits*6 + npoly, lower=-1e6, upper=1e6)

    

    m.solve(orbel)
   
    return m, flt


def tie_omega_function(tt, i):

    def calculate_omega(orbel):
        
        p = orbel[0+i*6]
        tp = orbel[1+i*6]
        ecc = orbel[2+i*6]

        theta_tt = calc_true_anomaly(p, tp, ecc, tt[i]) #in radians
        omega = (np.pi/2.0 - theta_tt)*180.0/np.pi
         
        if omega < 0:    
            omega += 360.0
        if omega > 360.0:
            omega -= 360.0
        
        return omega

    return calculate_omega

def rv_drive(orbel, t, norbits, npoly, telvec):
    #From rv_drive.pro in RVLIN by JTW

   
    if len(np.array(telvec).shape) > 0:
        ntel = np.unique(telvec).size

    rv = np.zeros_like(t)
   
    phase = np.zeros((rv.size,norbits))
    
    for i in range(norbits):  
        p = orbel[0+i*6]
        tp = orbel[1+i*6]
        ecc = orbel[2+i*6]
        om = orbel[3+i*6] *np.pi/180. #degrees to radians
        k = orbel[4+i*6]
        gamma = orbel[5+i*6]




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
            
        rv = rv + k*(np.cos(theta + om) + ecc*np.cos(om)) + gamma 

    #now add polynomial
    for i in range(npoly):
        rv = rv + orbel[i+norbits*6]*(t - epoch)**(i+1)
   
    #now add offsets
    tels = np.unique(telvec)
    for i in range(ntel-1):
        #print tels[i]
        a = np.squeeze(np.where(telvec == tels[i+1]))
        #print a.size
        rv[a] += orbel[i+norbits*6+npoly]

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
    

    nm = marr.size    #nm = nrv, ~100
    nec = ecc.size #nec = 1
 
    conv = 1e-12 #convergence criteria
    k = 0.85     #scale factor
    
    mphase = marr/(2.0*np.pi)
    #print mphase.size, mphase.shape

    
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


#to run as ...orbits.py
if __name__ == '__main__':
    orbits_test(webdat='yes',nboot=10)


#begin MCMC setup - following emcee line-fitting demo
def lnprior(theta, fullpars, flt, pnames, plo, phi):
    
    print fullpars

    #figure out which params are varied and the associated limits
    pfloat = pnames[flt.nonzero()]
    lfloat = plo[flt.nonzero()]
    hfloat = phi[flt.nonzero()]

    #flat priors for all
    if (theta >= lfloat).all() and (theta < hfloat).all():
        lnpri = 0.0

        #pper = norm.pdf(loc=)


        return lnpri 
    else:
        
        return -np.inf
    

def lnprob(theta, jdb, rv, srv, fullpars, flt, pnames, plo, phi, norbit, npoly, telvec):
    lp = lnprior(theta, fullpars, flt, pnames, plo, phi)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, jdb, rv, srv, fullpars, flt, norbit, npoly, telvec)

def lnlike(theta, jdb, rv, srv, fullpars, flt, norbit, npoly, telvec):
    
    newpars = np.copy(fullpars)
    
    newpars[flt.nonzero()] = theta
    
    #jitter = newpars[-1]

    model = rv_drive(newpars, jdb, norbit, npoly, telvec)

    stot = np.sqrt(srv**2 + newpars[-1]**2) #add floating jitter term

    l0 = np.sum(np.log(1/(np.sqrt(2*np.pi)*stot))) #penalize high jitter values
    chisq = -0.5*np.sum((rv - model)**2/stot**2)

    return l0 + chisq 

def lnlike_base(bestpars, jdb, rv, srv, norbit, npoly, telvec):
    
    model = rv_drive(bestpars, jdb, norbit, npoly, telvec)
    stot = np.sqrt(srv**2 + bestpars[-1]**2) #add floating jitter term

    l0 = np.sum(np.log(1/(np.sqrt(2*np.pi)*stot))) #penalize high jitter values
    chisq = -0.5*np.sum((rv - model)**2/stot**2)

    return l0 + chisq 

def setup_emcee(targname, m, jdb, rv, srv_in, nwalkers=200, circ=0, npoly=0, norbits=1, tt=np.zeros(1),jitter=0, pfix=1,nburn=200,telvec=-1,nsteps=1000): 

    bestpars = np.copy(m.params)
    bestpars = np.append(bestpars,0.5) #add placeholder for jitter
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

    if len(np.array(telvec).shape) > 0:
        ntel = np.unique(telvec).size
        tels = np.unique(telvec)
    else:
        ntel = 1
        

    #these will hold reasonable limits on priors
    plo = np.zeros(npars)
    phi = np.zeros(npars)

    #separate the params being varied from the full list
    flt = np.ones(npars) #all params float now, turn off individually  

    
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
        if circ[i] == 1:
            flt[3+i*6] = 0
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

    #need to limit offset terms if present!
    for i in range(ntel-1):
        
        plo[i+norbits*6+npoly] = -1e6
        phi[i+norbits*6+npoly] = 1e6

    #limit jitter
    plo[-1] = 0.001
    phi[-1] = 5.0


   
       
    
    f = np.squeeze(flt.nonzero())
    
    varpars = bestpars[f]
    ndim = varpars.size
    
    print 'MCMC params: ',pnames[f] 
    print varpars

    chain = run_emcee(targname, bestpars, varpars, flt, plo, phi, jdb, rv, srv, pnames, ndim, nwalkers=nwalkers, norbits=norbits, npoly=npoly, telvec=telvec,nsteps=nsteps)
  
    #It takes a number of iterations to spread walkers throughout param space
    #This is 'burning in'
    samples = chain[:, nburn:, :].reshape((-1, ndim))
    
    #combine samples with best-fit values
    mcpars = ut.fan(bestpars,samples.shape[0])
    
    for i in range(ndim):
        mcpars[:,f[i]] = samples[:,i]

    mcbest = np.percentile(mcpars,50, axis=0)

    #return m, flt, chain, samples, mcpars
    return mcbest, bestpars, pnames, flt, samples, mcpars, chain


def run_emcee(targname, bestpars, varpars, flt, plo, phi, jdb, rv, srv, pnames, ndim, nwalkers=200, nsteps=1000, norbits=1, npoly=0, telvec=-1):
    
    #Initialize walkers in tiny Gaussian ball around MLE results
    #number of params comes from varpars
    #pos = [varpars + 1e-3*np.random.randn(ndim) for i in range(nwalkers)] #??
    pos = [varpars + 1e-3*varpars*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(jdb,rv,srv,bestpars,flt, pnames, plo, phi, norbits, npoly, telvec))

    #Run MCMC
    sampler.run_mcmc(pos, nsteps)
    print("Mean acceptance fraction: {0:.3f}"
                .format(np.mean(sampler.acceptance_fraction)))

    return sampler.chain

