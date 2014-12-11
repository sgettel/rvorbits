#RV orbit modeling code based on RVLIN by Jason Wright (IDL) and orbits.f by Alex Wolszczan (FORTRAN77)

#
# Branch jit2
# PUT ALL TRANSITING PLANETS FIRST!
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

def orbits_test(targname='K00273',jitter=0.5,epoch=2.4568e6,circ=0,maxrv=1e6,minrv=-1e6,maxsrv=5, webdat='no', nwalkers=200, pfix=1,norbits=1,npoly=0,keck='no',outer_loop='no',nsteps=1000,nburn=300,ttfloat='no',fixjit='no',storeflat='no'):

    tag = ''
    if npoly > 4:
        print 'Must have <= 4th order polynomial'


    transit = np.zeros(1)
    telvec = np.zeros(1)
    circ = ut.arrayify(circ)
    pfix = ut.arrayify(pfix)
    jitter = ut.arrayify(jitter)

    for i in circ:
        if i == 1:
            tag += '_circ'
        else:
            tag += '_ecc'

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
            
            

#            sfile = open(home+'Dropbox/cfasgettel/research/keck/'+targname+'.dat')
#            kjdb, krv, ksrv = np.loadtxt(sfile,unpack=True,usecols=(2,3,4))
#            kjdb = kjdb + 2.45e6 

            sfile = open(home+'Dropbox/cfasgettel/research/keck/'+targname+'_full.dat')
            kjdb, krv, ksrv = np.loadtxt(sfile,unpack=True,usecols=(2,3,4))
            kjdb = kjdb + 2.44e6


            krvnorm = krv - np.median(krv)
            ktel = np.ones_like(kjdb)
                
            jdb = np.append(jdb,kjdb)
            rvnorm = rv - np.median(rv)
            rvnorm = np.append(rvnorm,krvnorm)
            srv = np.append(srv,ksrv)
            telvec = np.append(telvec,ktel)
            print kjdb.size,' Keck obs'
        else:
            rvnorm = rv - np.median(rv)

    #adjust values to be sensible
    #if jdb[0] > epoch:
    #    print 'truncating dates'
    tnorm = jdb - epoch #truncate
  
#    if jitter > 0.0: 
#        nsrv = np.sqrt(srv**2 + jitter**2)
#        print 'Adding ',str(jitter),' m/s fixed jitter'
#    else:
    nsrv = srv



    #process offsets  
    ntel = np.unique(telvec).size
    if ntel > 1:
        print 'Including ',ntel-1,' offset terms'
        
        offset = np.zeros(ntel-1)
    

    ip = np.arange(norbits)
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

    #jitter - final term/s, MCMC only
    rpl = 0.0

    if targname == 'HD209458':
        guesspars = np.array([3.524733, 2452826.628514, 0.0, 336.5415, 85.49157+10, -1.49-10])#HD209458
        transit = np.array([2452826.628514]) 
        mstar = 1.0
    
    if targname == 'K00273':

##################### UPDATE



        #guesspars = np.array([10.573769, 2455008.06601, 0.0, 90.0, 1.7358979, -3398.0498, 1.1889011, 0.010, 0.0])#, 0.0]) #K00273
        guesspars = np.array([10.573737, 2455008.06787, 0.0, 90.0, 1.7358979, -3398.0498, 1.1889011, 0.010, 0.0])#, 0.0])
        #transit = np.array([2455008.06601,0.0])
        transit = np.array([2455008.06787,0.0]) 
        psig = np.array([6.1e-6,0.0])
        ttsig = np.array([0.00061,0.0])
        porig = guesspars[0+ip*6] 
        ttorig = np.copy(transit)

        mstar = 1.069
        rpl = 1.854 #Me
        rple = 0.041
        inc = 89.9771
        ince = 0.30

        #2 planets...       
        if norbits > 1:
       # #guesspars = np.array([10.573769, 2455008.06601, 0.0, 90.0, 2.2, -16.0,  500.0, 2455041.9, 0.23, 40.0, 137.0,0.0])
            guesspars = np.array([10.573737, 2455008.06787, 0.0, 90.0, 2.2, -16.0,  500.0, 2455041.9, 0.23, 340.0, 137.0,0.0,0.0])


    
    if targname == 'K00069':
        guesspars = np.array([4.72673978, 2454944.29227, 0.0, 90.0, 1.733, -91.08, 0.0329])
        transit = np.array([2454944.29227])
        mstar = 0.887

    
    ip = np.arange(norbits)

    guesspars[1+ip*6] -= epoch

    #append offsets if needed
    if ntel > 1:
        guesspars = np.append(guesspars,offset)

    #append transit time if needed
    if ttfloat == 'yes':
        print 'Letting transit time/s float'
    for i in range(norbits):  
        if not transit[i] == 0.0:
            transit[i] -= epoch 

        #add all by default, turn them off in LM fitter
        if ttfloat == 'yes':
            guesspars = np.append(guesspars,transit[i]) 


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
        m, flt = rvfit_lsqmdl(guesspars, tnorm, rvnorm, nsrv, jitter=jitter,circ=circ, npoly=npoly,tt=transit,epoch=epoch,pfix=pfix,norbits=norbits,telvec=telvec,psig=psig,ttfloat=ttfloat,ttsig=ttsig)

    
        #display initial fit - want to show fixed params too
        m.print_soln()  
    
        #calc bayesian information criterion
        bic0 = calc_bic_lm(m,flt,srv)
        print 'BIC:          ',str(bic0)

        #mass estimate
        mpsini = mass_estimate(m, mstar, norbits=norbits)
        print 'mp*sin(i):         ',str(mpsini)

        #density estimate
        dpl = density_estimate(mpsini,rpl)
        print 'density:           ',str(dpl),'g/cc'

        #correct offset before plotting
        tels = np.unique(telvec)
        rvp = np.copy(rvnorm)
        
        par0 = np.copy(m.params)
        for i in range(ntel-1):
            a = np.squeeze(np.where(telvec == tels[i+1]))
            print 'offset: ', m.params[i+norbits*6+npoly]
            rvp[a] -= m.params[i+norbits*6+npoly]
            m.params[i+norbits*6+npoly] = 0
        
        #make plots & restore offset
        res0 = plot_rv(targname,tnorm,rvp,nsrv,guesspars,m,nmod=200,home=home,norbits=norbits,npoly=npoly,telvec=telvec,keck=keck)
        m.params = par0
       
        #call MCMC    
        if nwalkers > 0:
            mcbest, bestpars, pnames, flt, samples, mcpars, chain = setup_emcee(targname, m, tnorm, rvnorm, nsrv, circ=circ, npoly=npoly, tt=transit, jitter=jitter, nwalkers=nwalkers, pfix=pfix, telvec=telvec, norbits=norbits, nsteps=nsteps,psig=psig,porig=porig,nburn=nburn,ttfloat=ttfloat,ttsig=ttsig,fixjit=fixjit)

            #correct offset before plotting
            rvp = np.copy(rvnorm)
            print 'mcbest:', mcbest
            par1 = np.copy(mcbest)
            for i in range(ntel-1):
                a = np.squeeze(np.where(telvec == tels[i+1]))
                print 'MC offset: ', mcbest[i+norbits*6+npoly]
                rvp[a] -= mcbest[i+norbits*6+npoly]
                mcbest[i+norbits*6+npoly] = 0
            
            pres = plot_rv(targname,tnorm,rvp,nsrv,mcbest,m,nmod=200,home=home,norbits=norbits,npoly=npoly,telvec=telvec,mc=1,keck=keck)
            if storeflat == 'yes':
                f = open('/home/sgettel/Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+'_rvflat.dat','w')
                for i in range(tnorm.size):
                    f.write('#tnorm    res    srv    tel            \n') 

                    f.write(str(tnorm[i])+' '+str(pres[i])+' '+str(nsrv[i])+' '+str(telvec[i])+'\n')
                f.close()
            mcbest = par1
            
            mpsini, mparr_mc = mass_estimate(m, mstar, norbits=norbits, mcpar=mcpars)
            dpl = density_estimate(mpsini,rpl, mcmass=mparr_mc, rple=rple)

            #print output from mass_estimate for mc
            print_mc_errs(mcpars, mpsini, mparr_mc,norbits=norbits,npoly=npoly,telvec=telvec,ttfloat=ttfloat,ttsig=ttsig,mcdpl=dpl)
            bic = calc_bic_mc(mcbest, flt, tnorm, rvnorm, nsrv, norbits, npoly, telvec, transit, ttsig, ttfloat,circ)

            print 'MC BIC:          ',str(bic)

            #convergence testing...
            psrf = gelman_rubin(chain)

            #make a nice triangle plot
            f = np.squeeze(flt.nonzero())
            fig = triangle.corner(samples, labels=pnames[f], truths=bestpars[f]) 
            fig.savefig('/home/sgettel/Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+'_triangle.pdf')
            plt.close(fig)

            

        else:
            mcpars = -1
            mparr_mc = -1
            chain = -1

            
        
        write_full_soln(m, targname, mpsini, bic0, mcpars=mcpars, mparr_mc=mparr_mc, norbits=norbits, npoly=npoly, telvec=telvec, ttfloat=ttfloat, ttsig=ttsig,mbic=bic,psrf=psrf)

    return m,chain, mparr_mc

def calc_bic_lm(m, flt, srv):
    
    #calc BIC for LM fit, use original errors
    chisq = np.sum(m.resids**2/srv**2)
    
    k = np.sum(flt)
    n = srv.size
    print chisq, k, n
    bic = chisq + k*np.log(n) 

    return bic
    
def calc_bic_mc(mcbest,flt,jdb,rv,srv,norbit,npoly,telvec, tt, ttsig, ttfloat,circ):
  
    k = np.sum(flt)
    n = srv.size
  
    lnl = lnlike_base(mcbest, jdb, rv, srv, norbit, npoly, telvec, tt, ttsig, ttfloat,circ)
    bic = -2*lnl + k*np.log(n)

    return bic
     

def plot_rv(targname,jdb,rv,srv,gpars,m,nmod=1000,home='/home/sgettel/', norbits=1,npoly=0,telvec=-1,mc=0,keck='no'):

    if mc > 0:
        usepars = np.copy(gpars) #this is mcbest...
    else:
        usepars = np.copy(m.params)

    #print 'plot pars:', usepars
    tmod = np.linspace(np.min(jdb),np.max(jdb),nmod)

    model_final = rv_drive(usepars,tmod,norbits,npoly,telvec)
    
    if npoly > 0:
        parst =  np.copy(usepars) 
        parst[4] = 0.0
        poly = rv_drive(parst,tmod,norbits,npoly,telvec)

    k = np.squeeze(np.where(telvec == 1))

    #unphased data
    plt.figure(1)
    plt.errorbar(jdb,rv,yerr=srv,fmt='bo')
    if keck == 'yes':
        plt.errorbar(jdb[k],rv[k],yerr=srv[k],fmt='go') 
    plt.plot(tmod,model_final,'r-')
    plt.xlabel('Adjusted BJD')
    plt.ylabel('Normalized RV (m/s)')
    if npoly > 0:
        plt.plot(tmod,poly,'g-')
    
    if mc > 0:
       plt.savefig(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+'_autoplot_mc.pdf') 
    else:
        plt.savefig(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+'_autoplot.pdf')
    plt.close(1)

    #phase at 1st period
    pars = usepars[0:6]  #for model 
    pars[6:] = 0 #select first planet only
        
    parst = np.copy(usepars)
    parst[4] = 0.0 #other planets only
    parst[5] = 0.0
    phase = (jdb - pars[1])/pars[0] % 1.0

    rvt = rv_drive(parst,jdb,norbits,npoly,telvec)   
    pres = (rv-rvt) - pars[5] 
    telvec = np.zeros_like(tmod)

    plt.figure(2)
    plt.errorbar(phase, rv-rvt, yerr=srv,fmt='bo')
    if keck == 'yes':
        plt.errorbar(phase[k],rv[k]-rvt[k],yerr=srv[k],fmt='go') 
    plt.plot((tmod - pars[1])/pars[0] % 1.0, rv_drive(pars, tmod,1,0,telvec),'r.')
    plt.xlabel('Orbital Phase')
    plt.ylabel('Normalized RV (m/s)')
    #plt.plot((tmod - guess))
    if mc > 0:
        plt.savefig(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+'_phase_autoplot_mc.pdf')
    else:
        plt.savefig(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+'_phase_autoplot.pdf') 
    plt.close(2)
    
    return pres

def write_full_soln(m,targname,mpsini, bic, mcpars=-1, mparr_mc=-1,norbits=1,npoly=0,telvec=-1,tt=np.zeros(1),ttsig=-1,ttfloat='no',mbic=-1,psrf=-1):
    
    poly_names = ['dvdt:  ','quad:  ', 'cubic: ','quart: ']

    if len(np.array(telvec).shape) > 0:
        tels = np.unique(telvec)
        ntel = np.unique(telvec).size
    
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
    
        if ttfloat == 'yes' and ttsig[i] > 0:
             tind = i + norbits*6 + npoly + ntel-1
             f.write('TT: '+str(m.params[tind])+'\n')
             
        f.write('mp*sin(i): '+str(mpsini[i])+'\n')


    for i in range(npoly):
        f.write(str(poly_names[i])+str(m.params[i+norbits*6]) +'\n')
    f.write('r chi sq: '+str(m.rchisq)+'\n')
    f.write('BIC: '+str(bic)+'\n')
    
    for i in range(ntel-1):
        a = np.squeeze(np.where(telvec == tels[i+1]))
        f.write('offset: '+str(m.params[i+norbits*6+npoly])+'\n')
    
    if len(np.array(mcpars).shape) > 0: #print MCMC errs
        mcbest = np.percentile(mcpars,50, axis=0)
        mchi = np.percentile(mcpars,84, axis=0)
        mclo = np.percentile(mcpars,16, axis=0)

        for i in range(norbits):
            f.write('                                               \n')
            f.write('*****Planet '+str(i+1)+' MCMC Errors:***** \n')    
            f.write('Per: '+str(mcbest[0+i*6])+' +'+str(mchi[0+i*6] - mcbest[0+i*6])+' -'+str(mcbest[0+i*6] - mclo[0+i*6])+'\n')
            f.write('Tp: '+ str(mcbest[1+i*6])+' +'+str(mchi[1+i*6] - mcbest[1+i*6])+' -'+str(mcbest[1+i*6] - mclo[1+i*6])+'\n')
            f.write('ecc: '+ str(mcbest[2+i*6])+' +'+str(mchi[2+i*6] - mcbest[2+i*6])+' -'+str(mcbest[2+i*6] - mclo[2+i*6])+'\n')
            f.write('ecc1: '+ str(np.percentile(mcpars,0,axis=0)[2+i*6])+'-'+str(np.percentile(mcpars,68,axis=0)[2+i*6])+'\n')
            f.write('om: '+str(mcbest[3+i*6])+' +'+str(mchi[3+i*6] - mcbest[3+i*6])+' -'+str(mcbest[3+i*6] - mclo[3+i*6]) +'\n')
            f.write('K1: '+str(mcbest[4+i*6])+' +'+str(mchi[4+i*6] - mcbest[4+i*6])+' -'+str(mcbest[4+i*6] - mclo[4+i*6]) +'\n')
            f.write('gamma: '+ str(mcbest[5+i*6])+' +'+str(mchi[5+i*6] - mcbest[5+i*6])+' -'+str(mcbest[5+i*6] - mclo[5+i*6])+'\n')
            
            if ttfloat == 'yes' and ttsig[i] > 0:
                tind = i + norbits*6 + npoly + ntel-1
                f.write('TT: '+ str(mcbest[tind])+' +'+str(mchi[tind] - mcbest[tind])+' -'+str(mcbest[tind] - mclo[tind])+'\n')

            mpbest = np.percentile(mparr_mc[i,:], 50)
            mphi = np.percentile(mparr_mc[i,:], 84)
            mplo = np.percentile(mparr_mc[i,:], 16)
            f.write('nsamples: '+str(mparr_mc.shape[1])+'\n')
            f.write('convergence: '+str(psrf)+'\n')

            f.write('mp*sin(i): '+str(mpbest)+' +'+str(mphi-mpbest)+' -'+str(mpbest-mplo)+'\n')
            f.write('mass error:'+ str((mphi-mpbest)/mpbest*100)+','+str((mpbest-mplo)/mpbest*100)+ '%'+'\n') 


        for i in range(npoly):
            f.write(str(poly_names[i])+ str(mcbest[i+norbits*6])+' +'+str(mchi[i+norbits*6]-mcbest[i+norbits*6])+' -'+str(mcbest[i+norbits*6]-mclo[i+norbits*6]) +'\n')
        for i in range(ntel-1):
            a = np.squeeze(np.where(telvec == tels[i+1]))
            f.write('offset: '+str(mcbest[i+norbits*6+npoly])+' +'+str(mchi[i+norbits*6+npoly]-mcbest[i+norbits*6+npoly])+' -'+str(mcbest[i+norbits*6+npoly]-mclo[i+norbits*6+npoly])+'\n')

        for i in range(ntel):
            f.write('jitter: '+ str(mcbest[-ntel+i])+' +'+str(mchi[-ntel+i]-mcbest[-ntel+i])+' -'+str(mcbest[-ntel+i]-mclo[-ntel+i])+'\n')

        f.write('MC BIC: '+str(mbic)+'\n')
    f.close() 

def print_mc_errs(mcpars, mpsini, mparr_all,norbits=1,npoly=0,telvec=-1,tt=np.zeros(1),ttsig=-1,ttfloat='no',mcdpl=-1):
    poly_names = ['dvdt:  ','quad:  ', 'cubic: ','quart: ']

    if len(np.array(telvec).shape) > 0:
        tels = np.unique(telvec)
        ntel = np.unique(telvec).size

    mcbest = np.percentile(mcpars,50, axis=0)
    mchi = np.percentile(mcpars,84, axis=0)
    mclo = np.percentile(mcpars,16, axis=0) 

    for i in range(norbits):
        print '                                               '
        print '*****Planet ',str(i+1),' MCMC Errors:*****'  
        print 'Per: ', str(mcbest[0+i*6]),' +',str(mchi[0+i*6] - mcbest[0+i*6]),' -',str(mcbest[0+i*6] - mclo[0+i*6])
        print 'Tp: ', str(mcbest[1+i*6]),' +',str(mchi[1+i*6] - mcbest[1+i*6]),' -',str(mcbest[1+i*6] - mclo[1+i*6])
        print 'Ecc: ', str(mcbest[2+i*6]),' +',str(mchi[2+i*6] - mcbest[2+i*6]),' -',str(mcbest[2+i*6] - mclo[2+i*6])
        #f.write('ecc1: '+ str(np.percentile(mcpars,0,axis=0)[2+i*6])+'-'+str(np.percentile(mcpars,68,axis=0)[2+i*6])+'\n')
        print 'Om: ', str(mcbest[3+i*6]),' +',str(mchi[3+i*6] - mcbest[3+i*6]),' -',str(mcbest[3+i*6] - mclo[3+i*6])
        print 'K: ', str(mcbest[4+i*6]),' +',str(mchi[4+i*6] - mcbest[4+i*6]),' -',str(mcbest[4+i*6] - mclo[4+i*6])
        print 'gamma: ', str(mcbest[5+i*6]),' +',str(mchi[5+i*6] - mcbest[5+i*6]),' -',str(mcbest[5+i*6] - mclo[5+i*6])

        if ttfloat == 'yes' and ttsig[i] > 0:
            tind = i + norbits*6 + npoly + ntel-1
            print 'TT: ',str(mcbest[tind]),' +',str(mchi[tind] - mcbest[tind]),' -',str(mcbest[tind] - mclo[tind])

        mpbest = np.percentile(mparr_all[i,:], 50)
        mphi = np.percentile(mparr_all[i,:], 84)
        mplo = np.percentile(mparr_all[i,:], 16)

        print 'mp*sin(i): ',str(mpbest),' +',str(mphi-mpbest),' -',str(mpbest-mplo)
        print 'mass error:', str((mphi-mpbest)/mpbest*100),',',str((mpbest-mplo)/mpbest*100), '%'
        if i == 0:
            dpbest = np.percentile(mcdpl[:], 50)
            dphi = np.percentile(mcdpl[:], 84)
            dplo = np.percentile(mcdpl[:], 16)
            print 'density: ',str(dpbest),' +',str(dphi-dpbest),' -',str(dpbest-dplo)
    
    for i in range(npoly):
        print str(poly_names[i]), str(mcbest[i+norbits*6]),' +',str(mchi[i+norbits*6]-mcbest[i+norbits*6]),' -',str(mcbest[i+norbits*6]-mclo[i+norbits*6])

    
    for i in range(ntel-1):
        #a = np.squeeze(np.where(telvec == tels[i+1]))
        print 'offset: ',str(mcbest[i+norbits*6+npoly]),' +',str(mchi[i+norbits*6+npoly]-mcbest[i+norbits*6+npoly]),' -',str(mcbest[i+norbits*6+npoly]-mclo[i+norbits*6+npoly])
    
    for i in range(ntel):
        print 'jitter: ', str(mcbest[-ntel+i]),' +',str(mchi[-ntel+i]-mcbest[-ntel+i]),' -',str(mcbest[-ntel+i]-mclo[-ntel+i])
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

def density_estimate(mpsini,rpl,mcmass=-1,rple=-1):
    
    rpl = ut.arrayify(rpl)
    nplanets = rpl.size

    re2cm = 1./6.378e8
    vol = 4./3.*np.pi*(rpl/re2cm)**3 #cm**3
    #dpl = mpsini[0:nplanets]*5.973e27/vol
    dpl = mpsini[0]*5.973e27/vol

    if len(np.array(mcmass).shape) > 0:
        mcdpl = mcmass[0,:]*5.973e27/vol
        return mcdpl
    else:
        return dpl


   

def rvfit_lsqmdl(orbel,jdb,rv,srv,jitter=0, param_names=0,npoly=0,circ=0, tt=np.zeros(1),epoch=2.455e6,pfix=1,norbits=1,telvec=-1,psig=-1,ttfloat='no',ttsig=-1):
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
        
    #ttfloat = orbel[norbits*6 + npoly + ntel-1 + 0-norbits] norbits terms

   
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
        if ttfloat == 'yes':
            tr_names = ['transit']*norbits
            param_names.extend(tr_names)

    flt = np.ones_like(orbel) #Free param, Turn off individually

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
        elif psig[i] > 0: #limit with error from Kepler pipeline
            m.lm_prob.p_limit(0+i*6, lower=orbel[0+i*6]-5*psig[i], upper=orbel[0+i*6]+5*psig[i])
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

        #now include fixed transit time effects
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

    
    #transit time must be either fixed or constrained with ttsig
    #allow for variable transit time
    if ttfloat == 'yes':
        for i in range(norbits):
            ind = i + norbits*6 + npoly + ntel-1
            
            if ttsig[i] > 0:        #let float with limited range
                
                m.lm_prob.p_limit(ind, lower=orbel[ind]-5*ttsig[i], upper=orbel[ind]+5*ttsig[i])
                #now tie other parameters...
                if circ[i] == 1: #by convention tt=tp & omega=90
                    m.lm_prob.p_value(3+i*6, 90.0, fixed=True)
                    tiefunc2 = tie_tp_function(tt,i)
                    m.lm_prob.p_tie(1+i*6, tiefunc2)
                    flt[1+i*6] = 0
                    flt[3+i*6] = 0
                else: #tie omega to tt & tp
                    tiefunc  = tie_omega_function(tt,i)
                    m.lm_prob.p_tie(3+i*6, tiefunc)
                    flt[3+i*6] = 0
                
            else:  #fix at original value
                m.lm_prob.p_value(ind, tt[i], fixed=True)
    #else test that you have the right number of orbit params...
    

    m.solve(orbel)
   
    return m, flt

def tie_tp_function(tt,i):
    def calculate_tp(orbel):

        tp = tt[i]
        return tp

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

        #print tp, om
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

def rv_drive_mc(orbel, t, norbits, npoly, telvec, tt, ttsig, ttfloat,circ):
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
        
        #print tp, om

        #if transit time set
        if not tt[i] == 0:
            #if fixed, use tt[i] else use orbel[tind]
            if ttfloat == 'no': #use default value
                if circ[i] == 1:
                    tp = tt[i]
                    om = 90.0*np.pi/180.0
                else: 
                    theta_tt = calc_true_anomaly(p, tp, ecc, tt[i])
                    om = (np.pi/2.0 - theta_tt)*180.0/np.pi

                    om = om % 360.0
                    #if om < 0:    
                    #    om += 360.0
                    #if om > 360.0:  #do something smarter here...
                    #    om -= 360.0
            
            elif ttfloat == 'yes' and ttsig[i] > 0: #use value from orbel
                tind = i + norbits*6 + npoly + ntel-1
                tti = orbel[tind]
                if circ[i] == 1:
                    tp = tti
                    om = 90.0*np.pi/180.0
                else:
                    theta_tt = calc_true_anomaly(p, tp, ecc, tti)
                    om = (np.pi/2.0 - theta_tt)*180.0/np.pi
                    om = om % 360.0

                    #if om < 0:    
                    #    om += 360.0
                    #if om > 360.0:  #do something smarter here...
                    #    om -= 360.0
        #print tp, om
                        
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
def lnprior(theta, fullpars, flt, pnames, plo, phi, norbit, psig, porig, tt, ttsig, ttfloat,npoly,telvec):
    

    #print fullpars
    if len(np.array(telvec).shape) > 0:
        ntel = np.unique(telvec).size

    #figure out which params are varied and the associated limits
    f = np.squeeze(flt.nonzero())
    
    pfloat = pnames[f]
    lfloat = plo[f]
    hfloat = phi[f]

    
    pind = ut.arrayify(np.squeeze(np.where(pfloat == 'Per')))
    ttind = ut.arrayify(np.squeeze(np.where(pfloat == 'transit')))
 

    #flat priors for all
    if (theta >= lfloat).all() and (theta < hfloat).all():
        lnpri = 0.0
        
        for i in range(norbit):
            if psig[i] > 0:
                x = theta[pind[i]] #FIX THIS
                
                prper = norm.pdf(x,loc=porig[i],scale=psig[i])*psig[i]
                lnpri += np.log(prper)
                
            if ttfloat == 'yes':
                if ttsig[i] > 0:
                    x = theta[ttind[i]] #FIX THIS
                    
                    prtt = norm.pdf(x,loc=tt[i],scale=ttsig[i])*ttsig[i]
                    lnpri += np.log(prtt)
                    
        
        return lnpri
    else:
        
        return -np.inf
    

def lnprob(theta, jdb, rv, srv, fullpars, flt, pnames, plo, phi, norbit, npoly, telvec,psig,porig, tt, ttsig, ttfloat,circ):
    lp = lnprior(theta, fullpars, flt, pnames, plo, phi, norbit, psig, porig, tt, ttsig, ttfloat,npoly,telvec)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, jdb, rv, srv, fullpars, flt, norbit, npoly, telvec, tt, ttsig, ttfloat,circ)

def lnlike(theta, jdb, rv, srv, fullpars, flt, norbit, npoly, telvec, tt, ttsig, ttfloat,circ):
    
    newpars = np.copy(fullpars)
    
    newpars[flt.nonzero()] = theta
    
    #jitter = newpars[-ntel:]

    #model = rv_drive(newpars, jdb, norbit, npoly, telvec)
    model = rv_drive_mc(newpars, jdb, norbit, npoly, telvec, tt, ttsig, ttfloat,circ)

    #now add jitter
    tels = np.unique(telvec)
    ntel = tels.size
    stot = np.zeros(telvec.size)
    for i in range(ntel):
        a = np.squeeze(np.where(telvec == tels[i]))         
        stot[a] = np.sqrt(srv[a]**2 + newpars[-ntel+i]**2)

#    stot = np.sqrt(srv**2 + newpars[-1]**2) #add floating jitter term

    l0 = np.sum(np.log(1/(np.sqrt(2*np.pi)*stot))) #penalize high jitter values
    chisq = -0.5*np.sum((rv - model)**2/stot**2)

    return l0 + chisq 

def lnlike_base(bestpars, jdb, rv, srv, norbit, npoly, telvec, tt, ttsig, ttfloat,circ):
    
    #model = rv_drive(bestpars, jdb, norbit, npoly, telvec)
    model = rv_drive_mc(bestpars, jdb, norbit, npoly, telvec, tt, ttsig, ttfloat,circ)
    tels = np.unique(telvec)
    ntel = tels.size
    stot = np.zeros(telvec.size)
    for i in range(ntel):
        a = np.squeeze(np.where(telvec == tels[i]))         
        stot[a] = np.sqrt(srv[a]**2 + bestpars[-ntel+i]**2)
#    stot = np.sqrt(srv**2 + bestpars[-1]**2) #add floating jitter term

    l0 = np.sum(np.log(1/(np.sqrt(2*np.pi)*stot))) #penalize high jitter values
    chisq = -0.5*np.sum((rv - model)**2/stot**2)

    return l0 + chisq 

def setup_emcee(targname, m, jdb, rv, srv_in, nwalkers=200, circ=0, npoly=0, norbits=1, tt=np.zeros(1),jitter=0, pfix=1,nburn=300,telvec=-1,nsteps=1000,psig=-1,porig=-1,ttfloat='yes',ttsig=-1,fixjit='no'): 

    if len(np.array(telvec).shape) > 0:
        ntel = np.unique(telvec).size
        tels = np.unique(telvec)
    else:
        ntel = 1

    if jitter.size < ntel:
        jitter = np.ones(ntel)*jitter
    print 'jitter: ',jitter

    bestpars = np.copy(m.params)
    bestpars = np.append(bestpars,jitter) #add placeholder for jitter
    pnames = m.pnames
    jnames = ['jitter']*ntel
#    pnames.append(jnames)
    pnames = np.append(pnames,jnames)
    print 'smaller gaussian prior'

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
    #jitter = orbel[-ntel:]

    npars = bestpars.size

#    if jitter > 0:
#        print 'removing ',str(jitter),' m/s fixed jitter'
#        srv = np.sqrt(srv_in**2 - jitter**2)
#    else:
    srv = srv_in


    #these will hold reasonable limits on priors
    plo = np.zeros(npars)
    phi = np.zeros(npars)

    #separate the params being varied from the full list
    flt = np.ones(npars) #all params float now, turn off individually  

    
    for i in range(norbits):
       
        #fix normal orbit params first...
        
        #fix/limit period:
        if psig[i]> 0:
            plo[0+i*6] = porig[i] - 5*psig[i]#1000*psig[i]
            phi[0+i*6] = porig[i] + 5*psig[i]#1000*psig[i]  
        else:
            plo[0+i*6] = 0.1
            phi[0+i*6] = (np.max(jdb)-np.min(jdb))*20 
        if pfix[i] == 1:
            flt[0+i*6] = 0
        

        #limit Tp - within just over a full orbit
        plo[1+i*6] = bestpars[1+i*6]-bestpars[0+i*6]*0.6
        phi[1+i*6] = bestpars[1+i*6]+bestpars[0+i*6]*0.6

        #fix/limit ecc
        if circ[i] == 1:
            flt[2+i*6] = 0
        plo[2+i*6] = 0.0
        phi[2+i*6] = 0.99

        #limit omega
        if circ[i] == 1:
            flt[3+i*6] = 0
        else:
           bestpars[3+i*6] += 0.05 #poke it b/c it sticks at zero?
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


        #now consider fixed transit effects
        if not tt[i] == 0:
            if circ[i] == 1:
                flt[1+i*6] = 0
                bestpars[1+i*6] = tt[i] #this is not necessary?
                flt[3+i*6] = 0
                bestpars[3+i*6] = 90.0
            #else:
                #tie omega to ecc inside rv_drive_mc
                flt[3+i*6] = 0

    #limit polynomial terms
    for i in range(npoly):
        plo[i+norbits*6] = -1e6
        phi[i+norbits*6] = 1e6

    #need to limit offset terms if present!
    for i in range(ntel-1):
        plo[i+norbits*6+npoly] = -1e6
        phi[i+norbits*6+npoly] = 1e6

    #transit time must be either fixed or constrained with ttsig
    #now deal with floating transit effects
    if ttfloat == 'yes':
        for i in range(norbits):
            ind = i + norbits*6 + npoly + ntel-1
            
            if ttsig[i] > 0: #let float with limited range
                plo[ind] = bestpars[ind] - 5*ttsig[i]#1000*ttsig[i]
                phi[ind] = bestpars[ind] + 5*ttsig[i]#1000*ttsig[i]

                #now tie other params...
                if circ[i] == 1:
                    #tie tp to tt inside rv_drive_mc - needs tt,ttsig,ttfloat
                    bestpars[1+i*6] = tt[i]#this is not necessary?
                    flt[1+i*6] = 0

                    bestpars[3+i*6] = 90.0
                    flt[3+i*6] = 0
                else: #tie omega to tt & tp inside rv_drive
                    flt[3+i*6] = 0
                
            else: #dummy variable, fix at original value
                bestpars[ind] = tt[i]
                flt[ind] = 0
                plo[ind] = -1e8 #don't need this...
                phi[ind] = 1e8
    #else test that you have the right number of orbit params...

    #limit jitter
    if fixjit == 'yes':
        print 'adding ',str(jitter),'m/s fixed jitter '
        bestpars[-ntel:] = jitter
        flt[-ntel:] = 0
#    plo[-1] = 0.001
#    phi[-1] = 5.0
    plo[-ntel:] = 0.000
    phi[-ntel:] = 5.0

#    print pnames
#    print flt
#    print bestpars
#    print plo
#    print phi
       
    
    f = np.squeeze(flt.nonzero())
    
    varpars = bestpars[f]
    ndim = varpars.size
    
    print pnames

    print 'MCMC params: ',pnames[f] 
    print 'guesses from LM: ',varpars

    chain = run_emcee(targname, bestpars, varpars, flt, plo, phi, jdb, rv, srv, pnames, ndim, nwalkers=nwalkers, norbits=norbits, npoly=npoly, telvec=telvec,nsteps=nsteps,psig=psig,porig=porig,tt=tt,ttsig=ttsig,ttfloat=ttfloat,circ=circ)
  
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


def run_emcee(targname, bestpars, varpars, flt, plo, phi, jdb, rv, srv, pnames, ndim, nwalkers=200, nsteps=1000, norbits=1, npoly=0, telvec=-1,psig=-1,porig=-1,tt=np.zeros(1),ttsig=-1,ttfloat='no',circ=0):
    
    #Initialize walkers in tiny Gaussian ball around MLE results
    #number of params comes from varpars
    #pos = [varpars + 1e-3*np.random.randn(ndim) for i in range(nwalkers)] #??
    pos = [varpars + 1e-6*varpars*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(jdb,rv,srv,bestpars,flt, pnames, plo, phi, norbits, npoly, telvec,psig,porig, tt, ttsig, ttfloat,circ))

    #Run MCMC
    sampler.run_mcmc(pos, nsteps)
    print("Mean acceptance fraction: {0:.3f}"
                .format(np.mean(sampler.acceptance_fraction)))

    return sampler.chain


def gelman_rubin(chain):

    #(nchains, nsteps, npars)

    #discard first half of chain
    nchain = chain.shape[0] #m
    nsteps = chain.shape[1] #n
    npars = chain.shape[2]
    nkeep = nsteps/2

    last = chain[:,nkeep:,:]
 
    #calculate within-chain variance - this is an underestimate
    
    chvar = np.var(last,axis=1) #variance of each chain, shape(nchain,npars)
    W = 1.0/nchain*np.sum(chvar,axis=0) #shape(npars)
    print chvar.shape, W.shape

    #calculate between-chain variance
    chmean = np.mean(last,axis=1) #mean of each chain, shape(nchain,npars)
    meanmean = np.mean(chmean,axis=0) #shape(npars)

    B = np.zeros(npars)
    for i in range(npars):
        B[i] = np.float(nsteps)/(nchain-1)*np.sum((chmean[:,i] - meanmean[i])**2)

    #estimate variance of stationary dist
    esvar = (1 - 1.0/nsteps)*W + 1.0/nsteps*B

    R = np.sqrt(esvar/W)
    print 'Reduction factor: ', R
    return R
