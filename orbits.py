#RV orbit modeling code based on RVLIN by Jason Wright (IDL) and orbits.f by Alex Wolszczan (FORTRAN77)

#
# Branch ecoso
# PUT ALL TRANSITING PLANETS FIRST!
# Using time as generic, Tt or Tp as available
#


import emcee
import shelve
import socket
import triangle
import numpy as np
import read_rdb_harpsn as rr
import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
import utils as ut
from pwkit import lsqmdl, msmt
from scipy.stats import norm

def orbits_test(targname='K00273',jitter=0.5,epoch=2.4568478981528e6,circ=0,maxrv=1e6,minrv=-1e6,maxsrv=5, webdat='no', nwalkers=200, pfix=1,norbits=1,npoly=0,keck='no',outer_loop='no',nsteps=1000,nburn=500,fixjit='no',storeflat='yes',tfix=0,hd=0,machine='vonnegut0',thin=1,threads=1):

    tag = ''
    if npoly > 4:
        print 'Must have <= 4th order polynomial'


    #transit = np.zeros(1)
    telvec = np.zeros(1)
    circ = ut.arrayify(circ)
    pfix = ut.arrayify(pfix)
    jitter = ut.arrayify(jitter)
    tfix = ut.arrayify(tfix)

    for i in circ:
        if i == 1:
            tag += '_circ'
        else:
            tag += '_ecc'

    tag += '_'+str(npoly)
    tag += '_'+str(nsteps)

    host = socket.gethostname()
    if "macbook" in host:
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
        
    else:


        print 'Target: ',targname
        jdb, rv, srv, labels, fwhm, contrast, bis_span, rhk, sig_rhk, exptime = rr.process_all(targname,maxsrv=maxsrv,maxrv=maxrv,minrv=minrv,webdat=webdat)
        jdb += 2.4e6 #because process_all gives truncated JDBs
        telvec = np.zeros_like(jdb)
        print jdb.size,' HARPS-N obs'

        jdb0 = np.copy(jdb)
        rv0 = np.copy(rv)
        srv0 = np.copy(srv)

        if keck == 'yes':
            
            sfile = open(home+'Dropbox/cfasgettel/research/keck/'+targname+'_full.dat')
            kjdb, krv, ksrv, mdchi, kcts = np.loadtxt(sfile,unpack=True,usecols=(2,3,4,5,6),skiprows=3)
            sfile.close()
            
            kjdb = kjdb + 2.44e6


            krvnorm = krv - np.median(krv)
            ktel = np.ones_like(kjdb)
                
            jdb = np.append(jdb,kjdb)
            rvnorm = rv - np.median(rv)
            print 'median rv: ',np.median(rv)
            rvnorm = np.append(rvnorm,krvnorm)
            srv = np.append(srv,ksrv)
            telvec = np.append(telvec,ktel)
            print kjdb.size,' Keck obs'
            print np.mean(ksrv), ' typical Keck error '
        else:
            rvnorm = rv - np.median(rv)
            print 'median rv: ',np.median(rv)

    #adjust values to be sensible
    #if jdb[0] > epoch:
    #    print 'truncating dates'
    tnorm = jdb - epoch #truncate
    medrv = np.median(rv)

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
    #t0 = orbel[1+i*6] #either transit time or time of periastron
    #ecc = orbel[2+i*6]
    #om = orbel[3+i*6] *np.pi/180. #degrees to radians
    #k = orbel[4+i*6]
    #gamma = orbel[5+i*6]

    #dvdt = orbel[0+norbits*6] - optional polynomial fit
    #quad = orbel[1+norbits*6]
    #cubic = orbel[2+norbits*6]
    #quart = orbel[3+norbits*6]

    #offset = orbel[norbits*6 + npoly + (0-3)] #up to 4 offset terms

    #jitter - final term/s, MCMC only
    rpl = 0.0

    if targname == 'XO3':
        guesspars = np.array([3.19,2454024.728,0.29,346.3,1488,-1000.0])
        #transit = np.array([2454025.3967])
        ttime = np.array([1]) #yes, time is transit
        psig = np.array([1.4e-4])
        tsig = np.array([3.8e-3])
        porig = guesspars[0]*np.ones(1)
        torig = guesspars[1]*np.ones(1) - epoch
        mstar = 1.41
        rpl = 1.217
        rple = 0.073
        inc = 79.3
        ince = 1.36

    if targname == 'HD209458':
        guesspars = np.array([3.524733, 2452826.628514, 0.0, 336.5415, 85.49157+10, -1.49-10])#HD209458
        #transit = np.array([2452826.628514])
        ttime = np.array([1]) #time is transit
        mstar = 1.0
    
    if targname == 'K00273':

        #from literature
        mstar = 1.069
        rs = 1.081 #stellar radius
        ers = 0.019
        rs_dist = np.random.normal(loc=rs,scale=ers,size=10000) 

        #get transit params
        try:
            #read full transit posterior
            sfile = open(home+'Dropbox/cfasgettel/research/kepler/'+targname+'/fit_posteriors_koi273_shortcad.dat')
            p0_dist,tmod_dist,arstar_dist, rprs_dist, imp_dist= np.loadtxt(sfile,unpack=True,usecols=(0,1,2,3,4),skiprows=1)
            sfile.close()

            p0 = np.percentile(p0_dist,50)
            t0 = np.percentile(tmod_dist,50) + 2454833.0 + 16.0*p0
            p0_err = (np.percentile(p0_dist,84)-np.percentile(p0_dist,16))/2.
            t0_err = (np.percentile(tmod_dist,84)-np.percentile(tmod_dist,16))/2.
            print ' '
            print 'period: ',p0,' +',str(np.percentile(p0_dist,84)-p0),' -',str(p0-np.percentile(p0_dist,16))
            print 't0: ',t0,' +',str(np.percentile(tmod_dist,84)-np.percentile(tmod_dist,50)),' -',str(np.percentile(tmod_dist,50)-np.percentile(tmod_dist,16))
            print 'a/rstar: ',np.percentile(arstar_dist,50),' +',str(np.percentile(arstar_dist,84)-np.percentile(arstar_dist,50)),' -',str(np.percentile(arstar_dist,50)-np.percentile(arstar_dist,16))
            print 'rp/rstar: ',np.percentile(rprs_dist,50),' +',str(np.percentile(rprs_dist,84)-np.percentile(rprs_dist,50)),' -',str(np.percentile(rprs_dist,50)-np.percentile(rprs_dist,16))
            print 'impact: ',np.percentile(imp_dist,50),' +',str(np.percentile(imp_dist,84)-np.percentile(imp_dist,50)),' -',str(np.percentile(imp_dist,50)-np.percentile(imp_dist,16))

        except IOError:
            #reconstruct from error bars
            rprs = np.array([0.01596,0.00031,0.00085]) #rplanet/rstar, median, errlo, errhi
            #generate distribution of rprs - this isn't quite right...
            rprs_dist=msmt.sample_double_norm(rprs[0],rprs[2],rprs[1],4096)
            arstar = np.array([44.733919,8.3693767,3.2092898])
            arstar_dist = msmt.sample_double_norm(arstar[0],arstar[2],arstar[1],4096)
            imp = np.array([0.38190291,0.26185766,0.28049118]) #median, errlo, errhi
            imp_dist = msmt.sample_double_norm(imp[0],imp[2],imp[1],4096)

            p0 = 10.573763
            t0 = 2455008.0671344
            p0_err = 8.5e-6
            t0_err = 0.00078

        #calculate transit params
        rpl_dist = radius_estimate(rprs_dist,rs_dist)
        rpl = np.percentile(rpl_dist,50)
        print 'radius: ',str(rpl),' +',str(np.percentile(rpl_dist,84)-rpl),' -',str(rpl-np.percentile(rpl_dist,16))
        inc = inclination_estimate(arstar_dist,imp_dist)
        print 'inclination: ',np.percentile(inc,50),' +',str(np.percentile(inc,84)-np.percentile(inc,50)),' -',str(np.percentile(inc,50)-np.percentile(inc,16))
        print ' '
        #attempt to free some memory
        del p0_dist, tmod_dist, arstar_dist, rprs_dist, imp_dist

        guesspars = np.array([p0, t0, 0.0, 90.0, 2.2, -16.0,  500.0, 2455041.9, 0.23, 340.0, 137.0,0.0,0.0,0.0])
        ttime = np.array([1,0])

        psig = np.array([p0_err,0.0])
        tsig = np.array([t0_err,0.0])
        porig = guesspars[0+ip*6]
        torig = guesspars[1+ip*6] - epoch
        
    if targname == 'K00069':
    
        guesspars = np.array([4.72673978, 2454944.29227, 0.0, 90.0, 1.733, -91.08, 0.0329])
        ttime = np.array([1])
        #transit = np.array([2454944.29227])
        mstar = 0.887

    
    #trim JDB terms
    ip = np.arange(norbits)
    guesspars[1+ip*6] -= epoch

    #trim unused polynomial terms
    npars0 = norbits*6 + npoly
    guesspars = guesspars[0:npars0] 

    #append offsets if needed
    if ntel > 1:
        guesspars = np.append(guesspars,offset)



    if outer_loop == 'yes':
        print 'defunct!'

    else:
        m, flt = rvfit_lsqmdl(guesspars, tnorm, rvnorm, nsrv, jitter=jitter,circ=circ, npoly=npoly,ttime=ttime,epoch=epoch,pfix=pfix,norbits=norbits,telvec=telvec,psig=psig,tfix=tfix,tsig=tsig)

    
        #display initial fit - want to show fixed params too
        m.print_soln()  
    
        #calc bayesian information criterion
        bic0 = calc_bic_lm(m,flt,srv)
        print 'BIC:          ',str(bic0)

        #mass estimate*sin(i)
        mpsini, a2sini = mass_estimate(m, mstar, norbits=norbits, inc=inc)
        print 'mp*sin(i):         ',str(mpsini)
        #print 'mp:      ',str(mp)

        #density estimate
        dpl = density_estimate(mpsini,rpl)
        print 'density*sin(i):           ',str(dpl),'g/cc '

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
        res1, res0 = plot_rv(targname,tnorm,rvp,nsrv,guesspars,m,nmod=200,home=home,norbits=norbits,npoly=npoly,telvec=telvec,keck=keck,ttime=ttime,tag=tag)
        m.params = par0
       
        #call MCMC    
        if nwalkers > 0:
           

            mcbest, bestpars, mcnames, flt, mcpars, chain = setup_emcee_new(targname, m, tnorm, rvnorm, nsrv, tag, circ=circ, npoly=npoly, ttime=ttime, jitter=jitter, nwalkers=nwalkers, pfix=pfix, telvec=telvec, norbits=norbits, nsteps=nsteps,psig=psig,porig=porig,nburn=nburn,tfix=tfix,tsig=tsig,fixjit=fixjit,torig=torig,home=home,hd=hd,machine=machine,thin=thin,threads=threads)
            #mcbest = np.percentile(mcpars,50, axis=0) 

            #correct offset before plotting
            rvp = np.copy(rvnorm)
            print 'mcbest:', mcbest
            par1 = np.copy(mcbest)
            for i in range(ntel-1):
                a = np.squeeze(np.where(telvec == tels[i+1]))
                print 'MC offset: ', mcbest[i+norbits*6+npoly]
                rvp[a] -= mcbest[i+norbits*6+npoly]
                mcbest[i+norbits*6+npoly] = 0
            
            res1, pres = plot_rv(targname,tnorm,rvp,nsrv,mcbest,m,nmod=200,home=home,norbits=norbits,npoly=npoly,telvec=telvec,mc=1,keck=keck,ttime=ttime,tag=tag)
            
            if storeflat == 'yes':
                f = open(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+tag+'_rvflat.dat','w')
                f.write('#tnorm    res    srv    tel            \n') 
                for i in range(tnorm.size):
                    

                    f.write(str(tnorm[i])+' '+str(pres[i])+' '+str(nsrv[i])+' '+str(telvec[i])+'\n')
                f.close()
                f = open(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+tag+'_res.dat','w')
                for i in range(tnorm.size):
                    f.write(str(tnorm[i])+' '+str(res1[i])+' '+str(nsrv[i])+' '+str(telvec[i])+'\n') 
                f.close()
                
            mcbest = par1
            
            mpsini, a2sini, mparr_mc, a2arr_mc = mass_estimate(m, mstar, norbits=norbits, mcpar=mcpars)
            dpl = density_estimate(mpsini,rpl, mcmass=mparr_mc, rpl_dist=rpl_dist)
            
            #plot_hists(mcpars,mparr_mc,mcnames,flt,norbits=norbits)

            #print output from mass_estimate for mc
            print_mc_errs(mcpars, mpsini, a2sini, mparr_mc, a2arr_mc, norbits=norbits,npoly=npoly,telvec=telvec,tfix=tfix,tsig=tsig,mcdpl=dpl,inc=inc,ttime=ttime)
            bic = calc_bic_mc(mcbest, flt, tnorm, rvnorm, nsrv, norbits, npoly, telvec, ttime, tsig, tfix,circ)

            print 'MC BIC:          ',str(bic)

            #convergence testing...
            psrf = gelman_rubin(chain)            

        else:
            mcpars = -1
            mparr_mc = -1
            chain = -1
            bic = -1
            psrf = -1
            a2arr_mc = -1

            
        
        write_full_soln(m, targname, mpsini, a2sini, bic0, mcpars=mcpars, mparr_mc=mparr_mc, norbits=norbits, npoly=npoly, telvec=telvec, tfix=tfix, tsig=tsig,mbic=bic,psrf=psrf, a2arr_all=a2arr_mc,home=home,ttime=ttime,tag=tag,mcdpl=dpl,inc=inc)

        #store data as binary
        #np.save(home+'/Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+tag+'_chain.dat',chain)
        if home == '/home/sgettel/' and hd == 0:
            np.save('/pool/'+machine+'/harpsn/mass_estimate/'+targname+tag+'_mass.dat',mparr_mc)
            np.save('/pool/'+machine+'/harpsn/mass_estimate/'+targname+tag+'_density.dat',dpl)
            np.save('/pool/'+machine+'/harpsn/mass_estimate/'+targname+tag+'_mcpars.dat',mcpars)
            np.save('/pool/'+machine+'/harpsn/mass_estimate/'+targname+tag+'_mcnames.dat',mcnames)
            np.save('/pool/'+machine+'/harpsn/mass_estimate/'+targname+tag+'_flt.dat',flt)
            np.save('/pool/'+machine+'/harpsn/mass_estimate/'+targname+tag+'_bestpars.dat',bestpars)
        else:
            np.save(home+'/Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+tag+'_mass.dat',mparr_mc)
            np.save(home+'/Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+tag+'_density.dat',dpl)
            np.save(home+'/Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+tag+'_mcpars.dat',mcpars)
            np.save(home+'/Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+tag+'_mcnames.dat',mcnames)
            np.save(home+'/Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+tag+'_flt.dat',flt)
            np.save(home+'/Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+tag+'_bestpars.dat',bestpars)
        

    return m, chain, mparr_mc, mcpars, mcnames, flt, bestpars, tag, nwalkers, nsteps, nburn 

def fromeccom(ecc,om):
    #need om in radians
    sqesinom = np.sqrt(ecc)*np.sin(om) 
    sqecosom = np.sqrt(ecc)*np.cos(om)
    
    return sqesinom, sqecosom


def toeccom(sqesinom,sqecosom):
    #gives om in radians
    
    sqesinom = ut.arrayify(sqesinom)
    sqecosom = ut.arrayify(sqecosom)
    
    ecc = sqesinom**2 + sqecosom**2
    om = np.arctan2(sqesinom,sqecosom)
    bad = np.squeeze(np.where(sqecosom == 0))
    if bad.size > 0:
        om[bad] = 90.*np.pi/180.0*np.ones(bad.size)



    return ecc, om


def calc_bic_lm(m, flt, srv):
    
    #calc BIC for LM fit, use original errors
    chisq = np.sum(m.resids**2/srv**2)
    
    k = np.sum(flt)
    n = srv.size
    print chisq, k, n
    bic = chisq + k*np.log(n) 

    return bic
    
def calc_bic_mc(mcbest,flt,jdb,rv,srv,norbit,npoly,telvec, ttime, tsig, tfix,circ):
  
    k = np.sum(flt)
    n = srv.size
  
    lnl = lnlike_base(mcbest, jdb, rv, srv, norbit, npoly, telvec, ttime, tsig, tfix,circ)
    bic = -2*lnl + k*np.log(n)

    return bic
     

def plot_rv(targname,jdb,rv,srv,gpars,m,nmod=1000,home='/home/sgettel/', norbits=1,npoly=0,telvec=-1,mc=0,keck='no',ttime=0,tag=''):
    print 'ttime: ',ttime
    if mc > 0:
        usepars = np.copy(gpars) #this is mcbest...
        for i in range(norbits):
            ecc, om0 = toeccom(usepars[2+i*6],usepars[3+i*6])
            usepars[2+i*6] = ecc
            usepars[3+i*6] = om0*180./np.pi
    else:
        usepars = np.copy(m.params)

   
    tmod = np.linspace(np.min(jdb),np.max(jdb),nmod)

    model_final = rv_drive_new(usepars,tmod,norbits,npoly,telvec,ttime)
    res1 = rv - rv_drive_new(usepars,jdb,norbits,npoly,telvec,ttime)
    
    if npoly > 0:
        parst =  np.copy(usepars) 
        parst[4] = 0.0
        poly = rv_drive_new(parst,tmod,norbits,npoly,telvec,ttime)

    k = np.squeeze(np.where(telvec == 1))

    #unphased data, now with residuals!
    plt.figure(1)
    gs = grd.GridSpec(2,1, height_ratios=[3,1])
    ax1 = plt.subplot(gs[0])
    plt.errorbar(jdb,rv,yerr=srv,fmt='bo')
    if keck == 'yes':
        plt.errorbar(jdb[k],rv[k],yerr=srv[k],fmt='go') 
    plt.plot(tmod,model_final,'r-')
    #plt.xlabel('Adjusted BJD')
    plt.ylabel('Normalized RV (m/s)')

    ax2 = plt.subplot(gs[1])
    plt.errorbar(jdb,res1,yerr=srv,fmt='bo')
    if keck == 'yes':
        plt.errorbar(jdb[k],res1[k],yerr=srv[k],fmt='go') 
    plt.plot(tmod,np.zeros_like(tmod),'k-')
    plt.xlabel('Adjusted BJD')
    plt.ylabel('Residuals (m/s)')

    if mc > 0:
       plt.savefig(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+tag+'_autoplot_mc.pdf') 
       plt.savefig(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+tag+'_autoplot_mc.png')
    else:
        plt.savefig(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+tag+'_autoplot.pdf')
        plt.savefig(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+tag+'_autoplot.png')

    plt.close(1)

    #phase at 1st period
    pars = usepars[0:6]  #for model 
    pars[6:] = 0 #select first planet only
        
    parst = np.copy(usepars)
    parst[4] = 0.0 #other planets only
    parst[5] = 0.0
    phase = (jdb - pars[1])/pars[0] % 1.0

    rvt = rv_drive_new(parst,jdb,norbits,npoly,telvec,ttime)
    pres = (rv-rvt) - pars[5] 
    telvec = np.zeros_like(tmod)

    plt.figure(2)
    gs = grd.GridSpec(2,1, height_ratios=[3,1])
    ax1 = plt.subplot(gs[0])
    plt.errorbar(phase, rv-rvt, yerr=srv,fmt='bo')
    if keck == 'yes':
        plt.errorbar(phase[k],rv[k]-rvt[k],yerr=srv[k],fmt='go') 
    plt.plot((tmod - pars[1])/pars[0] % 1.0, rv_drive_new(pars, tmod,1,0,telvec,ttime),'r.')
    #plt.xlabel('Orbital Phase')
    plt.ylabel('Normalized RV (m/s)')
    #plt.plot((tmod - guess))

    ax2 = plt.subplot(gs[1])
    res2 = rv - rvt - rv_drive_new(pars,jdb,1,0,telvec,ttime)
    plt.errorbar(phase,res2,yerr=srv,fmt='bo')
    if keck == 'yes':
        plt.errorbar(phase[k],res2[k],yerr=srv[k],fmt='go') 
    plt.plot((tmod- pars[1])/pars[0] % 1.0,np.zeros_like(tmod),'k-')
    plt.xlabel('Orbital Phase')
    plt.ylabel('Residuals (m/s)')

    if mc > 0:
        plt.savefig(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+tag+'_phase_autoplot_mc.pdf')
        plt.savefig(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+tag+'_phase_autoplot_mc.png')
    else:
        plt.savefig(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+tag+'_phase_autoplot.pdf') 
        plt.savefig(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+tag+'_phase_autoplot.png')
    plt.close(2)
    
    return res1, pres

def write_full_soln(m,targname,mpsini, a2sini, bic, mcpars=-1, mparr_mc=-1,norbits=1,npoly=0,telvec=-1,tt=np.zeros(1),tsig=-1,tfix=0,mbic=-1,psrf=-1,a2arr_all=-1,home='/home/sgettel',ttime=0,tag='',mcdpl=-1,inc=-1):
    
    poly_names = ['dvdt:  ','quad:  ', 'cubic: ','quart: ']

    if len(np.array(telvec).shape) > 0:
        tels = np.unique(telvec)
        ntel = np.unique(telvec).size
    
    f = open(home+'/Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+tag+'_orbit.dat','w')
    for i in range(norbits):
        f.write('                                               \n') 
        f.write('*****Planet '+str(i+1)+' Solution:***** \n')
        f.write('Per: '+str(m.params[0+i*6])+'\n')
        if ttime[i] == 1:
            f.write('Tt: '+str(m.params[1+i*6])+'\n')
        else:
            f.write('Tp: '+str(m.params[1+i*6])+'\n')
        f.write('ecc: '+str(m.params[2+i*6])+'\n')
        f.write('om: '+str(m.params[3+i*6])+'\n')
        f.write('K1: '+str(m.params[4+i*6])+'\n')
        f.write('gamma: '+str(m.params[5+i*6])+'\n')
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
            eccs, oms = toeccom(mcpars[:,2+i*6],mcpars[:,3+i*6])

            f.write('                                               \n')
            f.write('*****Planet '+str(i+1)+' MCMC Errors:***** \n')    
            f.write('Per: '+str(mcbest[0+i*6])+' +'+str(mchi[0+i*6] - mcbest[0+i*6])+' -'+str(mcbest[0+i*6] - mclo[0+i*6])+'\n')
            if ttime[i] == 1:
                f.write('Tt: '+ str(mcbest[1+i*6])+' +'+str(mchi[1+i*6] - mcbest[1+i*6])+' -'+str(mcbest[1+i*6] - mclo[1+i*6])+'\n')
            else:
                f.write('Tp: '+ str(mcbest[1+i*6])+' +'+str(mchi[1+i*6] - mcbest[1+i*6])+' -'+str(mcbest[1+i*6] - mclo[1+i*6])+'\n')
            f.write('sqesinom: '+ str(mcbest[2+i*6])+' +'+str(mchi[2+i*6] - mcbest[2+i*6])+' -'+str(mcbest[2+i*6] - mclo[2+i*6])+'\n')
            f.write('sqecosom: '+str(mcbest[3+i*6])+' +'+str(mchi[3+i*6] - mcbest[3+i*6])+' -'+str(mcbest[3+i*6] - mclo[3+i*6]) +'\n')
            f.write('ecc: '+str(np.percentile(eccs,50))+' +'+str(np.percentile(eccs,84)-np.percentile(eccs,50))+' -'+str(np.percentile(eccs,50)-np.percentile(eccs,16))+'\n')
            f.write('ecc1: '+ str(np.percentile(eccs,0))+'-'+str(np.percentile(eccs,68))+'\n')
            f.write('om: '+str(np.percentile(oms,50))+' +'+str(np.percentile(oms,84)-np.percentile(oms,50))+' -'+str(np.percentile(oms,50)-np.percentile(oms,16))+'\n')
            f.write('K1: '+str(mcbest[4+i*6])+' +'+str(mchi[4+i*6] - mcbest[4+i*6])+' -'+str(mcbest[4+i*6] - mclo[4+i*6]) +'\n')
            f.write('gamma: '+ str(mcbest[5+i*6])+' +'+str(mchi[5+i*6] - mcbest[5+i*6])+' -'+str(mcbest[5+i*6] - mclo[5+i*6])+'\n')
            
            mpbest = np.percentile(mparr_mc[i,:], 50)
            mphi = np.percentile(mparr_mc[i,:], 84)
            mplo = np.percentile(mparr_mc[i,:], 16)
            a2best = np.percentile(a2arr_all[i,:], 50)
            a2hi = np.percentile(a2arr_all[i,:], 84)
            a2lo = np.percentile(a2arr_all[i,:], 16)

            f.write('nsamples: '+str(mparr_mc.shape[1])+'\n')
            f.write('convergence: '+str(psrf)+'\n')

            f.write('mp*sin(i): '+str(mpbest)+' +'+str(mphi-mpbest)+' -'+str(mpbest-mplo)+'\n')
            f.write('mass error:'+ str((mphi-mpbest)/mpbest*100)+','+str((mpbest-mplo)/mpbest*100)+ '%'+'\n') 
            f.write('a2*sin(i): '+str(a2best)+' +'+str(a2hi-a2best)+' -'+str(a2best-a2lo)+'\n')

            if i == 0:
                dpbest = np.percentile(mcdpl, 50)
                dplo = np.percentile(mcdpl, 16)
                dphi = np.percentile(mcdpl, 84)
                
                f.write('density*sin(i): '+str(dpbest)+' +'+str(dphi-dpbest)+' -'+str(dpbest-dplo)+'\n')

            if len(np.array(inc).shape) > 0:
                #draw from inclination distribution
                inc_dist = np.random.choice(inc,size=mparr_mc.shape[1])*np.pi/180.

                mpbest_cor = np.percentile(mparr_mc[i,:]/np.sin(inc_dist), 50)
                mphi_cor = np.percentile(mparr_mc[i,:]/np.sin(inc_dist), 84)
                mplo_cor = np.percentile(mparr_mc[i,:]/np.sin(inc_dist), 16)
                a2best_cor = np.percentile(a2arr_all[i,:]/np.sin(inc_dist), 50)
                a2hi_cor = np.percentile(a2arr_all[i,:]/np.sin(inc_dist), 84)
                a2lo_cor = np.percentile(a2arr_all[i,:]/np.sin(inc_dist), 16)

                f.write('mp: '+str(mpbest_cor)+' +'+str(mphi_cor-mpbest_cor)+' -'+str(mpbest_cor-mplo_cor)+'\n')
                f.write('mass error:'+ str((mphi_cor-mpbest_cor)/mpbest_cor*100)+','+str((mpbest_cor-mplo_cor)/mpbest_cor*100)+ '%'+'\n') 
                f.write('a2: '+str(a2best_cor)+' +'+str(a2hi_cor-a2best_cor)+' -'+str(a2best_cor-a2lo_cor)+'\n')

                if i == 0:
                    dpbest_cor = np.percentile(mcdpl/np.sin(inc_dist), 50)
                    dplo_cor = np.percentile(mcdpl/np.sin(inc_dist), 16)
                    dphi_cor = np.percentile(mcdpl/np.sin(inc_dist), 84)
                
                    f.write('density: '+str(dpbest_cor)+' +'+str(dphi_cor-dpbest_cor)+' -'+str(dpbest_cor-dplo_cor)+'\n')

        for i in range(npoly):
            f.write(str(poly_names[i])+ str(mcbest[i+norbits*6])+' +'+str(mchi[i+norbits*6]-mcbest[i+norbits*6])+' -'+str(mcbest[i+norbits*6]-mclo[i+norbits*6]) +'\n')
        for i in range(ntel-1):
            a = np.squeeze(np.where(telvec == tels[i+1]))
            f.write('offset: '+str(mcbest[i+norbits*6+npoly])+' +'+str(mchi[i+norbits*6+npoly]-mcbest[i+norbits*6+npoly])+' -'+str(mcbest[i+norbits*6+npoly]-mclo[i+norbits*6+npoly])+'\n')

        for i in range(ntel):
            f.write('jitter: '+ str(mcbest[-ntel+i])+' +'+str(mchi[-ntel+i]-mcbest[-ntel+i])+' -'+str(mcbest[-ntel+i]-mclo[-ntel+i])+'\n')

        f.write('MC BIC: '+str(mbic)+'\n')
    f.close() 

def print_mc_errs(mcpars, mpsini, a2sini, mparr_all, a2arr_all,norbits=1,npoly=0,telvec=-1,ttime=0,tsig=-1,tfix=0,mcdpl=-1,inc=-1):
    poly_names = ['dvdt:  ','quad:  ', 'cubic: ','quart: ']

    if len(np.array(telvec).shape) > 0:
        tels = np.unique(telvec)
        ntel = np.unique(telvec).size


    mcbest = np.percentile(mcpars,50, axis=0)
    mchi = np.percentile(mcpars,84, axis=0)
    mclo = np.percentile(mcpars,16, axis=0) 

    for i in range(norbits):
        
        eccs, oms = toeccom(mcpars[:,2+i*6],mcpars[:,3+i*6])  

        print '                                               '
        print '*****Planet ',str(i+1),' MCMC Errors:*****'  
        print 'Per: ', str(mcbest[0+i*6]),' +',str(mchi[0+i*6] - mcbest[0+i*6]),' -',str(mcbest[0+i*6] - mclo[0+i*6])
        if ttime[i] == 1:
            print 'Tt: ', str(mcbest[1+i*6]),' +',str(mchi[1+i*6] - mcbest[1+i*6]),' -',str(mcbest[1+i*6] - mclo[1+i*6])
        else:
           print 'Tp: ', str(mcbest[1+i*6]),' +',str(mchi[1+i*6] - mcbest[1+i*6]),' -',str(mcbest[1+i*6] - mclo[1+i*6]) 
        print 'sqesinom: ', str(mcbest[2+i*6]),' +',str(mchi[2+i*6] - mcbest[2+i*6]),' -',str(mcbest[2+i*6] - mclo[2+i*6])
        print 'sqecosom: ', str(mcbest[3+i*6]),' +',str(mchi[3+i*6] - mcbest[3+i*6]),' -',str(mcbest[3+i*6] - mclo[3+i*6])
        print 'ecc: ',str(np.percentile(eccs,50)),' +',str(np.percentile(eccs,84)-np.percentile(eccs,50)),' -',str(np.percentile(eccs,50)-np.percentile(eccs,16))
        print 'ecc1: ',str(np.percentile(eccs,0)),'-',str(np.percentile(eccs,68))
        print 'om: ',str(np.percentile(oms,50)),' +',str(np.percentile(oms,84)-np.percentile(oms,50)),' -',str(np.percentile(oms,50)-np.percentile(oms,16))
        print 'K1: ', str(mcbest[4+i*6]),' +',str(mchi[4+i*6] - mcbest[4+i*6]),' -',str(mcbest[4+i*6] - mclo[4+i*6])
        print 'gamma: ', str(mcbest[5+i*6]),' +',str(mchi[5+i*6] - mcbest[5+i*6]),' -',str(mcbest[5+i*6] - mclo[5+i*6])

        mpbest = np.percentile(mparr_all[i,:], 50)
        mphi = np.percentile(mparr_all[i,:], 84)
        mplo = np.percentile(mparr_all[i,:], 16)
        a2best = np.percentile(a2arr_all[i,:], 50)
        a2hi = np.percentile(a2arr_all[i,:], 84)
        a2lo = np.percentile(a2arr_all[i,:], 16)


        print 'mp*sin(i): ',str(mpbest),' +',str(mphi-mpbest),' -',str(mpbest-mplo)
        print 'mass error:', str((mphi-mpbest)/mpbest*100),',',str((mpbest-mplo)/mpbest*100), '%'
        print 'a2*sin(i): ',str(a2best),' +',str(a2hi-a2best),' -',str(a2best-a2lo)
        
        if i == 0: 
            dpbest = np.percentile(mcdpl, 50)
            dphi = np.percentile(mcdpl, 84)
            dplo = np.percentile(mcdpl, 16)
                
            print 'density*sin(i): '+str(dpbest)+' +'+str(dphi-dpbest)+' -'+str(dpbest-dplo)


        if len(np.array(inc).shape) > 0: 

            #draw from inclination distribution
            inc_dist = np.random.choice(inc,size=mparr_all.shape[1])*np.pi/180.

            mpbest_cor = np.percentile(mparr_all[i,:]/np.sin(inc_dist), 50)
            mphi_cor = np.percentile(mparr_all[i,:]/np.sin(inc_dist), 84)
            mplo_cor = np.percentile(mparr_all[i,:]/np.sin(inc_dist), 16)
            a2best_cor = np.percentile(a2arr_all[i,:]/np.sin(inc_dist), 50)
            a2hi_cor = np.percentile(a2arr_all[i,:]/np.sin(inc_dist), 84)
            a2lo_cor = np.percentile(a2arr_all[i,:]/np.sin(inc_dist), 16)

            print 'mp: ',str(mpbest_cor),' +',str(mphi_cor-mpbest_cor),' -',str(mpbest_cor-mplo_cor)
            print 'mass error:', str((mphi_cor-mpbest_cor)/mpbest_cor*100),',',str((mpbest_cor-mplo_cor)/mpbest_cor*100), '%'
            print 'a2: ',str(a2best_cor),' +',str(a2hi_cor-a2best_cor),' -',str(a2best_cor-a2lo_cor) 
            
            if i == 0:           
                dpbest_cor = np.percentile(mcdpl/np.sin(inc_dist), 50)
                dphi_cor = np.percentile(mcdpl/np.sin(inc_dist), 84)
                dplo_cor = np.percentile(mcdpl/np.sin(inc_dist), 16)
                
            print 'density: '+str(dpbest_cor)+' +'+str(dphi_cor-dpbest_cor)+' -'+str(dpbest_cor-dplo_cor)

    
    for i in range(npoly):
        print str(poly_names[i]), str(mcbest[i+norbits*6]),' +',str(mchi[i+norbits*6]-mcbest[i+norbits*6]),' -',str(mcbest[i+norbits*6]-mclo[i+norbits*6])

    
    for i in range(ntel-1):
        #a = np.squeeze(np.where(telvec == tels[i+1]))
        print 'offset: ',str(mcbest[i+norbits*6+npoly]),' +',str(mchi[i+norbits*6+npoly]-mcbest[i+norbits*6+npoly]),' -',str(mcbest[i+norbits*6+npoly]-mclo[i+norbits*6+npoly])
    
    for i in range(ntel):
        print 'jitter: ', str(mcbest[-ntel+i]),' +',str(mchi[-ntel+i]-mcbest[-ntel+i]),' -',str(mcbest[-ntel+i]-mclo[-ntel+i])
    return

def mass_estimate(m, mstar, norbits=1, bootpar=-1, mcpar=-1, inc=-1):
   
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
    #a1sini = np.sqrt(1. - eccs**2)/(2.*np.pi)*amps*(pers*86400.0)/1.496e11
    a2sini = 1.96e-2*mstar**(1./3.)*pers**(2./3.)

#    if len(np.array(inc).shape) > 0:
#        mp = mpsini/np.sin(np.percentile(inc,50)*np.pi/180.)
#        a2 = a2sini/np.sin(np.percentile(inc,50)*np.pi/180.)
   

   # mass estimate for mcmc
    if len(np.array(mcpar).shape) > 0:
        mparr_mc = np.zeros((norbits,mcpar.shape[0]))
        a2arr_mc = np.zeros((norbits,mcpar.shape[0]))

        for i in ip:
            eccs, oms = toeccom(mcpar[:,2+i*6],mcpar[:,3+i*6])
            #fmarr = (1 - mcpar[:,2+i*6]**2)**(1.5)*mcpar[:,i*6+4]**3*(mcpar[:,i*6]*86400.0)/(2.0*np.pi*G)
            fmarr = (1 - eccs**2)**(1.5)*mcpar[:,i*6+4]**3*(mcpar[:,i*6]*86400.0)/(2.0*np.pi*G)
            mparr = ((mstar*msun)**2*fmarr)**(1./3.)/mearth
            a2arr = 1.96e-2*mstar**(1./3.)*(mcpar[:,i*6])**(2./3.)
            mparr_mc[i,:] = mparr
            a2arr_mc[i,:] = a2arr      

        return mpsini, a2sini, mparr_mc, a2arr_mc
    else:
        return mpsini, a2sini

def density_estimate(mpsini,rpl,mcmass=-1,rpl_dist=-1):
    #first planet only...

    rpl = ut.arrayify(rpl)
    nplanets = rpl.size

    re2cm = 1./6.378e8
    vol = 4./3.*np.pi*(rpl/re2cm)**3 #cm**3
    #dpl = mpsini[0:nplanets]*5.973e27/vol
    dpl = mpsini[0]*5.973e27/vol

    if len(np.array(mcmass).shape) > 0:

        #sample radius distribution to match size of mass distribution
        rpl = np.random.choice(rpl_dist,size=mcmass.shape[1])
        vol = 4./3.*np.pi*(rpl/re2cm)**3 #cm**3
        mcdpl = mcmass[0]*5.973e27/vol
        
        #mpbest = np.percentile(mcmass[0,:], 50)
        #mphi = np.percentile(mcmass[0,:], 84)
        #mplo = np.percentile(mcmass[0,:], 16)
        #mperr = ((mphi-mpbest)+(mpbest-mplo))/2 #assume nearly Gaussian...
        #dpl = mpbest*5.973e27/vol

        #relerrsq = (mperr/mpbest)**2+(3*rple/rpl) #lo,hi
        #edpl = np.sqrt(dpl**2*relerrsq) 
        #mcdpl = np.concatenate([dpl, edpl])
        #print mcdpl.shape, mcdpl
        
        return mcdpl #this is now a distribution
    else:
        return dpl

def radius_estimate(rr,rs,size=10000):

    #draw samples from input distributions
    rr = np.random.choice(rr,size=size)
    rs = np.random.choice(rs,size=size)

    rp_dist = rr*rs * 6.9599e8/6.37814e6 #earth radii

#    rp = rprs * rs * 6.9599e8/6.37814e6 #earth radii
#    erp = np.sqrt((ers/rs)**2 + (erprs/rprs)**2)*rp
    
    return rp_dist

def inclination_estimate(arstar,b,size=10000):

    #draw samples from input distributions
    b = np.random.choice(b,size=size)
    arstar = np.random.choice(arstar,size=size)
    cosi = b/arstar
    inc_dist = np.arccos(cosi)*180.0/np.pi #degrees

    #berr = np.mean(b[1:3]) #assume nearly symmetric
    #ince = np.sqrt(berr**2*(1./(arstar[0]**2 - b[0]**2)) + arstar[1:3]**2*(1./(arstar[0]**4*(1-b[0]/arstar[0]))))

    return inc_dist


def rvfit_lsqmdl(orbel,jdb,rv,srv,jitter=0, param_names=0,npoly=0,circ=0, ttime=0,epoch=2.455e6,pfix=1,norbits=1,telvec=-1,psig=-1,tfix='no',tsig=-1):
    #p = orbel[0+i*6]
    #t0 = orbel[1+i*6] #periastron or transit
    #ecc = orbel[2+i*6]
    #om = orbel[3+i*6] *np.pi/180. #degrees to radians
    #k = orbel[4+i*6]
    #gamma = orbel[5+i*6]

    #dvdt = orbel[0+norbits*6] - optional polynomial fit
    #quad = orbel[1+norbits*6]
    #cubic = orbel[2+norbits*6]
    #quart = orbel[3+norbits*6]

    #offset = orbel[norbits*6 + npoly + (0-3)] #up to 4 offset terms
        
    ip = np.arange(norbits)

    if len(np.array(telvec).shape) > 0:
        ntel = np.unique(telvec).size
    
    if param_names == 0: 
        param_names = ['Per', 'T0', 'ecc', 'om', 'K1', 'gamma']*norbits
        poly_names = ['dvdt','quad','cubic','quart']
        param_names.extend(poly_names[:npoly]) 
        if ntel > 1:
            off_names = ['offset']*(ntel-1)
            param_names.extend(off_names)

    flt = np.ones_like(orbel) #Free param, Turn off individually

    m = lsqmdl.Model(None, rv, 1./srv) 
    m.set_func(rv_drive_new,param_names, args=(jdb,norbits,npoly,telvec,ttime) )

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
        
        #limit T0
        if tfix[i] == 1:
            m.lm_prob.p_value(1+i*6, orbel[1+i*6], fixed=True)
            flt[1+i*6] = 0
        elif tsig[i] > 0:
            m.lm_prob.p_limit(1+i*6, lower=orbel[1+i*6]-5*tsig[i], upper=orbel[1+i*6]+5*tsig[i])
        else:
            m.lm_prob.p_limit(1+i*6, lower=orbel[1+i*6]-orbel[0+i*6]/2., upper=orbel[1+i*6]+orbel[0+i*6]/2.) #T0 within one period guess


        #fix/limit ecc
        if circ[i] == 1:
            m.lm_prob.p_value(2+i*6, 0.0, fixed=True) 
            flt[2+i*6] = 0
        else:
            m.lm_prob.p_limit(2+i*6, lower=0.0, upper=1.0) #ecc must be physical
        #fix/limit omega
        if circ[i] == 1:
            m.lm_prob.p_value(3+i*6, 90.0, fixed=True) #convention
            flt[3+i*6] = 0
        else:
            m.lm_prob.p_limit(3+i*6, lower=0.0, upper=360.0)

        #limit K
        m.lm_prob.p_limit(4+i*6, lower=0.0, upper=1.0e5) #K must be physical
  
        #fix gamma except first
        if i > 0:
            m.lm_prob.p_value(5+i*6, 0.0, fixed=True)
            flt[5+i*6] = 0

 

    #limit polynomial terms
    for i in range(npoly):
        
        m.lm_prob.p_limit(i+norbits*6, lower=-1e6, upper=1e6) #dvdt and higher


    #limit offset terms
    for i in range(ntel-1):
        m.lm_prob.p_limit(i + norbits*6 + npoly, lower=-1e6, upper=1e6)

    


    m.solve(orbel)
   
    return m, flt



def rv_drive_new(orbel, t, norbits, npoly, telvec, ttime):
    #From rv_drive.pro in RVLIN by JTW
    
    
    if len(np.array(telvec).shape) > 0:
        ntel = np.unique(telvec).size
    
    rv = np.zeros_like(t)

#    phase = np.zeros((rv.size,norbits))

    for i in range(norbits):
        p = orbel[0+i*6]
        t0 = orbel[1+i*6]
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
    
    #if T0 = transit time, solve for time of periastron
        if ttime[i] == 1:
            theta_tt = np.pi/2.0 - om #true anomaly of transit time
            argt = np.tan(theta_tt/2.)*np.sqrt((1.-ecc)/(1.+ecc))
            eccanom_tt = 2.*np.arctan(argt) #eccentric anomaly of transit time
            meananom_tt = eccanom_tt - ecc*np.sin(eccanom_tt) #mean anomaly of transit time
            phase = meananom_tt/(2.*np.pi)
            tp = t0 - phase*p
        
        #else T0 = time of periastron
        else:
            tp = t0

        theta = calc_true_anomaly(p, tp, ecc, t)
    

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


###################################

def rv_drive_mc_new(orbel, t, norbits, npoly, telvec, ttime, tsig, tfix,circ):
    #From rv_drive.pro in RVLIN by JTW
    
    
    if len(np.array(telvec).shape) > 0:
        ntel = np.unique(telvec).size
    
    rv = np.zeros_like(t)

    phase = np.zeros((rv.size,norbits))
    
    for i in range(norbits):
        p = orbel[0+i*6]
        t0 = orbel[1+i*6] # Tp or Tt
        #        ecc = orbel[2+i*6]
        #        om = orbel[3+i*6] *np.pi/180. #degrees to radians
        
        sqesinom = orbel[2+i*6]
        sqecosom = orbel[3+i*6]
        
        k = orbel[4+i*6]
        gamma = orbel[5+i*6]
        
        #convert sqecoso & sqesino back into ecc and om
        ecc = sqesinom**2 + sqecosom**2
        
        if sqecosom == 0:
            om = 90.0*np.pi/180.0 #is this right?
        else:
            om = np.arctan2(sqesinom,sqecosom)
        
        
        #Error checking
        if p < 0 or ecc < 0 or ecc >= 1 or k < 0:
            print 'Bad inputs to rv_drive_mc_new'
            print p, ecc, k
            #print om, sqesinom, sqecosom
            if p < 0:
                p = 1e-2
            if ecc < 0:
                ecc = 0
            if ecc >= 1:
                ecc = 0.99
            if k < 0:
                k = 1e-2


        #if T0 = transit time, solve for time of periastron
        if ttime[i] == 1:
            theta_tt = np.pi/2.0 - om #true anomaly of transit time
            argt = np.tan(theta_tt/2.)*np.sqrt((1.-ecc)/(1.+ecc))
            eccanom_tt = 2.*np.arctan(argt) #eccentric anomaly of transit time
            meananom_tt = eccanom_tt - ecc*np.sin(eccanom_tt) #mean anomaly of transit time
            phase = meananom_tt/(2.*np.pi)
            tp = t0 - phase*p
        
        #else T0 = time of periastron
        else:
            tp = t0




        theta = calc_true_anomaly(p, tp, ecc, t)
    

        
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




#begin MCMC setup - following emcee line-fitting demo
def lnprior(theta, fullpars, flt, pnames, plo, phi, norbit, psig, porig, ttime, tsig, tfix,npoly,telvec,circ,termspp,torig):
    

    #print fullpars
    if len(np.array(telvec).shape) > 0:
        ntel = np.unique(telvec).size

    #figure out which params are varied and the associated limits
    f = np.squeeze(flt.nonzero())
    
    pfloat = pnames[f]
    lfloat = plo[f]
    hfloat = phi[f]

    
    pind = ut.arrayify(np.squeeze(np.where(pfloat == 'Per')))
    tind = ut.arrayify(np.squeeze(np.where(pfloat == 'T0')))
 

    #flat priors for most
  
    if (theta >= lfloat).all() and (theta <= hfloat).all():
        lnpri = 0.0
        
        ilo = 0
        ihi = 0
        for i in range(norbit):

            #test for bad sqesinom & sqecosom
            ihi += termspp[i] 
            
            qn = pfloat[ilo:ihi]
            q = theta[ilo:ihi]
            #print qn

            i2 = np.squeeze(np.where(qn == 'sqesin'))
            i3 = np.squeeze(np.where(qn == 'sqecos'))
            
            if i2.size == 1 and i3.size == 1:
                etest = q[i2]**2 + q[i3]**2
                if etest > 1:
                    #print 'bad ecc rejected'
                    #print q[i2],q[i3]
                    return -np.inf
                

            ilo += termspp[i]
            if psig[i] > 0:
                x = theta[pind[i]] #FIX THIS
                
                prper = norm.pdf(x,loc=porig[i],scale=psig[i])*psig[i]
                lnpri += np.log(prper)
                
            if tfix[i] == 0:
                if tsig[i] > 0:
                    x = theta[tind[i]] #FIX THIS
                    
                    prtt = norm.pdf(x,loc=torig[i],scale=tsig[i])*tsig[i]
                    lnpri += np.log(prtt)
                    
        
        return lnpri
    else:
        
        return -np.inf
    

def lnprob(theta, jdb, rv, srv, fullpars, flt, pnames, plo, phi, norbit, npoly, telvec,psig, porig, ttime,tsig, tfix,circ,termspp,torig):
    lp = lnprior(theta, fullpars, flt, pnames, plo, phi, norbit, psig, porig, ttime, tsig, tfix,npoly,telvec,circ,termspp,torig)
    #print lp
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, jdb, rv, srv, fullpars, flt, norbit, npoly, telvec, ttime, tsig, tfix,circ)

def lnlike(theta, jdb, rv, srv, fullpars, flt, norbit, npoly, telvec, ttime, tsig, tfix,circ):
    
    newpars = np.copy(fullpars)
    
    newpars[flt.nonzero()] = theta
    
    #jitter = newpars[-ntel:]

    #print newpars
    #model = rv_drive_mc(newpars, jdb, norbit, npoly, telvec, tt, ttsig, ttfloat,circ)
    model = rv_drive_mc_new(newpars, jdb, norbit, npoly, telvec, ttime, tsig, tfix,circ)

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

    #print l0, chisq
    return l0 + chisq 

def lnlike_base(bestpars, jdb, rv, srv, norbit, npoly, telvec, tt, ttsig, ttfloat,circ):
    
    
    model = rv_drive_mc_new(bestpars, jdb, norbit, npoly, telvec, tt, ttsig, ttfloat,circ)
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


##########################################################################
def setup_emcee_new(targname, m, jdb, rv, srv_in, tag, nwalkers=200, circ=0, npoly=0, norbits=1, ttime=0,jitter=0, pfix=1,nburn=300,telvec=-1,nsteps=1000,psig=-1,porig=-1,tfix=0,tsig=-1,fixjit='no',torig=-1,home='',hd=0,machine='vonnegut',thin=1,threads=1):
    
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
    pnames = np.append(pnames,jnames)
    
    #want to take steps in sqrt(e)*cos(omega) sqrt(e)*sin(omega)
    #setup - swap variables, change name array
    mcin = bestpars
    mcnames = pnames
    for i in range(norbits):
        ecc = bestpars[2+i*6]
        om = bestpars[3+i*6]*np.pi/180. #deg to rad
        mcin[2+i*6] = np.sqrt(ecc)*np.sin(om) + 0.001 #perturb in case = 0
        mcin[3+i*6] = np.sqrt(ecc)*np.cos(om) + 0.001
        mcnames[2+i*6] = 'sqesinom'
        mcnames[3+i*6] = 'sqecosom'
    
    
    
    #p = mcin[0+i*6]
    #t0 = mcin[1+i*6] - Tp or Tt
    #sqesinom = mcin[2+i*6] = sqrt(ecc)*sin(om)
    #sqecosom = mcin[3+i*6]  = sqrt(ecc)*cos(om)
    #k = mcin[4+i*6]
    #gamma = mcin[5+i*6]
    #dvdt = mcin[0+norbits*6] - optional polynomial fit
    #quad = mcin[1+norbits*6]
    #cubic = mcin[2+norbits*6]
    #quart = mcin[3+norbits*6]
    #jitter = mcin[-ntel:]
    
    npars = mcin.size

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
        mcin[0+i*6] = porig[i]
        if psig[i] > 0:
            plo[0+i*6] = porig[i] - 5*psig[i]#1000*psig[i]
            phi[0+i*6] = porig[i] + 5*psig[i]#1000*psig[i]
        else:
            plo[0+i*6] = 0.1
            phi[0+i*6] = (np.max(jdb)-np.min(jdb))*20
        if pfix[i] == 1:
            flt[0+i*6] = 0
        
        #limit T0 - within just over a full orbit
        mcin[1+i*6] = torig[i]
        if tsig[i] > 0:
            plo[1+i*6] = torig[i] - 5*tsig[i]
            phi[1+i*6] = torig[i] + 5*tsig[i]
        else:
            plo[1+i*6] = mcin[1+i*6]-mcin[0+i*6]*0.6
            phi[1+i*6] = mcin[1+i*6]+mcin[0+i*6]*0.6
            
        if tfix[i] == 1:
            flt[1+i*6] = 0
     
        #fix/limit sqesinom & sqcosom
        if circ[i] == 1:
            flt[2+i*6] = 0
            flt[3+i*6] = 0
            mcin[2+i*6] = 0
            mcin[3+i*6] = 0
        plo[2+i*6] = -1.0
        phi[2+i*6] = 1.0
        plo[3+i*6] = -1.0
        phi[3+i*6] = 1.0
        
        #limit K
        plo[4+i*6] = 0.0
        phi[4+i*6] = 1.0e5
        
        #fix gamma except first
        if i > 0:
            flt[5+i*6] = 0
        plo[5+i*6] = -1e8
        phi[5+i*6] = 1e8

    #limit polynomial terms
    for i in range(npoly):
        plo[i+norbits*6] = -1e6
        phi[i+norbits*6] = 1e6
    
    #need to limit offset terms if present!
    for i in range(ntel-1):
        plo[i+norbits*6+npoly] = -1e6
        phi[i+norbits*6+npoly] = 1e6




    #limit jitter
    if fixjit == 'yes':
        print 'adding ',str(jitter),'m/s fixed jitter '
        mcin[-ntel:] = jitter
        flt[-ntel:] = 0
        
    plo[-ntel:] = 0.001
    phi[-ntel:] = 10.0

#    print pnames
#    print flt
#print bestpars
#    print plo
#    print phi

#    for i in range(bestpars.size):
#        print bestpars[i], plo[i], phi[i]

    termspp = np.zeros(norbits) #count how many terms of 6 go with each planet
    for i in range(norbits):
        termspp[i] = np.sum(flt[0+i*6:6+i*6])
    print 'termspp: ',termspp

    f = np.squeeze(flt.nonzero())
    
    varpars = mcin[f]
    ndim = varpars.size
    
    print pnames
    
    print 'MCMC params: ',pnames[f]
    print 'guesses from LM: ',varpars
    
    chain = run_emcee(targname, mcin, varpars, flt, plo, phi, jdb, rv, srv, pnames, ndim, termspp, tag, nwalkers=nwalkers, norbits=norbits, npoly=npoly, telvec=telvec,nsteps=nsteps,psig=psig,porig=porig,ttime=ttime,tsig=tsig,tfix=tfix,circ=circ, torig=torig,home=home,hd=hd,machine=machine,thin=thin,threads=threads)
    print chain.shape
    #It takes a number of iterations to spread walkers throughout param space
    #This is 'burning in'
    samples = chain[:, nburn/thin:, :].reshape((-1, ndim))
   
    #combine samples with best-fit values
    mcpars = ut.fan(mcin,samples.shape[0])
    print mcpars.shape
    #mcpars.shape = [nwalkers*(nsteps-nburn)/thin,n_mcin]
    
    for i in range(ndim):
        mcpars[:,f[i]] = samples[:,i]
    mcbest = np.percentile(mcpars,50, axis=0)
    
    
#return m, flt, chain, samples, mcpars
    return mcbest, mcin, mcnames, flt, mcpars, chain



def run_emcee(targname, bestpars, varpars, flt, plo, phi, jdb, rv, srv, pnames, ndim, termspp, tag, nwalkers=200, nsteps=1000, norbits=1, npoly=0, telvec=-1,psig=-1,porig=-1,ttime=0,tsig=-1,tfix=0,circ=0,torig=-1,home='',hd=0,machine='vonnegut',thin=1,threads=1):
    
    #Initialize walkers in tiny Gaussian ball around MLE results
    #number of params comes from varpars
    #pos = [varpars + 1e-3*np.random.randn(ndim) for i in range(nwalkers)] #??
    pos = [varpars + 1e-6*varpars*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, threads=threads,args=(jdb,rv,srv,bestpars,flt, pnames, plo, phi, norbits, npoly, telvec,psig,porig, ttime, tsig, tfix,circ,termspp,torig))


    #Run MCMC
    sampler.run_mcmc(pos, nsteps,thin=thin)
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

def make_triangle_after(basename,path='/pool/vonnegut0/harpsn/mass_estimate/'):

    #restore variables
    mcpars = np.load(path+basename+'_mcpars.dat.npy')
    mcnames = np.load(path+basename+'_mcnames.dat.npy')
    flt = np.load(path+basename+'_flt.dat.npy')
    bestpars = np.load(path+basename+'_bestpars.dat.npy')

    f = np.squeeze(flt.nonzero())
    
    #use subsample if large number of steps - seems to help
    if mcpars.shape[0] > 3000000:
        sub = np.random.randint(0,mcpars.shape[0],size=3000000)
        mcpars = mcpars[:,f]
        fig = triangle.corner(mcpars[sub,:], labels=mcnames[f], truths=bestpars[f])
    else:
        fig = triangle.corner(mcpars[:,f], labels=mcnames[f], truths=bestpars[f])
    
    
    fig.savefig(path+basename+'_triangle.png')
    fig.savefig(path+basename+'_triangle.pdf')
    plt.close(fig)

def plot_hists_after(basename,path='/pool/vonnegut0/harpsn/mass_estimate/',norbits=1):
#mcpars.shape = [nwalkers*(nsteps-nburn),n_mcin]

    mcpars = np.load(path+basename+'_mcpars.dat.npy')
    mcnames = np.load(path+basename+'_mcnames.dat.npy')
    flt = np.load(path+basename+'_flt.dat.npy')
    bestpars = np.load(path+basename+'_bestpars.dat.npy')
    mass = np.load(path+basename+'_mass.dat.npy')

    f = np.squeeze(flt.nonzero())

    plt.figure(figsize=(18,15))
    
    nplots = len(f)+norbits+1
    print 'nplots ',nplots

   #auto-generate plot layout
    nrow = np.floor(np.sqrt(nplots)).astype(int)
    if nplots % nrow == 0:
        ncol = nrow
    else:
        ncol = nrow + 1 
    
    for i in range(f.size):
        #print mcnames[f[i]]
        ax = plt.subplot(nrow,ncol, i+1)
              
        t = ax.yaxis.get_offset_text()
        #t.set_size(40) #NOPE
        n, bins, patches = plt.hist(mcpars[:,f[i]],100)
        plt.axvline(np.median(mcpars[:,f[i]]),linewidth=2)
        plt.axvline(np.percentile(mcpars[:,f[i]],84),linewidth=2)
        plt.axvline(np.percentile(mcpars[:,f[i]],16),linewidth=2)
        ax.set_xlabel(mcnames[f[i]])
        plt.tight_layout()
        junk = ax.get_xticklabels() #why is this formatted like this?
        
        label = [junk[q].get_position()[0] for q in range(len(junk)) ]
        
        ax.set_xticklabels(label,rotation=-45)
    
#    for i in norbits:
#        ax = plt.subplot(nrow,ncol, f.size+i+1)
#        n, bins, patches = plt.hist(mcpars[:,f[i]],100)
    
    plt.savefig('/home/sgettel/Dropbox/cfasgettel/research/harpsn/mass_estimate/'+basename+'_hist.png')
    plt.savefig('/home/sgettel/Dropbox/cfasgettel/research/harpsn/mass_estimate/'+basename+'_hist.pdf')
    
    plt.close()
    
    return junk

def plot_rv_after(targname='K00273',nmod=1000, norbits=1,npoly=0,keck='no',epoch=2.4568478981528e6,home='/home/sgettel/',ttime = np.array([1,0])):
    
    jdb, rv, srv, labels, fwhm, contrast, bis_span, rhk, sig_rhk = rr.process_all(targname,maxrv=-50000,maxsrv=5)
    jdb += 2.4e6 #because process_all gives truncated JDBs
    telvec = np.zeros_like(jdb)
    print jdb.size,' HARPS-N obs'
    
    if keck == 'yes':
        
        sfile = open(home+'Dropbox/cfasgettel/research/keck/'+targname+'_full.dat')
        kjdb, krv, ksrv, mdchi, kcts = np.loadtxt(sfile,unpack=True,usecols=(2,3,4,5,6),skiprows=3)
        print np.mean(ksrv), ' typical Keck error '
            
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
    jdb = jdb - epoch #truncate
    medrv = np.median(rv)
    
    #path = '/pool/vonnegut0/harpsn/mass_estimate/K00273_circ_ecc_1_50002'
    path = '/pool/vonnegut0/harpsn/mass_estimate/K00273_ecc_ecc_1_80002'
    
    mcpars = np.load(path+'_mcpars.dat.npy')
    mcnames = np.load(path+'_mcnames.dat.npy')
    mcbest = np.percentile(mcpars,50, axis=0)
    
    #mcbest = np.copy(gpars) #this is mcbest...
    for i in range(norbits):
        ecc, om0 = toeccom(mcbest[2+i*6],mcbest[3+i*6])
        mcbest[2+i*6] = ecc
        mcbest[3+i*6] = om0*180/np.pi
    
    
    tmod = np.linspace(np.min(jdb),np.max(jdb),nmod)

    model_final = rv_drive_new(mcbest,tmod,norbits,npoly,telvec,ttime)
    res1 = rvnorm - rv_drive_new(mcbest,jdb,norbits,npoly,telvec,ttime)
    
    if npoly > 0:
        parst =  np.copy(mcbest)
        parst[4] = 0.0
        poly = rv_drive_new(parst,tmod,norbits,npoly,telvec,ttime)

    k = np.squeeze(np.where(telvec == 1))
    
    #unphased data, now with residuals!
    plt.figure(1)
    gs = grd.GridSpec(2,1, height_ratios=[3,1])
    ax1 = plt.subplot(gs[0])
    plt.errorbar(jdb,rvnorm,yerr=srv,fmt='bo')
    if keck == 'yes':
        plt.errorbar(jdb[k],rvnorm[k],yerr=srv[k],fmt='go')
    plt.plot(tmod,model_final,'r-')
    #plt.xlabel('Adjusted BJD')
    plt.ylabel('Normalized RV (m/s)')

    ax2 = plt.subplot(gs[1])
    plt.errorbar(jdb,res1,yerr=srv,fmt='bo')
    if keck == 'yes':
        plt.errorbar(jdb[k],res1[k],yerr=srv[k],fmt='go')
    plt.plot(tmod,np.zeros_like(tmod),'k-')
    plt.xlabel('Adjusted BJD')
    plt.ylabel('Residuals (m/s)')
    
    
    plt.savefig(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+'_afterplot.pdf')
    plt.savefig(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+'_afterplot.png')

    plt.close(1)
    
    #phase at 1st period
    pars = mcbest[0:6]  #for model
    pars[6:] = 0 #select first planet only
    
    parst = np.copy(mcbest)
    parst[4] = 0.0 #other planets only
    parst[5] = 0.0
    phase = (jdb - pars[1])/pars[0] % 1.0
    
    rvt = rv_drive_new(parst,jdb,norbits,npoly,telvec,ttime)
    pres = (rvnorm-rvt) - pars[5]
    telvec = np.zeros_like(tmod)
    
    plt.figure(2)
    gs = grd.GridSpec(2,1, height_ratios=[3,1])
    ax1 = plt.subplot(gs[0])
    plt.errorbar(phase, rvnorm-rvt, yerr=srv,fmt='bo')
    plt.errorbar(phase+1, rvnorm-rvt, yerr=srv,fmt='bo',mfc='none')
    plt.errorbar(phase-1, rvnorm-rvt, yerr=srv,fmt='bo',mfc='none')

    if keck == 'yes':
        plt.errorbar(phase[k],rvnorm[k]-rvt[k],yerr=srv[k],fmt='go',)
        plt.errorbar(phase[k]+1,rvnorm[k]-rvt[k],yerr=srv[k],fmt='go',mfc='none')
        plt.errorbar(phase[k]-1,rvnorm[k]-rvt[k],yerr=srv[k],fmt='go',mfc='none')

    plt.plot((tmod - pars[1])/pars[0] % 1.0, rv_drive_new(pars, tmod,1,0,telvec,ttime),'r.')
    plt.plot((tmod - pars[1])/pars[0] % 1.0 +1, rv_drive_new(pars, tmod,1,0,telvec,ttime),'r.')
    plt.plot((tmod - pars[1])/pars[0] % 1.0 -1, rv_drive_new(pars, tmod,1,0,telvec,ttime),'r.')

    #plt.xlabel('Orbital Phase')
    plt.ylabel('Normalized RV (m/s)')
    #plt.plot((tmod - guess))
    plt.xlim([-0.25,1.25])

    ax2 = plt.subplot(gs[1])
    res2 = rvnorm - rvt - rv_drive_new(pars,jdb,1,0,telvec,ttime)
    plt.errorbar(phase,res2,yerr=srv,fmt='bo')
    plt.errorbar(phase+1,res2,yerr=srv,fmt='bo',mfc='none')
    plt.errorbar(phase-1,res2,yerr=srv,fmt='bo',mfc='none')
    if keck == 'yes':
        plt.errorbar(phase[k],res2[k],yerr=srv[k],fmt='go')
        plt.errorbar(phase[k]+1,res2[k],yerr=srv[k],fmt='go',mfc='none')
        plt.errorbar(phase[k]-1,res2[k],yerr=srv[k],fmt='go',mfc='none')
    plt.plot(np.linspace(-0.25,1.25,num=tmod.size),np.zeros_like(tmod),'k-')
    #plt.plot((tmod- pars[1])/pars[0] % 1.0 + 1,np.zeros_like(tmod),'k-')
    #plt.plot((tmod- pars[1])/pars[0] % 1.0 - 1,np.zeros_like(tmod),'k-')
    plt.xlabel('Orbital Phase')
    plt.ylabel('Residuals (m/s)')
    plt.xlim([-0.25,1.25])
    
    
    plt.savefig(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+'_phase_afterplot.pdf')
    plt.savefig(home+'Dropbox/cfasgettel/research/harpsn/mass_estimate/'+targname+'_phase_afterplot.png')
    plt.close(2)
    
#return pres

#to run as ...orbits.py
if __name__ == '__main__':
    orbits_test(webdat='yes',nboot=10)
