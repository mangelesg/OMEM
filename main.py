import numpy as np
import numexpr as ne
from ecosys_constants import *

def eco_model(f_bio, f_phys, nvarstot = 0, ntracers =0, xprev = 0, xcur = 0, ncur=0,\
              tstart = 0, tend = 0, dt = 0, ndays = 0, ndayse=0,\
              x_westi =0, x_easti=0, x_northi=0, x_southi=0, tauini=0,\
              Ui = 0, Vi = 0, wi = 0, hordiff = 0, kvmix = 0,\
              rhokpi =0, rho_adiabi=0, convad=False,\
              no3res=0, sio3res=0, po4res=0, invtau=0,
              tempi = 0, ficei = 0, qswi = 0, \
              feflux_sedi=0, feflux_venti=0, \
              dusti=0, feflux_soli = 0, NHyi =0, NOxi = 0,\
              idx = 0, idy = 0, idz = 0, dzs = 0, z_t = 0, zw=0,latr = 0,\
              mask_zonal = 0, mask_merid = 0, maskt = 0, mask_kmt = 0,  mask_remin=0,\
              Navg = 5, OBC = False,  file_create_spinup = None,\
              restore=False, BIO=False, initial = False):
    
##**************************** Configuration of the Model ****************************

    # Create the variables needed for integration
    times = np.arange(tstart, tend, dt)
    Nt = len(times)
    print(Nt)
    Nz, Ny, Nx = np.shape(maskt)

  
    nprev, ncur, nnext = np.mod([ncur - 1, ncur, ncur + 1], 3)

    fcd = 86400

    deltat = np.mod(times/fcd, ndays) # BEC model is in seconds
    day = deltat.astype(np.int)
    dayp = np.mod(day + 1, ndays)
    deltat -= day
    
    ### extended times for the last period when we need to conect with the next year, january 01
    deltate = np.mod(times/fcd, ndayse) # BEC model is in seconds
    daye = deltate.astype(np.int)
    daype = np.mod(daye + 1, ndayse)
    deltate -= daye


    i = j = k =  np.s_[1:-1] # Index the grids. (m)inus is index-1 and (p)lus is index+1
    ip = jp = kp = np.s_[2:]
    im = jm = km = np.s_[:-2]



##*********************************** Start tracers ***************

    x = np.zeros([nvarstot, Nz, Ny, Nx, 3])
    dxdt_bio = np.zeros([nvarstot, Nz, Ny, Nx])
    
    x[:ntracers,:,:,:, nprev]  = xprev
    x[:ntracers,:,:,:, ncur]  = xcur
    
##**************************** Time-stepping the model****************************
    
    for n in range(Nt):
##**************************** Forcings **********************************************
        deltatn = deltat[n]
        deltatne = deltate[n]
        Uday   = Ui[daye[n],...]
        Udayp  = Ui[daype[n], ...]
#        UWday   = UWi[day[n],...]
#        UWdayp  = UWi[dayp[n], ...]
#        VNday   = VNi[day[n],...]
#        VNdayp  = VNi[dayp[n], ...]
        Vday   = Vi[daye[n],...]
        Vdayp  = Vi[daype[n], ...]
        Wday   = wi[daye[n], ...]
        Wdayp  = wi[daype[n], ...]

        U      = ne.evaluate('(1 - deltatne) * Uday  + deltatne * Udayp')
#        UW      = ne.evaluate('(1 - deltatn) * UWday  + deltatn * UWdayp')
#        VN      = ne.evaluate('(1 - deltatn) * VNday  + deltatn * VNdayp')
        V      = ne.evaluate('(1 - deltatne) * Vday  + deltatne * Vdayp')
        
        
        
        W       = ne.evaluate('(1 - deltatne) * Wday + deltatne * Wdayp')
        #kvmix   = ((1 - deltat[n]) * kvmixi[day[n], ...]    + deltat[n] * kvmixi[dayp[n], ...])
        
        surf_flux = np.zeros((ntracers,Ny,Nx))
        
#            ## time-dependant boundary conditions
        if OBC == True:
            x_east = (1-deltatn)*x_easti[:,day[n], :, :] + deltatn*x_easti[:,dayp[n], :, :]
            x_west = (1-deltatn)*x_westi[:,day[n], :, :] + deltatn*x_westi[:,dayp[n], :, :]

            x_north = (1-deltatn)*x_northi[:,day[n], :, :] + deltatn*x_northi[:,dayp[n], :, :]
            x_south = (1-deltatn)*x_southi[:,day[n], :, :] + deltatn*x_southi[:,dayp[n], :, :]

        if BIO == True:
            qswday = qswi[daye[n],...]
            qswdayp = qswi[daype[n],...]
            tempday = tempi[daye[n],...]
            tempdayp = tempi[daype[n],...]
            feflux_solday = feflux_soli[day[n],...]
            feflux_soldayp = feflux_soli[dayp[n],...]
            dustday = dusti[day[n],...]
            dustdayp = dusti[dayp[n],...]
            NOxday = NOxi[day[n],...]
            NOxdayp = NOxi[dayp[n],...]
            NHyday = NHyi[day[n],...]
            NHydayp = NHyi[dayp[n],...]

            
            qsw        = ne.evaluate('(1 - deltatne) * qswday + deltatne * qswdayp')
            temp       = ne.evaluate('(1 - deltatne) * tempday + deltatne * tempdayp')
            feflux_sol = ne.evaluate('(1 - deltatn) * feflux_solday + deltatn * feflux_soldayp')
            dust       = ne.evaluate('(1 - deltatn) * dustday + deltatn * dustdayp')
            NOx        = ne.evaluate('(1 - deltatn) * NOxday + deltatn * NOxdayp')
            NHy        = ne.evaluate('(1 - deltatn) * NHyday + deltatn * NHydayp')
            

            surf_flux[no3_ind,:,:]  = NOx
            surf_flux[nh4_ind,:,:]  = NHy
            surf_flux[po4_ind,:,:]  = (dust * (0.00105 *  0.15 / 30.974 * 1.0e9))
            surf_flux[sio3_ind,:,:] = (dust * (  0.308 * 0.075 / 28.085 * 1.0e9))
            surf_flux[fe_ind,:,:]   = feflux_sol
            
            
        if convad == True:
            rho_adiabday = rho_adiabi[day[n],...]
            rho_adiabdayp = rho_adiabi[dayp[n],...]
            rhokpday = rhokpi[day[n],...]
            rhokpdayp = rhokpi[dayp[n],...]
            rho_adiab  = ne.evaluate('(1 - deltatn) * rho_adiabday + deltatn * rho_adiabdayp')
            rhokp      = ne.evaluate('(1 - deltatn) * rhokpday + deltatn * rhokpdayp')


##**************************** First time-step Euler (to start the spinup) ****************************
        ## if file_spinup is differnt from 'none', it means we are creating a spinup
        xnvarsjinprev = x[:ntracers, :, j, i, nprev]
        masktji = maskt[:,j,i]
        
        
        if initial == True and n == 0:


            dxdt_phys = f_phys(x=x[:ntracers,:,:,:,:], U=U, V=V, w=W, hordiff=hordiff, kvmix=kvmix, \
                mask_zonal=mask_zonal,maskt=maskt,mask_merid=mask_merid,mask_kmt=mask_kmt, \
                idx = idx, idy = idy, idz=idz, dz=dzs, nprev=0, ncur=0,surf_flux=surf_flux[:], init=True)
                                 
            x[:ntracers,:,j, i,ncur]  =  ne.evaluate('xnvarsjinprev + dt*dxdt_phys')

##**************************** Integrate Physics **********************************************

        

        dxdt_phys = f_phys(x=x[:ntracers,:,:,:,:], U=U, V=V, w=W, hordiff=hordiff, kvmix=kvmix, \
                mask_zonal=mask_zonal,maskt=maskt,mask_merid=mask_merid,mask_kmt=mask_kmt, \
                idx = idx, idy = idy, idz=idz, dz=dzs, nprev=nprev, ncur=ncur,surf_flux=surf_flux[:], init=False)

        if BIO == False: # if BIO true, the physics are added below, with the biology
            x[:ntracers, :, j, i, nnext] = ne.evaluate('xnvarsjinprev + 2*dt*dxdt_phys')
    

        if convad == True:
            for tt in range(2):
                keven = np.asarray(np.arange(0,Nz-1,2))
                kodd = np.asarray(np.arange(1,Nz-1,2))
#                dztxcel = dzs
#                dzwxcel = dzs
                for kk in keven:

                    x[:ntracers, kk, j, i, nnext]  = np.ma.where(np.logical_and(rho_adiab[kk,j,i] > rhokp[kk,j,i], maskt[kk+1,j,i] !=0),  (1/(dzs[kk,j,i]+dzs[kk+1,j,i]))*(dzs[kk,j,i]*x[:ntracers, kk, j, i, nnext] + dzs[kk+1,j,i] * x[:ntracers, kk+1, j, i, nnext]), x[:ntracers, kk, j, i, nnext] )

                    x[:ntracers, kk + 1, j, i, nnext]  =  np.ma.where(np.logical_and(rho_adiab[kk,j,i] > rhokp[kk,j,i], maskt[kk+1,j,i] !=0),   x[:ntracers, kk, j, i, nnext],x[:ntracers, kk + 1, j, i, nnext] )

                for kk in kodd:

                    x[:ntracers, kk, j, i, nnext]  = np.ma.where(np.logical_and(rho_adiab[kk,j,i] > rhokp[kk,j,i], maskt[kk+1,j,i] !=0),  (1/(dzs[kk,j,i]+dzs[kk+1,j,i]))*(dzs[kk,j,i]*x[:ntracers, kk, j, i, nnext] + dzs[kk+1,j,i] * x[:ntracers, kk+1, j, i, nnext]), x[:ntracers, kk, j, i, nnext] )

                    x[:ntracers, kk + 1, j, i, nnext]  = np.ma.where(np.logical_and(rho_adiab[kk,j,i] > rhokp[kk,j,i], maskt[kk+1,j,i] !=0),  x[:ntracers, kk, j, i, nnext],x[:ntracers, kk + 1, j, i, nnext])


###**************************** Integrate Ecosystem **********************************************
    

        ## bec model has 22 variables that are ecosystem, adn from 22 to 33, is just fluxes of particulate material
        ## the dxdt_bio returns all the tendencies, plus the fluxes. For the fluxes there is no time derivative, is just
        ## "how much" we exported and imported.
        x[ntracers:,:,:,:,nnext] = 0
        if BIO == True:
            PAR_out = np.zeros((Ny,Nx))
            for kk in range(Nz):
                
                dxdt_bio[:,kk,j,i],PAR_out =  f_bio(kk=kk, nvarstot=nvarstot, ntracers=ntracers,  x = 0.5*(x[:ntracers,kk,j,i,ncur]+x[:ntracers,kk,j,i,nprev]), xprts = x[ntracers:,kk-1,j,i,nnext],
                       temp=temp[kk,j,i], qsw=qsw[j,i], PAR_in = PAR_out,\
                       no3res = no3res[kk,j,i], sio3res = sio3res[kk,j,i], po4res = po4res[kk,j,i], invtau = invtau[kk,j,i],\
                       feflux_sed = feflux_sedi[kk,j,i], dust = dust[j,i],\
                       dzs=dzs[kk,j,i], z_t=z_t[kk], zw=zw[:], idz =idz[kk,j,i], dt = dt, mask_remin= mask_remin[kk,j,i], maskt= maskt[kk,j,i], restore=restore)
                       
#                masktji = mask_kmt[kk,j,i]

#                xnvarskknnext = x[:ntracers, kk, j,i, nnext]
#                dxdt_bionvarskk = dxdt_bio[:ntracers,kk,j,i]
#                dxdt_particulateskk =dxdt_bio[ntracers:,kk,j,i]
#                x[:ntracers, kk, j,i, nnext] = ne.evaluate('xnvarskknnext + 2*dt*masktji*dxdt_bionvarskk')
#                x[ntracers:, kk, j,i, nnext] = ne.evaluate('dxdt_particulateskk*masktji') ## because in marbl code, there is no flux at bottom

            dxdt_bionvars = dxdt_bio[:ntracers,:,j,i]
            dxdt_particulates =dxdt_bio[ntracers:,:,j,i]
            masktji = mask_kmt[:,j,i]

            x[:ntracers, :, j,i, nnext] = ne.evaluate('xnvarsjinprev + 2*dt*masktji*(dxdt_bionvars + dxdt_phys)')
            x[ntracers:, :, j,i, nnext] = ne.evaluate('dxdt_particulates*masktji') ## because in marbl code, there is no flux at bottom

#
##### OBC nudging layer, we dont apply it to the particulated terms
###
#        x[:ntracers,:,:,0,nnext] = x[:ntracers,:,:,0,ncur] + dt*tauf*(x_west - x[:ntracers,:,:,0,ncur])
#        x[:ntracers,:,:,-1,nnext] = x[:ntracers,:,:,-1,ncur] + dt*tauf*(x_east - x[:ntracers,:,:,-1,ncur])
#        x[:ntracers,:,0,:,nnext] = x[:ntracers,:,0,:,ncur] + dt*tauf*(x_south - x[:ntracers,:,0,:,ncur])
#        x[:ntracers,:,-1,:,nnext] = x[:ntracers,:,-1,:,ncur] + dt*tauf*(x_north - x[:ntracers,:,-1,:,ncur])

        if OBC == True:
            tauin = tauini*86400
            tauout = 500*86400
            tauoutE = tauout
            tauoutW = tauout
            tauoutN = tauout
            tauoutS = tauout
####### ***************************************Calculate variable Tau_out as Chen et al., 2013***************************************
###IMPORTANT::::
###IMPORTANT::::
###IMPORTANT:::: For variable tau_out, we need to set a limit when to use this, because somettimes goes to zero. For example:
####            x[:ntracers,:,slj,-1,nnext] = np.ma.where(np.ma.logical_and(Usljim1 >0, tauoutE[:,:,slj]>dt), ne.evaluate('maskkmtsljim1*(xsljim1ncur - (Usljim1*dt*idxsljim1)*(xsljim1ncur-xsljim2ncur) - dt/tauoutEslj *(xsljim1ncur - x_eastslj))'),x[:ntracers,:,slj,-1,nnext])


#            xjim2ncur = x[:ntracers,:,:,-2,ncur]*maskt[:,:,-2]
#            xjim2nprev = x[:ntracers,:,:,-2,nprev]*maskt[:,:,-2]
#            xjim3ncur = x[:ntracers,:,:,-3,ncur]*maskt[:,:,-3]
#            idxjim2 = idx[:,-2]
#            Ujim2 = U[:,:,-2]*mask_zonal[:,:,-2]
#            #try:
#            tauoutE = ne.evaluate('-(xjim2ncur- x_east)/((xjim2ncur-xjim2nprev)/dt + Ujim2*(xjim2ncur-xjim3ncur)*idxjim2)')
#
#    #        except ZeroDivisionError:
#    #            tauoutE = 0
#            #tauoutE = np.nan_to_num(tauoutE*maskt[:,:,-2])
#            tauoutE = np.maximum(0,tauoutE)
#
#
#            xjip1ncur = x[:ntracers,:,:,1,ncur]
#            xjip2ncur = x[:ntracers,:,:,2,ncur]
#            xjip1nprev = x[:ntracers,:,:,1,nprev]
#            Uji0 = U[:,:,0]
#            idxjip1 = idx[:,1]
#
#            tauoutW = ne.evaluate('-(xjip1ncur - x_west)/((xjip1ncur-xjip1nprev)/dt + Uji0*(xjip2ncur-xjip1ncur)*idxjip1)')
#
#
#            tauoutW = np.maximum(0,tauoutW)
#
#            xjm2incur = x[:ntracers,:,-2,:,ncur]
#            xjm3incur = x[:ntracers,:,-3,:,ncur]
#            xjm2inprev = x[:ntracers,:,-2,:,nprev]
#            Vjm2i = V[:,-2,:]
#
#            tauoutN = ne.evaluate('-(xjm2incur - x_north)/((xjm2incur-xjm2inprev)/dt + Vjm2i*(xjm2incur-xjm3incur)*idy)')
#
#            tauoutN = np.maximum(0,tauoutN)
#
#            xjp1incur = x[:ntracers,:,1,:,ncur]
#            xjp1inprev = x[:ntracers,:,1,:,nprev]
#            xjp2incur = x[:ntracers,:,2,:,ncur]
#            Vj0i = V[:,0,:]
#
#            tauoutS = ne.evaluate('-(xjp1incur  - x_south)/((xjp1incur - xjp1inprev)/dt + Vj0i*(xjp2incur-xjp1incur )*idy)')
#            tauoutS = np.maximum(0,tauoutS)
#


##### *********************************************************************************************************************

        ### ***************************** eastern boundary **********************************



            #cx = np.zeros_like(x[:ntracers,:,:,-2,nprev])
    #            cx= np.ma.where(np.abs((x[:ntracers,:,:,-2,nprev]-x[:ntracers,:,:,-3,nprev]))>1e-10, -1/(dt*idx[:,-2])*(x[:ntracers,:,:,-2,ncur]-x[:ntracers,:,:,-2,nprev])/(x[:ntracers,:,:,-2,nprev]-x[:ntracers,:,:,-3,nprev]),cx)
    #            cx = np.maximum(0,cx)
    #            cx = np.minimum(cx, 1/(dt*idx[:,1]))


            # inward
            slj = np.s_[1:]
            xsljim1ncur = x[:ntracers,:,slj,-1,ncur]
            x_eastslj = x_east[:,:,slj]
            maskkmtsljim1 = mask_kmt[:,slj,-1]
            # inward
            x[:ntracers,:,slj,-1,nnext] =  ne.evaluate('maskkmtsljim1*(xsljim1ncur  - dt/tauin*(xsljim1ncur - x_eastslj))')

            # outward
            Usljim1 = U[:,slj,-1]
            idxsljim1 = idx[slj,-1]
            xsljim2ncur = x[:ntracers,:,slj,-2,ncur]

#            x[:ntracers,:,slj,-1,nnext] = np.ma.where(Usljim1 >0, ne.evaluate('maskkmtsljim1*(xsljim1ncur  - (Usljim1*dt*idxsljim1)*(xsljim1ncur - xsljim2ncur))'),x[:ntracers,:,slj,-1,nnext])

            tauoutEslj = tauoutE#[:,:,slj]

            x[:ntracers,:,slj,-1,nnext] = np.ma.where(Usljim1 >0, ne.evaluate('maskkmtsljim1*(xsljim1ncur - (Usljim1*dt*idxsljim1)*(xsljim1ncur-xsljim2ncur) - dt/tauoutEslj *(xsljim1ncur - x_eastslj))'),x[:ntracers,:,slj,-1,nnext])

     
     
     
    #### ***************************** western boundary **********************************
        #    cx = np.zeros_like(x[:ntracers,:,:,1,nprev])
    #            cx = np.ma.where(np.abs((x[:ntracers,:,:,2,nprev]-x[:ntracers,:,:,1,nprev]))>1e-10,-1/(dt*idx[:,1])*(x[:ntracers,:,:,1,ncur]-x[:ntracers,:,:,1,nprev])/(x[:ntracers,:,:,2,nprev]-x[:ntracers,:,:,1,nprev]),cx)
    #            cx = np.minimum(0,cx)
    #            cx = np.maximum(cx, -1/(dt*idx[:,-2]))

            # inward
            slj = np.s_[:-1]


            xslji0ncur = x[:ntracers,:,slj,0,ncur]
            xslji1ncur = x[:ntracers,:,slj,1,ncur]
            x_westslj = x_west[:,:,slj]
            maskkmtslji0 = mask_kmt[:,slj,0]

            x[:ntracers,:,slj,0,nnext] =  ne.evaluate('maskkmtslji0*(xslji0ncur - dt/tauin*(xslji0ncur - x_westslj))')
            idxslji0 = idx[slj,0]
            Uslji0 = U[:,slj,0]
            tauoutWslj = tauoutW#[:,:,slj]
            # outward
#            x[:ntracers,:,slj,0,nnext] = np.ma.where(U[:,slj,0]<0, ne.evaluate('maskkmtslji0*(xslji0ncur  - ((Uslji0)*dt*idxslji0)*(xslji1ncur -xslji0ncur ))') ,x[:ntracers,:,slj,0,nnext])

            x[:ntracers,:,slj,0,nnext] = np.ma.where(U[:,slj,0] < 0, ne.evaluate('maskkmtslji0*(xslji0ncur - (Uslji0*dt*idxslji0)*(xslji1ncur-xslji0ncur) - dt/tauoutWslj*(xslji0ncur - x_westslj))'),x[:ntracers,:,slj,0,nnext])#*maskt[:,slj,0]

          #  print('west',np.max(x[:ntracers,:,slj,0,nnext]))

     ##### ***************************** northern boundary **********************************


            #cy = np.zeros_like(x[:ntracers,:,-2,:,nprev])
    #            cy = np.ma.where(np.abs((x[:ntracers,:,-2,:,nprev]-x[:ntracers,:,-3,:,nprev]))>1e-10,-1/(dt*idy)*(x[:ntracers,:,-2,:,ncur]-x[:ntracers,:,-2,:,nprev])/(x[:ntracers,:,-2,:,nprev]-x[:ntracers,:,-3,:,nprev]),cy)
    #            cy = np.maximum(0,cy)
    #            cy = np.minimum(cy, 1/(dt*idy))


            # inward
            sli = np.s_[:-1]

            xjm1slincur = x[:ntracers,:,-1,sli,ncur]
            xjm2slincur = x[:ntracers,:,-2,sli,ncur]
            x_northsli= x_north[:,:,sli]
            maskkmtjm1sli = mask_kmt[:,-1,sli]

            x[:ntracers,:,-1,sli,nnext] = ne.evaluate('maskkmtjm1sli*(xjm1slincur  - dt/tauin*(xjm1slincur - x_northsli ))')
            Vjm1sli = V[:,-1,sli]
            # outward

#            x[:ntracers,:,-1,sli,nnext] = np.ma.where(V[:,-1,sli]>0, ne.evaluate('maskkmtjm1sli*(xjm1slincur - (Vjm1sli*dt*idy)*(xjm1slincur-xjm2slincur))'), x[:ntracers,:,-1,sli,nnext]  )

            tauoutNsli = tauoutN #[:,:,sli]
            Vjm1sli = V[:,-1,sli]
            x[:ntracers,:,-1,sli,nnext] = np.ma.where(V[:,-1,sli]>0, ne.evaluate('maskkmtjm1sli*(xjm1slincur - (Vjm1sli*dt*idy)*(xjm1slincur-xjm2slincur) - dt/tauoutNsli*(xjm1slincur- x_northsli))'),x[:ntracers,:,-1,sli,nnext] )

     ##### ***************************** southern boundary **********************************

           # cy = np.zeros_like(x[:ntracers,:,1,:,nprev])
    #            cy = np.ma.where(np.abs((x[:ntracers,:,2,:,nprev]-x[:ntracers,:,1,:,nprev]))>1e-10,-1/(dt*idy)*(x[:ntracers,:,1,:,ncur]-x[:ntracers,:,1,:,nprev])/(x[:ntracers,:,2,:,nprev]-x[:ntracers,:,1,:,nprev]),cy)
    #
    #            cy = np.minimum(0,cy)
    #            cy = np.maximum(cy, -1/(dt*idy))

            # inward
            sli = np.s_[1:]

            xj0slincur = x[:ntracers,:,0,sli,ncur]
            xj1slincur = x[:ntracers,:,1,sli,ncur]
            Vj0sli = V[:,0,sli]
            
            
            x_southsli = x_south[:,:,sli]
            maskkmtj0sli = mask_kmt[:,0,sli]
            tauoutSsli = tauoutS#[:,:,sli]
            
            x[:ntracers,:,0,sli,nnext] =  ne.evaluate('maskkmtj0sli*(xj0slincur - dt/tauin*(xj0slincur- x_southsli))')
            # outward
#            x[:ntracers,:,0,sli,nnext] = np.ma.where(V[:,0,sli]<0, ne.evaluate('maskkmtj0sli*(xj0slincur - (Vj0sli*dt*idy)*(xj1slincur -xj0slincur))') , x[:ntracers,:,0,sli,nnext] )

            x[:ntracers,:,0,sli,nnext] = np.ma.where(V[:,0,sli] < 0, ne.evaluate('maskkmtj0sli*(xj0slincur - (Vj0sli*dt*idy)*(xj1slincur -xj0slincur) - dt/tauoutSsli*(xj0slincur- x_southsli))'), x[:ntracers,:,0,sli,nnext] )
           # print('south',np.max(x[:ntracers,:,0,sli,nnext]))



    ####### ******************** just nudging

    #
    #        if OBC == 'nudging':
    #            x[:ntracers,:,:,0,nnext] = x[:ntracers,:,:,0,nnext] + dt*tauf*(x_west - x[:ntracers,:,:,0,ncur])
    #            x[:ntracers,:,:,-1,nnext] = x[:ntracers,:,:,-1,nnext] + dt*tauf*(x_east - x[:ntracers,:,:,-1,ncur])
    #            x[:ntracers,:,0,:,nnext] = x[:ntracers,:,0,:,nnext] + dt*tauf*(x_south - x[:ntracers,:,0,:,ncur])
    #            x[:ntracers,:,-1,:,nnext] = x[:ntracers,:,-1,:,nnext] + dt*tauf*(x_north - x[:ntracers,:,-1,:,ncur])
    #
    #
##**************************** END **********************************************
####**************************** Avergae time step for noise **********************************************
        if not np.mod(times[n], Navg*fcd):
  
            x[..., nnext] = 0.5 * (x[..., nnext] + x[..., ncur])
            x[..., ncur]  = 0.5 * (x[..., nprev] + x[..., ncur])



        
        yield times[n], x[:,:,:,:,nnext]                           # Yield back to the iterator

        nprev, ncur, nnext = ncur, nnext, nprev                     # Cycle the Records

    
    if file_create_spinup is not None:                                     # Save last time step
        np.savez(file_create_spinup, xprev = x[..., nprev], xcur = x[..., ncur], ncur=ncur, curtime=times[n])

