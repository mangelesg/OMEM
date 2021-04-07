### advection-diffusion code
import numpy as np
import numexpr as ne

#def velocities_wall(U=0, V=0, mask_zonal=0, mask_merid=0):
#    i = j = k =  np.s_[1:-1]
#    ip = jp = kp = np.s_[2:]
#    im = jm = km = np.s_[:-2]
#
#
#    Uji = U[:,:, j, i]
#    Ujmi = U[:,:,jm, i]
#    mask_zonalji = mask_zonal[:,j,i]
#    uE = ne.evaluate('(Uji+ Ujmi)*0.5*mask_zonalji')
#
#    Ujim = U[:,:, j, im]
#    Ujmim = U[:,:, jm, im]
#    mask_zonaljim = mask_zonal[:,j,im]
#    uW = ne.evaluate('(Ujim + Ujmim)*0.5*mask_zonaljim')
#
#
#    Vji = V[:,:, j, i]
#    Vjim = V[:,:, j, im]
#    mask_meridji = mask_merid[:,j,i]
#    vN = ne.evaluate('(Vji + Vjim)*0.5*mask_meridji')
#
#    Vjmi= V[:,:, jm, i]
#    Vjmim = V[:,:, jm, im]
#    mask_meridjmi = mask_merid[:,jm,i]
#    vS = ne.evaluate('( Vjmi + Vjmim)*0.5*mask_meridjmi')
#
#    return(uE, uW, vN, vS)


def adv_diff_Barakawa(x=0, U=0, V=0, w=0, mask_zonal = 0,  \
                      mask_merid = 0, maskt=0, mask_kmt=0, idx = 0, idy = 0, idz=0, dz=0,  nprev=0, ncur=0, surf_flux=0, hordiff=0, kvmix=0, init=False):
    
    i = j = k =  np.s_[1:-1]
    ip = jp = kp = np.s_[2:]
    im = jm = km = np.s_[:-2]
    

    zonal = 0
    meridional = 0
    diffx = 0
    diffy = 0
    diffz = 0
    upwell = 0
    

        
 # ********************** advection zonal and meridional

    Uji = U[:, j, i]
    Ujmi = U[:,jm, i]
    mask_zonalji = mask_zonal[:,j,i]
    uE = ne.evaluate('(Uji+ Ujmi)*0.5*mask_zonalji')
   
    Ujim = U[:, j, im]
    Ujmim = U[:, jm, im]
    mask_zonaljim = mask_zonal[:,j,im]
    uW = ne.evaluate('(Ujim + Ujmim)*0.5*mask_zonaljim')
    
    
    Vji = V[:, j, i]
    Vjim = V[:, j, im]
    mask_meridji = mask_merid[:,j,i]
    vN = ne.evaluate('(Vji + Vjim)*0.5*mask_meridji')
    
    Vjmi= V[:, jm, i]
    Vjmim = V[:, jm, im]
    mask_meridjmi = mask_merid[:,jm,i]
    vS = ne.evaluate('( Vjmi + Vjmim)*0.5*mask_meridjmi')
        
    
# ********************** advection zonal and meridional

    
    idxji = idx[j,i]
    
    xjipncur = x[:, :, j, ip, ncur]
    xjincur = x[:, :, j, i, ncur]
    xjimncur = x[:, :, j, im, ncur]
    
    #mask_kmtji = mask_kmt[:,j,i]
    zonal  = ne.evaluate('idxji *(uE*0.5*(xjipncur + xjincur) - uW*0.5*(xjimncur + xjincur))')

    xjpincur = x[:, :, jp, i, ncur]
    xjmincur = x[:, :, jm, i, ncur]
    meridional = ne.evaluate('idy*(vN*0.5*(xjpincur + xjincur) - vS*0.5*(xjmincur + xjincur))')

    
# ********************** diffussion zonal and meridional
    

    
    xjipnprev = x[:, :, j, ip, nprev]
    xjinprev = x[:, :, j, i, nprev]
    xjimnprev = x[:, :, j, im, nprev]
    
    idxjim = idx[j,im]
    
    diffx = ne.evaluate('idxji*hordiff*(mask_zonalji*(xjipnprev - xjinprev)*idxji -\
                       mask_zonaljim*( xjinprev - xjimnprev)*idxjim)')

    
    xjpinprev = x[:, :, jp, i, nprev]
    xjminprev = x[:, :, jm, i, nprev]
    diffy = ne.evaluate('idy*idy *hordiff*( mask_meridji*( xjpinprev - xjinprev) -\
                        mask_meridjmi*( xjinprev - xjminprev))')
    

# # ********************** upwelling (z positive upwards)
    wkmji = w[km,j,i]
    w_top = w[k,j,i]
    w_bottom = w[kp,j,i]
    w0ji = w[0,j,i]
    w1ji = w[1,j,i]
    w_topm1 = w[-1,j,i]
    
    idz0ji = idz[0,j,i]
    idzkji = idz[k,j,i]
    idzm1ji = idz[-1,j,i]
    
    #mask_kmtkji = mask_kmt[k,j,i]


    xkmjincur = x[:, km, j, i, ncur]
    xkjincur = x[:, k, j, i, ncur]
    xkpjincur = x[:, kp, j, i, ncur]
    xkm1jincur = x[:, -1, j, i, ncur]
    xkm2jincur = x[:, -2, j, i, ncur]
    x0jincur = x[:, 0, j, i, ncur]
    x1jincur = x[:, 1, j, i, ncur]
    #mask_kmt0ji = mask_kmt[0,j,i]

    
    upwell = np.ma.zeros(np.shape(x[:, :, j, i, ncur]))

    upwell[:,0,:,:] = ne.evaluate('idz0ji*(w0ji*x0jincur - w1ji*0.5*(x0jincur + x1jincur))')

    ## all levels up to kmt level (one before land in depth)
    upwell[:,k,:,:] = ne.evaluate('idzkji*(w_top*(xkmjincur + xkjincur)*0.5 - w_bottom*( xkpjincur + xkjincur)*0.5)')
    
    ## wbottom =0 at ocean floor
    upwell[:,-1,:,:] = ne.evaluate('idzm1ji*(w_topm1*(xkm1jincur + xkm2jincur)*0.5)')



    # # ####********************** vertical diffusvity (z positive upwards)

#
#    kvmixkmji = kvmix[km, j, i]
#    kvmixkji = kvmix[k, j, i]
#    kvmixkpji = kvmix[kp, j, i]
#    kvmix_top = ne.evaluate('(kvmixkmji + kvmixkji)*0.5')      #*mask_kmt[k,j,i]
#    kvmix_bottom = ne.evaluate('(kvmixkpji + kvmixkji)*0.5')    #*mask_kmt[k,j,i]


    dzkmji = dz[km,j, i]
    dzkji  = dz[k, j, i]
    dzkm1ji = dz[-1,j, i]
    dzkm2ji = dz[-2,j, i]
    dzkpji = dz[kp,j, i]
    dz0ji  = dz[0,j, i]
    dz1ji  = dz[1, j, i]
    
    idz_top    = ne.evaluate('1/((dzkmji + dzkji)*0.5)') # at top of cell
    idz_topm1    = ne.evaluate('1/((dzkm1ji + dzkm2ji)*0.5)') # at top of cell
    idz_bottom = ne.evaluate('1/((dzkpji + dzkji)*0.5)') # at bottom of cell
    idz_bottom0 = ne.evaluate('1/(0.5*(dz0ji + dz1ji))')



    #masktkji = maskt[k,j,i]

    xkmjinprev = x[:, km, j, i, nprev]
    xkjinprev = x[:, k, j, i, nprev]
    xkpjinprev = x[:, kp, j, i, nprev]
    


    x0jinprev = x[:, 0, j, i, nprev]
    x1jinprev = x[:, 1, j, i, nprev]

    diffz = np.zeros_like(x[:,:,j,i,ncur])
    surf_flux = surf_flux[:,j,i]
    
    diffz[:,0,:,:] = ne.evaluate('surf_flux - kvmix*idz_bottom0*(x0jinprev - x1jinprev )')

    diffz[:,k,:,:] = ne.evaluate('idzkji*kvmix*(idz_top*(xkmjinprev -  xkjinprev) - idz_bottom*(xkjinprev - xkpjinprev))')
    
    xkm2jinprev = x[:, -2, j, i, nprev]
    xkm1jinprev = x[:, -1, j, i, nprev]
    
    diffz[:,-1,:,:] = ne.evaluate('idz_topm1*kvmix*(idz_topm1*(xkm2jinprev -  xkm1jinprev))')


    mask_kmtji = mask_kmt[:,j,i]
    masktji = maskt[:,j,i]

    return(ne.evaluate('-mask_kmtji*zonal + mask_kmtji*diffx - mask_kmtji*meridional + mask_kmtji*diffy - mask_kmtji*upwell + masktji*diffz' ))
