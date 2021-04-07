import numpy as np
from scipy.interpolate import interp1d
import numexpr as ne
#import seapy ## only for smoothign functions

#### 1.2 Time-Interpolation
def intmonth(var12=0,nvrs=0,nt=0,nz=0,nla=0,nlo=0,ndim=0):
    # patch jan and dec at end and beginning of time series
    
    if ndim == 5:
        var14 = np.zeros((nvrs,14,nz,nla,nlo))
    elif ndim == 4 and nla !=0 and nlo!=0 :
        var14 = np.zeros((14,nz,nla,nlo))
    elif ndim == 4 and nla ==0 and nlo!=0 :
        var14 = np.zeros((nvrs,14,nz,nlo))
    elif ndim == 4 and nla !=0 and nlo ==0 :
        var14 = np.zeros((nvrs,14,nz,nla))
    elif ndim == 3:
        var14 = np.zeros((14,nla,nlo))
    if nvrs != 0:
        var14[:,0,...] = var12[:,-1,...]
        var14[:,1:13,...] = var12[:,:,...]
        var14[:,13,...] = var12[:,0,...]
    else:
        var14[0] = var12[-1]
        var14[1:13] = var12
        var14[13] = var12[0]
    # define initial monthly time and daily times
    tin = np.arange(-365/12,13*365/12,365/12) # 14 months
    tout=np.arange(0,365,1)
    if nvrs!=0:
        in1d = interp1d(tin, var14, axis =1)
    else:
        in1d = interp1d(tin, var14, axis =0)
    return(in1d(tout)[:])

### for input data in monthly resolution
#def intmonth(var12=0,ndim=0):
#    # patch jan and dec at end and beginning of time series
#    if ndim == 4:
#        nt,nd,nla,nlo = np.shape(var12)
#        var14 = np.zeros((14,nd,nla,nlo))
#    if ndim == 3:
#        nt,nla,nlo = np.shape(var12)
#        var14 = np.zeros((14,nla,nlo))
#    var14[0] = var12[-1]
#    var14[1:13] = var12
#    var14[13] = var12[0]
#    # define initial monthly time and daily times
#    tin = np.arange(-365/12,13*365/12,365/12) # 14 months
#    tout=np.arange(0,365,1)
#    in1d = interp1d(tin, var14, axis =0)
#    return(in1d(tout)[:])

## for input data in 5 days resolution
def int5days(var73=0,ndim=0):
    # patch jan and dec at end and beginning of time series
    if ndim == 4:
        nt,nd,nla,nlo = np.shape(var73)
        var77 = np.zeros((77,nd,nla,nlo))
    if ndim == 3:
        nt,nla,nlo = np.shape(var73)
        var77 = np.zeros((77,nla,nlo))
    var77[0:2] = var73[-2:]
    var77[2:75] = var73
    var77[75:] = var73[0:2]
    # define initial monthly time and daily times
    tin = np.arange(-10,365+10,5) # 14 months
    tout=np.arange(0,365,1)
    in1d = interp1d(tin, var77, axis =0)
    return(in1d(tout)[:])

## for input data in daily resolution
def int1day(var365=0,ndim=0):
    # patch jan and dec at end and beginning of time series
    if ndim == 4:
        nt,nd,nla,nlo = np.shape(var365)
        var385 = np.zeros((385,nd,nla,nlo))
    if ndim == 3:
        nt,nla,nlo = np.shape(var365)
        var385 = np.zeros((385,nla,nlo))
    var385[0:10] = var365[-10:]
    var385[10:375] = var365[:]
    var385[375:] = var365[0:10]
    # define initial monthly time and daily times
    tin = np.arange(-10,365+10,1) # 14 months
    tout=np.arange(0,365,1)
    in1d = interp1d(tin, var385, axis =0)
    return(in1d(tout)[:])
# mask everything where missing values are
#9.969209968386869e+36 
def mask_missing(var, zeroes=True):
    var = np.ma.masked_where(var>8e+20,var)
    if zeroes == True:
        var = var.filled(fill_value=0)
    return(var) # instead of nan, has 0
# to create latitudinal weights
def x_weights(latr, dy, Nx, Ny):
    latrad = np.deg2rad(latr) 
    weights = np.cos(latrad)
    dxs = np.zeros([Ny, Nx])
    for i in range(Ny):
        dxs[i,:] = dy*weights[i] # in cm
    return(dxs)
#def smooth3D(Var):
#    nt,nz,nla,nlo = np.shape(Var)
#    Varsm = np.zeros((nt,nz,nla,nlo))
#    for i in range(nt):
#        for k in range(nz):
#            Varsm[i,k,:,:] = seapy.smooth(Var[i,k,:,:], ksize=10)
#    ## now we put zeros in land, this step is necesary because the smooth function puts velocities
#    ## in the land points different from zero
#    #Varsm = np.ma.where(mask == 0, 0, Varsm[:,:,:,:])
#    return(Varsm)
#def smooth2D(Var):
#    nt,nla,nlo = np.shape(Var)
#    Varsm = np.zeros((nt,nla,nlo))
#    for i in range(nt):
#        Varsm[i,:,:] = seapy.smooth(Var[i,:,:], ksize=10)
#    ## now we put zeros in land, this step is necesary because the smooth function puts velocities
#    ## in the land points different from zero
#    #Varsm = np.ma.where(mask[0,:,:] == 0, 0, Varsm[:,:,:])
#    return(Varsm)
def mask_missing(var):
    var = np.ma.masked_where(var>8e+30,var)
    return(var.filled(fill_value=0)) # instead of nan, has 0
    
## to create third order upwind scheme factors based on
## https://www.cesm.ucar.edu/models/cesm1.0/cesm/cesmBbrowser/ ADVT_UPWIND3

def advt_upwind3_factors_vertical(dzs):
    c0=0; c1=1; c2=2
    dz = dzs.copy()
    nz,ny,nx = np.shape(dzs)
    dzc = np.zeros((nz+2,ny,nx))
    dzc[0] = dz[0]
    dzc[nz] = dz[nz-1]
    
    zeros = np.zeros((nz,ny,nx))
    talfzp = zeros.copy()
    tbetzp = zeros.copy()
    tgamzp = zeros.copy()
    talfzm = zeros.copy()
    tbetzm = zeros.copy()
    tdelzm = zeros.copy()

    for k in range(1,nz-1):
        talfzp[k-1] =  dz[k]*(c2*dz[k]+dzc[k-1])/ \
              ((dz[k]+dz[k+1])*          \
              (dzc[k-1]+c2*dz[k]+dz[k+1]))
        tbetzp[k-1] =  dz[k+1]*(c2*dz[k]+dzc[k-1])/ \
              ((dz[k]+dz[k+1])*            \
              (dz[k]+dzc[k-1]         ))
        tgamzp[k-1] =  -(dz[k]*dz[k+1])/ \
              ((dz[k]+dzc[k-1])*           \
              (dz[k+1]+dzc[k-1]+c2*dz[k]))

    tbetzp[0] = tbetzp[0] + tgamzp[0]
    tgamzp[0] = c0
    talfzp[nz-1] = c0
    tbetzp[nz-1] = c0
    tgamzp[nz-1] = c0
    
    
    for k in range(1,nz-1):
        talfzm[k-1] =  dz[k]*(c2*dz[k+1]+dzc[k+2])/ \
                  ((dz[k]+dz[k+1])*            \
                  (dz[k+1]+dzc[k+2]            ))
        tbetzm[k-1] =  dz[k+1]*(c2*dz[k+1]+dzc[k+2])/ \
                  ((dz[k]+dz[k+1])*              \
                  (dz[k]+dzc[k+2]+c2*dz[k+1]))
        tdelzm[k-1] =  -(dz[k]*dz[k+1])/ \
                  ((dz[k+1]+dzc[k+2])*         \
                  (dz[k]+dzc[k+2]+c2*dz[k+1]))

    talfzm[nz-2] = talfzm[nz-2] + tdelzm[nz-2]
    tdelzm[nz-2] = c0
    talfzm[nz-1] = c0
    tbetzm[nz-1] = c0
    tdelzm[nz-1] = c0

    return(talfzp, tbetzp, tgamzp, talfzm, tbetzm, tdelzm)


def advt_upwind3_factors_horizontal(nz, dxs, dy, KMT):
    c2=1; c1=1; c0=0
    ny,nx = np.shape(dxs)
    zeros = np.zeros_like(dxs)
    zerosk = np.zeros((nz,ny,nx))
    
    TALFXP = zeros.copy()
    TBETXP = zeros.copy()
    TGAMXP = zeros.copy()
    TALFXM = zeros.copy()
    TBETXM = zeros.copy()
    TDELXM = zeros.copy()
    TALFYP = zeros.copy()
    TBETYP = zeros.copy()
    TGAMYP = zeros.copy()
    TALFYM = zeros.copy()
    TBETYM = zeros.copy()
    TDELYM = zeros.copy()
    workx = zerosk.copy()
    worky = zerosk.copy()
    apx = zerosk.copy()
    bpx = zerosk.copy()
    gpx = zerosk.copy()
    amx = zerosk.copy()
    bmx = zerosk.copy()
    dmx = zerosk.copy()
    apy = zerosk.copy()
    bpy = zerosk.copy()
    gpy = zerosk.copy()
    amy = zerosk.copy()
    bmy = zerosk.copy()
    dmy = zerosk.copy()
    
    i = j =  np.s_[1:-1]
    ip = jp =  np.s_[2:]
    im = jm = np.s_[:-2]
    ipp = jpp = np.s_[3:]
    
    dxc   = dxs[j,i]#DXT(i  ,j,iblock)
    dxcw  = dxs[j,im]#DXT(i-1,j,iblock)
    dxce  = dxs[j,ip]#DXT(i+1,j,iblock)
    
    dxs_exti = np.zeros((ny,nx+1))
    dxs_exti[:,:-1] = dxs[:,:]
    
    dxce2 = dxs_exti[j,ipp]#DXT(i+2,j,iblock)

    TALFXP[j,i] = dxc*(c2*dxc+dxcw)/ \
                     ((dxc+dxce)*(dxcw+c2*dxc+dxce))
    TBETXP[j,i] = dxce*(c2*dxc+dxcw)/ \
                     ((dxc+dxcw)*(        dxc+dxce))
    TGAMXP[j,i] =     -(   dxc*dxce)/ \
                     ((dxc+dxcw)*(dxcw+c2*dxc+dxce))

    TALFXM[j,i] = dxc *(c2*dxce+dxce2)/ \
                     ((dxc  +dxce)*(       dxce+dxce2))
    TBETXM[j,i] = dxce*(c2*dxce+dxce2)/ \
                     ((dxc  +dxce)*(dxc+c2*dxce+dxce2))
    TDELXM[j,i] =     -(   dxc *dxce )/ \
                     ((dxce2+dxce)*(dxc+c2*dxce+dxce2))


    dyc   = dy
    dycs  = dy
    dycn  = dy
    dycn2 = dy

    TALFYP[:,:] = dyc *(c2*dyc+dycs)/ \
                     ((dyc+dycn)*(dycs+c2*dyc+dycn))
    TBETYP[:,:] = dycn*(c2*dyc+dycs)/ \
                     ((dyc+dycn)*(dycs+   dyc     ))
    TGAMYP[:,:] =     -(   dyc*dycn)/ \
                     ((dyc+dycs)*(dycs+c2*dyc+dycn))

    TALFYM[:,:] = dyc *(c2*dycn+dycn2)/ \
                     ((dyc+dycn)*(       dycn+dycn2))
    TBETYM[:,:] = dycn*(c2*dycn+dycn2)/ \
                     ((dyc+dycn)*(dyc+c2*dycn+dycn2))
    TDELYM[:,:] =     -(   dyc *dycn )/ \
                     ((dycn2+dycn)*(dyc+c2*dycn+dycn2))


    KMT_exti = np.zeros((ny,nx+1))
    KMT_extj = np.zeros((ny+1,nx))
    KMT_exti[:,:-1] = KMT[:,:].copy()
    KMT_extj[:-1,:] = KMT[:,:].copy()
    
    KMTE = KMT[j,ip]
    KMTW = KMT[j,im]
    KMTEE = KMT_exti[j,ipp]
    KMTN = KMT[jp,i]
    KMTS = KMT[jm,i]
    KMTNN = KMT_extj[jpp,i]



    for kk in range(nz):
    
        
        workx[kk,:,:]  = ne.evaluate('TBETXP + TALFXP')
        apx[kk,:,:]     = c0
        
        workx[kk,j,i]  = np.ma.where( kk <= KMTE, TBETXP[j,i], workx[kk,j,i])
        apx[kk,j,i]  = np.ma.where( kk <= KMTE, TALFXP[j,i], apx[kk,j,i])
     
        bpx[kk,:,:]    = workx[kk,:,:] + TGAMXP
        gpx[kk,:,:]    = c0
     
        bpx[kk,j,i]  = np.ma.where(kk <= KMTW, workx[kk,j,i], bpx[kk,j,i])
        gpx[kk,j,i]  = np.ma.where(kk <= KMTW, TGAMXP[j,i], gpx[kk,j,i])

        amx[kk,:,:]    = ne.evaluate('TALFXM + TDELXM')
        dmx[kk,:,:]    = c0
        
        amx[kk,j,i]    = np.ma.where(kk <= KMTEE, TALFXM[j,i], amx[kk,j,i])
        dmx[kk,j,i]    = np.ma.where(kk <= KMTEE, TDELXM[j,i], dmx[kk,j,i])

        bmx[kk,:,:]    = TBETXM.copy()
    
        worky[kk,:,:]   = ne.evaluate('TBETYP + TALFYP')
        apy[kk,:,:]     = c0
        
        worky[kk,j,i]   = np.ma.where(kk <= KMTN, TBETYP[j,i], worky[kk,j,i])
        apy[kk,j,i]     = np.ma.where(kk <= KMTN, TALFYP[j,i], apy[kk,j,i])

        bpy[kk,:,:]    = worky[kk,:,:] + TGAMYP
        gpy[kk,:,:]    = c0
        
        bpy[kk,j,i]    = np.ma.where(kk <= KMTS, worky[kk,j,i], bpy[kk,j,i])
        gpy[kk,j,i]    = np.ma.where(kk <= KMTS, TGAMYP[j,i], gpy[kk,j,i])


        amy[kk,:,:]    = ne.evaluate('TALFYM + TDELYM')
        dmy[kk,:,:]    = c0
        
        amy[kk,j,i]    = np.ma.where(kk <= KMTNN, TALFYM[j,i], amy[kk,j,i])
        dmy[kk,j,i]    = np.ma.where(kk <= KMTNN, TDELYM[j,i], dmy[kk,j,i])

        bmy[kk,:,:]    = TBETYM
        
        
        
        
        
        
    return(apx,bpx,gpx,amx,bmx,dmx,apy,bpy,gpy,amy,bmy,dmy)
