### ecosystem model from cesm
import numpy as np
from ecosys_constants import *
import numexpr as ne
from decimal import *
getcontext().prec = 14

# input z_t and dz in cm, temp in C, par in? feflux in?
def dxdt_ecosys(kk=0, ntracers=0,  nvarstot=0, x=0, xprts=0, z_t=0, zw =0, dzs=0, idz=0, dt=0, temp=0, qsw=0, PAR_in=0, dust=0,\
 feflux_sed =0, no3res=0, sio3res=0, po4res=0,invtau=0, mask_remin = 0, maskt=0, restore=False):
    
    c0 = 0
    frem = 8e9


    Ny, Nx = np.shape(x[po4_ind])
   
    dxdt = np.zeros((nvarstot,Ny,Nx))
    #zeros = np.zeros((Ny,Nx))


#    z_t is zt from marbl, dzs is delta_z from marbl, but zw is the bottom of the cell
#   while z_w_top in pop is the top of the cell

    zw_all = zw[:]
    zw = zw[kk]
    
    
    p_remin_scalef = c1 ## UNLESS given byinput file
    grazing_type = 1 #  1 for michaelis-menten, 2 for sigmoidal
#!-----------------------------------------------------------------------
#!  create local copies of model tracers
#!  treat negative values as zero
#! in BEC they use x = 0.5(x_old + x_current)
#!-----------------------------------------------------------------------

    x = np.maximum(x,0)
    xprts = np.maximum(xprts,0)

    #O2 = x[o2_ind]
    # DIC = x[dic_ind]
    # ALK = x[alk_ind]
    
    PO4 = x[po4_ind]
    NO3 = x[no3_ind]
    SiO3 = x[sio3_ind]
    NH4 = x[nh4_ind]
    Fe = x[fe_ind]
    Lig = x[lig_ind]
    
    DOC = x[doc_ind]
    DON = x[don_ind]
    DOP = x[dop_ind]
    DOCr = x[docr_ind]
    DONr = x[donr_ind]
    DOPr = x[dopr_ind]
    
    zooC = x[zooC_ind]
        
    spC = x[spC_ind]
    spChl = x[spChl_ind]
    spCaCO3 = x[spCaCO3_ind]
    spFe = x[spFe_ind]
    spP = x[spP_ind]
    
    diatC = x[diatC_ind]
    diatChl = x[diatChl_ind]
    diatSi = x[diatSi_ind]
    diatFe = x[diatFe_ind]
    diatP = x[diatP_ind]
    
    diazC = x[diazC_ind]
    diazChl = x[diazChl_ind]
    diazFe = x[diazFe_ind]
    diazP = x[diazP_ind]
            
##*************************** Calculate restoring HERE, but add at the end ************************
    restore_no3 = 0
    restore_po4 = 0
    restore_sio3 = 0
    if restore == True:
        restore_no3 = ne.evaluate('(no3res - NO3) * invtau')
    
        restore_po4 = ne.evaluate('(po4res - PO4) * invtau')
    
        restore_sio3 = ne.evaluate('(sio3res - SiO3) * invtau')

##*************************** zero consistency for auothrophs *************************************
        

    # If any phyto box are zero, set others to zeros.

    spC = np.ma.where(np.logical_or(np.logical_or(np.logical_or(spC == c0, spChl == c0) , spFe == c0), spP == c0), c0, spC)
    spChl = np.ma.where(np.logical_or(np.logical_or(np.logical_or(spC == c0, spChl == c0) , spFe == c0), spP == c0), c0, spChl)
    spFe = np.ma.where(np.logical_or(np.logical_or(np.logical_or(spC == c0, spChl == c0) , spFe == c0), spP == c0), c0, spFe)
    spP = np.ma.where(np.logical_or(np.logical_or(np.logical_or(spC == c0, spChl == c0) , spFe == c0), spP == c0), c0, spP)
    spCaCO3 = np.ma.where(np.logical_or(np.logical_or(np.logical_or(spC == c0, spChl == c0) , spFe == c0), spP == c0), c0, spCaCO3)

    diatC = np.ma.where(np.logical_or(np.logical_or(np.logical_or(np.logical_or(diatC == c0, diatChl == c0) , diatFe == c0), diatSi==0), diatP==c0), c0, diatC)
    diatChl = np.ma.where(np.logical_or(np.logical_or(np.logical_or(np.logical_or(diatC == c0, diatChl == c0) , diatFe == c0), diatSi==0), diatP==c0), c0,  diatChl)
    diatFe = np.ma.where(np.logical_or(np.logical_or(np.logical_or(np.logical_or(diatC == c0, diatChl == c0) , diatFe == c0), diatSi==0), diatP==c0), c0,  diatFe)
    diatP = np.ma.where(np.logical_or(np.logical_or(np.logical_or(np.logical_or(diatC == c0, diatChl == c0) , diatFe == c0), diatSi==0), diatP==c0), c0, diatP)
    diatSi = np.ma.where(np.logical_or(np.logical_or(np.logical_or(np.logical_or(diatC == c0, diatChl == c0) , diatFe == c0), diatSi==0), diatP==c0), c0, diatSi)


    diazC = np.ma.where(np.logical_or(np.logical_or(np.logical_or(diazC == c0, diazChl == c0) , diazFe == c0), diazP==c0), c0, diazC)
    diazChl = np.ma.where(np.logical_or(np.logical_or(np.logical_or(diazC == c0, diazChl == c0) , diazFe == c0), diazP==c0), c0,  diazChl)
    diazFe = np.ma.where(np.logical_or(np.logical_or(np.logical_or(diazC == c0, diazChl == c0) , diazFe == c0), diazP==c0), c0,  diazFe)
    diazP = np.ma.where(np.logical_or(np.logical_or(np.logical_or(diazC == c0, diazChl == c0) , diazFe == c0), diazP==c0), c0,  diazP)

    
    

##************************************** particulates
#!-----------------------------------------------------------------------
#!  Hard POC is QA flux and soft POC is excess POC.
#!-----------------------------------------------------------------------
    ## -------------------- INITIALIZE PARTICULATES
    # !-----------------------------------------------------------------------
    # !  parameters, from Armstrong et al. 2000
    # !
    # !  July 2002, length scale for excess POC and bSI modified by temperature
    # !  Value given here is at Tref of 30 deg. C, JKM
    # !-----------------------------------------------------------------------
    # !-----------------------------------------------------------------------
    # !  Initial flux_out is the out from previous run
    # !-----------------------------------------------------------------------

    P_CaCO3_sflux_out = c0
    P_CaCO3_hflux_out = c0
    P_SiO2_sflux_out = c0
    P_SiO2_hflux_out = c0
    
    dust_sflux_out = c0
    dust_hflux_out = c0

    POC_sflux_out = c0
    POC_hflux_out = c0
    
    POP_sflux_out = c0
    POP_hflux_out = c0
    
    P_iron_sflux_out = c0
    P_iron_hflux_out = c0
    
    if kk == 0:
        P_CaCO3_sflux_in = c0
        P_CaCO3_hflux_in = c0
        P_SiO2_sflux_in = c0
        P_SiO2_hflux_in = c0


        POC_sflux_in = c0
        POC_hflux_in = c0
        POP_sflux_in = c0
        POP_hflux_in = c0
        P_iron_sflux_in = c0
        P_iron_hflux_in = c0

        
        dust_sflux_in = ne.evaluate('(c1 - dust_gamma) * dust')
        dust_hflux_in = ne.evaluate('dust_gamma * dust') # dust_sflux_out
        # !-----------------------------------------------------------------------
        # !  Compute initial QA(dust) POC flux deficit.
        # !-----------------------------------------------------------------------
        QA_dust_def = ne.evaluate('dust_rho *(dust_sflux_in + dust_hflux_in)')
        
    else:
        
        
        P_CaCO3_sflux_in = xprts[pic_sind]   # P_CaCO3_sflux_out
        P_CaCO3_hflux_in = xprts[pic_hind]   # P_CaCO3_hflux_out
        P_SiO2_sflux_in = xprts[psio2_sind]   # P_SiO2_sflux_out
        P_SiO2_hflux_in = xprts[psio2_hind]  # P_SiO2_hflux_out

        dust_sflux_in = xprts[pdust_sind]   # dust_hflux_out
        dust_hflux_in = xprts[pdust_hind]   # dust_hflux_out

        POC_sflux_in = xprts[poc_sind]   # POC_sflux_out
        POC_hflux_in = xprts[poc_hind]   # POC_hflux_out
        
        POP_sflux_in = xprts[pop_sind]   # POC_sflux_out
        POP_hflux_in = xprts[pop_hind]   # POC_hflux_out
        
        P_iron_sflux_in = xprts[pfe_sind]   # P_iron_sflux_out
        P_iron_hflux_in = xprts[pfe_hind]   # P_iron_hflux_out
        
        QA_dust_def = xprts[qadustdef_ind]   # P_iron_hflux_out
        
### ************************************** surface fluxs

### ------------------------------- MARBL PAR calculation from QSW ------------------------------------
#
#    !-----------------------------------------------------------------------
#    !  compute PAR related quantities
#    !  Morel, Maritorena, JGR, Vol 106, No. C4, pp 7163--7180, 2001
#    !  0.45   fraction of incoming SW -> PAR (non-dim)
#    !-----------------------------------------------------------------------

    f_qsw_par = 0.45 # from marbl_setting_mod
#
#    ! PAR below this threshold is changed to 0.0
#    ! tendencies from PAR below this threshold are small enough to not affect tracer values
    PAR_threshold = 1.0e-19

#      !-----------------------------------------------------------------------
#      ! compute attenuation coefficient over column
#      !-----------------------------------------------------------------------
    KPARdz = c0
    WORK1 = np.maximum(ne.evaluate('spChl + diatChl + diazChl'), 0.02)
    
    KPARdz = np.ma.where(WORK1 < 0.13224, ne.evaluate('0.000919*(WORK1**0.3536)*dzs'), KPARdz)
    KPARdz= np.ma.where(WORK1 >= 0.13224, ne.evaluate('0.001131*(WORK1**0.4562)*dzs'), KPARdz)

#      !-----------------------------------------------------------------------
#      ! propagate PAR values through column, only on subcolumns with PAR>0
#      ! note that if col_frac is 0, then so is PAR
#      !-----------------------------------------------------------------------

    WORK1 = np.exp(-KPARdz)
    PAR_out = c0
    PAR_avg = c0
    
#      !-----------------------------------------------------------------------
#      ! set depth independent quantities, sub-column fractions and PAR at surface
#      ! ignore provided shortwave where col_frac == 0
#      !-----------------------------------------------------------------------
    if kk == 0:

        PAR_in_top = ne.evaluate('f_qsw_par*qsw') # at level k=0, par is a fraction of qsw
     
        PAR_in_top = np.ma.where(PAR_in_top  < PAR_threshold, 0, PAR_in_top )
        PAR_in = PAR_in_top
        PAR_out = ne.evaluate('PAR_in_top*WORK1')
        PAR_avg = ne.evaluate('PAR_in_top*(1-WORK1)/KPARdz')
        
        
    else:
        PAR_out =  ne.evaluate('PAR_in * WORK1')
        PAR_out =  np.ma.where(PAR_out <PAR_threshold,0,PAR_out)
   

        PAR_avg = ne.evaluate('PAR_in*(1-WORK1)/KPARdz')
        
        

##******************** set elemental ratios *****************************

    # set local variables, with incoming ratios

    thetaC_sp   = ne.evaluate('spChl   / (spC + epsC)')
    thetaC_diat = ne.evaluate('diatChl / (diatC + epsC)')
    thetaC_diaz = ne.evaluate('diazChl / (diazC + epsC)')


    
    Qsi         = ne.evaluate('diatSi  / (diatC + epsC)')
    Qsi = np.minimum(gQsi_max, Qsi)
    
    Qfe_diat    = ne.evaluate('diatFe  / (diatC + epsC)')
    Qfe_sp      = ne.evaluate('spFe    / (spC + epsC)')
    Qfe_diaz    = ne.evaluate('diazFe / (diazC + epsC)')
    
    

    # DETERMINE NEW ELEMENTAL RATIOS FOR GROWTH (NEW BIOMASS)

    gQsi      = gQsi_0
    
    gQfe_diat = gQfe_diat_0
    gQfe_sp   = gQfe_sp_0
    gQfe_diaz = gQfe_diaz_0
    
    #Qp_sp = Qp_sp_fixed
#    Qp_diat = Qp_diat_fixed
#    Qp_diaz = Qp_diaz_fixed


    #gQp_sp = gQp_sp_fixed
#    gQp_diat = gQp_diat_fixed
#    gQp_diaz = gQp_diaz_fixed
#
    ##tthis is actually more complex in MARBL, so let's keep it fixed
   
    Qp_sp      = ne.evaluate('spP    / (spC + epsC)')
    gQp_sp = np.minimum(ne.evaluate('(((PquotaSlope * PO4) + PquotaIntercept) * 0.001)'), PquotaMinNP)
    Qp_sp = np.ma.where(spP <= c0, 0, Qp_sp)
    gQp_sp = np.ma.where(Qp_sp <= c0, 0, gQp_sp)

    Qp_diat     = ne.evaluate('diatP    / (diatC + epsC)')
    gQp_diat = np.minimum(ne.evaluate('(((PquotaSlope * PO4) + PquotaIntercept) * 0.001)'), PquotaMinNP)
    Qp_diat = np.ma.where(diatP <= c0, 0, Qp_diat)
    gQp_diat = np.ma.where(Qp_diat <= c0, 0, gQp_diat)
    
    Qp_diaz     = ne.evaluate('diazP    / (diazC + epsC)')
    gQp_diaz = np.minimum(ne.evaluate('(((PquotaSlope * PO4) + PquotaIntercept) * 0.001)'), PquotaMinNP)
    Qp_diaz = np.ma.where(diazP <= c0, 0, Qp_diaz)
    gQp_diaz = np.ma.where(Qp_diaz <= c0, 0, gQp_diaz)
    
    
##########  VARIABLE RATiOS

    gQfe_sp= np.ma.where(Fe < ne.evaluate('gQ_Fe_kFe_thres * parm_sp_kFe'), np.maximum(ne.evaluate('(gQfe_sp*Fe/ (gQ_Fe_kFe_thres * parm_sp_kFe))'), sp_gQfe_min), gQfe_sp)
    
    gQfe_diat= np.ma.where(Fe < ne.evaluate('gQ_Fe_kFe_thres * parm_diat_kFe'), np.maximum(ne.evaluate('(gQfe_diat*Fe/ (gQ_Fe_kFe_thres * parm_diat_kFe))'), diat_gQfe_min), gQfe_diat)
    
    gQfe_diaz = np.ma.where(Fe < ne.evaluate('gQ_Fe_kFe_thres * parm_diaz_kFe'), np.maximum(ne.evaluate('(gQfe_diaz*Fe/ (gQ_Fe_kFe_thres * parm_diaz_kFe))'), diaz_gQfe_min), gQfe_diaz)
    
    


    gQsi = np.ma.where(np.logical_and(np.logical_and(Fe < ne.evaluate('gQ_Fe_kFe_thres * parm_diat_kFe'), Fe > c0), SiO3 > ne.evaluate('gQ_Si_kSi_thres * parm_diat_kSiO3')), np.minimum(ne.evaluate('(gQsi * gQ_Fe_kFe_thres * parm_diat_kFe / Fe)'), gQsi_max),gQsi)
    
    
    gQsi = np.ma.where(Fe == c0, gQsi_max, gQsi)

    #  !  Modify the initial si/C ratio under low ambient Si conditions
    gQsi = np.ma.where(SiO3 < ne.evaluate('(gQ_Si_kSi_thres * parm_diat_kSiO3)'), np.maximum(ne.evaluate('(gQsi * SiO3/ (gQ_Si_kSi_thres * parm_diat_kSiO3))'), gQsi_min),gQsi)


#        !-----------------------------------------------------------------------
#        !  QCaCO3 is the percentage of sp organic matter which is associated
#        !  with coccolithophores
#        !-----------------------------------------------------------------------

    QCaCO3 = np.minimum(ne.evaluate('spCaCO3 / (spC + epsC)'), QCaCO3_max)
   

###******************************* temperature functional form


    # Tref = 30.0 reference temperature (deg. C)
    # Using q10 formulation with Q10 value of 2.0 (Doney et al., 1996).
    Q_10 = 1.7
    Tfunc = ne.evaluate('Q_10**(((temp + T0_Kelvin)-(Tref + T0_Kelvin))/10)')



####*********************** Pprime for all phyto groups **********
## ************************* compute sp Pprime *************************

   
    f_loss_thres = np.minimum( np.maximum(ne.evaluate('(thres_z2_phyto - z_t)/(thres_z2_phyto - thres_z1_phyto)'), c0)  , c1)
    
    C_loss_thres = ne.evaluate('f_loss_thres * loss_thres_sp')

    C_loss_thres = np.ma.where(temp < temp_thres_sp, ne.evaluate('f_loss_thres * loss_thres2_sp'),C_loss_thres)
    
    Pprime_sp = np.maximum(ne.evaluate('spC - C_loss_thres'), c0)
   
### ********************* diat Pprime
    
    C_loss_thres = ne.evaluate('f_loss_thres * loss_thres_diat')

    C_loss_thres = np.ma.where(temp < temp_thres_diat, ne.evaluate('f_loss_thres * loss_thres2_diat'), C_loss_thres)
    
    Pprime_diat = np.maximum(ne.evaluate('diatC - C_loss_thres'), c0)

### ********************** diaz Pprime ***************************

    C_loss_diaz = ne.evaluate('f_loss_thres * loss_thres_diaz')
    
    C_loss_diaz = np.ma.where(temp < temp_thres_diaz, ne.evaluate('f_loss_thres * loss_thres_diaz2'), C_loss_diaz)
    
    Pprime_diaz = np.maximum(ne.evaluate('diazC - C_loss_diaz'), c0)
   
 
    workPprime = ne.evaluate('Pprime_sp + Pprime_diat + Pprime_diaz')
 
 ######### ********************* SMALL PHYTOPLANKTON ******************************************
######### ************************************************************************

## ************************* compute sp uptake *************************

    # !-----------------------------------------------------------------------
    # !  Get relative nutrient uptake rates for phytoplankton,
    # !  min. relative uptake rate modifies C fixation in the manner
    # !  that the min. cell quota does in GD98.

    VNO3_sp = ne.evaluate('(NO3 / parm_sp_kNO3) / (c1 + (NO3 / parm_sp_kNO3) + (NH4 / parm_sp_kNH4))')
    VNH4_sp = ne.evaluate('(NH4 / parm_sp_kNH4) / (c1 + (NO3 / parm_sp_kNO3) + (NH4 / parm_sp_kNH4))')
    VNtot_sp = ne.evaluate('VNO3_sp + VNH4_sp')

    # !-----------------------------------------------------------------------
    # !  get relative Fe uptake by phytoplankton
    # !  get relative P uptake rates for phytoplankton
    # !-----------------------------------------------------------------------
    VFe_sp = ne.evaluate('Fe / (Fe + parm_sp_kFe)')
 
    # !-----------------------------------------------------------------------
    # !  Partition diaz P uptake in same manner to no3/nh4 partition
    # !-----------------------------------------------------------------------

    VPO4_sp = ne.evaluate('(PO4 /parm_sp_kPO4)/(c1 + (PO4/parm_sp_kPO4) + (DOP/parm_sp_kDOP))')
    VDOP_sp = ne.evaluate('(DOP /parm_sp_kDOP)/(c1 + (PO4/parm_sp_kPO4) + (DOP/parm_sp_kDOP))')

    VPtot_sp = ne.evaluate('VPO4_sp + VDOP_sp')
    # !-----------------------------------------------------------------------
    # !  Small Phytoplankton C-fixation - given light and relative uptake rates
    # !  determine small phyto nutrient limitation factor for carbon fixation
    # !-----------------------------------------------------------------------
    
    f_nut = np.minimum(VNtot_sp, VFe_sp)
    f_nut = np.minimum(f_nut, VPtot_sp)

### ************************* sp photosytnthesis *****************************

    PCmax = ne.evaluate('PCref_sp * f_nut * Tfunc')
    PCmax = np.ma.where(temp<temp_thres_sp, c0, PCmax)
    
    photoacc_sp = c0
    light_lim = c0
    PCphoto_sp = c0
    photoC_sp = c0
    pChl = c0
    
    WORK1 = ne.evaluate('parm_alphaChlsp * thetaC_sp * PAR_avg')
        
    light_lim = np.ma.where(np.logical_and(WORK1>c0, thetaC_sp>c0), (c1 - np.exp((-c1 * WORK1) / (PCmax + epsTinv))), light_lim)
    
    PCphoto_sp = np.ma.where(np.logical_and(WORK1>c0, thetaC_sp>c0), ne.evaluate('PCmax * light_lim'), PCphoto_sp)
       
    # !  calculate pChl, (used in photoadapt., GD98)
    # !  2.3   max value of thetaN (Chl/N ratio) (mg Chl/mmol N)
   
    pChl = np.ma.where(np.logical_and(WORK1>c0, thetaC_sp>c0), ne.evaluate('thetaN_max_sp * PCphoto_sp / (WORK1 )'),pChl)
    
    photoacc_sp = np.ma.where(np.logical_and(WORK1>c0, thetaC_sp>c0), ne.evaluate('(pChl * PCphoto_sp * Q / (thetaC_sp)) * spChl'), photoacc_sp)
  
    photoC_sp = np.ma.where(thetaC_sp>c0, ne.evaluate('PCphoto_sp * spC'), photoC_sp)
    
##***************************** nutrient uptake *****************************
    # !-----------------------------------------------------------------------
    # !  Get nutrient uptakes by small phyto based on calculated C fixation
    # !  total N uptake (VNC_sp) is used in photoadaption
    # !-----------------------------------------------------------------------
    NO3_V_sp = c0
    NH4_V_sp = c0
    PO4_V_sp = c0
    DOP_V_sp = c0

    NO3_V_sp = np.ma.where(VNtot_sp > c0, ne.evaluate('(VNO3_sp / VNtot_sp ) * photoC_sp * Q'),  NO3_V_sp)
    NH4_V_sp = np.ma.where(VNtot_sp > c0, ne.evaluate('(VNH4_sp / VNtot_sp ) * photoC_sp * Q'),  NH4_V_sp)
    PO4_V_sp = np.ma.where(VPtot_sp > c0, ne.evaluate('(VPO4_sp / VPtot_sp ) *photoC_sp * gQp_sp'),PO4_V_sp)
    DOP_V_sp = np.ma.where(VPtot_sp > c0, ne.evaluate('(VDOP_sp / VPtot_sp ) *photoC_sp * gQp_sp'),DOP_V_sp)
    
    
    photoFe_sp = ne.evaluate('photoC_sp * gQfe_sp')

### ************************* sp calcification

    # !-----------------------------------------------------------------------
    # !  CaCO3 Production, parameterized as function of small phyto production
    # !  decrease CaCO3 as function of nutrient limitation decrease CaCO3 prod
    # !  at low temperatures increase CaCO3 prod under bloom conditions
    # !  maximum calcification rate is 40% of primary production
    # !-----------------------------------------------------------------------

    CaCO3_PROD = ne.evaluate('f_prod_sp_CaCO3 * photoC_sp')
    CaCO3_PROD = ne.evaluate('CaCO3_PROD * f_nut * f_nut')
    CaCO3_PROD = np.ma.where(temp < CaCO3_temp_thres1, CaCO3_PROD * \
                np.maximum(temp-CaCO3_temp_thres2, c0)/(CaCO3_temp_thres1-CaCO3_temp_thres2), CaCO3_PROD)
    CaCO3_PROD = np.ma.where(spC > CaCO3_sp_thres, \
            np.minimum(CaCO3_PROD*spC/CaCO3_sp_thres,ne.evaluate('f_photosp_CaCO3*photoC_sp')), CaCO3_PROD)


##************************ sp losses **********************************

    sp_loss = ne.evaluate('sp_mort * Pprime_sp * Tfunc')
    
    sp_agg = np.minimum(ne.evaluate('(sp_agg_rate_max * dps) * Pprime_sp'), ne.evaluate('sp_mort2 * Pprime_sp **auto_mort2_exp'))
    
    sp_agg = np.maximum(ne.evaluate('(sp_agg_rate_min * dps) * Pprime_sp'), sp_agg)
     
    #routing of graze_sp & sp_loss
    # sp_agg all goes to POC
    # min.%C routed from sp_loss = 0.59 * QCaCO3, or P_CaCO3%rho

    sp_loss_poc = ne.evaluate('QCaCO3 * sp_loss')
    
    sp_loss_doc = ne.evaluate('(c1 - parm_labile_ratio) * (sp_loss - sp_loss_poc)')
    sp_loss_dic = ne.evaluate('parm_labile_ratio * (sp_loss - sp_loss_poc)')
    
##*************************** sp grazing **********************************
    ## sigmoidal grazing
    #called auto_graze in marbl
    
    
    graze_sp = c0
    graze_sp = np.ma.where(Pprime_sp>c0,ne.evaluate('parm_sp_z_umax_0*Tfunc* zooC * (Pprime_sp**grazing_type / (Pprime_sp**grazing_type + parm_z_grz**grazing_type))'),graze_sp)

    #graze_sp is named graze_rate in marbl, z_ingest is graze_zoo

    graze_sp_zoo = ne.evaluate('f_sp_graze_zoo * graze_sp')
    
    #graze_sp_zootot = graze_sp_zoo
    
    graze_sp_poc = graze_sp * np.maximum(ne.evaluate('caco3_poc_min * QCaCO3'), \
                np.minimum(ne.evaluate('spc_poc_fac * (Pprime_sp + 0.6)**1.6'),f_graze_sp_poc_lim))
                
    
    graze_sp_doc = ne.evaluate('f_sp_graze_doc * graze_sp')
    
    
    graze_sp_dic = ne.evaluate('graze_sp  - graze_sp_zoo - graze_sp_poc - graze_sp_doc')
    

##************************* sp routing ***********************
    remaining_spP_pop =ne.evaluate('(graze_sp_poc + sp_loss_poc + sp_agg)*Qp_sp')
    
    remaining_spP =  ne.evaluate('((graze_sp + sp_loss + sp_agg) * Qp_sp) - \
    ( graze_sp_zoo * Qp_sp - remaining_spP_pop)')
    
    remaining_spP_pop = np.ma.where(remaining_spP<c0, ne.evaluate('remaining_spP_pop + remaining_spP'), remaining_spP_pop)
    
    remaining_spP = np.maximum(remaining_spP,0)
    
    remaining_spP_dop = ne.evaluate('f_toDOP* remaining_spP')
    
    remaining_spP_dip = ne.evaluate('(c1 - f_toDOP) * remaining_spP')


######### ********************* DIATOMS ******************************************
######### ************************************************************************
    # !-----------------------------------------------------------------------
    # !  0.01 small diatom threshold C concentration (mmol C/m^3)
    # !  get diatom loss(in C units)
    # !  Diatom agg loss, min. 1%/day
    # !  get grazing rate (graze_diat) on diatoms  (in C units)
    # !----------------------------------------------------------------------
    

## ************************* diatoms uptake

    VNO3_diat = ne.evaluate('(NO3 / parm_diat_kNO3) / (c1 + (NO3 / parm_diat_kNO3) + (NH4 / parm_diat_kNH4))')
    VNH4_diat = ne.evaluate('(NH4 / parm_diat_kNH4) / (c1 + (NO3 / parm_diat_kNO3) + (NH4 / parm_diat_kNH4))')
    VNtot_diat = ne.evaluate('VNO3_diat + VNH4_diat')

    # !-----------------------------------------------------------------------
    # !  get relative Fe uptake by diatoms
    # !  get relative P uptake rates for diatoms
    # !  get relative SiO3 uptake rate for diatoms
    # !-----------------------------------------------------------------------

    VFe_diat = ne.evaluate('Fe / (Fe + parm_diat_kFe)')
    
    VPO4_diat = ne.evaluate('(PO4 /parm_diat_kPO4)/(c1 + (PO4/parm_diat_kPO4) + (DOP/parm_diat_kDOP))')
    VDOP_diat = ne.evaluate('(DOP /parm_diat_kDOP)/(c1 + (PO4/parm_diat_kPO4) + (DOP/parm_diat_kDOP))')

    VPtot_diat = ne.evaluate('VPO4_diat + VDOP_diat')
    
    
    VSiO3_diat = ne.evaluate('SiO3 / (SiO3 + parm_diat_kSiO3)')
    
    # !-----------------------------------------------------------------------
    # !  Diatom carbon fixation and photoadapt.
    # !  determine diatom nutrient limitation factor for carbon fixation
    # !-----------------------------------------------------------------------

    f_nut = np.minimum(VNtot_diat, VFe_diat)
    f_nut = np.minimum(f_nut, VPtot_diat)
    f_nut = np.minimum(f_nut, VSiO3_diat)
    
###************************** diatom photosynthesis
    # !-----------------------------------------------------------------------
    # !  get diatom photosynth. rate, phyto C biomass change, photoadapt
    # !-----------------------------------------------------------------------

    PCmax = ne.evaluate('PCref_diat * f_nut * Tfunc')
    PCmax = np.ma.where(temp<temp_thres_diat, c0, PCmax)
    
    photoacc_diat = c0
    light_lim = c0
    PCphoto_diat = c0
    pChl = c0
    photoC_diat = c0
    
    WORK1 = ne.evaluate('parm_alphaChldiat * thetaC_diat * PAR_avg')
    
    
    light_lim = np.ma.where(np.logical_and(WORK1>c0, thetaC_diat>c0), (c1 - np.exp((-c1 * WORK1) / (PCmax + epsTinv))), light_lim)
    
    PCphoto_diat =  np.ma.where(np.logical_and(WORK1>c0, thetaC_diat>c0), ne.evaluate('PCmax * light_lim'), PCphoto_diat)
   
    # !  calculate pChl, (used in photoadapt., GD98)
    # !  2.3   max value of thetaN (Chl/N ratio) (mg Chl/mmol N)
   
    pChl = np.ma.where(np.logical_and(WORK1>c0, thetaC_diat>c0), ne.evaluate('thetaN_max_diat * PCphoto_diat / (WORK1 )'),pChl)
  
    photoacc_diat = np.ma.where(np.logical_and(WORK1>c0, thetaC_diat>c0), ne.evaluate('(pChl * PCphoto_diat * Q / (thetaC_diat)) * diatChl'), photoacc_diat)
  
    photoC_diat = np.ma.where(thetaC_diat>c0, ne.evaluate('PCphoto_diat * diatC'), photoC_diat)
    

##************************* diat nutrient uptake

    # !-----------------------------------------------------------------------
    # !  Get nutrient uptake by diatoms based on C fixation
    # !-----------------------------------------------------------------------

    NO3_V_diat = c0
    NH4_V_diat = c0
    PO4_V_diat = c0
    DOP_V_diat = c0

    NO3_V_diat = np.ma.where(VNtot_diat > c0, ne.evaluate('(VNO3_diat / VNtot_diat ) * photoC_diat * Q'),  NO3_V_diat)
    NH4_V_diat = np.ma.where(VNtot_diat > c0, ne.evaluate('(VNH4_diat / VNtot_diat ) * photoC_diat * Q'),  NH4_V_diat)
    PO4_V_diat = np.ma.where(VPtot_diat > c0, ne.evaluate('(VPO4_diat / VPtot_diat ) *photoC_diat * Qp_diat'), PO4_V_diat)
    DOP_V_diat = np.ma.where(VPtot_diat > c0, ne.evaluate('(VDOP_diat / VPtot_diat ) *photoC_diat * Qp_diat'), DOP_V_diat)
    
    
    photoFe_diat = ne.evaluate('photoC_diat * gQfe_diat')
    photoSi_diat = ne.evaluate('photoC_diat * gQsi')
    

    
### ********************* diatoms losses


    diat_loss = ne.evaluate('diat_mort * Pprime_diat * Tfunc')
    
    diat_agg = np.minimum(ne.evaluate('(diat_agg_rate_max * dps) * Pprime_diat'), ne.evaluate('diat_mort2 * Pprime_diat **auto_mort2_exp'))
    
    diat_agg = np.maximum(ne.evaluate('(diat_agg_rate_min * dps) * Pprime_diat'), diat_agg)


    diat_loss_poc = ne.evaluate('f_diat_loss_poc * diat_loss') # f_diat_loss_poc =0
    
    
    diat_loss_doc = ne.evaluate('(c1-parm_labile_ratio) * (diat_loss - diat_loss_poc)')
    diat_loss_dic = ne.evaluate('parm_labile_ratio * (diat_loss - diat_loss_poc)')

##************************* diat grazing

    graze_diat =  ne.evaluate('parm_diat_z_umax_0*Tfunc* zooC * (Pprime_diat**grazing_type / (Pprime_diat**grazing_type + parm_z_grz**grazing_type))')
    
    graze_diat_zoo = ne.evaluate('z_ingest_diat * graze_diat')
    
    graze_diat_poc = ne.evaluate('f_graze_diat_poc * graze_diat')
    
    graze_diat_doc = ne.evaluate('graze_doc * graze_diat')
    
    graze_diat_dic = ne.evaluate('graze_diat - graze_diat_zoo - graze_diat_poc - graze_diat_doc')

##************************* diat routing

    remaining_diatP_pop =  ne.evaluate('(graze_diat_poc + diat_loss_poc + diat_agg) * Qp_diat')
    
    remaining_diatP = ne.evaluate('(graze_diat + diat_loss + diat_agg)*Qp_diat - graze_diat_zoo*Qp_diat - remaining_diatP_pop')
    
    remaining_diatP_pop = np.ma.where(remaining_diatP<c0, ne.evaluate('remaining_diatP_pop + remaining_diatP'), remaining_diatP_pop)
    
    remaining_diatP = np.maximum(remaining_diatP,0)

    remaining_diatP_dop = ne.evaluate('f_toDOP* remaining_diatP')
    
    remaining_diatP_dip = ne.evaluate('(c1 - f_toDOP) * remaining_diatP')


######### ********************* DIAZOTHROPS ******************************************
######### ************************************************************************
    # !-----------------------------------------------------------------------
    # !  0.03 small diazotroph threshold C concentration (mmol C/m^3)
    # !  Lower value used at temperatures < 16 deg. C, negligible biomass
    # !  get diazotroph loss(in C units)
    # !  get grazing rate (graze_diaz) on diazotrophs  (in C units)
    # !  no aggregation loss for diazotrophs
    # !-----------------------------------------------------------------------

##********************** diazothrops uptake
    # !-----------------------------------------------------------------------
    # !  Relative uptake rates for diazotrophs nitrate is VNO3, ammonium is VNH4
    # !-----------------------------------------------------------------------

    VNO3_diaz =ne.evaluate('(NO3 / parm_diaz_kNO3) / (c1 + (NO3 / parm_diaz_kNO3) + (NH4 / parm_diaz_kNH4))')
    VNH4_diaz = ne.evaluate('(NH4 / parm_diaz_kNH4) / (c1 + (NO3 / parm_diaz_kNO3) + (NH4 / parm_diaz_kNH4))')
   
    VNtot_diaz = c1 # nitrogen fixers

    # !-----------------------------------------------------------------------
    # !  get relative Fe uptake by diazotrophs
    # !  get relative P uptake rates for diazotrophs
    # !-----------------------------------------------------------------------

    VFe_diaz = ne.evaluate('Fe/(Fe + parm_diaz_kFe)')
    
    VPO4_diaz = ne.evaluate('(PO4 /parm_diaz_kPO4)/(c1 + (PO4/parm_diaz_kPO4) + (DOP/parm_diaz_kDOP))')
    VDOP_diaz = ne.evaluate('(DOP/parm_diaz_kDOP)/(c1 + (PO4/parm_diaz_kPO4) + (DOP/parm_diaz_kDOP))')

    VPtot_diaz = ne.evaluate('VPO4_diaz + VDOP_diaz')
    
    f_nut = np.minimum(VNtot_diaz, VFe_diaz)
    f_nut = np.minimum(VPtot_diaz, f_nut)

##***************************** diaz photosynthesis

    # !-----------------------------------------------------------------------
    # !  get diazotroph photosynth. rate, phyto C biomass change
    # !-----------------------------------------------------------------------

    PCmax = ne.evaluate('PCref_diaz * f_nut * Tfunc')
    PCmax = np.ma.where(temp<temp_thres_diaz, c0, PCmax)

    photoacc_diaz = c0
    light_lim = c0
    PCphoto_diaz = c0
    pChl = c0
    photoC_diaz = c0

    WORK1 = ne.evaluate('parm_alphaChldiaz * thetaC_diaz * PAR_avg')

    light_lim = np.ma.where(np.logical_and(WORK1>c0, thetaC_diaz >c0) ,c1 - np.exp((-c1 * WORK1) / (PCmax + epsTinv)), light_lim)
    
    PCphoto_diaz =  np.ma.where(np.logical_and(WORK1>c0, thetaC_diaz >c0),  ne.evaluate('PCmax * light_lim'), PCphoto_diaz)
   
    # !  calculate pChl, (used in photoadapt., GD98)
    # !  2.3   max value of thetaN (Chl/N ratio) (mg Chl/mmol N)

    pChl = np.ma.where(np.logical_and(WORK1>c0, thetaC_diaz >c0),  ne.evaluate('thetaN_max_diaz * PCphoto_diaz / (WORK1 )'),pChl)
  
    photoacc_diaz = np.ma.where(np.logical_and(WORK1>c0, thetaC_diaz >c0),  ne.evaluate('(pChl * PCphoto_diaz * Q / (thetaC_diaz)) * diazChl'), photoacc_diaz)
  
    photoC_diaz = np.ma.where(thetaC_diaz>c0, ne.evaluate('PCphoto_diaz * diazC'), photoC_diaz)

###************************* diaz nutrient uptake

    NO3_V_diaz = c0
    NH4_V_diaz = c0
    PO4_V_diaz = c0
    DOP_V_diaz = c0

    NO3_V_diaz = np.ma.where(VNtot_diaz > c0, ne.evaluate('(VNO3_diaz / VNtot_diaz ) * photoC_diaz * Q'),  NO3_V_diaz)
    NH4_V_diaz = np.ma.where(VNtot_diaz > c0, ne.evaluate('(VNH4_diaz / VNtot_diaz ) * photoC_diaz * Q'),  NH4_V_diaz)
    PO4_V_diaz = np.ma.where(VPtot_diaz > c0, ne.evaluate('(VPO4_diaz / VPtot_diaz ) *photoC_diaz * gQp_diaz'),PO4_V_diaz)
    DOP_V_diaz = np.ma.where(VPtot_diaz > c0, ne.evaluate('(VDOP_diaz / VPtot_diaz ) *photoC_diaz * gQp_diaz'),DOP_V_diaz)
    
    
    photoFe_diaz = ne.evaluate('photoC_diaz * gQfe_diaz')

##*********************** diaz nitrogen fixation
    # !-----------------------------------------------------------------------
    # !  Get N fixation by diazotrophs based on C fixation,
    # !  Diazotrophs fix more than they need then 30% is excreted
    # !-----------------------------------------------------------------------

    Nfix     = ne.evaluate('(photoC_diaz * Q  * r_Nfix_photo) - NO3_V_diaz - NH4_V_diaz')
    Nexcrete = ne.evaluate('Nfix + NO3_V_diaz  + NH4_V_diaz - photoC_diaz * Q')


    
### ********************** diazothrops losses ***************************

    diaz_loss = ne.evaluate('diaz_mort * Pprime_diaz * Tfunc')
    
    diaz_agg = np.minimum(ne.evaluate('(diaz_agg_rate_max * dps) * Pprime_diaz'),ne.evaluate('diaz_mort2 * Pprime_diaz **auto_mort2_exp'))
    
    diaz_agg = np.maximum(ne.evaluate('(diaz_agg_rate_min * dps) * Pprime_diaz'), diaz_agg)

   
    diaz_loss_poc = ne.evaluate('f_diaz_loss_poc * diaz_loss') # = 0
    
    diaz_loss_doc = ne.evaluate('(c1 - parm_labile_ratio) * (diaz_loss - diaz_loss_poc)')
    diaz_loss_dic = ne.evaluate('parm_labile_ratio * (diaz_loss - diaz_loss_poc)')

##************************* diaz grazing

    graze_diaz = ne.evaluate('parm_diaz_z_umax_0* Tfunc * zooC * (Pprime_diaz**grazing_type / (Pprime_diaz**grazing_type + parm_z_grz**grazing_type))')

    # !-----------------------------------------------------------------------
    # !  routing of graze_diaz & diaz_loss
    # !-----------------------------------------------------------------------

    graze_diaz_zoo = ne.evaluate('z_ingest_diaz * graze_diaz')
    
    graze_diaz_poc = ne.evaluate('f_graze_diaz_poc * graze_diaz')
    
    graze_diaz_doc = ne.evaluate('graze_doc * graze_diaz')
    
    graze_diaz_dic = ne.evaluate('graze_diaz - graze_diaz_zoo - graze_diaz_poc - graze_diaz_doc')

    # !-----------------------------------------------------------------------
    # !  Note as diazotrophs have different Qp, we must route enough P into zoopl
    # !  and sinking detritus pools to fill their fixed p/C ratios.  The left over
    # !  P (remaining_diazP) is split between DOP and DIP pools
    # !-----------------------------------------------------------------------

    remaining_diazP_pop =  ne.evaluate('(graze_diaz_poc + diaz_loss_poc + diaz_agg) * Qp_diaz')
    
    remaining_diazP = ne.evaluate('(graze_diaz + diaz_loss + diaz_agg)*Qp_diaz - graze_diaz_zoo*Qp_diaz - remaining_diazP_pop')
    
    remaining_diazP_pop = np.ma.where(remaining_diazP<c0, ne.evaluate('remaining_diazP_pop + remaining_diazP'), remaining_diazP_pop)
    
    remaining_diazP = np.maximum(remaining_diazP,0)
    
    remaining_diazP_dop = ne.evaluate('f_toDOP* remaining_diazP')
    
    remaining_diazP_dip = ne.evaluate('(c1 - f_toDOP) * remaining_diazP')
    

    

######### ********************* ZOOPLANKTON ******************************************
######### ************************************************************************

    # !-----------------------------------------------------------------------
    # !  get fractional factor for routing of zoo losses, based on food supply
    # !  more material is routed to large detrital pool when diatoms eaten
    # !-----------------------------------------------------------------------

    f_zoo_detr = ne.evaluate('(f_diat_zoo_detr * (graze_diat + epsC * epsTinv) + \
                  f_sp_zoo_detr * (graze_sp + epsC * epsTinv) + \
                  f_diaz_zoo_detr * (graze_diaz + epsC * epsTinv)) / \
                (graze_diat + graze_sp + graze_diaz + 3 * epsC * epsTinv)')

    # !-----------------------------------------------------------------------
    # !  0.01 small zoo threshold C concentration (mmol C/m^3)
    # !  zoo losses
    # !-----------------------------------------------------------------------
    
    f_loss_thres = np.minimum( np.maximum(ne.evaluate('(thres_z2_zoo - z_t)/(thres_z2_zoo - thres_z1_zoo)'), c0)  , c1)

    C_loss_thres = ne.evaluate('f_loss_thres * loss_thres_zoo')
    
    Zprime = np.maximum(ne.evaluate('zooC - C_loss_thres'), c0)
    
    # zoo ortality not modified by Tfunc like in bec
    
    zoo_loss = ne.evaluate('(parm_z_mort2_0 * Zprime**zoo_mort2_exp + parm_z_mort_0 * Zprime)*Tfunc')
    
    zoo_loss_poc = ne.evaluate('f_zoo_detr*zoo_loss')
    
    zoo_loss_doc = ne.evaluate('(c1 - parm_labile_ratio) * (c1 - f_zoo_detr) * zoo_loss')
    
    zoo_loss_dic = ne.evaluate('parm_labile_ratio * (c1 - f_zoo_detr) * zoo_loss')
##**************************** zoo grazing (not implemented)

    zoo_graze = c0
    zoo_graze_zoo = c0
    zoo_graze_zootot = c0
    x_graze_zoo = c0
    zoo_graze_poc = c0
    zoo_graze_doc = c0
    zoo_graze_dic = c0

######### ********************* ORGANIC MATTER ******************************************
######### ************************************************************************

##*************************** dissolved organic matter *************************
    
    DOC_prod = ne.evaluate('zoo_loss_doc + sp_loss_doc +  diat_loss_doc + diaz_loss_doc  +\
                  graze_sp_doc + graze_diat_doc + graze_diaz_doc + zoo_graze_doc')
    DON_prod = ne.evaluate('DOC_prod * Q*f_toDON')
    
    DOP_prod = ne.evaluate('(zoo_loss_doc + zoo_graze_doc) * Qp_zoo + remaining_spP_dop +remaining_diatP_dop +remaining_diazP_dop ')
               
    DOC_reminR = DOC_reminR_dark
    DON_reminR = DON_reminR_dark
    DOP_reminR = DOP_reminR_dark
    
    DOC_reminR = np.ma.where(PAR_avg > 1, DOC_reminR_light, DOC_reminR)
    DON_reminR = np.ma.where(PAR_avg > 1, DON_reminR_light, DON_reminR)
    DOP_reminR = np.ma.where(PAR_avg > 1, DOP_reminR_light, DOP_reminR)

    DOCr_reminR = DOCr_reminR0
    DONr_reminR = DONr_reminR0
    DOPr_reminR = DOPr_reminR0

    if kk == 0:
        work = np.log(PAR_in)*0.4373*(1000*idz)
        DOCr_reminR = np.ma.where(PAR_in > 1, ne.evaluate('work*DOMr_reminR_photo'), DOCr_reminR)
        DONr_reminR = np.ma.where(PAR_in > 1, ne.evaluate('work*DOMr_reminR_photo'), DONr_reminR)
        DOPr_reminR = np.ma.where(PAR_in > 1, ne.evaluate('work*DOMr_reminR_photo'), DOPr_reminR)


    
    
    DOC_remin = ne.evaluate('DOC * DOC_reminR')
    DON_remin = ne.evaluate('DON * DON_reminR')
    DOP_remin = ne.evaluate('DOP * DOP_reminR')
    
    DOCr_remin = ne.evaluate('DOCr * DOCr_reminR')
    DONr_remin = ne.evaluate('DONr * DONr_reminR')
    DOPr_remin = ne.evaluate('DOPr * DOPr_reminR')
    
    
# from BEC model
#    DOFe_prod = ne.evaluate('(zoo_loss_doc * Qfe_zoo)  + (Qfe_sp * (graze_sp_doc + sp_loss_doc)) \
#               + (Qfe_diat * (graze_diat_doc + diat_loss_doc)) \
#               + (Qfe_diaz * (graze_diaz_doc + diaz_loss_doc))')
 

##******************************* scavenge ***************************
    # !-----------------------------------------------------------------------
    # !  Compute iron scavenging :
    # !  1) compute in terms of loss per year per unit iron (%/year/fe)
    # !  2) scale by sinking POMx6 + Dust + bSi + CaCO3 flux
    # !  3) increase scavenging at higher iron (>0.6nM)
    # !  4) convert to net loss per second
    # !-----------------------------------------------------------------------
    
    
    p0 = c0
    p1 = c0
    p2 = c0
    Fefree = c0
    FeLig1 = c0
    FeLig2 = c0 ## allways 0 for us.

    p0 = np.ma.where(Fe > c0, -Fe, p0)
    p1 = np.ma.where(Fe > c0, ne.evaluate('c1 + KFeLig1*(Lig-Fe)'), p1)
    p2 = np.ma.where(Fe > c0, KFeLig1, p2)
    Fefree = np.ma.where(Fe > c0, (-p1 + np.sqrt(p1**2 - 4*p2*p0))/(2*p2), Fefree)
    FeLig1 = np.ma.where(Fe > c0, ne.evaluate('KFeLig1*Fefree*Lig/(c1 + KFeLig1 * Fefree)'),FeLig1)

    sinking_mass = ne.evaluate('(POC_sflux_in + POC_hflux_in) * 3*12.01 + \
           (P_CaCO3_sflux_in + P_CaCO3_hflux_in) * P_CaCO3_mass + \
           (P_SiO2_sflux_in+P_SiO2_hflux_in)*P_SiO2_mass + (dust_sflux_in + dust_hflux_in) * dust_fescav_scale')
   
   
    Fe_scavenge_rate = ne.evaluate('parm_Fe_scavenge_rate0 * sinking_mass')
    Lig_scavenge_rate = ne.evaluate('parm_Lig_scavenge_rate0 * sinking_mass')
    FeLig_scavenge_rate = ne.evaluate('parm_FeLig_scavenge_rate0 * sinking_mass')
    
    Lig_scavenge = ne.evaluate('yps * FeLig1 * Lig_scavenge_rate')
    Fe_scavenge = ne.evaluate('yps * Fefree * Fe_scavenge_rate + yps*FeLig1*FeLig_scavenge_rate')

##***************************** large detritus production **********************

    POC_prod = ne.evaluate('zoo_loss_poc + sp_agg + diat_agg + diaz_agg + sp_loss_poc + diat_loss_poc+ diaz_loss_poc   + graze_sp_poc + graze_diat_poc + graze_diaz_poc ')
    
    POP_prod = ne.evaluate('Qp_zoo*(zoo_loss_poc + zoo_graze_poc) + remaining_spP_pop + remaining_diatP_pop + remaining_diazP_pop')
    
    DOP_loss_P_bal = c0
    DOP_loss_P_bal = np.ma.where(POP_prod < c0, -POP_prod, DOP_loss_P_bal)
  
    POP_prod = np.maximum(0,POP_prod)
    # !-----------------------------------------------------------------------
    # !  large detrital CaCO3
    # !  33% of CaCO3 is remin when phyto are grazed
    # !-----------------------------------------------------------------------
    
    P_CaCO3_prod = ne.evaluate('((c1 - f_graze_CaCO3_remin) * graze_sp +  sp_loss + sp_agg) * QCaCO3')

    # !-----------------------------------------------------------------------
    # !  large detritus SiO2
    # !  grazed diatom SiO2, 60% is remineralized
    # !-----------------------------------------------------------------------

    P_SiO2_prod = ne.evaluate('((c1 - f_graze_si_remin) * graze_diat + diat_agg  \
              + f_diat_loss_poc * diat_loss) * Qsi') # only diat has Si
    
    dust_prod = c0

    P_iron_prod =  ne.evaluate('(sp_agg + graze_sp_poc + sp_loss_poc) * Qfe_sp +\
                                (diat_agg + graze_diat_poc + diat_loss_poc) * Qfe_diat +\
                                (diaz_agg + graze_diaz_poc + diaz_loss_poc) * Qfe_diaz+\
                                (zoo_loss_poc + zoo_graze_poc)* Qfe_zoo +\
                                Fe_scavenge')


######### ********************* PARTICULATES ******************************************
######### ************************************************************************

    # subroutine compute_particulate_terms(k, POC, P_CaCO3, P_SiO2, dust, P_iron, & 1,30
    #        QA_dust_def, temp, O2_loc, this_block)

    # ! !DESCRIPTION:
    # !  Compute outgoing fluxes and remineralization terms. Assumes that
    # !  production terms have been set. Incoming fluxes are assumed to be the
    # !  outgoing fluxes from the previous level.
    # !
    # !  It is assumed that there is no production of dust.
    # !
    # !  Instantaneous remineralization in the bottom cell is implemented by
    # !  setting the outgoing flux to zero.
    # !
    # !  For POC, the hard subclass is the POC flux qualitatively associated
    # !  with the ballast flux. The soft subclass is the excess POC flux.
    # !
    # !  Remineralization for the non-iron particulate pools is computing
    # !  by first computing the outgoing flux and then computing the
    # !  remineralization from conservation, i.e.
    # !     flux_in - flux_out + prod * dz - remin * dz == 0.
    # !
    # !  For iron, remineralization is first computed from POC remineralization
    # !  and then flux_out is computed from conservation. If the resulting
    # !  flux_out is negative or should be zero because of the sea floor, the
    # !  remineralization is adjusted.
    # !  Note: all the sinking iron is in the P_iron_sflux pool, hflux Fe not
    # !        explicitly tracked, it is assumed that total iron remin is
    # !        proportional to total POC remin.
    # !
    # !  Based upon Armstrong et al. 2000
    # !
    # !  July 2002, added temperature effect on remin length scale of
    # !  excess POC (all soft POM& Iron) and on SiO2.
    # !  new variable passed into ballast, Tfunc, main temperature function
    # !  computed in ecosystem routine.  scaling factor for dissolution
    # !  of excess POC, Fe, and Bsi now varies with location (f(temperature)).
    # !
    # !  Added diffusive iron flux from sediments at depths < 1100m,
    # !  based on Johnson et al., 1999, value of 5 umolFe/m2/day,
    # !      this value too high, using 2 umolFe/m2/day here
    # !
    # !  Allow hard fraction of ballast to remin with long length scale 40,000m
    # !     thus ~ 10_ of hard ballast remins over 4000m water column.
    # !
    # !  Sinking dust flux is decreased by assumed instant solubility/dissolution
    # !     at ocean surface from the parm_Fe_bioavail.
    # !
    # !  Modified to allow different Q10 factors for soft POM and bSI remin,
    # !  water temp is now passed in instead of Tfunc (1/2005, JKM)


    # POC,          & ! base units = nmol C
    # P_CaCO3,      & ! base units = nmol CaCO3
    # P_SiO2,       & ! base units = nmol SiO2
    # dust,         & ! base units = g
    # P_iron          ! base units = nmol Fe

    # diss,        & ! dissolution length for soft subclass
    # gamma,       & ! fraction of production -> hard subclass
    # mass,        & ! mass of 1e9 base units in g
    # rho            ! QA mass ratio of POC to this particle class
    # sflux_in,    & ! incoming flux of soft subclass (base units/cm^2/sec)
    # hflux_in,    & ! incoming flux of hard subclass (base units/cm^2/sec)
    # prod,        & ! production term (base units/cm^3/sec)
    # sflux_out,   & ! outgoing flux of soft subclass (base units/cm^2/sec)
    # hflux_out,   & ! outgoing flux of hard subclass (base units/cm^2/sec)
    # remin          ! remineralization term (base units/cm^3/sec)

    sed_denitrif = c0
    other_remin = c0
##********************** compute scalelength and decay factors ***************

    scalelength = c1
    
    if zw < parm_scalelen_z1:
        scalelength =  parm_scalelen_vals1
    elif zw >= parm_scalelen_z4:
        scalelength =  parm_scalelen_vals4
    else:
        if zw < parm_scalelen_z2:
            scalelength = ne.evaluate('parm_scalelen_vals1 + (parm_scalelen_vals2 - parm_scalelen_vals1)*(zw - parm_scalelen_z1)/(parm_scalelen_z2 - parm_scalelen_z1)')
        
        
        elif zw < parm_scalelen_z3:
            scalelength = ne.evaluate('parm_scalelen_vals2 + (parm_scalelen_vals3 - parm_scalelen_vals2)*(zw - parm_scalelen_z2)/(parm_scalelen_z3 - parm_scalelen_z2)')
    if p_remin_scalef != c1:
        scalelength = scalelength/p_remin_scalef

    
    DECAY_Hard     = np.exp(-dzs*p_remin_scalef  / 4.0e6)
    DECAY_HardDust = np.exp(-dzs*p_remin_scalef  / 1.2e8)

    # !-----------------------------------------------------------------------
    # !  increase POC diss where O2 concentrations are low ## not implementes
    # !-----------------------------------------------------------------------
    poc_diss = POC_diss
    sio2_diss = P_SiO2_diss
    caco3_diss = P_CaCO3_diss
    dust_diss = parm_dust_diss


#    o2_scalefactor = c0
#
#    o2_scalefactor = np.ma.where(O2 < o2_sf_o2_range_hi, c1 + (o2_sf_val_lo_o2 - c1) * np.minimum(c1, ne.evaluate('(o2_sf_o2_range_hi - O2)/(o2_sf_o2_range_hi - o2_sf_o2_range_lo)')), o2_scalefactor)
#
#    poc_diss   = np.ma.where(O2 < o2_sf_o2_range_hi,poc_diss   * o2_scalefactor, poc_diss)
#    sio2_diss  = np.ma.where(O2 < o2_sf_o2_range_hi,sio2_diss  * o2_scalefactor, sio2_diss)
#    caco3_diss = np.ma.where(O2 < o2_sf_o2_range_hi,caco3_diss * o2_scalefactor, caco3_diss)
#    dust_diss  = np.ma.where(O2 < o2_sf_o2_range_hi,dust_diss  * o2_scalefactor, dust_diss)
#



    poc_diss = ne.evaluate('poc_diss*scalelength')
    sio2_diss = ne.evaluate('sio2_diss*scalelength')
    caco3_diss = ne.evaluate('caco3_diss*scalelength')
    dust_diss = ne.evaluate('dust_diss*scalelength')

    DECAY_POC_E    = np.exp(-dzs  / poc_diss)
    DECAY_SiO2     = np.exp(-dzs  / sio2_diss)
    DECAY_CaCO3    = np.exp(-dzs  / caco3_diss)
    DECAY_dust     = np.exp(-dzs  / dust_diss)


    # !----------------------------------------------------------------------
    # !   Tref = 30.0 reference temperature (deg. C)
    # !
    # !   Using q10 formulation with Q10 value of 1.1 soft POM (TfuncP) and
    # !       a Q10 value of 2.5 soft bSi (TfuncS)
    # !-----------------------------------------------------------------------

    # !
    # !  NOTE: Turning off temperature effect on POM lengthscale, three instances
    # !    of TfuncP below have been removed, see comment lines.
    # !-------------------------------------------------------------------------
    #TfuncP = 1.1**(((temp + T0_Kelvin) - (Tref + T0_Kelvin)) / 10)
    TfuncS = ne.evaluate('2.5**(((temp + T0_Kelvin) - (Tref + T0_Kelvin)) / 10)')

    # !-----------------------------------------------------------------------
    # !  Set outgoing fluxes for non-iron pools.
    # !  The outoing fluxes for ballast materials are from the
    # !  solution of the coresponding continuous ODE across the model
    # !  level. The ODE has a constant source term and linear decay.
    # !  It is assumed that there is no sub-surface dust production.
    # !-----------------------------------------------------------------------

    P_CaCO3_sflux_out = ne.evaluate('P_CaCO3_sflux_in * DECAY_CaCO3 + \
       P_CaCO3_prod * ((c1 - P_CaCO3_gamma) * (c1 - DECAY_CaCO3) * caco3_diss)')

    P_CaCO3_hflux_out = ne.evaluate('P_CaCO3_hflux_in * DECAY_Hard + P_CaCO3_prod * (P_CaCO3_gamma * dzs )')

    P_SiO2_sflux_out = ne.evaluate('P_SiO2_sflux_in* DECAY_SiO2 + P_SiO2_prod * ((c1 - P_SiO2_gamma) * (c1 - DECAY_SiO2)*sio2_diss)')

    P_SiO2_hflux_out = ne.evaluate('P_SiO2_hflux_in * DECAY_Hard + P_SiO2_prod * (P_SiO2_gamma * dzs )')

    dust_sflux_out = ne.evaluate('dust_sflux_in * DECAY_dust')

    dust_hflux_out = ne.evaluate('dust_hflux_in * DECAY_HardDust')

    # !-----------------------------------------------------------------------
    # !  Compute how much POC_PROD is available for deficit reduction
    # !  and excess POC flux after subtracting off fraction of non-dust
    # !  ballast production from net POC_PROD.
    # !-----------------------------------------------------------------------

    POC_PROD_avail = ne.evaluate('POC_prod - P_CaCO3_rho * P_CaCO3_prod -  P_SiO2_rho * P_SiO2_prod')

    # !-----------------------------------------------------------------------
    # !  Check for POC production bounds violations
    # !-----------------------------------------------------------------------

    #             if (POC_PROD_avail < c0) then
    #                poc_error = .true.
    #             endif

    # !-----------------------------------------------------------------------
    # !  Compute 1st approximation to new QA_dust_def, the QA_dust
    # !  deficit leaving the cell. Ignore POC_PROD_avail at this stage.
    # !-----------------------------------------------------------------------
    new_QA_dust_def = c0
    # QA_dust def is defined as bla*(dust_sflux_in + dust_hflux_in), so no division by 0
    new_QA_dust_def = np.ma.where(QA_dust_def > 0, ne.evaluate('QA_dust_def *\
            (dust_sflux_out + dust_hflux_out) / (dust_sflux_in + dust_hflux_in)'), new_QA_dust_def)
    
    # !-----------------------------------------------------------------------
    # !  Use POC_PROD_avail to reduce new_QA_dust_def.
    # !-----------------------------------------------------------------------
    ### BEC has only new_QA_dust_def > c0 MARBL has np.logical_and(new_QA_dust_def > c0,POC_prod>0),
    new_QA_dust_def = np.ma.where(np.logical_and(new_QA_dust_def > c0, POC_prod >c0),ne.evaluate('new_QA_dust_def - POC_PROD_avail * dzs') , new_QA_dust_def)

    POC_PROD_avail = np.ma.where(new_QA_dust_def < c0, ne.evaluate('-new_QA_dust_def * idz') , POC_PROD_avail)

    POC_PROD_avail = np.ma.where(new_QA_dust_def >= c0, 0, POC_PROD_avail)

    new_QA_dust_def = np.maximum(new_QA_dust_def , c0)
    QA_dust_def = new_QA_dust_def

    # !-----------------------------------------------------------------------
    # !  Compute outgoing POC fluxes. QA POC flux is computing using
    # !  ballast fluxes and new_QA_dust_def. If no QA POC flux came in
    # !  and no production occured, then no QA POC flux goes out. This
    # !  shortcut is present to avoid roundoff cancellation errors from
    # !  the dust_rho * dust_flux_out - QA_dust_def computation.
    # !  Any POC_PROD_avail still remaining goes into excess POC flux.
    # !-----------------------------------------------------------------------


    POC_hflux_out  = ne.evaluate('P_CaCO3_rho * (P_CaCO3_sflux_out  + P_CaCO3_hflux_out ) + \
            P_SiO2_rho * (P_SiO2_sflux_out  + P_SiO2_hflux_out ) +  dust_rho * \
            (dust_sflux_out  + dust_hflux_out ) - new_QA_dust_def')
            
            
            
    POC_hflux_out  = np.maximum(POC_hflux_out , c0)

    POC_hflux_out  = np.ma.where(np.logical_and(POC_hflux_in  == c0,  POC_prod  == c0), c0, POC_hflux_out)
    
    POC_sflux_out  = ne.evaluate('POC_sflux_in  * DECAY_POC_E + POC_PROD_avail *((c1 - DECAY_POC_E) * poc_diss)')



    # !-----------------------------------------------------------------------
    # !  Compute remineralization terms. It is assumed that there is no
    # !  sub-surface dust production.
    # !-----------------------------------------------------------------------

    P_CaCO3_remin  = ne.evaluate('P_CaCO3_prod  + ((P_CaCO3_sflux_in  - P_CaCO3_sflux_out ) + \
       (P_CaCO3_hflux_in  - P_CaCO3_hflux_out )) * idz')

    P_SiO2_remin  = ne.evaluate('P_SiO2_prod  + ((P_SiO2_sflux_in  - P_SiO2_sflux_out ) + \
       (P_SiO2_hflux_in  - P_SiO2_hflux_out )) * idz')

    POC_remin  = ne.evaluate('POC_prod  +  ((POC_sflux_in  - POC_sflux_out ) + \
       (POC_hflux_in  - POC_hflux_out )) * idz')
       
    PON_remin = ne.evaluate('Q*POC_remin')

    dust_remin  =  ne.evaluate('((dust_sflux_in  - dust_sflux_out ) + \
       (dust_hflux_in  - dust_hflux_out )) * idz')

    # !-----------------------------------------------------------------------
    # !  Compute iron remineralization and flux out.
    # !-----------------------------------------------------------------------
    P_iron_remin = ne.evaluate('(POC_remin  *  (P_iron_sflux_in  + P_iron_hflux_in ) / \
    (POC_sflux_in  + POC_hflux_in ))')
 
    P_iron_remin  = np.ma.where(POC_sflux_in  + POC_hflux_in  == c0, ne.evaluate('(POC_remin  * parm_Red_Fe_C)'), P_iron_remin)
    
    P_iron_remin = ne.evaluate('P_iron_remin + P_iron_sflux_in*parm_Fe_desorption_rate0')
    
    
    P_iron_sflux_out  = ne.evaluate('P_iron_sflux_in  + dzs  *  ((c1 - P_iron_gamma) * P_iron_prod  - P_iron_remin )')

    P_iron_remin  = np.ma.where(P_iron_sflux_out  < c0, ne.evaluate('P_iron_sflux_in  * idz  + \
      (c1 - P_iron_gamma) * P_iron_prod'),  P_iron_remin)
    
    P_iron_sflux_out  =  np.maximum(P_iron_sflux_out  , c0)

    # !-----------------------------------------------------------------------
    # !  Compute iron release from dust remin/dissolution
    # !
    # !  dust remin gDust = 0.035 / 55.847 * 1.0e9 = 626712.0 nmolFe
    # !                      gFe     molFe     nmolFe
    # !  Also add in Fe source from sediments if applicable to this cell.
    # !-----------------------------------------------------------------------


    P_iron_remin  = ne.evaluate('P_iron_remin  + dust_remin  * dust_to_Fe  + (feflux_sed * idz )')

    #!maltrud what is up here--jkm does not have this, in MARBL this is not commented
    P_iron_hflux_out  = P_iron_hflux_in

    POP_remin = c0
    POP_remin =  np.ma.where(POC_sflux_in + POC_hflux_in > c0, ne.evaluate('POC_remin*((POP_sflux_in + POP_hflux_in)/(POC_sflux_in  + POC_hflux_in))'), POP_remin)
    POP_remin = np.ma.where(POC_prod>c0,ne.evaluate('POC_remin*(POP_prod/POC_prod)'),POP_remin)
   

    POP_sflux_out = ne.evaluate('POP_sflux_in + dzs*((1 - POP_gamma)*POP_prod - POP_remin)')
    
    
    
    POP_remin = np.ma.where(POP_sflux_out <c0, ne.evaluate('POP_sflux_in*idz + (1 - POP_gamma)*POP_prod'), POP_remin)
    
    
    POP_sflux_out = np.maximum(POP_sflux_out, c0)
    
    
    POP_hflux_out = POP_hflux_in
    
    POC_remin =  np.ma.where(mask_remin==frem, ne.evaluate('POC_remin  + (POC_sflux_out  + POC_hflux_out) * idz'),POC_remin)
    POC_sflux_out = np.ma.where(mask_remin==frem, c0,POC_sflux_out)
    
    POC_hflux_out =  np.ma.where(mask_remin==frem, c0,  POC_hflux_out)
    
    
    
#    # !-----------------------------------------------------------------------
#    # !  Remineralize everything in bottom cell.
    # !-----------------------------------------------------------------------


#! extract particulate fluxes at particulate_flux_ref_depth, if this layer contains that depth
    if (kk > 1):
        ztop = zw_all[kk-1]
    else:
       ztop = c0

     

     
#    if (np.logical_and(ztop <= particulate_flux_ref_depth_cm , particulate_flux_ref_depth_cm < zw)):
#        if (particulate_flux_ref_depth_cm == ztop):
#        # ! expressions are simplified if particulate_flux_ref_depth is exactly the top of the layer
#            POC_flux_at_ref_depth     = POC_sflux_in     + POC_hflux_in
#            POP_flux_at_ref_depth     = POP_sflux_in     + POP_hflux_in
#            P_CaCO3_flux_at_ref_depth = P_CaCO3_sflux_in + P_CaCO3_hflux_in
#            P_SiO2_flux_at_ref_depth  = P_SiO2_sflux_in  + P_SiO2_hflux_in
#            P_iron_flux_at_ref_depth  = P_iron_sflux_in  + P_iron_hflux_in
#        else:
#           wbot = (particulate_flux_ref_depth_cm - ztop)*idz
#
#           flux_top = POC_sflux_in + POC_hflux_in
#           flux_bot = POC_sflux_out + POC_hflux_out
#           POC_flux_at_ref_depth = flux_top + wbot * (flux_bot - flux_top)
#
#           flux_top = POP_sflux_in + POP_hflux_in
#           flux_bot = POP_sflux_out + POP_hflux_out
#           POP_flux_at_ref_depth = flux_top + wbot * (flux_bot - flux_top)
#
#           flux_top = P_CaCO3_sflux_in + P_CaCO3_hflux_in
#           flux_bot = P_CaCO3_sflux_out + P_CaCO3_hflux_out
#           P_CaCO3_flux_at_ref_depth = flux_top + wbot * (flux_bot - flux_top)
#
#           flux_top = P_SiO2_sflux_in + P_SiO2_hflux_in
#           flux_bot = P_SiO2_sflux_out + P_SiO2_hflux_out
#           P_SiO2_flux_at_ref_depth = flux_top + wbot * (flux_bot - flux_top)
#
#           flux_top = P_iron_sflux_in + P_iron_hflux_in
#           flux_bot = P_iron_sflux_out + P_iron_hflux_out
#           P_iron_flux_at_ref_depth = flux_top + wbot * (flux_bot - flux_top)

##*********************************  ligand terms ********************

    Lig_prod = ne.evaluate('POC_remin*remin_to_Lig + DOC_prod*remin_to_Lig')
    
    #!----------------------------------------------------------------------
    #!  Ligand losses due to photochemistry in first ocean layer
    #!  ligand photo-oxidation a function of PAR (2/5/2015)
    #!----------------------------------------------------------------------
    Lig_photochem = c0
    
    rate_per_sec = c0
    PAR_in_top0 = f_qsw_par*qsw[:,:]
    
    rate_per_sec = np.ma.where(PAR_in_top0 > 1, np.log(PAR_in_top0)*0.4373*(10e2/1000)*(4*yps), rate_per_sec)
    
    if kk == 0:
        Lig_photochem = ne.evaluate('Lig*rate_per_sec')

    Lig_deg = ne.evaluate('POC_remin*parm_Lig_degrade_rate0')
    
    Lig_loss = ne.evaluate('Lig_scavenge + 0.2*(photoFe_sp + photoFe_diat + photoFe_diaz) + Lig_photochem + Lig_deg')
    



##********************************** nitrificaton *********************

    
    NITRIF = 0

    NITRIF = np.ma.where(np.logical_and(PAR_out < parm_nitrif_par_lim,NH4 !=0), ne.evaluate('parm_kappa_nitrif * NH4'), NITRIF)


    NITRIF = np.ma.where(np.logical_and(PAR_in > parm_nitrif_par_lim,NH4 !=0), NITRIF*np.log(PAR_out/parm_nitrif_par_lim)/(-KPARdz), NITRIF)

   
    DENITRIF = c0
    
    sed_denitrif = c0
    


######### ********************* TIME DERIVATIVES ******************************************
######### ************************************************************************


    # !-----------------------------------------------------------------------
    # !  nitrate \ ammonium
    # !-----------------------------------------------------------------------

    dxdt[no3_ind] = ne.evaluate('restore_no3 + NITRIF - DENITRIF  - sed_denitrif - NO3_V_diat -  NO3_V_sp  - NO3_V_diaz ')

    dxdt[nh4_ind] = ne.evaluate('- NH4_V_diat - NH4_V_sp -NH4_V_diaz - NITRIF + \
    DON_remin + DONr_remin +\
    Q * (zoo_loss_dic + zoo_graze_dic + sp_loss_dic + diat_loss_dic + diaz_loss_dic + \
    graze_sp_dic + graze_diat_dic + graze_diaz_dic + DOC_prod*(c1 - f_toDON)) + PON_remin*(1 - PONremin_refract) + Nexcrete')

    # !-----------------------------------------------------------------------
    # !  dissolved iron
    # !-----------------------------------------------------------------------

    dxdt[fe_ind] =ne.evaluate('P_iron_remin - photoFe_sp - photoFe_diat - photoFe_diaz  - Fe_scavenge +\
        Qfe_zoo*(zoo_loss_dic + zoo_loss_doc + zoo_graze_dic + zoo_graze_doc) +\
        Qfe_sp*(sp_loss_dic + graze_sp_dic) + graze_sp_zoo * (Qfe_sp-Qfe_zoo) +\
        Qfe_sp*(sp_loss_doc + graze_sp_doc) + \
        Qfe_diat*(diat_loss_dic + graze_diat_dic) + graze_diat_zoo * (Qfe_diat-Qfe_zoo) +\
        Qfe_diat*(diat_loss_doc + graze_diat_doc) +\
        Qfe_diaz*(diaz_loss_dic + graze_diaz_dic) + graze_diaz_zoo * (Qfe_diaz-Qfe_zoo) +\
        Qfe_diaz*(diaz_loss_doc + graze_diaz_doc)')
    
    
    # !-----------------------------------------------------------------------
    # !  iron binding ligand
    # !-----------------------------------------------------------------------

    dxdt[lig_ind] = ne.evaluate('Lig_prod - Lig_loss')
    
    # !-----------------------------------------------------------------------
    # !  dissolved SiO3
    # !-----------------------------------------------------------------------

    dxdt[sio3_ind] = ne.evaluate('restore_sio3 + P_SiO2_remin  - photoSi_diat + \
            Qsi * (f_graze_si_remin * graze_diat + (1 - f_diat_loss_poc) * diat_loss)')

    # !-----------------------------------------------------------------------
    # !  phosphate
    # !-----------------------------------------------------------------------

    dxdt[po4_ind] = ne.evaluate('restore_po4 + DOPr_remin + \
            (c1 - POPremin_refract)*POP_remin +\
            remaining_spP_dip + remaining_diatP_dip + remaining_diazP_dip +\
            Qp_zoo * (zoo_loss_dic + zoo_graze_dic ) \
              - PO4_V_sp - PO4_V_diaz - PO4_V_diat')
    # !-----------------------------------------------------------------------
    # !  zoo Carbon
    # !-----------------------------------------------------------------------

    dxdt[zooC_ind] = ne.evaluate('graze_sp_zoo + graze_diat_zoo + graze_diaz_zoo -zoo_graze - zoo_loss')


    # !-----------------------------------------------------------------------
    # !  small phyto Carbon
    # !-----------------------------------------------------------------------

    dxdt[spC_ind] = ne.evaluate('photoC_sp - graze_sp - sp_loss - sp_agg')
    
    # !-----------------------------------------------------------------------
    # !  small phyto P
    # !-----------------------------------------------------------------------

    dxdt[spP_ind] = ne.evaluate('PO4_V_sp + DOP_V_sp - Qp_sp * (graze_sp + sp_loss + sp_agg)')
    
    # !-----------------------------------------------------------------------
    # !  small phyto Chlorophyll
    # !-----------------------------------------------------------------------

    dxdt[spChl_ind] = ne.evaluate('photoacc_sp - thetaC_sp * (graze_sp + sp_loss + sp_agg)')

    # !-----------------------------------------------------------------------
    # !  small phytoplankton CaCO3
    # !-----------------------------------------------------------------------

    dxdt[spCaCO3_ind] = ne.evaluate('CaCO3_PROD - QCaCO3*(graze_sp + sp_loss + sp_agg)')


    # !-----------------------------------------------------------------------
    # !  small phyto Fe
    # !-----------------------------------------------------------------------

    dxdt[spFe_ind] = ne.evaluate('photoFe_sp  - Qfe_sp * (graze_sp + sp_loss + sp_agg)')





    # !-----------------------------------------------------------------------
    # !  diatom Carbon
    # !-----------------------------------------------------------------------

    dxdt[diatC_ind] = ne.evaluate('photoC_diat - graze_diat -  diat_loss - diat_agg')
    
    # !-----------------------------------------------------------------------
    # !  Diatom P
    # !-----------------------------------------------------------------------

    dxdt[diatP_ind] = ne.evaluate('PO4_V_diat + DOP_V_diat - Qp_diat * (graze_diat + diat_loss + diat_agg)')

    # !-----------------------------------------------------------------------
    # !  diatom Chlorophyll
    # !-----------------------------------------------------------------------

    dxdt[diatChl_ind] = ne.evaluate('photoacc_diat - thetaC_diat * (graze_diat + diat_loss + diat_agg)')

    # !-----------------------------------------------------------------------
    # !  Diatom Fe
    # !-----------------------------------------------------------------------

    dxdt[diatFe_ind] =  ne.evaluate('photoFe_diat - Qfe_diat * (graze_diat + diat_loss + diat_agg)')
    
    # !-----------------------------------------------------------------------
    # !  Diatom Si
    # !-----------------------------------------------------------------------

    dxdt[diatSi_ind] =  ne.evaluate('photoSi_diat  - Qsi * (graze_diat + diat_loss + diat_agg)')



    # !-----------------------------------------------------------------------
    # !  Diazotroph C
    # !-----------------------------------------------------------------------

    dxdt[diazC_ind] = ne.evaluate('photoC_diaz - graze_diaz - diaz_loss - diaz_agg')
    
    # !-----------------------------------------------------------------------
    # !  Diazotroph P
    # !-----------------------------------------------------------------------

    dxdt[diazP_ind] = ne.evaluate('PO4_V_diaz + DOP_V_diaz - Qp_diaz * (graze_diaz + diaz_loss + diaz_agg)')
    
    # !-----------------------------------------------------------------------
    # !  diazotroph Chlorophyll
    # !-----------------------------------------------------------------------

    dxdt[diazChl_ind] =ne.evaluate('photoacc_diaz - thetaC_diaz * (graze_diaz + diaz_loss + diaz_agg)')

    # !-----------------------------------------------------------------------
    # !  Diazotroph Fe
    # !-----------------------------------------------------------------------

    dxdt[diazFe_ind] = ne.evaluate('photoFe_diaz - Qfe_diaz * (graze_diaz + diaz_loss + diaz_agg)')




    # !-----------------------------------------------------------------------
    # !  dissolved organic Matter
    # !-----------------------------------------------------------------------

    dxdt[doc_ind] = ne.evaluate('DOC_prod*(c1 - DOCprod_refract) - DOC_remin')
    dxdt[docr_ind] = ne.evaluate('DOC_prod*DOCprod_refract - DOCr_remin + POC_remin*POCremin_refract')
    
    dxdt[don_ind] = ne.evaluate('DON_prod*(c1 - DONprod_refract) - DON_remin')
    
    dxdt[donr_ind] = ne.evaluate('DON_prod * DONprod_refract - DONr_remin + PON_remin * PONremin_refract')
    
    dxdt[dop_ind] = ne.evaluate('DOP_prod * (c1 - DOPprod_refract) - DOP_remin - DOP_V_sp - DOP_V_diat - DOP_V_diaz   - DOP_loss_P_bal')
    
    dxdt[dopr_ind] = ne.evaluate('DOP_prod * DOPprod_refract - DOPr_remin + POP_remin * POPremin_refract')






    # # !-----------------------------------------------------------------------
    # # !   dissolved inorganic Carbon
    # # !-----------------------------------------------------------------------

    # dxdt[dic_ind] =  DOC_remin + POC_remin  + P_CaCO3_remin  + \
    #         f_graze_CaCO3_remin * graze_sp * QCaCO3 + zoo_loss_dic + sp_loss_dic + \
    #         graze_sp_dic + diat_loss_dic + graze_diat_dic - photoC_sp - photoC_diat \
    #         - CaCO3_PROD + graze_diaz_dic + diaz_loss_dic - photoC_diaz

    # !-----------------------------------------------------------------------
    # !  alkalinity
    # !-----------------------------------------------------------------------

    #    DTRACER_MODULE(:,:,alk_ind) = -DTRACER_MODULE(:,:,no3_ind) + \
    #       DTRACER_MODULE(:,:,nh4_ind) + \
    #       c2 * (P_CaCO3_remin  + f_graze_CaCO3_REMIN * graze_sp * QCaCO3 \
    #             - CaCO3_PROD)

    # !-----------------------------------------------------------------------
    # !  oxygen
    # !-----------------------------------------------------------------------

    # !  DTRACER_MODULE(:,:,o2_ind) = (photoC_sp + photoC_diat) / \
    # !     parm_Red_D_C_O2 + photoC_diaz/parm_Red_D_C_O2_diaz

    # !  WORK1 = (O2  - parm_o2_min) / parm_o2_min_delta
    # !  WORK1 = min(max(WORK1,c0),c1)
    # !  DTRACER_MODULE(:,:,o2_ind) = DTRACER_MODULE(:,:,o2_ind) + WORK1 \
    # !     * ((- POC_remin  - DOC_remin - zoo_loss_dic - sp_loss_dic  \
    # !        - graze_sp_dic - diat_loss_dic - graze_diat_dic   \
    # !        - graze_diaz_dic - diaz_loss_dic) / parm_Red_P_C_O2)

#    O2_PRODUCTION = zeros
#
#    O2_PRODUCTION = np.ma.where (photoC_diat > c0, O2_PRODUCTION + photoC_diat * \
#             ((NO3_V_diat/(NO3_V_diat+NH4_V_diat)) / parm_Red_D_C_O2 + \
#              (NH4_V_diat/(NO3_V_diat+NH4_V_diat)) / parm_Remin_D_C_O2), O2_PRODUCTION)
#
#    O2_PRODUCTION = np.ma.where(photoC_sp > c0, O2_PRODUCTION + photoC_sp * \
#             ((NO3_V_sp/(NO3_V_sp+NH4_V_sp)) / parm_Red_D_C_O2 + \
#              (NH4_V_sp/(NO3_V_sp+NH4_V_sp)) / parm_Remin_D_C_O2), O2_PRODUCTION)
#
#    O2_PRODUCTION = np.ma.where(photoC_diaz > c0, O2_PRODUCTION + photoC_diaz * \
#             ((photoNO3_diaz/(photoNO3_diaz+photoNH4_diaz+diaz_Nfix)) / parm_Red_D_C_O2 + \
#              (photoNH4_diaz/(photoNO3_diaz+photoNH4_diaz+diaz_Nfix)) / parm_Remin_D_C_O2 + \
#              (diaz_Nfix/(photoNO3_diaz+photoNH4_diaz+diaz_Nfix)) / parm_Red_D_C_O2_diaz), O2_PRODUCTION)
#
#    WORK1 = (O2  - parm_o2_min) / parm_o2_min_delta
#    WORK1 = np.minimum(np.maximum(WORK1,c0),c1)
#    O2_CONSUMPTION = WORK1 * ((POC_remin  + DOC_remin + zoo_loss_dic + \
#            sp_loss_dic + graze_sp_dic + diat_loss_dic + graze_diat_dic + \
#            graze_diaz_dic + diaz_loss_dic)/ parm_Remin_D_C_O2 + (2*NITRIF))
#
#    dxdt[o2_ind] = O2_PRODUCTION - O2_CONSUMPTION


    # Now assign particulates
   
    
    
    
    
    dxdt[ntracers+pic_sind] = P_CaCO3_sflux_out
    dxdt[ntracers+pic_hind] = P_CaCO3_hflux_out
    dxdt[ntracers+psio2_sind] = P_SiO2_sflux_out
    dxdt[ntracers+psio2_hind] = P_SiO2_hflux_out
    dxdt[ntracers+pdust_sind] = dust_sflux_out
    dxdt[ntracers+pdust_hind] = dust_hflux_out
    dxdt[ntracers+poc_sind] = POC_sflux_out
    dxdt[ntracers+poc_hind] = POC_hflux_out
    dxdt[ntracers+pop_sind] = POP_sflux_out
    dxdt[ntracers+pop_hind] = POP_hflux_out
    dxdt[ntracers+pfe_sind] = P_iron_sflux_out
    dxdt[ntracers+pfe_hind] = P_iron_hflux_out
    dxdt[ntracers+qadustdef_ind] = QA_dust_def
    
    
#    if kk == 0:
#        dxdt[ntracers+pic_sind] = c0
#        dxdt[ntracers+pic_hind] = c0
#        dxdt[ntracers+psio2_sind] = c0
#        dxdt[ntracers+psio2_hind] = c0
#        dxdt[ntracers+poc_sind] = c0
#        dxdt[ntracers+poc_hind] = c0
#        dxdt[ntracers+pop_sind] = c0
#        dxdt[ntracers+pop_hind] = c0
#        dxdt[ntracers+pfe_sind] = c0
#        dxdt[ntracers+pfe_hind] = c0
#    if kk < 10:
    
    #print('my code',  graze_diaz[0,0], graze_sp_zoo[0,0])
    #print('my code',  DOP[0,0],  PO4[0,0])

    return(np.nan_to_num(dxdt[:,:,:]), np.nan_to_num(PAR_out[:,:]))
