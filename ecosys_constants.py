from decimal import *
getcontext().prec = 14
## Constants used in ecosyscesm.py file
## time conversions


fcd = 86400 # factor to convert from secods to days
spd = 86400.0   # number of seconds in a day
dps = 1 / spd            # number of days in a second
yps = 1 / (365.0*spd) # number of years in a second
c0=0
frem = 8e9 ## value to indicate the last vertical element before land
c1 = 1
T0_Kelvin = 273.15
## Redfield ratios
parm_Red_D_C_P  = 117.0 # carbon:phosphorus
parm_Red_D_N_P  =  16.0 # nitrogen:phosphorus
parm_Red_D_O2_P = 170.0 # oxygen:phosphorus
parm_Remin_D_O2_P = 138.0 # oxygen:phosphorus
parm_Red_P_C_P  = parm_Red_D_C_P # carbon:phosphorus
parm_Red_D_C_N  = parm_Red_D_C_P/parm_Red_D_N_P # carbon:nitrogen
parm_Red_P_C_N  = parm_Red_D_C_N # carbon:nitrogen
parm_Red_D_C_O2 = parm_Red_D_C_P/parm_Red_D_O2_P # carbon:oxygen
parm_Remin_D_C_O2 = parm_Red_D_C_P/parm_Remin_D_O2_P # carbon:oxygen
parm_Red_P_C_O2 = parm_Red_D_C_O2 # carbon:oxygen
parm_Red_Fe_C   = 3.0e-6 # iron:carbon
parm_Red_D_C_O2_diaz = parm_Red_D_C_P/150.0 # carbon:oxygen for diazotrophs

# ecosystem parameters accessible via input file
parm_Fe_bioavail    = 0.02
parm_o2_min         = 4.0
parm_o2_min_delta   = 2.0
parm_no3_min        = 110.0
parm_kappa_nitrif   = 0.06 * dps       # (= 1/( days) to 1/s)
parm_nitrif_par_lim = 1.0

Q             = 0.137 # N/C ratio (mmol/mmol) of phyto & zoo
Qp            = 8.5470085470085479e-03 # P/C ratio (mmol/mmol) sphyto,diat,zoo


parm_z_mort_0       = 0.1* dps
parm_z_mort2_0      = 0.4* dps
parm_z_grz          = 1.2 # same for all phyto groups



PquotaSlope = 7
PquotaIntercept = 5.571
PquotaMinNP = 0.00854701
Qp_zoo = 1/117
Qfe_zoo       = 3.0e-6 # zooplankton fe/C ratio
thres_z1_zoo          = 110.0e2 # threshold = C_loss_thres for z shallower than this (cm)
thres_z2_zoo          = 150.0e2 # threshold = 0 for z deeper than this (cm)
zoo_mort2_exp = 1.5
loss_thres_zoo    = 0.075 # zoo conc. where losses go to zero
loss_thres2_zoo    = 0.001 # zoo conc. where losses go to zero

gQ_Fe_kFe_thres = 10
gQ_Si_kSi_thres  = 6


parm_sp_kNO3        = 0.25
parm_sp_kNH4        = 0.01
parm_sp_kFe         = 0.03e-3
parm_sp_kPO4        = 0.01
parm_sp_kDOP        = 0.3
parm_sp_CO2         = 0
parm_sp_SiO3        = 0
parm_alphaChlsp     = 0.39 * dps
parm_sp_z_umax_0    = 3.3 * dps
#f_sp_z_grz          = 1.2 ## same as parm_z_grz for all phyto groups
f_sp_graze_zoo      = 0.3
f_sp_graze_poc      = 0.0
f_sp_graze_doc      = 0.06
f_sp_zoo_detr       = 0.12
PCref_sp = 5* dps # max phyto C-spec. grth rate at tref (1/sec)
sp_mort    = 0.1  * dps # sphyto mort rate (1/sec)
sp_mort2   = 0.01 * dps # sphyto quad. mort rate (1/sec/((mmol C/m3))
sp_agg_rate_max   = 0.5 # max agg. rate for small phyto (1/d)
sp_agg_rate_min   = 0.01 # min agg. rate for small phyto (1/d)
Qp_sp_fixed = 8.5470085470085479e-03 #1/117
gQp_sp_fixed = Qp_sp_fixed
gQfe_sp_0     = 3.0000000000000001e-05 # initial sphyto fe/C ratio
sp_gQfe_min   = 2.5e-6 # min sphyto fe/C ratio
sp_alphaPI_per_day = 0.39
thetaN_max_sp   = 2.5 # sp max thetaN (Chl/N) (mg Chl/mmol N)
QCaCO3_max    = 0.4 # max QCaCO3
loss_thres_sp     = 0.01 # small phyto conc. where losses go to zero
loss_thres2_sp     = 0 # small phyto conc. where losses go to zero
temp_thres_sp     = -10 # temp where losses go to zero
spc_poc_fac      = 0.13 # small phyto grazing factor (1/mmolC)
f_graze_sp_poc_lim = 0.36
f_prod_sp_CaCO3  = 0.07 # fraction of sp prod. as CaCO3 prod.
f_photosp_CaCO3  = 0.4 # proportionality between small phyto production and CaCO3 pro
f_sp_loss_poc = 0

#f_graze_sp_dic   = 1 - z_ingest_sp - f_graze_sp_doc # fraction to DIC
z_ingest_sp         = 0.3 # graze_zoo in marbl


parm_diat_kFe       = 7e-5
parm_diat_kNO3      = 0.5
parm_diat_kNH4      = 0.05
parm_diat_kSiO3     = 0.7
parm_diat_kDOP      = 0.5
parm_diat_kPO4      = 0.05
parm_diat_kDOP      = 0.5
parm_diat_CO2       = 0
parm_alphaChldiat     = 0.28 * dps
parm_diat_z_umax_0    = 3.15 * dps
PCref_diat = 5* dps # max phyto C-spec. grth rate at tref (1/sec)
diat_mort  = 0.1   * dps # diatom mort rate (1/sec)
diat_mort2 = 0.01 * dps # diatom quad mort rate (1/sec/((mmol C/m3))
diat_agg_rate_max = 0.5 # max agg. rate for diatoms (1/d)
diat_agg_rate_min = 0.02 # min agg. rate for diatoms (1/d)
Qp_diat_fixed = 8.5470085470085479e-03# 1/117
gQp_diat_fixed = Qp_diat_fixed
gQsi_0        = 0.137 # initial diatom Si/C ratio
gQfe_diat_0   = 3.0000000000000001e-05 # initial diatom fe/C ratio
diat_gQfe_min = 2.5e-6 # min diatom fe/C ratio
gQsi_max      = 0.822 # max diatom Si/C ratio
gQsi_min      = 0.0457 # min diatom Si/C ratio
thetaN_max_diat = 4.0 # diat max thetaN (Chl/N) (mg Chl/mmol N)
diat_alphaPI_per_day = 0.29
f_graze_diat_poc = 0.39 # fraction diatom grazing to POC
f_graze_diat_doc = 0.06
f_diat_loss_poc  = 0 # fraction diatom loss to POC
f_diat_loss_dc   = 1-f_diat_loss_poc # fraction diatom loss to DOC
f_diat_zoo_detr = 0.24 # fraction of zoo losses to detrital pool when eating diatoms
#f_graze_diat_dic = 1 - z_ingest_diat - f_graze_diat_poc  - f_graze_diat_doc # fraction diatom grazing to DIC
loss_thres_diat   = 0.02 # diat conc. where losses go to zero
loss_thres2_diat   = 0 # diat conc. where losses go to zero
temp_thres_diat   = -10 #
z_ingest_diat       = 0.25 # graze_zoo in marbel (diat graze by zoo)


parm_diaz_kNO3      = 2
parm_diaz_kNH4      = 0.2
parm_diaz_kPO4      = 0.015 # diaz half-sat. const. for P (diatom value)
parm_diaz_kDOP      = 0.075 # diaz half-sat. const. for P (diatom value)
parm_diaz_kFe       = 0.045e-3
parm_diaz_kCO2      = 0
parm_diaz_kSiO3     = 0
parm_alphaChldiaz     = 0.39 * dps
parm_diaz_z_umax_0    = 3.33 * dps
PCref_diaz = 2.5* dps # max phyto C-spec. grth rate at tref (1/sec)
diaz_mort  = 0.1* dps # diaz mort rate (1/sec)
diaz_mort2  = 0.01* dps # diaz mort rate (1/sec)
diaz_agg_rate_max = 0.5 # max agg. rate for diaz (1/d)
diaz_agg_rate_min = 0.01 # min agg. rate for diaz(1/d)
Qp_diaz_fixed = 2.7350427350427355e-03#0.32/117
gQp_diaz_fixed = Qp_diaz_fixed
gQfe_diaz_0   = 6.0000000000000002e-05#from init.mod, 70.0e-6 from marbl pft.mod # initial diaz. fe/C ratio
diaz_gQfe_min = 2.5e-06 #5.4e-6 from pft.mod # min diaz fe/C ratio
thetaN_max_diaz = 2.5 # diaz max thetaN (Chl/N) (mg Chl/mmol N) carbon:nitrogen ratio
diaz_alphaPI_per_day = 0.39
#for denitrification net removal of 120 mols NO3 for 117 mols C (136 = 120 + 16)
f_diaz_loss_poc = 0 # fraction diaz loss to sinking pool
f_graze_diaz_zoo = 0.3 # fraction diaz. grazing to zoo
f_graze_diaz_poc = 0.1 # fraction diaz grazing to POC
f_graze_diaz_doc = 0.06
f_diaz_zoo_detr = 0.12 # fraction of zoo losses to detrital pool when eating diaz
#f_graze_diaz_dic = 1-f_graze_diaz_zoo-f_graze_diaz_poc - f_graze_diaz_doc # fraction diaz grazing to DIC
loss_thres_diaz   = 0.02 # diaz conc. where losses go to zero
loss_thres_diaz2  = 0.001 # diaz conc. thres at low temp
temp_thres_diaz   = 15.0 # Temp. where diaz conc thres drops
z_ingest_diaz       = 0.3 # graze_zoo in marbl




parm_sd_remin_0     = 0.006667 * dps     # (= 1/(100 days))



parm_alphaChl       = 0.3 * dps

parm_labile_ratio   = 0.94


auto_mort2_exp = 1.75 # autothrophs

KFeLig1 = 10e13*1e-6
parm_Fe_scavenge_rate0 = 18 # base scavenging rate (1/y)
parm_Lig_scavenge_rate0 = 0.015
parm_FeLig_scavenge_rate0 = 1.4
parm_Lig_degrade_rate0 = 0.000094
remin_to_Lig = 0.0001

fe_scavenge_rate0 = 18 # base scavenging rate (1/y)
fe_scavenge_thres1 = 0.6e-3 # upper thres. for Fe scavenging
dust_fescav_scale  = 1.0e9 # dust scavenging scale factor
fe_max_scale2      = 1000.0 # unitless scaling coeff.
f_fescav_P_iron    = 0.9 # fraction of Fe scavengingto particulate Fe


#  ! dust_to_Fe: conversion of dust to iron (nmol Fe/g Dust)
#  ! dust remin gDust = 0.035 gFe       mol Fe     1e9 nmolFe
#  !                    --------- *  ----------- * ----------
#  !                      gDust      molw_Fe gFe      molFe
#

dust_to_Fe =0.035/55.847*1.0e9 # dust to iron conversion


caco3_poc_min    = 0.4 # minimum proportionality between QCaCO3 and grazing losses to POC (mmol C/mmol CaCO3)




graze_doc   = 0.06 # fraction  phyto. grazing to DOC, same for all species





f_graze_CaCO3_remin = 0.33 # fraction of spCaCO3 grazing which is remin
f_graze_si_remin    = 0.5 # fraction of diatom Si grazing which is remin
r_Nfix_photo=1.25 # N fix relative to C fix (non-dim)

#    SET FIXED RATIOS for N/C, P/C, SiO3/C, Fe/C
#    assumes C/N/P of 117/16/1 based on Anderson and Sarmiento, 1994
#    for diazotrophs a N/P of 45 is assumed based on Letelier & Karl, 1998
f_toDON               = 0.70 #  & ! fraction DON relative to DOC
f_toDOP               = 0.15 #    ! fraction of remaining_P to DOP






denitrif_C_N  = parm_Red_D_C_P/136.0

## loss term threshold parameters, chl:c ratios

thres_z1_phyto          = 80.0e2 # threshold = C_loss_thres for z shallower than this (cm), for autothrophs
thres_z2_phyto          = 120.0e2 # threshold = 0 for z deeper than this (cm)



CaCO3_temp_thres1 = 4.0 # upper temp threshold for CaCO3 prod
CaCO3_temp_thres2 = -2.0 # lower temp threshold
CaCO3_sp_thres    = 2.5 # bloom condition thres (mmolC/m3)
diaz_kNO3         = 1.0 # diazotroph Ks for nitrate (mmolN/m3)
diaz_kNH4         = 0.1 # diazotroph Ks for ammonimum (mmolN/m3)

## attenuation coefficients for PAR and related parameters

k_chl = 0.03e-2 # Chl atten. coeff. (1/cm/(mg Chl/m^3))
k_h2o = 0.04e-2 # water atten. coeff (1/cm)
f_qsw_par = 0.45 # PAR fraction

## Temperature parameters

Tref = 30.0 # reference temperature (C)
Q_10 = 2.0 # factor for temperature dependence (non-dim)

phlo_init = 7.0 # low bound for ph for no prev soln
phhi_init = 9.0 # high bound for ph for no prev soln
del_ph = 0.20 # delta-ph for prev soln

surf_avg_dic_const = 1944.0
surf_avg_alk_const = 2225.0
atm_co2_const = 280.0

epsC      = 1.00e-8 # small C concentration (mmol C/m^3)
epsTinv   = 3.17e-8 # small inverse time scale (1/year) (1/sec)
epsnondim = 1.00e-6 # small non-dimensional number (non-dim)





po4_ind          =  0 #  dissolved inorganic phosphate
no3_ind          =  1 #  dissolved inorganic nitrate
sio3_ind         =  2 #  dissolved inorganic silicate
nh4_ind          =  3 #  dissolved ammonia
fe_ind           =  4 #  dissolved inorganic iron
lig_ind          = 5
#o2_ind           =  6 #  dissolved oxygen
#dic_ind          =  6 #  dissolved inorganic carbon
#alk_ind          =  7 #  alkalinity
doc_ind          =  6 #  dissolved organic carbon
don_ind          = 7 #  dissolved organic nitrogen
dop_ind          = 8 #  dissolved organic phosphorus
dopr_ind         = 9
donr_ind         = 10
docr_ind         = 11

zooC_ind         = 12 #  zooplankton carbon

spC_ind          =  13 #  small phytoplankton carbon
spP_ind          = 14
spChl_ind        = 15 #  small phytoplankton chlorophyll
spFe_ind         = 16 #  small phytoplankton iron
spCaCO3_ind      = 17 #  small phytoplankton caco3

diatC_ind        = 18 #  diatom carbon
diatChl_ind      = 19 #  diatom chlorophyll
diatSi_ind       = 20 #  diatom silicon
diatFe_ind       = 21 #  diatom iron
diatP_ind        = 22 # added for marbl

diazC_ind        = 23 #  diazotroph carbon
diazChl_ind      = 24 #  diazotroph Chlorophyll
diazFe_ind       = 25 #  diazotroph iron
diazP_ind        = 26 # added for marbl

pic_sind = 0#27
pic_hind = 1#28
psio2_sind = 2#29
psio2_hind = 3#30
pdust_sind = 4#31
pdust_hind = 5#32
poc_sind = 6#33
poc_hind = 7#34
pfe_sind = 8#35
pfe_hind = 9#36
pop_sind = 10#37 # added for marbl
pop_hind = 11#38 # added for marbl
qadustdef_ind = 12#39


#sflux_in,    & ! incoming flux of soft subclass (base units/cm^2/sec)
#hflux_in,    & ! incoming flux of hard subclass (base units/cm^2/sec)


#nutr_rest_time_inv  #! inverse restoring time scale for nutrients (1/secs)







    # !-----------------------------------------------------------------------
    # !  parameters, from Armstrong et al. 2000
    # !
    # !  July 2002, length scale for excess POC and bSI modified by temperature
    # !  Value given here is at Tref of 30 deg. C, JKM
    # !-----------------------------------------------------------------------

parm_hPOC_dust_ratio = 0.01
parm_POC_diss = 10000


POC_diss      =  parm_POC_diss # diss. length (cm), modified by temp
POC_gamma     = c0      #  not used
POC_mass      = 12.01#  molecular weight of POC
POC_rho       = c0     # not used

POP_diss = parm_POC_diss
POP_gamma = c0
POP_mass = 30.974
POP_rho = c0

parm_CaCO3_diss = 50000.0
parm_hPOC_CaCO3_ratio = 0.01
P_CaCO3_diss  =  parm_CaCO3_diss # diss. length (cm)
P_CaCO3_gamma =  0.02 # prod frac -> hard subclass
P_CaCO3_mass  = 100.09 # molecular weight of CaCO
P_CaCO3_rho   = parm_hPOC_CaCO3_ratio* P_CaCO3_mass / POC_mass # QA mass ratio for CaCO3, This ratio is used in ecos_set_interior

P_SiO2_diss   = 65000 # diss. length (cm), modified by temp
P_SiO2_gamma  = 0.00 # prod frac -> hard subclass
P_SiO2_mass   =  60.08 # molecular weight of SiO2
parm_hPOC_SiO2_ratio = 0.01
P_SiO2_rho    =  parm_hPOC_SiO2_ratio * P_SiO2_mass / POC_mass # QA mass ratio for SiO2

parm_dust_diss     =  40000.0 # diss. length (cm)
dust_gamma    =  0.98 # prod frac -> hard subclass
dust_mass     =  1.0e9 # base units are already grams
dust_rho      =  parm_hPOC_dust_ratio * dust_mass / POC_mass # QA mass ratio for dust

P_iron_diss   = 60000.0 # diss. length (cm) - not used
P_iron_gamma  = c0        # prod frac -> hard subclass - not used
P_iron_mass   = c0        # not used
P_iron_rho    =  c0        # not used

parm_Fe_desorption_rate0 = 1e-6

parm_f_prod_sp_CaCO3 =   0.07
parm_sed_denitrif_coeff =   1
bury_coeff_rmean_timescale_years =   10
#caco3_bury_thres_opt = 'omega_calc'
caco3_bury_thres_depth =   300000
caco3_bury_thres_omega_calc = 1
PON_bury_coeff =   0.5

parm_scalelen_z1 = 10000
parm_scalelen_z2 = 25000
parm_scalelen_z3 = 50000
parm_scalelen_z4 = 100000

parm_scalelen_vals1 = 1
parm_scalelen_vals2 = 3.6
parm_scalelen_vals3 = 4.7
parm_scalelen_vals4 = 4.8

particulate_flux_ref_depth_cm = 100*100
#     particulate_flux_ref_depth_cm = cmperm * particulate_flux_ref_depth, cmperm = 100


DOCprod_refract  = 0.01       # & ! fraction of DOCprod to refractory pool
DONprod_refract  = 0.0115      # & ! fraction of DONprod to refractory pool
DOPprod_refract  = 0.003      # & ! fraction of DOPprod to refractory pool
POCremin_refract = DOCprod_refract * 0.06      # & ! fraction of POCremin to refractory pool
PONremin_refract = DONprod_refract * 0.03      # & ! fraction of POCremin to refractory pool
POPremin_refract = DOPprod_refract * 0.06      #  ! fraction of POCremin to refractory pool



DOC_reminR_light  = (c1/(365.0*15.0)) * dps #, & ! remin rate for semi-labile DOC, 1/15yr
DON_reminR_light  = (c1/(365.0*15.0)) * dps #, & ! remin rate for semi-labile DON, 1/15yr
DOP_reminR_light  = (c1/(365.0*60.0)) * dps #, & ! remin rate for semi-labile DOP, 1/60yr
DOC_reminR_dark   = (c1/(365.0*6.0)) * dps #,  & ! remin rate in the dark, 1/6yr
DON_reminR_dark   = (c1/(365.0*5.5)) * dps #  & ! remin rate in the dark, 1/5.5yr
DOP_reminR_dark   = (c1/(365.0*4.5)) * dps  #    ! remin rate in the dark, 1/4.5yr


DOCr_reminR0      = (c1/(365.0*16000.0)) * dps # & ! remin rate for refractory DOC, 1/16000yr
DONr_reminR0      = (c1/(365.0*9500.0)) * dps #  & ! remin rate for refractory DON, 1/9500yr
DOPr_reminR0      = (c1/(365.0*5500.0)) * dps #  & ! remin rate for refractory DOP, 1/5500yr
DOMr_reminR_photo = (c1/(365.0*18.0)) * dps   #    ! additional remin from photochemistry, 1/18yrs over top 10m


o2_sf_o2_range_hi = 45
o2_sf_o2_range_lo = 5
o2_sf_val_lo_o2 = 2.6
