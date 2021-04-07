
%matplotlib inline

import numpy as np
import matplotlib.pyplot as plt
import warnings
import netCDF4
import gc
import sys
import importlib # to reload our modules if changes are made
# give the path where the model modules are stored
sys.path.insert(1, '/Users/angeles/Desktop/python/python2/Barakawa/modules')
import continuity
import adv_diff
import main
import ecosystem
import utilities
import density
import saver
import ecosys_constants
from ecosys_constants import *

# directory where data for running model is stored
datadir = '/Volumes/WD_BLACK/data/final_run_chapter3/'


# IMPORTANT: If the data files are already processed and in correct input form, you can skip the processing parts and go to line 280 od the code, or where it says "INPUT DATA LIST"




# slice of depths you want to use
zli = np.s_[:]









# ************************ Load U and V ************************
# load U
fileyr = '0120'
fileext = '.day.'+fileyr+'.15N_45N.115E_145E_interp_depth.nc'
fil = netCDF4.Dataset(datadir+'velocities/UVEL/processed/UVEL'+fileext)
uvelg = fil.variables['UVEL'][:,zli,1:,1:] #cm/s
time = fil.variables['time'][:]
ulat = fil.variables['lat'][1:]
ulon = fil.variables['lon'][1:]
z_t = fil.variables['z_t'][zli]
fil.close()


# load V
fileyr = '0120'
fileext = '.day.'+fileyr+'.15N_45N.115E_145E_interp_depth.nc'
fil = netCDF4.Dataset(datadir+'velocities/VVEL/processed/VVEL'+fileext)
vvelg = fil.variables['VVEL'][:,zli,1:,1:] #cm/s
fil.close()

# ************************ Create masks ************************
## create mask, u and v have the same mask
nlac = len(ulat[:])
nloc = len(ulon[:])
nzc = len(z_t)

# obtain the mask of the velocities
mask1 = np.ones((nzc,nlac,nloc))
m_mask = np.ma.getmask(uvelg[3,:nzc,:,:])
masku = np.ma.masked_array(mask1, mask=m_mask)
masku0 = np.ma.filled(masku ,0)

# create mask for tracers based on u
maskT = masku0[:,:-1,:-1]*masku0[:,1:,1:] # this is how the mask for tracers should be
maskU = masku0[:,1:,1:] # now U has to start one lat and lon to the east-north of T

# now create no-flux mask based on tracer, at the zonal ocean-land boundaries
mask_zonal = np.zeros(np.shape(maskT))
mask_zonal[:,:,:-1] = maskT[:,:,:-1]*maskT[:,:,1:]
mask_zonal[:,:,-1] = 0
# and at the meridional ocean-land boundaries
mask_merid = np.zeros(np.shape(maskT))
mask_merid[:,:-1,:] = maskT[:,:-1,:]*maskT[:,1:,:]
mask_merid[:,-1,:] = 0

# now we need to remove one lat and lon for U, because of how the mask was created
uvelg = (uvelg[:,:,1:,1:]).filled(fill_value=0)*maskU
vvelg = (vvelg[:,:,1:,1:]).filled(fill_value=0)*maskU

# ************************ Load forcings ************************

# Load Temp
fileyr = '0120'
fileext = '.day.'+fileyr+'.15N_45N.115E_145E_interp_depth_setmisstonn.nc'
fil = netCDF4.Dataset(datadir+'forcing/TEMP/processed/TEMP'+fileext)
tempg = fil.variables['TEMP'][:,zli,1:,1:]*maskT
tlat = fil.variables['lat'][1:]
tlon = fil.variables['lon'][1:]
z_t2 = fil.variables['z_t'][zli]
latre = np.abs(tlat[:])
fil.close()

# load short wave radiation
fileyr = '0120'
fileext = '.day.'+fileyr+'.15N_45N.115E_145E_processed.nc'
fil = netCDF4.Dataset(datadir+'forcing/FSDS/processed/FSDS'+fileext)
qswg = fil.variables['FSDS'][:,1:,1:]*maskT[0,:,:]
fil.close()



# ************************ Load Interior forcings ************************

diric = datadir + 'surface_flux/'
file_ext = '_feventflux_gx1v6_5gmol_cesm1_97_2017_setmisstonn.nc'
filr = netCDF4.Dataset(diric+'FESEDFLUXIN'+file_ext)
## we remove one lat and lon, because we lost one when we created the mask
feflux_ventg = filr.variables['FESEDFLUXIN'][:,1:,1:]*maskT*1.1574e-6 # from mumol/m2/day to  nmol/cm2/s # z, lat, lon (no units)
filr.close()

file_ext = '_fesedfluxTot_gx1v6_cesm2_2018_c180618_setmisstonn.nc'
filr = netCDF4.Dataset(diric+'FESEDFLUXIN'+file_ext)
feflux_sedg = filr.variables['FESEDFLUXIN'][:,1:,1:]*maskT*1.1574e-6 #from mumol/m2/day to  nmol/cm2/s, z, lat, lon
#in CESM site the file is multiplied by this factor,
#the field does not contain the units in it.
#### from https://github.com/ESCOMP/POP2-CESM/blob/master/bld/namelist_files/namelist_defaults_pop.xml:
# <fesedflux_input%scale_factor>1.1574e-6</fesedflux_input%scale_factor>
filr.close()


# ************************ Load surface climatological forcings ************************

file_ext = '_dst79gnx_gx1v6_090416_processed.nc'
filr = netCDF4.Dataset(diric+'DSTSF'+file_ext)
dustg = filr.variables['DSTSF'][:,1:,1:]*maskT[0,:,:]*0.1 #kg/m2/s to g/cm2/s, time, lat, lon
filr.close()

file_ext = '_solFe_scenario4_current_gx1v6_8gmol_cesm1_93_20161114_processed.nc'
filr = netCDF4.Dataset(diric+'DSTSF'+file_ext)
feflux_solg = filr.variables['DSTSF'][:,1:,1:]*maskT[0,:,:]*6.2668e4 #kg/m2/s to nmol/cm2/s , time, lat, lon
filr.close()


#### CLIMATOLOGIES for years 1950-1960
file_ext = '_ndep_ocn_1850-2000_w_nhx_emis_gx1v6_c180926_processed_1950_1960.nc'
filr = netCDF4.Dataset(diric+'NHx_deposition'+file_ext)
NHyg = filr.variables['NHx_deposition'][:,1:,1:]*maskT[0,:,:]*7.1429e6 #  kg/m2/s to nmol/cm2/s # time, lat, lon

filr = netCDF4.Dataset(diric+'NOy_deposition'+file_ext)
NOxg = filr.variables['NOy_deposition'][:,1:,1:]*maskT[0,:,:]*7.1429e6  #  kg/m2/s to nmol/cm2/s, time, lat, lon
timen = filr.variables['time'][:]
filr.close()


# ************************ Interpolate climatological forcings to daily *********************

importlib.reload(utilities)
nt,nla,nlo = np.shape(dustg)
dustg  = utilities.intmonth(dustg,nt=nt,nla=nla,nlo=nlo,ndim=3)
feflux_solg = utilities.intmonth(feflux_solg,nt=nt,nla=nla,nlo=nlo,ndim=3)
NOxg = utilities.intmonth(NOxg,nt=nt,nla=nla,nlo=nlo,ndim=3)
NHyg = utilities.intmonth(NHyg,nt=nt,nla=nla,nlo=nlo,ndim=3)

#************************** Initial Conditions *************************

diric = datadir + 'ic_ecosystem/'
file_ext = '_20N_36N_120E_1240E_interp_depth_setmisstonn.nc'
# diric = '/Volumes/WD_BLACK/data/CESM2/model_inputs/cesm2_original/'
# grid = 'g01x01'
# diric = diric+grid+'/ic_ecosys/ecosys_2020/'
#file_ext = '_init_Barakawa_setmisstonn_25N_50N.120E_160E_latlon.nc'

var_name = ['PO4', 'NO3', 'SiO3', 'NH4', 'Fe', 'Lig', 'DOC', 'DON', \
            'DOP', 'DOPr', 'DONr', 'DOCr', 'zooC', 'spC', 'spP', 'spChl',\
            'spFe', 'spCaCO3', 'diatC', 'diatChl', 'diatSi', 'diatFe',\
            'diatP', 'diazC', 'diazChl', 'diazFe', 'diazP']

var_index = [po4_ind, no3_ind, sio3_ind, nh4_ind, fe_ind, lig_ind,\
             doc_ind, don_ind, dop_ind, dopr_ind, donr_ind, docr_ind,\
             zooC_ind, spC_ind, spP_ind, spChl_ind, spFe_ind, spCaCO3_ind,\
            diatC_ind, diatChl_ind, diatSi_ind, diatFe_ind, diatP_ind,\
            diazC_ind, diazChl_ind, diazFe_ind, diazP_ind]

ntracers = 27

## get size of array

filr = netCDF4.Dataset(diric+'Fe'+file_ext)
Nz, Ny, Nx = np.shape(filr.variables['Fe'][:,1:,1:]) # remember to remove one lat lon
filr.close()

x0csm = np.ma.zeros((ntracers,Nz,Ny,Nx))

for var,ind in zip(var_name,var_index):
    filr = netCDF4.Dataset(diric+var+file_ext)
    x0csm[ind,:,:,:] = filr.variables[var][:,1:,1:]*maskT[:,:,:]
    filr.close()

# ### Important: there are  0.0 in  places that are not land
x0csm = np.maximum(0,x0csm)

#************************** Open Boundary Conditions *************************

diric = datadir + 'obc_ecosystem/climatology/'
file_ext = '.15N_45N.115E_145E_clim_interp_depth_setmisstonn.nc'

var_name = ['PO4', 'NO3', 'SiO3', 'NH4', 'Fe', 'Lig', 'DOC', 'DON', \
            'DOP', 'DOPr', 'DONr', 'DOCr', 'zooC', 'spC', 'spP', 'spChl',\
            'spFe', 'spCaCO3', 'diatC', 'diatChl', 'diatSi', 'diatFe',\
            'diatP', 'diazC', 'diazChl', 'diazFe', 'diazP']

var_index = [po4_ind, no3_ind, sio3_ind, nh4_ind, fe_ind, lig_ind,\
             doc_ind, don_ind, dop_ind, dopr_ind, donr_ind, docr_ind,\
             zooC_ind, spC_ind, spP_ind, spChl_ind, spFe_ind, spCaCO3_ind,\
            diatC_ind, diatChl_ind, diatSi_ind, diatFe_ind, diatP_ind,\
            diazC_ind, diazChl_ind, diazFe_ind, diazP_ind]

ntracers = 27

## get size of array

filr = netCDF4.Dataset(diric+'Fe'+file_ext)
Nt, Nz, Ny, Nx = np.shape(filr.variables['Fe'][:,zli,1:,1:])
filr.close()

xclim = np.ma.zeros((ntracers,12,Nz,Ny,Nx))

for var,ind in zip(var_name,var_index):
    filr = netCDF4.Dataset(diric+var+file_ext)
    Nz = np.shape(filr.variables[var][:,zli,:,:])[1]
    if ind < 12:
        xclim[ind,:,:,:,:] = filr.variables[var][:,zli,1:,1:]*maskT[:,:,:]
    else:
        xclim[ind,:,:14,:,:] = filr.variables[var][:,:14,1:,1:]*maskT[:14,:,:]
    filr.close()

# obcs do not contain any values below zero (done with cdo)


## put monthly data into daily
importlib.reload(utilities)
nvrs,nt,nz,nla,nlo = np.shape(xclim[:,:,:,:,:])

x_westg  = utilities.intmonth(xclim[:,:,:,:,0],nvrs=nvrs, nt=nt, nz=nz,nla=nla,nlo=0,ndim=4)
x_eastg  = utilities.intmonth(xclim[:,:,:,:,-1],nvrs=nvrs, nt=nt, nz=nz,nla=nla,nlo=0,ndim=4)
x_northg = utilities.intmonth(xclim[:,:,:,-1,:],nvrs=nvrs, nt=nt, nz=nz,nla=0,nlo=nlo,ndim=4)
x_southg = utilities.intmonth(xclim[:,:,:,0,:],nvrs=nvrs, nt=nt, nz=nz,nla=0,nlo=nlo,ndim=4)


# remove large values
x_westg = np.ma.where(x_westg > 1000000,0,x_westg)
x_eastg = np.ma.where(x_eastg > 1000000,0,x_eastg)
x_northg = np.ma.where(x_northg > 1000000,0,x_northg)
x_southg = np.ma.where(x_southg > 1000000,0,x_southg)


# ******************* Create Depth coordinate *****************************
# delta_z
dz = np.zeros(len(z_t))
dz[0] = 1000
dz[1:-1] = (z_t[2:] - z_t[1:-1])/2 + (z_t[1:-1] - z_t[:-2])/2
dz[-1] = (z_t[-1] - z_t[-2])

# z depth
zw = z_t + dz/2 # depth of bottom of cell

# *********************** Grid spacing and inverse grid spacing **********************
### create spatial differentials
# dx, dy and dz
importlib.reload(utilities)
nt, nz, Ny, Nx = np.shape(uvelg[:,:,:,:])
dy = 11100*100; #cm
dxs = utilities.x_weights(latre, dy, Nx, Ny)
dzs2D =np.repeat(dz[:,np.newaxis], Ny, axis=1)
dzs = np.repeat(dzs2D[:,:,np.newaxis],Nx,axis=2)
dxs = dxs[:,:]
dzs = dzs[:,:,:]
idx = 1.0 / dxs
idy = 1.0 / dy
idz = 1.0/ dzs


# ************************* Calculate vertical velocity ****************************

importlib.reload(continuity)
ndaysv = 365
wcontg =  continuity.wcont(uvelg[:,...], vvelg[:,...], idy, idx, dzs, maskT, mask_zonal, mask_merid, ndaysv)


# ************************* Calculate density if convective adjustment will be used ****************************


# import density
# nt,nz,ny,nx = np.shape(saltg)
# pressy =np.repeat(z_t[:,np.newaxis], ny, axis=1)
# pressxy =np.repeat(pressy[:,:,np.newaxis], nx, axis=2)
# rhok_adiab = np.zeros_like(saltg)
# rhok_adiab[:,:-1,:,:] = density.calculate_rho(tempg[:,:-1,:,:],saltg[:,:-1,:,:],pressxy[1:,:,:]) # temp, sal, press
# rhokp = np.zeros_like(saltg)
# rhokp[:,:-1,:,:] = density.calculate_rho(tempg[:,1:,:,:],saltg[:,1:,:,:],pressxy[1:,:,:]) # temp, sal, press

# rhok_adiab = rhok_adiab[:,:-1,:,:]*maskT[1:,:,:]
# rhokp = rhokp[:,:-1,:,:]*maskT[1:,:,:]





### *********************** INPUT DATA LIST *************************

# ******************** Remove lat and lon edges *********************************
# this step is necesray since when calculating W we lose one lat and lon

# slices for lat, lon and depth
sla = np.s_[1:-1]; slo = np.s_[1:-1]; zs = np.s_[:-1]

x_eastr = x_eastg[:,:,zs,sla]
x_westr = x_westg[:,:,zs,sla]
x_northr = x_northg[:,:,zs,slo]
x_southr = x_southg[:,:,zs,slo]
x0csmr = x0csm[:,zs,sla,slo]
uvelr = uvelg[:,zs,sla,slo]
vvelr  = vvelg[:,zs,sla,slo]

qswg = qswg[:,sla,slo]
tempr = tempg[:,zs,sla,slo]
feflux_sedg = feflux_sedg[zs,sla,slo]
feflux_ventg = feflux_ventg[zs,sla,slo]
dustg = dustg[:,sla,slo]
feflux_solg = feflux_solg[:,sla,slo]
NOxg = NOxg[:,sla,slo]
NHyg = NHyg[:,sla,slo]

idxr = idx[sla,slo]
idyr = idy
idzr = idz[zs,sla,slo]
dzsr = dzs[zs,sla,slo]
z_tr = z_t[zs]
zwr = zw[zs]
mask_zonalr = mask_zonal[zs,sla,slo]
mask_meridr = mask_merid[zs,sla,slo]
maskTr = maskT[zs,sla,slo]
mask_kmtr = mask_kmt[zs,sla,slo]
mask_reminr = mask_remin[zs,sla,slo]
wcontr = wcontg[:,zs,sla,slo]

# if density is used
# rhokp_adiabr = rhok_adiab[:,zs,sla,slo]
# rhokpr = rhokp[:,zs,sla,slo]

### *********************** CALL MAIN *************************

%%time

# reload if changes where made after importing
#importlib.reload(main)
#importlib.reload(adv_diff)
#importlib.reload(ecosystem)
#importlib.reload(run)
#importlib.reload(utilities)


test    = 'test2'
dirsave = '/Users/angeles/Desktop/python/python2/Barakawa/moore/'+test+'/'
caract  ='.bio.chen.2e7.025.noconvad.'

vlevsi = 37
ntracers = 27
nvarstot = 40
nn = 1
ndaysv = 365

#for nn in range(3,5):
dti = 800
dtiprev = 1200

if nn == 0:
    filespin = None
    tstarti = 0
    xprevi = x0csmr[:ntracers,:,:,:]
    xcuri = 0
    ncuri = 1
    initi = True


elif nn != 0 :
    nin = '0'+str(nn)
    filespin = dirsave+'spinup.'+test+caract+nin+'.npz'

    data = np.load(filespin)
    ncuri = data['ncur']


    xprevi  = data['xprev'][:ntracers,:,:,:]
    xcuri  = data['xcur'][:ntracers,:,:,:]

    tstarti = np.load(filespin)['curtime']+dtiprev

    print('ble',tstarti)
    initi = False
#     if nn > 11:
#         vlevsi = 60
#     else:
#         vlevsi = 10

nout = '0'+str(nn+1)

tendi = tstarti+100*fcd

print(nn, 'input = ', filespin, tstarti, tendi)

nvars,Nz,Ny,Nx = np.shape(x0csmr[:,:,:,:])

run = main.eco_model(f_bio = ecosystem.dxdt_ecosys, \
      f_phys = adv_diff.adv_diff_Barakawa,\
      nvarstot = nvarstot, ntracers = ntracers, xprev = xprevi, xcur = xcuri, ncur = ncuri,\
      dt = dti, tstart = tstarti, tend = tendi, Navg=5,  ndays=365, ndayse=ndaysv,\
      x_westi = x_westr[:ntracers,:,:,:].copy(), x_easti = x_eastr[:ntracers,:,:,:].copy(),\
      x_northi = x_northr[:ntracers,:,:,:].copy(), x_southi = x_southr[:ntracers,:,:,:].copy(), tauini=5,\
      no3res = np.zeros_like(x0csm[0,:,:,:]), po4res = np.zeros_like(x0csm[0,:,:,:]),\
      sio3res = np.zeros_like(x0csm[0,:,:,:]), invtau = np.zeros_like(x0csm[0,:,:,:]),\
      Ui = uvelr, Vi = vvelr, wi = wcontr,\
      tempi=tempr[:,:,:,:], ficei=0, qswi=qswg[:,:,:], \
      feflux_sedi=feflux_sedg[:,:,:], feflux_venti=feflux_ventg[:,:,:],\
      dusti=dustg[:,:,:], feflux_soli = feflux_solg[:,:,:],\
      NOxi = NOxg[:,:,:], NHyi = NHyg[:,:,:],\
      hordiff = 2e7, kvmix = 0.25,\
      rhokpi = 0,rho_adiabi = 0,\
      idx = idxr, idy = idyr, idz = idzr, dzs = dzsr, z_t = z_tr, zw = zwr,\
      mask_zonal = mask_zonalr, mask_merid = mask_meridr, \
      maskt = maskTr, mask_kmt = mask_kmtr,mask_remin = mask_reminr,\
      BIO = True, convad = False, restore = False, OBC = True,\
      file_create_spinup=dirsave+'spinup.'+test+caract+nout+'.npz', initial=initi)


run_saver.saver(run, vlevs=vlevsi, nvars=nvarstot, dt=dti, Nx=Nx, Ny=Ny, Nz=Nz, tstart=tstarti, tend=tendi,
          time_step = 1*fcd, filename=dirsave+test+caract+nout+'.npz')
print('output = '+ dirsave+test+caract+nout+'.npz')
