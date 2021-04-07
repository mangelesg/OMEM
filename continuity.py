import numpy as np
def wcont(U, V, idy, idx, dzs, maskT, mask_zonal, mask_merid, ndays):

    Nz,Ny,Nx = np.shape(U[0,:,:,:])
    wcont = np.ma.zeros((ndays,Nz,Ny,Nx))



    # Index the grids. (m)inus is index-1 and (p)lus is index+1
    i = j  =  np.s_[1:-1]
    ip = jp = kp = np.s_[2:]
    im = jm = km = np.s_[:-2]
    nz = Nz - 1

    U_east = (U[:,:,j,i] + U[:,:,jm,i])*0.5#*mask_zonal[:,j,i]
    U_west = (U[:,:,j,im]+ U[:,:,jm,im])*0.5#*mask_zonal[:,j,im]
    V_north = (V[:,:,j,i] + V[:,:,j,im])*0.5#*mask_merid[:,j,i]
    V_south = (V[:,:,jm,i] + V[:,:,jm,im])*0.5#*mask_merid[:,jm,i]



    for nt in range(ndays):
       # wcont[nt,0,j,i] = wvelTd[nt,0,j,i]*maskT[0,j,i]

        wcont[nt,1,j,i] = (wcont[nt,0,j,i]\
        + dzs[0,j,i]*idx[j,i]*(U_east[nt,0,:,:] - U_west[nt,0,:,:] )\
        + dzs[0,j,i]*idy*( V_north[nt,0,:,:] - V_south[nt,0,:,:]))*maskT[1,j,i]

        for k in range(1,nz):
#            dzk = (dzs[k,j,i] +dzs[k-1,j,i]*0.5 + dzs[k+1,j,i]*0.5)
#
#            wcont[nt,k+1,j,i] = (wcont[nt,k-1,j,i]\
#            + dzk*idx[j,i]*( U_east[nt,k,:,:] - U_west[nt,k,:,:])\
#                        + dzk*idy*( V_north[nt,k,:,:]- V_south[nt,k,:,:]))*maskT[k+1,j,i]
   
            dzk = dzs[k,j,i] #+dzs[k-1,j,i]*0.5 + dzs[k+1,j,i]*0.5)

            wcont[nt,k+1,j,i] = (wcont[nt,k,j,i]\
            + dzk*idx[j,i]*( U_east[nt,k,:,:] - U_west[nt,k,:,:])\
                        + dzk*idy*( V_north[nt,k,:,:]- V_south[nt,k,:,:]))*maskT[k+1,j,i]
   
   
   
    return(wcont)
