import numpy as np

## to save the data in files every tme_step points
def saver(run, dt=0, tstart=0, tend=0, time_step = 7,\
              Nx=0, Ny=0, Nz=0, vlevs=0, nvars = 5, filename='filename'):
    times = np.arange(tstart, tend, dt)
    Nt = len(times)
    
    Nts = int(Nt*dt/time_step)     # Number of t points to save
    x_array = np.zeros((nvars,Nts,vlevs,Ny,Nx))
    t_array = np.zeros(Nts)

    count = 0
    for t, x in run:
        if not np.mod(t, time_step):
        
            x_array[:,count,:,:,:] = x[:,:vlevs,:,:]
            t_array[count] = t
            count = count + 1
  
    np.savez(filename, x=x_array,  t=t_array)
    del  t_array, x_array
