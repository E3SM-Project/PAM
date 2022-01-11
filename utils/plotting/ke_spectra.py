import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

def rm(x,N):
    ret = x
    h = int(np.floor(N/2))
    nx = x.shape[0]
    for i in range(nx) :
        if (i-h<0) :
            ret[i] = np.average(x[0:i+1])
        elif (i-h+N>nx-1) :
            ret[i] = np.average(x[i:nx])
        else :
            ret[i] = np.average(x[i-h:i-h+N])
    return ret


def ke_spectra(uin,vin,win) :
    if (len(uin.shape) == 4) :
        u = np.mean( np.mean( np.mean(u4d , axis=0) , axis=0 ) , axis=0 )
        v = np.mean( np.mean( np.mean(u4d , axis=0) , axis=0 ) , axis=0 )
        w = np.mean( np.mean( np.mean(u4d , axis=0) , axis=0 ) , axis=0 )
    if (len(uin.shape) == 3) :
        u = np.mean( np.mean(u4d , axis=0) , axis=0 )
        v = np.mean( np.mean(u4d , axis=0) , axis=0 )
        w = np.mean( np.mean(u4d , axis=0) , axis=0 )
    if (len(uin.shape) == 2) :
        u = np.mean(u4d , axis=0)
        v = np.mean(u4d , axis=0)
        w = np.mean(u4d , axis=0)
    if (len(uin.shape) == 1) :
        u = u4d
        v = u4d
        w = u4d

    ke = (u*u + v*v + w*w) / 2
    spd = np.abs( np.fft.rfft(ke) )**2
    freq = np.fft.rfftfreq(len(u), d=500)
    return freq,spd
    

ncid = nc.Dataset("supercell_500m_ord3.nc_0.nc", "r")
nx = ncid.dimensions['x'].size
ny = ncid.dimensions['y'].size
nz = ncid.dimensions['z'].size
nt = ncid.dimensions['t'].size
t  = ncid.variables['t'][:]
u4d = ncid.variables['u'][0:nt-1,9,:,:]
v4d = ncid.variables['v'][0:nt-1,9,:,:]
w4d = ncid.variables['w'][0:nt-1,9,:,:]
freq,spd = ke_spectra(u4d,v4d,w4d)
plt.loglog( freq , rm(spd,3) )

ncid = nc.Dataset("supercell_500m_ord5.nc_0.nc", "r")
nx = ncid.dimensions['x'].size
ny = ncid.dimensions['y'].size
nz = ncid.dimensions['z'].size
nt = ncid.dimensions['t'].size
t  = ncid.variables['t'][:]
u4d = ncid.variables['u'][0:nt-1,:,:,:]
v4d = ncid.variables['v'][0:nt-1,:,:,:]
w4d = ncid.variables['w'][0:nt-1,:,:,:]
freq,spd = ke_spectra(u4d,v4d,w4d)
plt.loglog( freq , rm(spd,3) )

ncid = nc.Dataset("supercell_500m_ord7.nc_0.nc", "r")
nx = ncid.dimensions['x'].size
ny = ncid.dimensions['y'].size
nz = ncid.dimensions['z'].size
nt = ncid.dimensions['t'].size
t  = ncid.variables['t'][:]
u4d = ncid.variables['u'][0:nt-1,:,:,:]
v4d = ncid.variables['v'][0:nt-1,:,:,:]
w4d = ncid.variables['w'][0:nt-1,:,:,:]
freq,spd = ke_spectra(u4d,v4d,w4d)
plt.loglog( freq , rm(spd,3) )

ncid = nc.Dataset("supercell_500m_ord9.nc_0.nc", "r")
nx = ncid.dimensions['x'].size
ny = ncid.dimensions['y'].size
nz = ncid.dimensions['z'].size
nt = ncid.dimensions['t'].size
t  = ncid.variables['t'][:]
u4d = ncid.variables['u'][0:nt-1,:,:,:]
v4d = ncid.variables['v'][0:nt-1,:,:,:]
w4d = ncid.variables['w'][0:nt-1,:,:,:]
freq,spd = ke_spectra(u4d,v4d,w4d)
plt.loglog( freq , rm(spd,3) )

plt.legend(["3rd","5th","7th","9th"])
plt.show()

