
import numpy as np
import sys
import matplotlib.pyplot as plt
from netCDF4 import Dataset

nc = Dataset(sys.argv[1],"r")
nt = nc.dimensions["t"].size
nx = nc.dimensions["x"].size
ny = nc.dimensions["y"].size
nz = nc.dimensions["z"].size
x        = nc.variables["x"][:]
y        = nc.variables["y"][:]
z        = nc.variables["z"][:]
rho      = nc.variables["dens_pert"]    [nt-1,:,ny/2-1,:]
u        = nc.variables["u"]            [nt-1,:,ny/2-1,:]
v        = nc.variables["v"]            [nt-1,:,ny/2-1,:]
w        = nc.variables["w"]            [nt-1,:,ny/2-1,:]
theta    = nc.variables["pot_temp_pert"][nt-1,:,ny/2-1,:]

xx,zz = np.meshgrid(x,z)

levs = np.arange(-2,14,0.5)
plt.contourf(xx,zz,theta,cmap="jet",levels=levs)
plt.colorbar(orientation="horizontal")
ax = plt.gca()
ax.set_aspect("equal")
plt.savefig("theta.png",dpi=600)




