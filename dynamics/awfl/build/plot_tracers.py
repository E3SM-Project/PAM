from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

nc = Dataset("dry_tracers.nc","r")
uniform = nc.variables["tracer_uniform"][:,:,0,:]
theta   = nc.variables["tracer_theta"  ][:,:,0,:]
block   = nc.variables["tracer_block"  ][:,:,0,:]
nt = nc.variables["t"][:].shape[0]
nz = nc.variables["z"][:].shape[0]

X,Z=np.meshgrid(nc.variables["x"][:],nc.variables["z"][:])
X = X / 1000;
Z = Z / 1000;

plt.contourf(X,Z,uniform[nt-1,:,:],cmap="jet")
plt.colorbar(shrink=0.5)
plt.axis('scaled')
plt.tight_layout()
plt.xlabel("x-coordinate (km)")
plt.ylabel("z-coordinate (km)")
plt.savefig("uniform_final.png", bbox_inches = 'tight', pad_inches=0.05)
plt.close()

plt.contourf(X,Z,theta[0,:,:],cmap="jet")
plt.colorbar(shrink=0.5)
plt.axis('scaled')
plt.tight_layout()
plt.xlabel("x-coordinate (km)")
plt.ylabel("z-coordinate (km)")
plt.savefig("theta_init.png", bbox_inches = 'tight', pad_inches=0.05)
plt.close()

plt.contourf(X,Z,theta[nt-1,:,:],cmap="jet")
plt.colorbar(shrink=0.5)
plt.axis('scaled')
plt.tight_layout()
plt.xlabel("x-coordinate (km)")
plt.ylabel("z-coordinate (km)")
plt.savefig("theta_final.png", bbox_inches = 'tight', pad_inches=0.05)
plt.close()

plt.contourf(X,Z,block[0,:,:],cmap="jet")
plt.colorbar(shrink=0.5)
plt.axis('scaled')
plt.tight_layout()
plt.xlabel("x-coordinate (km)")
plt.ylabel("z-coordinate (km)")
plt.savefig("block_init.png", bbox_inches = 'tight', pad_inches=0.05)
plt.close()

plt.contourf(X,Z,block[nt-1,:,:],cmap="jet")
plt.colorbar(shrink=0.5)
plt.axis('scaled')
plt.tight_layout()
plt.xlabel("x-coordinate (km)")
plt.ylabel("z-coordinate (km)")
plt.savefig("block_final.png", bbox_inches = 'tight', pad_inches=0.05)
plt.close()

plt.plot(nc.variables["x"][:]/1000,block[nt-1,nz/2,:])
plt.plot(nc.variables["x"][:]/1000,block[nt-1,nz/2,:],'bo')
plt.tight_layout()
plt.xlabel("x-coordinate (km)")
plt.ylabel("Block tracer value at z==10km")
plt.savefig("block_slice.png", bbox_inches = 'tight', pad_inches=0.05)
plt.close()

