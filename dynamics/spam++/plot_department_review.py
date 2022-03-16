import xarray as xr
from plot_helpers import plotvar_scalar2D, plotvar_vector2D, plot_stat, plot_rawstat
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib

DS_weno = xr.open_dataset('WENOFUNC9-DualWENOFUNC9-100/output.nc')
DS_cfv = xr.open_dataset('CFV10-DualCFV10-100/output.nc')
DS_weno.load()
DS_cfv.load()


energy_weno = DS_weno.energy
energy_cfv = DS_cfv.energy


matplotlib.rcParams.update({'font.size': 30})

plt.figure(figsize=(10,8))
plt.plot( (energy_cfv.isel(energy_ndofs=2) - energy_cfv.isel(energy_ndofs=2)[0])/energy_cfv.isel(energy_ndofs=2)[0]*100., label='cfv')
plt.plot( (energy_weno.isel(energy_ndofs=2) - energy_weno.isel(energy_ndofs=2)[0])/energy_weno.isel(energy_ndofs=2)[0]*100., label='weno')
plt.legend()
plt.xlabel('Nsteps')
plt.ylabel('Fractional Change in Energy')
plt.tight_layout()
plt.savefig('energy.png')



# Probably control what to plot based on the test case run...

densl_weno = DS_weno.densl
densl_cfv = DS_cfv.densl


plt.figure(figsize=(10,8))
plt.contourf(densl_cfv.isel(t=20,densl_ndofs=1,ncells_z=0))
plt.colorbar()
plt.contour(densl_cfv.isel(t=20,densl_ndofs=1,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('sl_cfv.png')
plt.close('all')

plt.figure(figsize=(10,8))
plt.contourf(densl_weno.isel(t=20,densl_ndofs=1,ncells_z=0))
plt.colorbar()
plt.contour(densl_weno.isel(t=20,densl_ndofs=1,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('sl_weno.png')
plt.close('all')
