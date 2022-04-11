import xarray as xr
from plot_helpers import plotvar_scalar2D, plotvar_vector2D, plot_stat, plot_rawstat
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Circle

DS_weno = xr.open_dataset('WENO9-DualWENO9-100/output.nc')
DS_wenofunc = xr.open_dataset('WENOFUNC9-DualWENOFUNC9-100/output.nc')
DS_cfv = xr.open_dataset('CFV10-DualCFV10-100/output.nc')
DS_weno.load()
DS_wenofunc.load()
DS_cfv.load()

matplotlib.rcParams.update({'font.size': 30})

# Probably control what to plot based on the test case run...

densl_weno = DS_weno.densl
densl_wenofunc = DS_wenofunc.densl
densl_cfv = DS_cfv.densl


plt.figure(figsize=(30,8))

plt.subplot(1,3,1)

plt.contourf(densl_cfv.isel(t=20,densl_ndofs=1,ncells_z=0))
plt.colorbar()
#plt.contour(densl_cfv.isel(t=20,densl_ndofs=1,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.text(60,20,'spurious\noscillations', fontdict={'color':'r'})
plt.title('CFV', fontdict={'color':'r', 'size':40})
c = plt.Circle((50, 50), 15, facecolor='none', edgecolor='r', linewidth=5)
plt.gca().add_artist(c)


plt.subplot(1,3,2)
plt.contourf(densl_weno.isel(t=20,densl_ndofs=1,ncells_z=0))
plt.colorbar()
#plt.contour(densl_weno.isel(t=20,densl_ndofs=1,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.text(60,20,'oscillation-\nlimited', fontdict={'color':'r'})
plt.title('WENO', fontdict={'color':'r', 'size':40})
c = plt.Circle((50, 50), 15, facecolor='none', edgecolor='r', linewidth=5)
plt.gca().add_artist(c)

plt.subplot(1,3,3)
plt.contourf(densl_wenofunc.isel(t=20,densl_ndofs=1,ncells_z=0))
plt.colorbar()
#plt.contour(densl_wenofunc.isel(t=20,densl_ndofs=1,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.text(60,20,'oscillation-\nlimited', fontdict={'color':'r'})
plt.title('WENOFUNC', fontdict={'color':'r', 'size':40})
c = plt.Circle((50, 50), 15, facecolor='none', edgecolor='r', linewidth=5)
plt.gca().add_artist(c)

plt.tight_layout()
plt.savefig('esmd.png',transparent=True)
plt.close('all')

energy_weno = DS_weno.energy
energy_wenofunc = DS_wenofunc.energy
energy_cfv = DS_cfv.energy


matplotlib.rcParams.update({'font.size': 30})

plt.figure(figsize=(10,8))
plt.plot( (energy_cfv.isel(energy_ndofs=2) - energy_cfv.isel(energy_ndofs=2)[0])/energy_cfv.isel(energy_ndofs=2)[0]*100., label='cfv')
plt.plot( (energy_weno.isel(energy_ndofs=2) - energy_weno.isel(energy_ndofs=2)[0])/energy_weno.isel(energy_ndofs=2)[0]*100., label='weno')
plt.plot( (energy_wenofunc.isel(energy_ndofs=2) - energy_wenofunc.isel(energy_ndofs=2)[0])/energy_wenofunc.isel(energy_ndofs=2)[0]*100., label='wenofunc')
plt.legend()
plt.xlabel('Nsteps')
plt.ylabel('Fractional Change in Energy')
plt.tight_layout()
plt.savefig('energy.png')

