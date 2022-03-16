import xarray as xr
from plot_helpers import plotvar_scalar2D, plotvar_vector2D, plot_stat, plot_rawstat
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Circle

DS_weno = xr.open_dataset('WENOFUNC9-DualWENOFUNC9-100/output.nc')
DS_cfv = xr.open_dataset('CFV10-DualCFV10-100/output.nc')
DS_weno.load()
DS_cfv.load()

matplotlib.rcParams.update({'font.size': 30})

# Probably control what to plot based on the test case run...

densl_weno = DS_weno.densl
densl_cfv = DS_cfv.densl


plt.figure(figsize=(20,8))

plt.subplot(1,2,1)

plt.contourf(densl_cfv.isel(t=20,densl_ndofs=1,ncells_z=0),  np.linspace(9.76,10.32,15))
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.text(60,20,'spurious\noscillations', fontdict={'color':'r'})
plt.title('CFV', fontdict={'color':'r', 'size':40})
c = plt.Circle((50, 50), 15, facecolor='none', edgecolor='r', linewidth=5)
plt.gca().add_artist(c)


plt.subplot(1,2,2)
plt.contourf(densl_weno.isel(t=20,densl_ndofs=1,ncells_z=0),  np.linspace(9.76,10.32,15))
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.text(60,20,'oscillation-\nlimited', fontdict={'color':'r'})
plt.title('WENO', fontdict={'color':'r', 'size':40})
c = plt.Circle((50, 50), 15, facecolor='none', edgecolor='r', linewidth=5)
plt.gca().add_artist(c)

plt.tight_layout()
plt.savefig('ECP-tswe.png',transparent=True)
plt.close('all')





DS_weno = xr.open_dataset('MRB/MRB100-WENOFUNC9-HODGE2/output.nc')
DS_cfv = xr.open_dataset('MRB/MRB100-CFV9-HODGE2/output.nc')
DS_weno.load()
DS_cfv.load()

matplotlib.rcParams.update({'font.size': 30})

# Probably control what to plot based on the test case run...

densl_weno = DS_weno.densl
densl_cfv = DS_cfv.densl


plt.figure(figsize=(20,8))
plt.subplot(1,2,1)
plt.contourf(densl_cfv.isel(t=9,densl_ndofs=1,primal_ncells_z=0)- 300, np.linspace(-0.08,0.64,19))
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.text(35,65,'spurious\noscillations', fontdict={'color':'r'})
plt.title('CFV', fontdict={'color':'r', 'size':40})
c = plt.Circle((25, 100), 15, facecolor='none', edgecolor='r', linewidth=5)
plt.gca().add_artist(c)


plt.subplot(1,2,2)
plt.contourf(densl_weno.isel(t=9,densl_ndofs=1,primal_ncells_z=0) - 300, np.linspace(-0.08,0.64,19))
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.text(35,65,'oscillation-\nlimited', fontdict={'color':'r'})
plt.title('WENO', fontdict={'color':'r', 'size':40})
c = plt.Circle((25, 100), 15, facecolor='none', edgecolor='r', linewidth=5)
plt.gca().add_artist(c)

plt.tight_layout()
plt.savefig('ECP-mrb-s.png',transparent=True)
plt.close('all')


densfctl_weno = DS_weno.densfctl
densfctl_cfv = DS_cfv.densfctl

plt.figure(figsize=(20,8))
plt.subplot(1,2,1)
plt.contourf(densfctl_cfv.isel(t=9,densfctl_ndofs=0,primal_ncells_z=0), np.linspace(0,0.028,15))
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.text(35,65,'spurious\noscillations', fontdict={'color':'r'})
plt.title('CFV', fontdict={'color':'r', 'size':40})
c = plt.Circle((25, 100), 15, facecolor='none', edgecolor='r', linewidth=5)
plt.gca().add_artist(c)


plt.subplot(1,2,2)
plt.contourf(densfctl_weno.isel(t=9,densfctl_ndofs=0,primal_ncells_z=0), np.linspace(0,0.028,15))
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.text(35,65,'oscillation-\nlimited', fontdict={'color':'r'})
plt.title('WENO', fontdict={'color':'r', 'size':40})
c = plt.Circle((25, 100), 15, facecolor='none', edgecolor='r', linewidth=5)
plt.gca().add_artist(c)

plt.tight_layout()
plt.savefig('ECP-mrb-qv.png',transparent=True)
plt.close('all')

