import xarray as xr
from plot_helpers import plotvar_scalar2D, plotvar_vector2D, plot_stat, plot_rawstat
import numpy as np

DS = xr.open_dataset('output.nc')
DS.load()

M = DS.mass
KE = DS.kinetic_energy
PE = DS.potential_energy
TE = DS.total_energy
PV = DS.potential_vorticity
PENS = DS.potential_enstrophy

plot_stat('mass', M)
plot_stat('total_energy', TE)
plot_rawstat('kinetic_energy', KE)
plot_rawstat('potential_energy', PE)
plot_stat('pv', PV)
plot_stat('pens', PENS)


h = DS.h
v = DS.v
q = DS.q

Nlist = np.arange(21)

# FIX THIS STUFF UP...
# Probably control what to plot based on the test case run...

for i in Nlist:
    plotvar_scalar2D('q', q.isel(t=i,q_ndofs=0,ncells_z=0),i)
    plotvar_scalar2D('h', h.isel(t=i,h_ndofs=0,ncells_z=0),i)
    plotvar_vector2D('v', v.isel(t=i,v_ndofs=0,ncells_z=0), v.isel(t=i,v_ndofs=1,ncells_z=0),i)
