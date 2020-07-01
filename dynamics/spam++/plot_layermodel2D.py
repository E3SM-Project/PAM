import xarray as xr
from plot_helpers import plotvar_scalar2D, plotvar_vector2D, plot_stat, plot_rawstat
import numpy as np
import sys

DS = xr.open_dataset('output.nc')
DS.load()

Nlist = np.arange(21)
ntracers = int(sys.argv[2])
model = sys.argv[1]

mass = DS.mass
energy = DS.energy
pens = DS.pens
pv = DS.pv
if ntracers > 0:
    tracer = DS.massfct
    tracermin = DS.fctmin
    tracermax = DS.fctmax


plot_stat('mass', mass.isel(mass_ndofs=0))
if (model == 'tswe'):
    plot_stat('bouyancy', mass.isel(mass_ndofs=1))

plot_stat('total_energy', energy.isel(energy_ndofs=2))
plot_rawstat('kinetic_energy', energy.isel(energy_ndofs=1))
plot_rawstat('potential_energy', energy.isel(energy_ndofs=0))

plot_stat('pens', pens.isel(pens_ndofs=0))
plot_stat('pv', pv.isel(pv_ndofs=0))

for l in range(ntracers):
    plot_stat('tracer' + str(l), tracer.isel(massfct_ndofs=l))
    plot_rawstat('tracer' + str(l) + 'min', tracermin.isel(fctmin_ndofs=l))
    plot_rawstat('tracer' + str(l) + 'max', tracermax.isel(fctmax_ndofs=l))


# Probably control what to plot based on the test case run...

v = DS.v
dens = DS.dens
q = DS.q
densl = DS.densl
if (ntracers>0):
    densfct = DS.densfct
    densfctl = DS.densfctl

for i in Nlist:
    plotvar_scalar2D('q', q.isel(t=i,q_ndofs=0,ncells_z=0),i)
    plotvar_scalar2D('h', dens.isel(t=i,dens_ndofs=0,ncells_z=0),i)
    if (model == 'tswe'):    
        plotvar_scalar2D('S', dens.isel(t=i,dens_ndofs=1,ncells_z=0),i)
        plotvar_scalar2D('sl', densl.isel(t=i,densl_ndofs=1,ncells_z=0),i)
    for l in range(ntracers):
        plotvar_scalar2D('tr' + str(l), densfct.isel(t=i,densfct_ndofs=l,ncells_z=0),i)
        plotvar_scalar2D('trl' + str(l), densfctl.isel(t=i,densfctl_ndofs=l,ncells_z=0),i)
    plotvar_vector2D('v', v.isel(t=i,v_ndofs=0,ncells_z=0), v.isel(t=i,v_ndofs=1,ncells_z=0),i)

