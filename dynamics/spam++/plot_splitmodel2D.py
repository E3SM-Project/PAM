import xarray as xr
from plot_helpers import plotvar_scalar2D, plotvar_vector2D, plot_stat, plot_rawstat
import numpy as np
import sys

DS = xr.open_dataset('output.nc')
DS.load()

# SHOULD AUTO LOAD FROM input-layermodel2D.txt
Nlist = np.arange(21)

model = sys.argv[1]
ntracers = int(sys.argv[2])
ntracers_fct = int(sys.argv[3])
if (model == 'swe'): ndensity = 1 + ntracers
if (model == 'tswe'): ndensity = 2 + ntracers

mass = DS.mass
energy = DS.energy
pens = DS.pens
pv = DS.pv

densmax = DS.densmax
densmin = DS.densmin

if ntracers_fct > 0:
    massfct = DS.massfct
    densfctmax = DS.densfctmax
    densfctmin = DS.densfctmin


plot_stat('mass', mass.isel(mass_ndofs=0))
if (model == 'tswe'):
    plot_stat('bouyancy', mass.isel(mass_ndofs=1))

plot_stat('total_energy', energy.isel(energy_ndofs=2))
plot_rawstat('kinetic_energy', energy.isel(energy_ndofs=1))
plot_rawstat('potential_energy', energy.isel(energy_ndofs=0))

plot_stat('pens', pens.isel(pens_ndofs=0))
plot_stat('pv', pv.isel(pv_ndofs=0))

for l in range(ndensity-ntracers,ndensity):
    plot_stat('tracer' + str(l-ntracers), mass.isel(mass_ndofs=l))
    plot_rawstat('tracer' + str(l-ntracers) + 'min', densmin.isel(densmin_ndofs=l))
    plot_rawstat('tracer' + str(l-ntracers) + 'max', densmax.isel(densmax_ndofs=l))

if ntracers_fct > 0:
    for l in range(ntracers):
        plot_stat('tracerfct' + str(l), massfct.isel(massfct_ndofs=l))
        plot_rawstat('tracerfct' + str(l) + 'min', densfctmin.isel(densfctmin_ndofs=l))
        plot_rawstat('tracerfct' + str(l) + 'max', densfctmax.isel(densfctmax_ndofs=l))


# Probably control what to plot based on the test case run...

v = DS.v
dens = DS.dens
q = DS.q
densl = DS.densl

if (ntracers_fct>0):
    densfct = DS.densfct
    densfctl = DS.densfctl


for i in Nlist:
    plotvar_scalar2D('q', q.isel(t=i,q_ndofs=0,ncells_z=0),i)
    plotvar_scalar2D('h', dens.isel(t=i,dens_ndofs=0,ncells_z=0),i)
    if (model == 'tswe'):    
        plotvar_scalar2D('S', dens.isel(t=i,dens_ndofs=1,ncells_z=0),i)
        plotvar_scalar2D('sl', densl.isel(t=i,densl_ndofs=1,ncells_z=0),i)
    for l in range(ndensity-ntracers,ndensity):
        plotvar_scalar2D('tr' + str(l-ntracers), dens.isel(t=i,dens_ndofs=l,ncells_z=0),i)
        plotvar_scalar2D('trl' + str(l-ntracers), densl.isel(t=i,densl_ndofs=l,ncells_z=0),i)
    if ntracers_fct > 0:
        for l in range(ntracers_fct):
            plotvar_scalar2D('trfct' + str(l), densfct.isel(t=i,densfct_ndofs=l,ncells_z=0),i)
            plotvar_scalar2D('trfctl' + str(l), densfctl.isel(t=i,densfctl_ndofs=l,ncells_z=0),i)
    plotvar_vector2D('v', v.isel(t=i,v_ndofs=0,ncells_z=0), v.isel(t=i,v_ndofs=1,ncells_z=0),i)

