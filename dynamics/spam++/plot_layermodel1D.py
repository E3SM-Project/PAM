import xarray as xr
from plot_helpers import plotvar_scalar1D, plot_stat, plot_rawstat
import numpy as np
import sys

DS = xr.open_dataset('output.nc')
DS.load()

ndensity = DS.dims['dens_ndofs']
ndensity_fct = DS.dims.get('densfct_ndofs',0)
nt = DS.dims['t']
model = sys.argv[1]

mass = DS.mass
energy = DS.energy
densmax = DS.densmax
densmin = DS.densmin
if ndensity_fct > 0:
    massfct = DS.massfct
    densfctmax = DS.densfctmax
    densfctmin = DS.densfctmin

#maybe mass/tracer stuff gets a list of names of size ndensity/ndensityfct?
#YES DO IT LIKE THIS!

if model == 'swe':
    dens_names = ['h',]
    densfct_names = []
    dens_stat_names = ['mass',]
    densfct_stat_names = []
    nprogdens = 1
    nprogdensfct = 0  
if model == 'tswe':
    dens_names = ['h','S',]
    densfct_names = []
    dens_stat_names = ['mass','bouyancy',]
    densfct_stat_names = []
    nprogdens = 2
    nprogdensfct = 0  

for k in range(ndensity-nprogdens):
    dens_names.append('T'+str(k))
    dens_stat_names.append('tracer'+str(k))
for k in range(ndensity_fct-nprogdensfct):
    densfct_names.append('Tfct'+str(k))        
    densfct_stat_names.append('tracerfct'+str(k))  

for l,name in zip(range(ndensity), dens_stat_names):
    plot_stat('total_'+name, mass.isel(mass_ndofs=l))
    plot_rawstat('min_' +name, densmin.isel(densmin_ndofs=l))
    plot_rawstat('max_' +name, densmax.isel(densmax_ndofs=l))
    
if ndensity_fct > 0:
    for l,name in zip(range(ndensity_fct), densfct_stat_names):
        plot_stat('total_'+name, massfct.isel(massfct_ndofs=l))
        plot_rawstat('min_' +name, densfctmin.isel(densfctmin_ndofs=l))
        plot_rawstat('max_' +name, densfctmax.isel(densfctmax_ndofs=l))
                
plot_rawstat('internal_energy', energy.isel(energy_ndofs=3))
plot_rawstat('potential_energy', energy.isel(energy_ndofs=2))
plot_rawstat('kinetic_energy', energy.isel(energy_ndofs=1))
plot_stat('total_energy', energy.isel(energy_ndofs=0))

Nlist = np.arange(0,nt)

v = DS.v
dens = DS.dens
densl = DS.densl
hs = DS.hs

if (ndensity_fct>0):
    densfct = DS.densfct
    densfctl = DS.densfctl

plotvar_scalar1D('hs', hs.isel(hs_ndofs=0,dual_nlayers=0),0)


# WHERE EXACTLY SHOULD sl/trl/trfctl live? are they straight 0-forms? twisted n-forms?
for i in Nlist:
    plotvar_scalar1D('v', v.isel(t=i,v_ndofs=0,primal_ncells_y=0,primal_nlayers=0),i)
    for l,name in zip(range(ndensity), dens_names):
            plotvar_scalar1D(name, dens.isel(t=i,dens_ndofs=l,dual_ncells_y=0,dual_nlayers=0),i)
            plotvar_scalar1D(name+'l', densl.isel(t=i,densl_ndofs=l,primal_ncells_y=0,primal_nlayers=0),i)
    if ndensity_fct > 0:
        for l,name in zip(range(ndensity_fct), densfct_names):
            plotvar_scalar1D(name, densfct.isel(t=i,densfct_ndofs=l,dual_ncells_y=0,dual_nlayers=0),i)
            plotvar_scalar1D(name+'l', densfctl.isel(t=i,densfctl_ndofs=l,primal_ncells_y=0,primal_nlayers=0),i)


