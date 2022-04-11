import xarray as xr
from plot_helpers import plotvar_scalar2D, plotvar_vector2D, plot_stat, plot_rawstat
import numpy as np
import sys

DS = xr.open_dataset('test_dycore.nc')
DS.load()

ndensity = DS.dims['dens_ndofs']
nt = DS.dims['t']
model = sys.argv[1]
nens = DS.dims['nens']

mass = DS.mass
energy = DS.energy
pens = DS.pens
pv = DS.pv

densmax = DS.densmax
densmin = DS.densmin

#maybe mass/tracer stuff gets a list of names of size ndensity/ndensityfct?
#YES DO IT LIKE THIS!

if model == 'swe':
    dens_names = ['h',]
    densfct_names = []
    dens_stat_names = ['mass',]
    nprogdens = 1
if model == 'tswe':
    dens_names = ['h','S',]
    dens_stat_names = ['mass','bouyancy',]
    nprogdens = 2

if model == 'ce':
    dens_names = ['rho','Theta',]
    dens_stat_names = ['mass','entropic_var_density',]
    nprogdens = 2

#THIS BREAKS FOR RHOD VARIANTS...
#probably ok, this is just a quick and dirty plotting script...
# MAYBE MODEL ARGUMENT BECOMES MCED THOUGH?
if model == 'mce':
    dens_names = ['rho','Theta','rho_v', 'rho_l', 'rho_i']
    dens_stat_names = ['mass','entropic_var_density','vapor', 'liquid', 'ice']
    nprogdens = 5
    
for k in range(ndensity-nprogdens):
    dens_names.append('T'+str(k))
    dens_stat_names.append('tracer'+str(k))

for n in range(nens):
    for l,name in zip(range(ndensity), dens_stat_names):
        plot_stat('total_'+name + '.' + str(n), mass.isel(mass_ndofs=l, nens=n))
        plot_rawstat('min_' +name + '.' + str(n), densmin.isel(densmin_ndofs=l, nens=n))
        plot_rawstat('max_' +name + '.' + str(n), densmax.isel(densmax_ndofs=l, nens=n))
		        
    plot_rawstat('internal_energy.' + str(n), energy.isel(energy_ndofs=3, nens=n))
    plot_rawstat('potential_energy.'+ str(n), energy.isel(energy_ndofs=2, nens=n))
    plot_rawstat('kinetic_energy.'+ str(n), energy.isel(energy_ndofs=1, nens=n))
    plot_stat('total_energy.'+ str(n), energy.isel(energy_ndofs=0, nens=n))

    plot_stat('total_pens.'+ str(n), pens.isel(pens_ndofs=0, nens=n))
    plot_stat('total_pv.'+ str(n), pv.isel(pv_ndofs=0, nens=n))


Nlist = np.arange(0,nt)

v = DS.v
dens = DS.dens
q = DS.q
densl = DS.densl
hs = DS.hs
coriolis = DS.coriolis

for n in range(nens):

    plotvar_scalar2D('hs.' + str(n), hs.isel(hs_ndofs=0,dual_nlayers=0,nens=n),0)
    plotvar_scalar2D('coriolis.' + str(n), coriolis.isel(coriolis_ndofs=0,primal_nlayers=0,nens=n),0)

    for i in Nlist:
        plotvar_scalar2D('q.' + str(n), q.isel(t=i,q_ndofs=0,dual_nlayers=0,nens=n),i)
        plotvar_vector2D('v.' + str(n), v.isel(t=i,v_ndofs=0,primal_nlayers=0,nens=n), v.isel(t=i,v_ndofs=1,primal_nlayers=0,nens=n),i)
        for l,name in zip(range(ndensity), dens_names):
            plotvar_scalar2D(name +'.' + str(n), dens.isel(t=i,dens_ndofs=l,dual_nlayers=0,nens=n),i)
            plotvar_scalar2D(name+'l.' + str(n), densl.isel(t=i,densl_ndofs=l,primal_nlayers=0,nens=n),i)
#THIS ASSUMES TOTAL DENSITY IS IN DENS(0)
            plotvar_scalar2D(name+'c.' + str(n), dens.isel(t=i,dens_ndofs=l,dual_nlayers=0,nens=n) / dens.isel(t=i,dens_ndofs=0,dual_nlayers=0,nens=n),i)

