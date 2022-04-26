import xarray as xr
from plot_helpers import plotvar_scalar2D, plotvar_vector2D, plot_stat, plot_rawstat
import numpy as np
import sys
import yaml


model = sys.argv[1]
nprocs = int(sys.argv[3])
with open(sys.argv[2], 'r') as file:
    config = yaml.safe_load(file)
    ncbase = config['out_prefix']
    crm_nx = config['crm_nx']
    crm_ny = config['crm_ny']
    nprocx = config['nprocx']
    nprocy = config['nprocy']

DSarr = []
varr = []
warr = []
densarr = []
QXZlarr = []
denslarr = []
hsarr = []
coriolisxzarr = []

for i in range(nprocs):
    DS = xr.open_dataset(ncbase + str(i) + '.nc')
    DS.load()
    ndensity = DS.dims['dens_ndofs']
    nt = DS.dims['t']
    nens = DS.dims['nens']
    if (i==0):
        mass = DS.mass
        energy = DS.energy
        pens = DS.pens
        pv = DS.pv
        densmax = DS.densmax
        densmin = DS.densmin
    varr.append(DS.v)
    warr.append(DS.w)
    densarr.append(DS.dens)
    QXZlarr.append(DS.QXZl)
    denslarr.append(DS.densl)
    hsarr.append(DS.hs)
    coriolisxzarr.append(DS.coriolisxz)

if model == 'swe':
    dens_names = ['h',]
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

#THIS IS A LITLTE BROKEN FOR RHOD VARIANTS...
#probably ok, this is just a quick and dirty plotting script...

#SUPER WRONG NOW FOR NAMES
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

def build_data(

Nlist = np.arange(0,nt)



for n in range(nens):

    plotvar_scalar2D('hs.'+ str(n), hs.isel(hs_ndofs=0, dual_ncells_y=0,nens=n),0)
    plotvar_scalar2D('coriolisxz.'+ str(n), coriolisxz.isel(coriolisxz_ndofs=0, primal_ncells_y=0,nens=n),0)

    for i in Nlist:
        plotvar_scalar2D('qxz.'+ str(n), QXZl.isel(t=i,QXZl_ndofs=0, dual_ncells_y=0,nens=n),i)

        plotvar_scalar2D('v.'+ str(n), v.isel(t=i,v_ndofs=0, primal_ncells_y=0,nens=n),i)
        plotvar_scalar2D('w.'+ str(n), w.isel(t=i,w_ndofs=0, primal_ncells_y=0,nens=n),i)


        for l,name in zip(range(ndensity), dens_names):
            plotvar_scalar2D(name + '.' + str(n), dens.isel(t=i,dens_ndofs=l, dual_ncells_y=0,nens=n),i)
            plotvar_scalar2D(name+'l.'+ str(n), densl.isel(t=i,densl_ndofs=l, primal_ncells_y=0,nens=n),i)
#THIS ASSUMES TOTAL DENSITY IS IN DENS(0)
            plotvar_scalar2D(name+'c.'+ str(n), dens.isel(t=i,dens_ndofs=l, dual_ncells_y=0,nens=n) / dens.isel(t=i,dens_ndofs=0, dual_ncells_y=0,nens=n),i)
    #if model in ['tswe','ce','mce']:
    #        plotvar_scalar2D('thetal', dens.isel(t=i,dens_ndofs=1, dual_ncells_y=0,nens=n) / dens.isel(t=i,dens_ndofs=0, dual_ncells_y=0,nens=n),i)




