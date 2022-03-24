import xarray as xr
from plot_helpers import plotvar_scalar2D, plotvar_vector2D, plot_stat, plot_rawstat
import numpy as np
import sys

DS = xr.open_dataset('output.nc')
DS.load()

ndensity = DS.dims['dens_ndofs']
nt = DS.dims['t']
model = sys.argv[1]

mass = DS.mass
energy = DS.energy
pens = DS.pens
pv = DS.pv

densmax = DS.densmax
densmin = DS.densmin




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
if model == 'mce':
    dens_names = ['rho','Theta','rho_v', 'rho_l', 'rho_i']    
    dens_stat_names = ['mass','entropic_var_density','vapor', 'liquid', 'ice']    
    nprogdens = 5
    
for k in range(ndensity-nprogdens):
    dens_names.append('T'+str(k))
    dens_stat_names.append('tracer'+str(k))

for l,name in zip(range(ndensity), dens_stat_names):
    plot_stat('total_'+name, mass.isel(mass_ndofs=l))
    plot_rawstat('min_' +name, densmin.isel(densmin_ndofs=l))
    plot_rawstat('max_' +name, densmax.isel(densmax_ndofs=l))
                
plot_rawstat('internal_energy', energy.isel(energy_ndofs=3))
plot_rawstat('potential_energy', energy.isel(energy_ndofs=2))
plot_rawstat('kinetic_energy', energy.isel(energy_ndofs=1))
plot_stat('total_energy', energy.isel(energy_ndofs=0))

plot_stat('total_pens', pens.isel(pens_ndofs=0))
plot_stat('total_pv', pv.isel(pv_ndofs=0))


Nlist = np.arange(0,nt)

v = DS.v
w = DS.w
dens = DS.dens
QXZl = DS.QXZl
densl = DS.densl
hs = DS.hs
coriolisxz = DS.coriolisxz

K = DS.K
F = DS.F
Fw = DS.Fw
he = DS.he
hew = DS.hew

plotvar_scalar2D('hs', hs.isel(hs_ndofs=0, dual_ncells_y=0),0)
plotvar_scalar2D('coriolisxz', coriolisxz.isel(coriolisxz_ndofs=0, primal_ncells_y=0),0)

for i in Nlist:
    plotvar_scalar2D('qxz', QXZl.isel(t=i,QXZl_ndofs=0, dual_ncells_y=0),i)

    plotvar_scalar2D('v', v.isel(t=i,v_ndofs=0, primal_ncells_y=0),i)
    plotvar_scalar2D('w', w.isel(t=i,w_ndofs=0, primal_ncells_y=0),i)

    #plotvar_scalar2D('F', F.isel(t=i,F_ndofs=0, dual_ncells_y=0),i)
    #plotvar_scalar2D('FW', Fw.isel(t=i,Fw_ndofs=0, dual_ncells_y=0),i)
    #plotvar_scalar2D('he', he.isel(t=i,he_ndofs=0, dual_ncells_y=0),i)
    #plotvar_scalar2D('hew', hew.isel(t=i,hew_ndofs=0, dual_ncells_y=0),i)
    #plotvar_scalar2D('K', K.isel(t=i,K_ndofs=0, dual_ncells_y=0),i)

#    plotvar_vector2D('v', v.isel(t=i,v_ndofs=0), w.isel(t=i,w_ndofs=0),i)
    for l,name in zip(range(ndensity), dens_names):
            plotvar_scalar2D(name, dens.isel(t=i,dens_ndofs=l, dual_ncells_y=0),i)
            plotvar_scalar2D(name+'l', densl.isel(t=i,densl_ndofs=l, primal_ncells_y=0),i)
#THIS ASSUMES TOTAL DENSITY IS IN DENS(0)
            plotvar_scalar2D(name+'c', dens.isel(t=i,dens_ndofs=l, dual_ncells_y=0) / dens.isel(t=i,dens_ndofs=0, dual_ncells_y=0),i)
    #if model in ['tswe','ce','mce']:
    #        plotvar_scalar2D('thetal', dens.isel(t=i,dens_ndofs=1, dual_ncells_y=0) / dens.isel(t=i,dens_ndofs=0, dual_ncells_y=0),i)




