import xarray as xr
from plot_helpers import plotvar_scalar2D, plotvar_vector2D, plot_stat, plot_rawstat
import numpy as np
import sys

DS = xr.open_dataset('test_dycore0.nc')
DS.load()

ndensity = DS.dims['dens_ndofs']
nt = DS.dims['t']
model = sys.argv[1]
nens = DS.dims['nens']

plane = 'xz'
if len(sys.argv) > 2:
    plane = sys.argv[2]
print(f"plane = {plane}")

mass = DS.mass
energy = DS.energy
pens = DS.pens
pv = DS.pv

densmax = DS.densmax
densmin = DS.densmin

plane_idx = DS.dims['dual_ncells_y'] // 2
ndims = DS.dims['v_ndofs']
print(f"plane_idx = {plane_idx}")
if len(sys.argv) > 3:
    plane_idx = int(sys.argv[3])


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

if model == 'an':
    dens_names = ['Theta',]
    dens_stat_names = ['entropic_var_density',]
    nprogdens = 1

#THIS IS A LITLTE BROKEN FOR RHOD VARIANTS...
#probably ok, this is just a quick and dirty plotting script...
if model == 'mce':
    dens_names = ['rho','Theta','rho_v', 'rho_l', 'rho_i']    
    dens_stat_names = ['mass','entropic_var_density','vapor', 'liquid', 'ice']    
    nprogdens = 5

if model == 'man':
    dens_names = ['Theta','rho_v', 'rho_l', 'rho_i']    
    dens_stat_names = ['entropic_var_density','vapor', 'liquid', 'ice']    
    nprogdens = 4
    
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
w = DS.w
dens = DS.dens
total_dens = DS.total_dens
QXZl = DS.QXZl
densl = DS.densl
hs = DS.hs
coriolisxz = DS.coriolisxz

for n in range(nens):

    if plane == 'xz':
        plotvar_scalar2D('.'.join(['hs', plane, str(n)]), hs.isel(hs_ndofs=0, dual_ncells_y=plane_idx,nens=n),0)
        plotvar_scalar2D('.'.join(['coriolisxz', plane, str(n)]), coriolisxz.isel(coriolisxz_ndofs=0, primal_ncells_y=plane_idx,nens=n),0)
    if plane == 'yz':
        plotvar_scalar2D('.'.join(['hs', plane, str(n)]), hs.isel(hs_ndofs=0, dual_ncells_x=plane_idx,nens=n),0)
        plotvar_scalar2D('.'.join(['coriolisxz', plane, str(n)]), coriolisxz.isel(coriolisxz_ndofs=0, primal_ncells_x=plane_idx,nens=n),0)

    for i in Nlist:
        if plane == 'xz':
            plotvar_scalar2D('.'.join(['qxz', plane, str(n)]), QXZl.isel(t=i,QXZl_ndofs=0, dual_ncells_y=plane_idx,nens=n),i)
            #if ndims > 1:
            #    plotvar_scalar2D('.'.join(['qyz', plane, str(n)]), QXZl.isel(t=i,QXZl_ndofs=1, dual_ncells_y=plane_idx,nens=n),i)
            plotvar_scalar2D('.'.join(['vx', plane, str(n)]), v.isel(t=i,v_ndofs=0, primal_ncells_y=plane_idx,nens=n),i)
            if ndims > 1:
                plotvar_scalar2D('.'.join(['vy', plane, str(n)]), v.isel(t=i,v_ndofs=1, primal_ncells_y=plane_idx,nens=n),i)
            plotvar_scalar2D('.'.join(['w', plane, str(n)]), w.isel(t=i,w_ndofs=0, primal_ncells_y=plane_idx,nens=n),i)
        if plane == 'yz':
            plotvar_scalar2D('.'.join(['qxz', plane, str(n)]), QXZl.isel(t=i,QXZl_ndofs=0, dual_ncells_x=plane_idx,nens=n),i)
            #if ndims > 1:
            #    plotvar_scalar2D('.'.join(['qyz', plane, str(n)]), QXZl.isel(t=i,QXZl_ndofs=1, dual_ncells_x=plane_idx,nens=n),i)
            plotvar_scalar2D('.'.join(['vx', plane, str(n)]), v.isel(t=i,v_ndofs=0, primal_ncells_x=plane_idx,nens=n),i)
            if ndims > 1:
                plotvar_scalar2D('.'.join(['vy', plane, str(n)]), v.isel(t=i,v_ndofs=1, primal_ncells_x=plane_idx,nens=n),i)
            plotvar_scalar2D('.'.join(['w', plane, str(n)]), w.isel(t=i,w_ndofs=0, primal_ncells_x=plane_idx,nens=n),i)


        for l,name in zip(range(ndensity), dens_names):
            if plane == 'xz':
                plotvar_scalar2D('.'.join([name, plane, str(n)]), dens.isel(t=i,dens_ndofs=l, dual_ncells_y=plane_idx,nens=n),i)
                plotvar_scalar2D('.'.join([name, 'l', plane, str(n)]), densl.isel(t=i,densl_ndofs=l, primal_ncells_y=plane_idx,nens=n),i)
                plotvar_scalar2D('.'.join([name, 'c', plane, str(n)]), dens.isel(t=i,dens_ndofs=l, dual_ncells_y=plane_idx,nens=n) / total_dens.isel(t=i,total_dens_ndofs=0, dual_ncells_y=plane_idx,nens=n),i)

            if plane == 'yz':
                plotvar_scalar2D('.'.join([name, plane, str(n)]), dens.isel(t=i,dens_ndofs=l, dual_ncells_x=plane_idx,nens=n),i)
                plotvar_scalar2D('.'.join([name, 'l', plane, str(n)]), densl.isel(t=i,densl_ndofs=l, primal_ncells_x=plane_idx,nens=n),i)
                plotvar_scalar2D('.'.join([name, 'c', plane, str(n)]), dens.isel(t=i,dens_ndofs=l, dual_ncells_x=plane_idx,nens=n) / total_dens.isel(t=i,total_dens_ndofs=0, dual_ncells_x=plane_idx,nens=n),i)




