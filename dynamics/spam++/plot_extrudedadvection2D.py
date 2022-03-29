import xarray as xr
from plot_helpers import plotvar_scalar2D, plotvar_vector2D, plot_stat, plot_rawstat
import numpy as np
import sys

DS = xr.open_dataset('output.nc')
DS.load()

ntdofs = DS.dims['dens_ndofs']
nQdofs = DS.dims['QXZ_ndofs']
nt = DS.dims['t']

mass = DS.mass
densmax = DS.densmax
densmin = DS.densmin
qmass = DS.qmass

dens_names = []
dens_stat_names = []
for k in range(ntdofs):
    dens_names.append('T'+str(k))
    dens_stat_names.append('tracer'+str(k))

for l,name in zip(range(ntdofs), dens_stat_names):
    plot_stat('total_'+name, mass.isel(mass_ndofs=l))
    plot_rawstat('min_' +name, densmin.isel(densmin_ndofs=l))
    plot_rawstat('max_' +name, densmax.isel(densmax_ndofs=l))

for l in range(nQdofs):
    plot_stat('qmass' + str(l), qmass.isel(qmass_ndofs=l))
        
Nlist = np.arange(0,nt)

v = DS.v
w = DS.w
U = DS.U
Uw = DS.Uw
vt = DS.vt
wt = DS.wt
dens = DS.dens
densl = DS.densl
QXZ = DS.QXZ
QXZl = DS.QXZl

qxzvertrecon = DS.qxzvertrecon
qxzrecon = DS.qxzrecon
qxzflux = DS.qxzflux 
qxzvertflux = DS.qxzvertflux 
D2 = DS.D2
D2Q = DS.D2Q
               
plotvar_scalar2D('v', v.isel(v_ndofs=0, primal_ncells_y=0),0)
plotvar_scalar2D('w', w.isel(w_ndofs=0, primal_ncells_y=0),0)
plotvar_scalar2D('U', U.isel(U_ndofs=0, dual_ncells_y=0),0)
plotvar_scalar2D('Uw', Uw.isel(Uw_ndofs=0, dual_ncells_y=0),0)
plotvar_scalar2D('vt', vt.isel(vt_ndofs=0, primal_ncells_y=0),0)
plotvar_scalar2D('wt', wt.isel(wt_ndofs=0, primal_ncells_y=0),0)

#plotvar_scalar2D('D2', D2.isel(D2_ndofs=0, primal_ncells_y=0),0)

for i in Nlist:
    for l,name in zip(range(ntdofs), dens_names):
            plotvar_scalar2D(name, dens.isel(t=i,dens_ndofs=l, dual_ncells_y=0),i)
            plotvar_scalar2D(name+'l', densl.isel(t=i,densl_ndofs=l, primal_ncells_y=0),i)
    for l in range(nQdofs):
        plotvar_scalar2D('QXZ' + str(l), QXZ.isel(t=i,QXZ_ndofs=l,primal_ncells_y=0),i)
        plotvar_scalar2D('QXZl' + str(l), QXZl.isel(t=i,QXZl_ndofs=l,dual_ncells_y=0),i)

    #for l in range(nQdofs):
    #    plotvar_scalar2D('qxzrecon' + str(l), qxzrecon.isel(t=i,qxzrecon_ndofs=l,primal_ncells_y=0),i)
    #    plotvar_scalar2D('qxzvertrecon' + str(l), qxzvertrecon.isel(t=i,qxzvertrecon_ndofs=l,primal_ncells_y=0),i)
    #    plotvar_scalar2D('qxzflux' + str(l), qxzflux.isel(t=i,qxzflux_ndofs=l,primal_ncells_y=0),i)
    #    plotvar_scalar2D('qxzvertflux' + str(l), qxzvertflux.isel(t=i,qxzvertflux_ndofs=l,primal_ncells_y=0),i)
    #    plotvar_scalar2D('D2Q' + str(l), D2Q.isel(t=i,D2Q_ndofs=l,primal_ncells_y=0),i)
