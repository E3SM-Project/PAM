import xarray as xr
from plot_helpers import plotvar_scalar2D, plotvar_vector2D, plot_stat, plot_rawstat
import numpy as np
import sys

DS = xr.open_dataset('output.nc')
DS.load()

# SHOULD AUTO LOAD FROM input-layermodel2D.txt
Nlist = np.arange(21)

# ADJUST THIS LOGIC A LITTLE FOR SWE/TSWE VS. CE
model = sys.argv[1]
ntracers = int(sys.argv[2])
ntracers_fct = int(sys.argv[3])
if (model == 'swe'): 
    ndensity = 1 + ntracers
if (model == 'tswe'): 
    ndensity = 2 + ntracers
if (model == 'ce'): 
    ndensity = 2
if (model == 'mcerho'): 
    ndensity = 2
if (model == 'mcerhod'): 
    ndensity = 2

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

#maybe mass/tracer stuff gets a list of names of size ndensity/ndensityfct?
#YES DO IT LIKE THIS!

mass_names = []
massfct_names = []
dens_names = []
densfct_names = []

plot_stat('total_mass', mass.isel(mass_ndofs=0))
if (model == 'tswe'):
    plot_stat('total_bouyancy', mass.isel(mass_ndofs=1))
if model in ['ce','mcerho','mcerhod']:    
    plot_stat('total_entropic_density', mass.isel(mass_ndofs=1))

plot_rawstat('internal_energy', energy.isel(energy_ndofs=3))
plot_rawstat('potential_energy', energy.isel(energy_ndofs=2))
plot_rawstat('kinetic_energy', energy.isel(energy_ndofs=1))
plot_stat('total_energy', energy.isel(energy_ndofs=0))

plot_stat('total_pens', pens.isel(pens_ndofs=0))
plot_stat('total_pv', pv.isel(pv_ndofs=0))

# ADJUST THIS LOGIC A LITTLE FOR SWE/TSWE VS. CE
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
# actually, no, this is just a quick output plotting script
# add more detailed test case specific plotting stuff for a paper later on...

v = DS.v
dens = DS.dens
q = DS.q
densl = DS.densl
hs = DS.hs
coriolis = DS.coriolis

# ADJUST THIS LOGIC A LITTLE FOR SWE/TSWE VS. CE
if (ntracers_fct>0):
    densfct = DS.densfct
    densfctl = DS.densfctl

plotvar_scalar2D('hs', hs.isel(hs_ndofs=0,dual_ncells_z=0),0)
plotvar_scalar2D('coriolis', coriolis.isel(coriolis_ndofs=0,primal_ncells_z=0),0)

# ADJUST THIS LOGIC A LITTLE FOR SWE/TSWE VS. CE
# WHERE EXACTLY SHOULD sl/trl/trfctl live? are they straight 0-forms? twisted n-forms?
for i in Nlist:
    plotvar_scalar2D('q', q.isel(t=i,q_ndofs=0,dual_ncells_z=0),i)
    plotvar_scalar2D('h', dens.isel(t=i,dens_ndofs=0,dual_ncells_z=0),i)
    if model in ['tswe','ce','mcerho','mcerhod']:    
        plotvar_scalar2D('S', dens.isel(t=i,dens_ndofs=1,dual_ncells_z=0),i)
        plotvar_scalar2D('sl', densl.isel(t=i,densl_ndofs=1,primal_ncells_z=0),i)
    for l in range(ndensity-ntracers,ndensity):
        plotvar_scalar2D('tr' + str(l-ntracers), dens.isel(t=i,dens_ndofs=l,dual_ncells_z=0),i)
        plotvar_scalar2D('trl' + str(l-ntracers), densl.isel(t=i,densl_ndofs=l,primal_ncells_z=0),i)
    if ntracers_fct > 0:
        for l in range(ntracers_fct):
            plotvar_scalar2D('trfct' + str(l), densfct.isel(t=i,densfct_ndofs=l,dual_ncells_z=0),i)
            plotvar_scalar2D('trfctl' + str(l), densfctl.isel(t=i,densfctl_ndofs=l,primal_ncells_z=0),i)
    plotvar_vector2D('v', v.isel(t=i,v_ndofs=0,primal_ncells_z=0), v.isel(t=i,v_ndofs=1,primal_ncells_z=0),i)

