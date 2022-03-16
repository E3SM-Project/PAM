import xarray as xr
from plot_helpers import plotvar_scalar2D, plotvar_vector2D, plot_stat, plot_rawstat
import numpy as np
import sys

DS = xr.open_dataset('output.nc')
DS.load()

#Nlist = [0,1,2,3,4,5,6,7,8,9,10]
#Nlist = [0,2,4,6,8,10,20,40,60,80,100,120,140,160,180,200]
#Nlist = [0,20,40,60,80,100,120,140,160,180,200]
Nlist = [0,40,80,120,160,200]
ntdofs = int(sys.argv[1])
ntfctdofs = int(sys.argv[2])
nQdofs = int(sys.argv[3])

Tmass = DS.Tmass
Tmin = DS.Tmin
Tmax = DS.Tmax
Tfctmass = DS.Tfctmass
Tfctmin = DS.Tfctmin
Tfctmax = DS.Tfctmax
Qmass = DS.Qmass

for l in range(ntdofs):
    plot_stat('Tmass' + str(l), Tmass.isel(Tmass_ndofs=l))
    plot_rawstat('Tmin' + str(l), Tmin.isel(Tmin_ndofs=l))
    plot_rawstat('Tmax' + str(l), Tmax.isel(Tmax_ndofs=l))
for l in range(ntfctdofs):
    plot_stat('Tfctmass' + str(l), Tfctmass.isel(Tfctmass_ndofs=l))
    plot_rawstat('Tfctmin' + str(l), Tfctmin.isel(Tfctmin_ndofs=l))
    plot_rawstat('Tfctmax' + str(l), Tfctmax.isel(Tfctmax_ndofs=l))
for l in range(nQdofs):
    plot_stat('Qmass' + str(l), Qmass.isel(Qmass_ndofs=l))

v = DS.v
U = DS.U
UT = DS.UT
plotvar_vector2D('v', v.isel(v_ndofs=0,primal_ncells_z=0), v.isel(v_ndofs=1,primal_ncells_z=0),0)
plotvar_vector2D('U', U.isel(U_ndofs=0,dual_ncells_z=0), U.isel(U_ndofs=1,dual_ncells_z=0),0)
plotvar_vector2D('UT', UT.isel(UT_ndofs=0,primal_ncells_z=0), UT.isel(UT_ndofs=1,primal_ncells_z=0),0)

plotvar_scalar2D('UT00', UT.isel(UT_ndofs=0,primal_ncells_z=0),0)
plotvar_scalar2D('UT01', UT.isel(UT_ndofs=1,primal_ncells_z=0),0)

print(UT.isel(UT_ndofs=0,primal_ncells_z=0))
print(np.max(UT.isel(UT_ndofs=0,primal_ncells_z=0)), np.min(UT.isel(UT_ndofs=0,primal_ncells_z=0)))

print(UT.isel(UT_ndofs=1,primal_ncells_z=0))
print(np.max(UT.isel(UT_ndofs=1,primal_ncells_z=0)), np.min(UT.isel(UT_ndofs=1,primal_ncells_z=0)))

T = DS.T
T0 = DS.T0
Tfct = DS.Tfct
Tfct0 = DS.Tfct0
Q = DS.Q
Q0 = DS.Q0

for i in Nlist:

    for l in range(ntdofs):
        plotvar_scalar2D('T' + str(l), T.isel(t=i,T_ndofs=l,dual_ncells_z=0),i)
        plotvar_scalar2D('Tl' + str(l), T0.isel(t=i,T0_ndofs=l,primal_ncells_z=0),i)

    for l in range(ntfctdofs):
        plotvar_scalar2D('Tfct' + str(l), Tfct.isel(t=i,Tfct_ndofs=l,dual_ncells_z=0),i)
        plotvar_scalar2D('Tfctl' + str(l), Tfct0.isel(t=i,Tfct0_ndofs=l,primal_ncells_z=0),i)

    for l in range(nQdofs):
        plotvar_scalar2D('Q' + str(l), Q.isel(t=i,Q_ndofs=l,primal_ncells_z=0),i)
        plotvar_scalar2D('Ql' + str(l), Q0.isel(t=i,Q0_ndofs=l,dual_ncells_z=0),i)

# def update(i):
#      ax.clear()
#      ax.set_xlabel('x')
#      ax.set_ylabel('y')
#      ax.contourf(q.isel(t=i%40,q_ndofs=0,ncells_z=0))
#      ax.contour(q.isel(t=i%40,q_ndofs=0,ncells_z=0))

# Construct the animation, using the update function as the animation director.
#animation = FuncAnimation(fig, update, interval=80, blit=False)
#animation.save('test.mp4')

#plt.show()
