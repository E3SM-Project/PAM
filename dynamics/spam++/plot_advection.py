import xarray as xr
from plot_helpers import plotvar_scalar2D, plotvar_vector2D, plot_stat, plot_rawstat
import numpy as np

DS = xr.open_dataset('output.nc')
DS.load()



qmass = DS.qmass
qmin = DS.qmin
qmax = DS.qmax

plot_stat('mass', qmass)
plot_rawstat('min', qmin)
plot_rawstat('max', qmax)


v = DS.v
q = DS.q
q0 = DS.q0

#Nlist = [0,1,2,3,4,5,6,7,8,9,10]
Nlist = [0,1,10,40,80,120,160,200]

plotvar_vector2D('v', v.isel(v_ndofs=0,ncells_z=0), v.isel(v_ndofs=1,ncells_z=0),0)
for i in Nlist:
    plotvar_scalar2D('q', q.isel(t=i,q_ndofs=0,ncells_z=0),i)
    plotvar_scalar2D('q0', q0.isel(t=i,q0_ndofs=0,ncells_z=0),i)

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
