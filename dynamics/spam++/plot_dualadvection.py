import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

DS = xr.open_dataset('output.nc')
DS.load()

u = DS.u
q = DS.q


qmass = DS.qmass
qmin = DS.qmin
qmax = DS.qmax

def plot_stat(statname, data):
    plt.figure(figsize=(10,8))
    plt.plot( (data - data[0])/data[0]*100. )
    plt.xlabel('Nsteps')
    plt.ylabel('Fractional Change in ' + statname)
    plt.savefig(statname + '.png')

def plot_rawstat(statname, data):
    plt.figure(figsize=(10,8))
    plt.plot(data)
    plt.xlabel('Nsteps')
    plt.ylabel(statname)
    plt.savefig(statname + 'raw.png')

plot_stat('mass', qmass)
plot_rawstat('min', qmin)
plot_rawstat('max', qmax)

Nlist = [0,1,10,40]

plt.figure(figsize=(10,8))
plt.quiver(u.isel(u_ndofs=0,ncells_z=0), u.isel(u_ndofs=1,ncells_z=0))
# FIX THIS TO BE BASED ON WIND MAGNITUDE
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('vvec.png')


plt.figure(figsize=(10,8))
plt.contourf(u.isel(u_ndofs=0,ncells_z=0))
plt.colorbar()
plt.contour(u.isel(u_ndofs=0,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('u.png')

plt.figure(figsize=(10,8))
plt.contourf(u.isel(u_ndofs=1,ncells_z=0))
plt.colorbar()
plt.contour(u.isel(u_ndofs=1,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('v.png')

for i in Nlist:
    plt.figure(figsize=(10,8))
    plt.contourf(q.isel(t=i,q_ndofs=0,ncells_z=0))
    plt.colorbar()
    plt.contour(q.isel(t=i,q_ndofs=0,ncells_z=0))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig('q' + str(i) + '.png')




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
