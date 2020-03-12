import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

DS = xr.open_dataset('output.nc')
DS.load()

u = DS.u
q = DS.q


mass = DS.mass
min = DS.min
max = DS.max
print(mass)
print(min)
print(max)

Nlist = [0,1,10,40]

plt.figure(figsize=(10,8))
plt.quiver(u.isel(u_ndofs=0,ncells_z=0), u.isel(u_ndofs=1,ncells_z=0))
# FIX THIS TO BE BASED ON WIND MAGNITUDE
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('u.png')

for i in Nlist:
    plt.figure(figsize=(10,8))
    plt.contourf(q.isel(t=i,q_ndofs=0,ncells_z=0))
    plt.colorbar()
    plt.contour(q.isel(t=i,q_ndofs=0,ncells_z=0))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig('q' + str(i) + '.png')



fig = plt.figure(figsize=(10,8))
ax = fig.add_axes([0.05, 0.05, 0.8, 0.8])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.contourf(q.isel(t=0,q_ndofs=0,ncells_z=0))
#fig.colorbar()
ax.contour(q.isel(t=0,q_ndofs=0,ncells_z=0))
plt.savefig('test.png')

def update(i):
     ax.clear()
     ax.set_xlabel('x')
     ax.set_ylabel('y')
     ax.contourf(q.isel(t=i%40,q_ndofs=0,ncells_z=0))
     ax.contour(q.isel(t=i%40,q_ndofs=0,ncells_z=0))

# Construct the animation, using the update function as the animation director.
animation = FuncAnimation(fig, update, interval=80, blit=False)

animation.save('test.mp4')

plt.show()
