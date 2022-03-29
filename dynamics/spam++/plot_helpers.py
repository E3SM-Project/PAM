import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

def plot_stat(statname, data):
    plt.figure(figsize=(10,8))
    plt.plot( (data - data[0])/data[0]*100. )
    plt.xlabel('Nsteps')
    plt.ylabel('Fractional Change in ' + statname)
    plt.tight_layout()
    plt.savefig(statname + '.png')

def plot_rawstat(statname, data):
    plt.figure(figsize=(10,8))
    plt.plot(data)
    plt.xlabel('Nsteps')
    plt.ylabel(statname)
    plt.tight_layout()
    plt.savefig(statname + 'raw.png')

def plotvar_scalar1D(plotname,vardat,i):
    plt.figure(figsize=(10,8))
    plt.plot(vardat)
    plt.xlabel('x')
    plt.ylabel(plotname)
    plt.tight_layout()
    plt.savefig(plotname + '.' + str(i) + '.png')
    plt.close('all')

def plotvar_scalar2D(plotname,vardat,i):
    plt.figure(figsize=(10,8))
    plt.contourf(vardat)
    plt.colorbar()
    plt.contour(vardat)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.tight_layout()
    plt.savefig(plotname + '.' + str(i) + '.png')
    plt.close('all')


def plotvar_vector2D(plotname,vardat0, vardat1,i):
    plt.figure(figsize=(10,8))
    plt.contourf(vardat0)
    plt.colorbar()
    plt.contour(vardat0)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.tight_layout()
    plt.savefig(plotname + '0.' + str(i) + '.png')
    plt.close('all')

    plt.figure(figsize=(10,8))
    plt.contourf(vardat1)
    plt.colorbar()
    plt.contour(vardat1)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.tight_layout()
    plt.savefig(plotname + '1.' + str(i) + '.png')
    plt.close('all')
