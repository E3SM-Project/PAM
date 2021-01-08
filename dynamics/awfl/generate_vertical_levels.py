
import numpy as np
import sys, getopt
import matplotlib.pyplot as plt
from netCDF4 import Dataset

avail = ["help","function=","exp-base=","tanh-inflect=","tanh-steep=","tanh-scale=","z0=","ztop=","nlev=","no-plots","output="]

function="tanh"
tanh_inflect=0.4
tanh_steep=5.
tanh_scale=10.
exp_base=10.
z0=0.
ztop=10000.
nlev=128
showplots=True
output="vcoords.nc"

def print_help() :
  print("Usage: python generate_vertical_levels.py")
  print("Options:")
  print("--function=[exp|tanh]")
  print("  equal: self-explanatory")
  print("  exp:")
  print("    * This gives you a consistent coarsening as you increase in height")
  print("    --exp-base=10")
  print("      * Top cell grid spacing will be 10x larger than bottom cell")
  print("  tanh:")
  print("    * This allows you to concentrate layers at the surface, then keep grid spacing fairly constant after that")
  print("    --tanh-inflect=0.4")
  print("      * tanh inflection point is 40% through the vertical indices")
  print("    --tanh-steep=5")
  print("      * How quickly grid spacing changes around tanh inflection point")
  print("        (generally between 2 and 10)")
  print("    --tanh-scale=10")
  print("      * The maximum bound on the grid spacing ratio of top cell to bottom cell")
  print("--help")
  print("  * Print this message")
  print("--z0=[Bottom of the vertical grid]")
  print("--ztop=[Top of the vertical grid]")
  print("--nlev=[Number of cells]")
  print("--no-plots")
  print("  * Don't show matplotlib plots, e.g., over a remote server w/ a laggy connection")
  print("--output=[vcoords.nc]")


try:
  opts , args = getopt.getopt(sys.argv[1:] , "h" , avail)
except:
  print("ERROR: Invalid argument")
  print_help()
  sys.exit(-1)
for opt,arg in opts:
  if (opt == "-h" or opt == "--help") :
    print_help()
    sys.exit()
  elif opt == "--function" :
    function = arg
  elif opt == "--exp-base" :
    exp_base = float(arg)
  elif opt == "--tanh-inflect" :
    tanh_inflect = float(arg)
  elif opt == "--tanh-steep" :
    tanh_steep = float(arg)
  elif opt == "--tanh-scale" :
    tanh_scale = float(arg)
  elif opt == "--z0" :
    z0 = float(arg)
  elif opt == "--ztop" :
    ztop = float(arg)
  elif opt == "--nlev" :
    nlev = int(arg)
  elif opt == "--no-plots" :
    showplots = False
  elif opt == "--output" :
    output = arg
  else :
    print("ERROR: Invalid argument")
    # print_help()
    sys.exit(-1)
    

print("nlev: "+str(nlev))
print("z0: "+str(z0))
print("ztop: "+str(ztop))
print("function: "+function)
if function == "exp" :
  print("  exp_base: "+str(exp_base))
if function == "tanh" :
  print("  tanh_inflect: "+str(tanh_inflect))
  print("  tanh_steep: "+str(tanh_steep))
  print("  tanh_scale: "+str(tanh_scale))

zthick = ztop - z0

dz   = np.array([0. for i in range(nlev)])
zint = np.array([0. for i in range(nlev+1)])
zmid = np.array([0. for i in range(nlev)])

zint[0] = 0;
  
for i in range(nlev) :
  zloc = float(i)/float(nlev-1)
  if (function == "exp") :
    dz[i] = pow(exp_base,zloc)
  elif (function == "tanh") :
    dz[i] = ((np.tanh((zloc-tanh_inflect)*tanh_steep)+1)/2*(tanh_scale-1))+1
  elif (function == "equal") :
    dz[i] = 1
  zint[i+1] = zint[i] + dz[i]

scale = zthick / zint[nlev]
for i in range(nlev+1) :
  zint[i] = zint[i] * scale + z0

for i in range(nlev) :
  dz[i] = zint[i+1] - zint[i]
  zmid[i] = ( zint[i+1] + zint[i] ) / 2

if showplots :
  plt.plot(range(nlev),dz)
  plt.xlabel("vertical level index")
  plt.ylabel("grid spacing")
  plt.show()

  plt.close()

  skip = int( nlev / 50 )
  plt.hlines(zint[::skip],0,1,linewidth=0.2)
  plt.ylabel("height")
  plt.show()

print("\n**************** grid spacing ****************")
print(dz)
print("\n**************** interface locations ****************")
print(zint)


nc = Dataset(output,"w")
nc.createDimension("num_interfaces",nlev+1)
nc_zint = nc.createVariable("vertical_interfaces","f8",("num_interfaces",))
nc_zint[:] = zint
nc.close()


