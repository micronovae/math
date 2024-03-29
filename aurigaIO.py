from env import Auriga 
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl

basePATH = "/Simulations/Auriga/level4_DM/halo_6/output"
dm_pos = Auriga.snapshot.loadSubset(basePATH,127,'dm',['Coordinates'])
pos = dm_pos['Coordinates']

x = pos[:,0]
y = pos[:,1]
z = pos[:,2]

plt.figure(0)
plt.hist2d(x,y,norm=mpl.colors.LogNorm(),bins=64)
plt.savefig(fname='x-y.png')

plt.figure(1)
plt.hist2d(x,z,norm=mpl.colors.LogNorm(),bins=64)
plt.savefig(fname='x-z.png')

plt.figure(2)
plt.hist2d(y,z,norm=mpl.colors.LogNorm(),bins=64)
plt.savefig(fname='y-z.png')

print(dm_pos['Coordinates'])
print(type(dm_pos))
