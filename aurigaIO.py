from env import Auriga 
import sys
import matplotlib as plt

basePATH = "/Simulations/Auriga/level4_DM/halo_6/output"
dm_pos = Auriga.snapshot.loadSubset(basePATH,127,'dm',['Coordinates'])
x = dm_pos[:,0]
y = dm_pos[:,1]

plt.plot
print(type(dm_pos))
