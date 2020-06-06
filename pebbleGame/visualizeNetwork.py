import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline

pebbleGameOutput = np.loadtxt('./pebbleGameOutput.dat',dtype = int)
coordinates = np.loadtxt('coordinates.dat',dtype = float)

N = len(coordinates) # number of sites
cutoff = 0.5

# Find the short bonds
points = pebbleGameOutput[:,0:2] - 1 # minus one is because the indexing in python starts from 0
d = coordinates[points]
lengths = np.linalg.norm(d[:,0]-d[:,1],axis=-1)
short_idx = lengths < cutoff
pebbleGameOutput = pebbleGameOutput[short_idx]

# Find stressed and non-stressed bonds
stressed = pebbleGameOutput[pebbleGameOutput[:,3] == 1]
non_stressed = pebbleGameOutput[pebbleGameOutput[:,3] == 0]

# Find hinge and non-hinge
clusters = {i:set() for i in range(1,N+1)}
for item in pebbleGameOutput:
    site1,site2,cluster,stress = item
    clusters[site1].add(cluster)
    clusters[site2].add(cluster)
hinges =np.array([item for item in clusters if len(clusters[item]) > 1])
non_hinges = np.setdiff1d(np.arange(N)+1,hinges)

# Now the plotting begins
fig,ax = plt.subplots(1,1,figsize = (15,15))
ax.set_axis_off()
# First plot the bonds
def bond_plot(coordinates,bondlist):
    x,y = [], []
    for item in bondlist:
        xs = coordinates[item][:,0]
        ys = coordinates[item][:,1]
        x.extend(xs)
        x.append(None)
        y.extend(ys)
        y.append(None)
    return x,y
x_stressed, y_stressed = bond_plot(coordinates,stressed[:,0:2]-1)
x_non_stressed, y_non_stressed = bond_plot(coordinates,non_stressed[:,0:2]-1)
ax.plot(x_stressed, y_stressed, color = 'k', lw = 4.5,zorder=1)
ax.plot(x_non_stressed, y_non_stressed, color = 'r', lw = 4.5,zorder=1)

# Then plot the sites, which include hinges and non-hinges
try:
    ax.scatter(coordinates[hinges-1,0],coordinates[hinges-1,1], color = 'g', s = 100,zorder=2)
except:
    pass
try:
    ax.scatter(coordinates[non_hinges-1,0],coordinates[non_hinges-1,1], color = 'k', s =100,zorder=2)
except:
    pass

for i in range(1,N+1):
   ax.text(coordinates[i-1,0],coordinates[i-1,1],i, fontsize = 14, color = 'r',zorder = 3)

plt.tight_layout()
plt.savefig('./network.pdf')
