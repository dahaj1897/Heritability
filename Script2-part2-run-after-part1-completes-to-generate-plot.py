import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import numpy as np

levels1=np.arange(.85, 1.65,.01)
levels2=np.arange(.85, 1.65,.01)
edge_removal=0

fig = plt.figure()
data = np.genfromtxt("plot_me.csv", dtype=float, delimiter=',') 
x = data[:,0]
y = data[:,1]
z = data[:,2]
xi = np.linspace(min(x)-edge_removal, max(x)+edge_removal)
yi = np.linspace(min(y)-edge_removal, max(y)+edge_removal)
X, Y = np.meshgrid(xi, yi)
Z = griddata(x, y, z, xi, yi, interp='linear')
plt.contourf(X,Y,Z, levels1, inline=1, fontsize=10, extend='both')
plt.colorbar()
plt.contour(X, Y, Z, (0,), colors = 'g', linewidths = 2, hold='on')    
plt.xlabel('Developmental instability affecting just cells (Sigma)')
plt.ylabel("Stdv of the clusters as a fraction of the mean cluster size ")    
#plt.title('Relative heritability')
plt.savefig("Relative heritability.png", dpi=300) 
plt.savefig("Relative heritability.pdf") 


