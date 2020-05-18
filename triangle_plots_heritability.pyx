#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 09:43:50 2020

@author: seyed
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ternary
from scipy import stats
#%%
#Importing bootsrapped data
ar_raw = np.genfromtxt("/home/seyed/GT/Will/Matt/experiment/raw_aspect_ratio.csv", dtype=float, delimiter=',') 

fs_raw = np.genfromtxt("/home/seyed/GT/Will/Matt/experiment/raw_fracture_size.csv", dtype=float, delimiter=',') 

"""
names=['Clb2','Akr1','Arp8','Ace2']
ar=np.zeros(4)
ar_sd=np.zeros(4)
fs=np.zeros(4)
fs_sd=np.zeros(4)
for i in range(len(ar_raw[0,:])):
    ar[i]=np.mean(ar_raw[:,i])
    ar_sd[i]=np.std(ar_raw[:,i])
    fs[i]=np.mean(fs_raw[:,i])
    fs_sd[i]=np.std(fs_raw[:,i])
"""

#%%
#Function to create calulate heritabilities and create the data for ternary plots 
import ternary

def heat_map_gen(scale=100):
    #importing files
    ar = np.genfromtxt("/home/seyed/GT/Will/Matt/experiment/aspect_ratio.csv", dtype=float, delimiter=',') 
    fs = np.genfromtxt("/home/seyed/GT/Will/Matt/experiment/fracture_size.csv", dtype=float, delimiter=',')
    from ternary.helpers import simplex_iterator
    d=dict()
    dcel=dict()
    dclu=dict()
    #this loops creates all possible combination of three genotypes 
    for (i,j,k) in simplex_iterator(100):
        if i > 3 and j > 3 and k > 3:
            x=1000 #number of uni cells divided by 100 
            y=10 #number of clusters divided by 100
            """
            genotype codes based on the input file:
            0 ~ Clb2, 1 ~ Akr1, 2 ~ Arp8 3 ~ Y55
            """
            n=1 #first genotype
            m=2 #second genotype
            p=3 #third genotype
            #unicells
            ar1=np.random.choice(ar[:,n],i*x) 
            ar2=np.random.choice(ar[:,m],j*x)
            ar3=np.random.choice(ar[:,p],k*x)
            #clusters
            fs1=np.random.choice(fs[:,n],i*y)
            fs2=np.random.choice(fs[:,m],j*y)
            fs3=np.random.choice(fs[:,p],k*y)
            
            """
            Calculating the heritability of the unicells (See the math appendix for detail)
            """
            arave=np.mean(np.concatenate((ar1,ar2,ar3)))
            
            ssf=len(ar1)*np.square(np.mean(ar1)-arave)+len(ar2)*np.square(np.mean(ar2)-arave)+len(ar3)*np.square(np.mean(ar3)-arave)
            sss=np.sum(np.square(ar1-np.mean(ar1)))+np.sum(np.square(ar2-np.mean(ar2)))+np.sum(np.square(ar3-np.mean(ar3)))
            
            numcell=np.sum([len(ar1),len(ar2),len(ar3)])
            MSf=ssf/2
            MSw=sss/(numcell-3)

            varf=(MSf-MSw)*2/(numcell-np.sum([np.square(len(ar1)),np.square(len(ar2)),np.square(len(ar3))])/numcell)
            varsa=MSw
            h2ce=varf/(varf+varsa) #Cell level heritability
            """
            Calculating the heritability of the Clusters (See the math appendix for detail)
            """
            numclus=np.sum([len(fs1),len(fs2),len(fs3)])
            fsave=np.mean(np.concatenate((fs1,fs2,fs3)))
            SSf=len(fs1)*np.square(np.mean(fs1)-fsave)+len(fs2)*np.square(np.mean(fs2)-fsave)+len(fs3)*np.square(np.mean(fs3)-fsave)
            SSs=np.sum(np.square(fs1-np.mean(fs1)))+np.sum(np.square(fs2-np.mean(fs2)))+np.sum(np.square(fs3-np.mean(fs3)))
            
            msf=SSf/2
            msw=SSs/(numclus-3)
            VARf=(msf-msw)*2/(numclus-np.sum([np.square(len(fs1)),np.square(len(fs2)),np.square(len(fs3))])/numclus)
            VARsa=msw
            h2cl=VARf/(VARf+VARsa) #Cluster level heritability
                
            r=h2cl/h2ce #Ratio of heritabilities
            d[(i,j)]=r 
            dclu[(i,j)]=h2cl   
            dcel[(i,j)]=h2ce
            
    return [d, dcel, dclu, n, m, p]

#%%
scale=100
heatmap=heat_map_gen(scale)

#%%
"""
Plot the ration of heritabilities
"""
names=[r"$\Delta  Clb2$, $\Delta Ace2$",r"$\Delta Arp8$, $\Delta Ace2$",r"$\Delta  Akr1$, $\Delta Ace2$",r"$\Delta Ace2$"]
name=['c2','a8','a1','a2']
title=name[heatmap[3]]+name[heatmap[4]]+name[heatmap[5]]
print(title)
fontsize = 11
offset = 0.15
scale=100

import matplotlib.colors as colors
d=heatmap[0]
cmap = colors.ListedColormap(['black','red','orange', 'gold','c', 'lime', 'blue','navy', 'purple'])
figure, tax = ternary.figure(scale=scale)
tax.heatmap(d, style="triangular", cmap=cmap, colorbar=False, vmin=0.5, vmax=5)
#tax.set_title("Ratio of heritabilities")

tax.bottom_axis_label(names[heatmap[3]], fontsize=fontsize, offset=offset)

tax.right_axis_label(names[heatmap[4]], fontsize=fontsize, offset=offset)

tax.left_axis_label(names[heatmap[5]], fontsize=fontsize, offset=offset)

tax.ticks(axis='lbr', linewidth=2, multiple=25, offset=0.02)
tax.clear_matplotlib_ticks()
tax.get_axes().axis('off')
tax.boundary()
tax.gridlines(color="black", multiple=25)
tax.show()

plt.savefig("/home/seyed/GT/Will/Matt/triangle_plots/ratio_%s.svg"%(title))
plt.savefig("/home/seyed/GT/Will/Matt/triangle_plots/ratio_%s.png"%(title), dpi=300)

#%%
"""
Plot cell level heritability
"""
dcel=heatmap[1]
figure, tax = ternary.figure(scale=scale)

tax.heatmap(dcel, style="triangular", colorbar=False, vmin=0.2, vmax=1)
#tax.set_title("Cell level heritability")

tax.bottom_axis_label(names[heatmap[3]], fontsize=fontsize, offset=offset)

tax.right_axis_label(names[heatmap[4]], fontsize=fontsize, offset=offset)

tax.left_axis_label(names[heatmap[5]], fontsize=fontsize, offset=offset)

tax.ticks(axis='lbr', linewidth=2, multiple=25, offset=0.02)
tax.clear_matplotlib_ticks()
tax.get_axes().axis('off')
tax.boundary()
tax.gridlines(color="black", multiple=25)
tax.show()
plt.savefig("/home/seyed/GT/Will/Matt/triangle_plots/uni_%s.svg"%(title))
plt.savefig("/home/seyed/GT/Will/Matt/triangle_plots/uni_%s.png"%(title), dpi=300)
#%%
"""
Plot cluster level heritability
"""
dclu=heatmap[2]
figure, tax = ternary.figure(scale=scale)
tax.heatmap(dclu, style="triangular", colorbar=False, vmin=0.2, vmax=1)
#tax.set_title("Collective level heritability")

tax.bottom_axis_label(names[heatmap[3]], fontsize=fontsize, offset=offset)

tax.right_axis_label(names[heatmap[4]], fontsize=fontsize, offset=offset)

tax.left_axis_label(names[heatmap[5]], fontsize=fontsize, offset=offset)

tax.ticks(axis='lbr', linewidth=2, multiple=25, offset=0.02)

tax.clear_matplotlib_ticks()
tax.get_axes().axis('off')
tax.boundary()
tax.gridlines(color="black", multiple=25)

tax.show()
plt.savefig("/home/seyed/GT/Will/Matt/triangle_plots/col_%s.svg"%(title))
plt.savefig("/home/seyed/GT/Will/Matt/triangle_plots/col_%s.png"%(title), dpi=300)











