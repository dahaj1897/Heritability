# -*- coding: utf-8 -*-
"""
Created on Wed Jan 06 10:26:07 2016

@author: wratcliff3
"""

import numpy as np
np.set_printoptions(threshold=np.nan)
import math
from scipy.stats import linregress
import matplotlib.pyplot as plt
import pandas as pd



cell_genos=np.linspace(1,2,10)
cluster_size=32
st_dev_cell=0.25
st_dev_group=0.25
timesteps_to_run=4 #this is equivalent to the number of generations the seed population of 10 clusters (each with a different genetic mean size) goes through.
replicates=1 #keep to 1 for now until an averaging function is implemented, otherwise it won't plot in 3d
granularity=32 #number of x and y values to calculate (higher=more fine grained). This will generate an n x n grid for plotting relative heritability (in the paper we use 32 x 32 = 1024 points).



group_sizes=[]
st_dev_cells=[]
st_dev_groups=[]
slope_cell=[]
slope_group_volume=[]
slope_group_radius=[]
slope_group_settling=[]
cell_stdev_plotting_list=[]
group_stdev_plotting_list=[]
relative_heritability=[]

"""So, here's how the population is initialized:
Col 1: Cell ID (#)
Col 2: Cell genetic size
Col 3: Cellular phenotype
Col 4: Cell parent
Col 5: Cluster ID (#)
Col 6: Cluster parent
Col 7: empty
Col 8: empty
Col 9: empty
Col 10: Cell parental phenotype"""
for ii in np.linspace(0.0001,st_dev_cell,granularity):
    for a in np.linspace(0.0001,st_dev_group,granularity):
        for b in range(0,replicates):
            pop=np.zeros((len(cell_genos)*cluster_size,10))
            st_dev_cell=ii
            st_dev_group=a
            #initialize the population
            for i in range(0,np.shape(pop)[0]):
                pop[i][0]=i #Number each cell
                pop[i][1]=cell_genos[math.floor(i/cluster_size)]
                pop[i][2]=np.random.normal(pop[i][1],st_dev_cell)
                pop[i][4]=math.floor(i/cluster_size)  
            timestep=1
            #run through a round of reproduction
            for j in range(0,timesteps_to_run):
                #print pop
                cell_max=int(max(pop[:,0]))+1
                cluster_max=int(max(pop[:,4]))+1
                cells_added=len(cell_genos)*cluster_size*2**(timestep-1)*2
                cells_added_first=len(cell_genos)*cluster_size*2**(timestep-1) #this counts up the first reproductive event within the timepoint, total cells added is for both offspring clusters
                #print"cell stdev", ii
               # print "group stdev", a
                #print "generation number", timestep
                #first cluster produced
                cluster_variance_factor=np.random.normal(1,st_dev_group)
                for i in range(0,cells_added_first): #this loops through every additional cell for the first cluster offspring
                    if (cluster_max+math.floor(i/cluster_size))!=(cluster_max+math.floor((i-1)/cluster_size)): #if your cluster number does not equal the one lower down from you, then you get a new cluster-level variance factor. 
                        cluster_variance_factor=np.random.normal(1,st_dev_group)
                    #print "cluster_variance_factor and cluster number is", cluster_variance_factor, (cluster_max+math.floor(i/cluster_size))    #I've confirmed that every cell in the group is getting the same environmental variance               
                    pop=np.vstack([pop,[cell_max+i,pop[(cell_max+i)-cells_added_first][1],np.random.normal(pop[(cell_max+i)-cells_added_first][1],st_dev_cell)*cluster_variance_factor,pop[(cell_max+i)-cells_added_first][0],cluster_max+math.floor(i/cluster_size),pop[(cell_max+i)-cells_added_first][4],0,0,0,pop[(cell_max+i)-cells_added_first][2]]])
                cell_max=int(max(pop[:,0]))+1
                cluster_max=int(max(pop[:,4]))+1
                #second cluster produced
                for i in range(0,cells_added_first):
                    pop=np.vstack([pop,[cell_max+i,pop[(cell_max+i)-cells_added][1],np.random.normal(pop[(cell_max+i)-cells_added][1],st_dev_cell)*cluster_variance_factor,pop[(cell_max+i)-cells_added][0],cluster_max+math.floor(i/cluster_size),pop[(cell_max+i)-cells_added][4],0,0,0,pop[(cell_max+i)-cells_added][2]]])
                timestep+=1
            #np.savetxt("full population.csv", pop, delimiter=",") #this will save a CSV of the whole run, useful for statistics or error checking
            cell_x=pop[:,9]
            cell_y=pop[:,2]
            cell_x=cell_x[len(cell_genos)*cluster_size:]
            cell_y=cell_y[len(cell_genos)*cluster_size:]
            #linear regression of parent on offspring phenotype
            
            #Pandas dataframe work isolating groups
            df=pd.DataFrame(pop)
            size_by_ID=df.groupby(4)[2].sum()
            parent_by_ID=df.groupby(4)[5].mean()
            joined=pd.concat([size_by_ID,parent_by_ID], axis=1, ignore_index=True)
            parent_size=[]
            for i in range(0,len(joined[0])):
                j=joined[1][i]
                parent_size.append(joined[0][j])  
            offspring_size=joined[0]
            parent_size_cleaned=list(parent_size[len(cell_genos):])
            offspring_size_cleaned=list(offspring_size[len(cell_genos):])
            tempratio=(linregress(parent_size_cleaned,offspring_size_cleaned)[0]) / (linregress(cell_x,cell_y)[0])   
            cell_stdev_plotting_list.append(ii)
            group_stdev_plotting_list.append(a)
            relative_heritability.append(tempratio)

#print "cell_stdev_plotting_list", cell_stdev_plotting_list
#print "group_stdev_plotting_list", group_stdev_plotting_list
#print "relative_heritability", relative_heritability
np.savetxt("plot_me!.csv", zip(cell_stdev_plotting_list, group_stdev_plotting_list, relative_heritability), delimiter=",")

