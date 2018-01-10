# -*- coding: utf-8 -*-
"""
Created on Wed Jan 06 10:26:07 2016

@author: Will Ratcliff
"""

import numpy as np
np.set_printoptions(threshold=np.nan)
import math
from scipy.stats import linregress
import matplotlib.pyplot as plt
import pandas as pd


cell_genos=np.linspace(1,2,10)
st_dev_cell=0.25
st_dev_group=0.25
timesteps_to_run=4 #this is equivalent to the number of generations the seed population of 10 clusters (each with a different genetic mean size) goes through.
number_of_stdevs=5 #this name is a bit confusing, but this is the number of different st_devs to be run in the iteration
replicates=5

#use these lists for running multiple iterations of cell and group stdevs
group_stdev_iterator=np.linspace(.0001,st_dev_cell,number_of_stdevs)
cell_stdev_iterator=np.linspace(.0001,st_dev_group,number_of_stdevs)


group_sizes=[]
st_dev_cells=[]
st_dev_groups=[]
slope_cell=[]
slope_group_volume=[]
slope_group_radius=[]
slope_group_settling=[]
stdev_list = [ [] for i in range(number_of_stdevs) ]	
group_sizes_list = [ [] for i in range(len(cell_genos)) ]

"""So, here's how the popu lation is initialized:
Col A: Cell ID (#)
Col B: Cell genetic size
Col C: Cellular phenotype
Col D: Cell parent
Col E: Cluster ID (#)
Col F: Cluster parent
Col G: empty
Col H: empty
Col I: empty
Col J: Cell parental phenotype"""
for ii in range(0,len(stdev_list)):
    for a in np.arange(2,12,2):	#'a' is the number of cells per group from 2 to 32
        for b in range(0,replicates):
            #print "cluster size %d, replicate %d", %(a,b)
            cluster_size=a
            pop=np.zeros((len(cell_genos)*cluster_size,10))
            st_dev_cell=cell_stdev_iterator[ii] #change this to st_dev_group=group_stdev_iterator[ii] to iterate across group-level variance
            #st_dev_group=group_stdev_iterator[ii]
            #initialize the population
            for i in range(0,np.shape(pop)[0]):
                pop[i][0]=i #Number each cell
                pop[i][1]=cell_genos[math.floor(i/cluster_size)]
                pop[i][2]=np.random.normal(pop[i][1],st_dev_cell)
                pop[i][4]=math.floor(i/cluster_size)  
            timestep=1
            #run through a round of reproduction
            for j in range(0,timesteps_to_run):
                cell_max=int(max(pop[:,0]))+1
                cluster_max=int(max(pop[:,4]))+1
                cells_added=len(cell_genos)*cluster_size*2**(timestep-1)*2
                cells_added_first=len(cell_genos)*cluster_size*2**(timestep-1) #this counts up the first reproductive event within the timepoint, total cells added is for both offspring clusters
                print("st_dev value %d of %d" %(ii, number_of_stdevs))
                print("cluster size", a)
                print("replicate", b)
                print("generation number", timestep)
                print("population size", len(cell_genos)*cluster_size*2**(timestep-1))
                print("number of cells added this timestep", cells_added)
                #first cluster produced
                cluster_variance_factor=np.random.normal(1,st_dev_group)
                for i in range(0,cells_added_first): #this loops through every additional cell for the first cluster offspring
                    if (cluster_max+math.floor(i/cluster_size))!=(cluster_max+math.floor((i-1)/cluster_size)): #if your cluster number does not equal the one lower down from you, then you get a new cluster-level variance factor. 
                        cluster_variance_factor=np.random.normal(1,st_dev_group)
                    pop=np.vstack([pop,[cell_max+i,pop[(cell_max+i)-cells_added_first][1],np.random.normal(pop[(cell_max+i)-cells_added_first][1],st_dev_cell)*cluster_variance_factor,pop[(cell_max+i)-cells_added_first][0],cluster_max+math.floor(i/cluster_size),pop[(cell_max+i)-cells_added_first][4],0,0,0,pop[(cell_max+i)-cells_added_first][2]]])
                cell_max=int(max(pop[:,0]))+1
                cluster_max=int(max(pop[:,4]))+1
                #second cluster produced
                for i in range(0,cells_added_first):
                    pop=np.vstack([pop,[cell_max+i,pop[(cell_max+i)-cells_added][1],np.random.normal(pop[(cell_max+i)-cells_added][1],st_dev_cell)*cluster_variance_factor,pop[(cell_max+i)-cells_added][0],cluster_max+math.floor(i/cluster_size),pop[(cell_max+i)-cells_added][4],0,0,0,pop[(cell_max+i)-cells_added][2]]])
                timestep+=1
            np.savetxt("full-population.csv", pop, delimiter=",") #this will save a CSV of the whole run, use for statistics or error-checking
            cell_x=pop[:,9]
            cell_y=pop[:,2]
            cell_x=cell_x[len(cell_genos)*cluster_size:]
            cell_y=cell_y[len(cell_genos)*cluster_size:]
            #linear regression of parent on offspring phenotype
            print("slope of parent-offspring regression for CELL size is", linregress(cell_x,cell_y)[0])
           
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
            parent_radius=[]
            offspring_radius=[]
            for i in range(0,len(parent_size_cleaned)):
                parent_radius.append((3.*parent_size_cleaned[i]/(4.*math.pi))**(1./3.)) #manual check of this calculation confirmed it is correct
            for i in range(0,len(offspring_size_cleaned)):
                offspring_radius.append((3.*offspring_size_cleaned[i]/(4.*math.pi))**(1./3.))
            parent_stokes=[]
            offspring_stokes=[]
            for i in range(0,len(parent_size_cleaned)):
                parent_stokes.append((9.81*(2*parent_radius[i])**2*(.1)) / (18.*1.002)) #Manual check of this calculation confirms it is correct. #9.81 is m/s gravity, then diameter (in Meters, which we might want to change!), then difference in density of particles from fluid, dividied by 18*the dynamic viscosity of water, which I chose 20deg C as the temp. http://www.calculatoredge.com/new/stroke.htm and http://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html
            for i in range(0,len(offspring_size_cleaned)):
                offspring_stokes.append((9.81*(2*offspring_radius[i])**2*(.1)) / (18.*1.002)) #9.81 is m/s gravity, then diameter, then difference in density of particles from fluid, dividied by 18*the dynamic viscosity of water, which I chose 20deg C as the temp. http://www.calculatoredge.com/new/stroke.htm and http://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html        
            print("slope of parent-offspring regression for GROUP volume is", linregress(parent_size_cleaned,offspring_size_cleaned)[0])
            print("slope of parent-offspring regression for GROUP radius is", linregress(parent_radius,offspring_radius)[0])
            print("slope of parent-offspring regression for GROUP settling speed is", linregress(parent_stokes,offspring_stokes)[0])
            print("size", parent_size_cleaned[1], len(parent_size_cleaned))
            print("radius", parent_radius[1], len(parent_radius))
            print("stokes", parent_stokes[1], len(parent_stokes))
            #group_sizes.append(a)
            group_sizes_list[ii].append(a)
            #slope_cell.append(linregress(cell_x,cell_y)[0])
            #slope_group_volume.append(linregress(parent_size_cleaned,offspring_size_cleaned)[0])
            #slope_group_radius.append(linregress(parent_radius,offspring_radius)[0])
            #slope_group_settling.append(linregress(parent_stokes,offspring_stokes)[0])
            print("heritability groups", (linregress(parent_size_cleaned,offspring_size_cleaned)[0]))
            print("heritability cells", (linregress(cell_x,cell_y)[0]))
            tempratio=(linregress(parent_size_cleaned,offspring_size_cleaned)[0]) / (linregress(cell_x,cell_y)[0])   
            stdev_list[ii].append(tempratio)


cmap = plt.get_cmap('rainbow')
colors = [cmap(i) for i in np.linspace(0, 1, len(stdev_list))]

#print "group_sizes", group_sizes, len(group_sizes)
#print "stdev_list[4]", stdev_list[0], len(stdev_list[0])
#print "group_sizes_list[4]", group_sizes_list[0], len(group_sizes_list[0])

for i, color in enumerate(colors, start=0):
    plt.scatter(group_sizes_list[i],stdev_list[i], color=color, alpha=.5)
plt.xlabel('Number of cells per group')
plt.ylabel('Ratio of group to cellular-level heritability for size')
plt.savefig("Ratio of group to cell heritability iterator=cell variance, group sd=%s.png" %(st_dev_group), dpi=300)
plt.savefig("Ratio of group to cell heritability iterator=cell variance, group sd=%s.pdf" %(st_dev_group))
