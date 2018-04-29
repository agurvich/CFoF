
# coding: utf-8

# In[27]:


import h5py
import sys
import getopt
import os
import random

from all_utils import *

import matplotlib.pyplot as plt 
from distinct_colours import get_distinct#,cm_linear,cm_plusmin

from readsnap import readsnap

import time


# In[28]:


NSNe=10
np.random.seed(516)

## set-up some random positions
x,y,z = ((np.random.rand(3,NSNe)-0.5)*15).astype('f')


## fixed linking length of .1, for now
linkingLengths = np.ones(NSNe,dtype='f')*50.0/NSNe*10*2

## launch times, fixed at 1 for now
launchTimes = np.ones(NSNe,dtype='f')

## cooling times 
coolingTimes = np.ones(NSNe,dtype='f')

## ids, just integers now 
ids = np.arange(NSNe,dtype='f')


# In[21]:


import ctypes

class SupernovaCluster(ctypes.Structure):
    pass

SupernovaCluster._fields_ = [
                ("xs", ctypes.POINTER(ctypes.c_float)),
                ("ys", ctypes.POINTER(ctypes.c_float)),
                ("zs", ctypes.POINTER(ctypes.c_float)),
                ("ids",ctypes.POINTER(ctypes.c_float)),
                ("launchTimes", ctypes.POINTER(ctypes.c_float)),
                ("coolingTimes", ctypes.POINTER(ctypes.c_float)),
                ("linkingLengths", ctypes.POINTER(ctypes.c_float)),
                ("numNGB",ctypes.c_int),
                ("cluster_id",ctypes.c_int),
                ("NextCluster",ctypes.POINTER(SupernovaCluster))
            ]

head = SupernovaCluster()


# In[26]:


exec_call = "/home/abg6257/CFoF/fof_sne.so"
c_obj = ctypes.CDLL(exec_call)

h_out_cast=ctypes.c_int
H_OUT=h_out_cast()

print "Executing c code"
init_time=time.time()
numClusters = c_obj.FoFNGB(
    ctypes.c_int(NSNe),
    x.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
    y.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
    z.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),

    launchTimes.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
    coolingTimes.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
    linkingLengths.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
    
    ids.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),

    ctypes.byref(head),
    ctypes.byref(H_OUT))

h=np.ctypeslib.as_array(H_OUT)
print numClusters,'many clusters found'
print time.time()-init_time,'s elapsed'
## skip the empty head node
head = head.NextCluster.contents


# In[23]:


def extractLinkedListValues(numClusters,head,write_to_file=0):    
    ## initialize arrays
    numNGBs=[]
    masterListIndices=[0]
    clusterIDs=[]
    keys = np.array(head._fields_)
    ## assumes that the final 3 in this list are numNGB,cluster_id,NextCluster (which they are, since I define _fields_ above)
    valss = [[] for i in xrange(len(keys)-3)]
    
    ## loop through links of the linked list and add their values to a single flattened array
    for i in xrange(numClusters):
        extractSNClusterObjValues(head,keys,valss,numNGBs,masterListIndices,clusterIDs)
        ## iterate the linked list
        try:
            if i < (numClusters-1):
                head = head.NextCluster.contents
        except:
            print i,numClusters
            raise
            
    ## unpack the flattened arrays
    flat_xss,flat_yss,flat_zss,flat_idss,flat_ltss,flat_ctss,flat_llss=valss
    if write_to_file:
        raise Exception("Unimplemented")
        with h5py.File("file.hdf5",'w') as handle:
            pass
        
    ##split the flattened arrays
    xss = splitFlattenedArray(flat_xss,masterListIndices)
    yss = splitFlattenedArray(flat_yss,masterListIndices)
    zss = splitFlattenedArray(flat_zss,masterListIndices)
    idss = splitFlattenedArray(flat_idss,masterListIndices)
    ltss = splitFlattenedArray(flat_ltss,masterListIndices)
    ctss = splitFlattenedArray(flat_ctss,masterListIndices)
    llss = splitFlattenedArray(flat_llss,masterListIndices)
    return xss,yss,zss,idss,ltss,ctss,llss,np.array(numNGBs),np.array(clusterIDs)
    
    
def splitFlattenedArray(flat_arr,masterListIndices):
    return np.split(flat_arr,masterListIndices[1:-1])

        
def extractSNClusterObjValues(head,keys,valss,numNGBs,masterListIndices,clusterIDs):
    for i,(key,val) in enumerate(keys):
        if key == 'numNGB':
            numNGBs+=[head.numNGB]
            masterListIndices+=[masterListIndices[-1]+head.numNGB]
        elif key == 'cluster_id':
            clusterIDs+=[head.cluster_id]
        elif key =='NextCluster':
            pass
        else:
            valss[i]=np.append(valss[i],[np.ctypeslib.as_array(getattr(head,key),shape=(head.numNGB,))])
            


# In[24]:



print dir(head)
print 

xss,yss,zss,idss,ltss,ctss,llss,numNGBs,clusterIDs=extractLinkedListValues(numClusters,head)

    #print np.ctypeslib.as_array(getattr(head,key),shape=head.numNGB)

    
"""
print "Last one:",
print np.ctypeslib.as_array(head.ids,shape=(head.numNGB,))
print i,numClusters
"""
print "",


# In[25]:


plt.hist(numNGBs)
index = list(np.logical_and(numNGBs>5,numNGBs<10)).index(1)

my_coords = np.array([xss[index],yss[index],zss[index]]).T


plt.figure()
plt.scatter(my_coords[:,0],my_coords[:,1],marker='x')
for (x,y,r) in zip(my_coords[:,0],my_coords[:,1],llss[index]):
    print x,y,r
    plt.gca().add_artist(plt.Circle((x,y),r,ls='--',lw=3,fill=None))
plt.axes().set_aspect('equal', 'datalim')


for coord in my_coords:
    dist = np.sum((my_coords-coord)**2,axis=1)/llss[index]**2
    print dist
print 'done'
print

