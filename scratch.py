
# coding: utf-8

# In[4]:



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


# In[5]:


NSNe=10000

## set-up some random positions
x,y,z = ((np.random.rand(3,NSNe)-0.5)*15).astype('f')


## fixed linking length of .1, for now
linkingLengths = np.ones(NSNe,dtype='f')*50.0/NSNe*10.

## launch times, fixed at 1 for now
launchTimes = np.ones(NSNe,dtype='f')

## cooling times 
coolingTimes = np.ones(NSNe,dtype='f')

## ids, just integers now 
ids = np.arange(NSNe,dtype='f')


# In[6]:


import ctypes

class SupernovaCluster(ctypes.Structure):
    pass

SupernovaCluster._fields_ = [
                ("xs", ctypes.POINTER(ctypes.c_float)),
                ("ys", ctypes.POINTER(ctypes.c_float)),
                ("zs", ctypes.POINTER(ctypes.c_float)),
                ("ids",ctypes.POINTER(ctypes.c_float)),
                ("linkingLengths", ctypes.POINTER(ctypes.c_float)),
                ("numNGB",ctypes.c_int),
                ("cluster_id",ctypes.c_int),
                ("NextCluster",ctypes.POINTER(SupernovaCluster))
            ]

head = SupernovaCluster()


# In[ ]:


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


# In[13]:


## skip the empty head node
head = head.NextCluster.contents
for i in xrange(numClusters):
    #print head.numNGB,'|',
    #print np.ctypeslib.as_array(head.ids,shape=(head.numNGB,))
    try:
        if i < (numClusters-1):
            head = head.NextCluster.contents
    except:
        print i,numClusters
        raise
    
print "Last one:",
print np.ctypeslib.as_array(head.ids,shape=(head.numNGB,))
print i,numClusters

