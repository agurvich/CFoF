
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

NSNe=10000

a = np.linspace(1,19,NSNe)
b = np.linspace(11,21,NSNe)

a = np.array(a,dtype='f',ndmin=1)
b = np.array(b,dtype='f',ndmin=1)




## set-up some random positions
x,y,z = ((np.random.rand(3,NSNe)-0.5)*15).astype('f')

"""
print 'python distances'
print x
print y
print z
print '---------------------'

print (x-x[0])**2+(y-y[0])**2+(z-z[0])**2
"""

## fixed linking length of .1, for now
linkingLengths = np.ones(NSNe,dtype='f')*50.0/NSNe*10.

## launch times, fixed at 1 for now
launchTimes = np.ones(NSNe,dtype='f')

## cooling times 
coolingTimes = np.ones(NSNe,dtype='f')

## ids, just integers now 
ids = np.arange(NSNe,dtype='f')


import ctypes
class Supernova(ctypes.Structure):
    _fields_ = [
                ("x", ctypes.c_float),
                ("y", ctypes.c_float)]

class LLSupernova(ctypes.Structure):
    pass
LLSupernova._fields_ = [
                ("x", ctypes.c_float),
                ("y", ctypes.c_float),
                ("next_LLSN",ctypes.POINTER(LLSupernova))]

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



"""
c = (Supernova*5)()
c.x = 2.0
c.y = 10.0


"""
head = SupernovaCluster()

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
    a.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
    b.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
    #ctypes.byref(c),

    ctypes.byref(head),
    ctypes.byref(H_OUT))

h=np.ctypeslib.as_array(H_OUT)
print h,'many clusters found (h)'
print numClusters,'many clusters found'
print time.time()-init_time,'s elapsed'


## skip the empty head node
head = head.NextCluster.contents
for i in xrange(numClusters):
    #print head.numNGB,'|',
    #print np.ctypeslib.as_array(head.ids,shape=(head.numNGB,))
    if i == (numClusters-2):
        print np.ctypeslib.as_array(head.ids,shape=(head.numNGB,))
    try:
        if i < (numClusters-1):
            head = head.NextCluster.contents
    except:
        print i,numClusters
        raise
    
print "Last one:",
print np.ctypeslib.as_array(head.ids,shape=(head.numNGB,))
print i,numClusters

"""
print c[1].x
print c[1].y
"""

## step through the linked list
"""
e = d.next_LLSN

print d.x
print d.y

print dir(e)
print e.contents.x
print e.contents.y

print '29 and 13'
print e.contents.next_LLSN.contents.x
print e.contents.next_LLSN.contents.y
"""

