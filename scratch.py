
import h5py
import sys
import getopt
import os
import random

from all_utils import *

import matplotlib.pyplot as plt 
from distinct_colours import get_distinct#,cm_linear,cm_plusmin

from readsnap import readsnap

a = np.linspace(0,10,11)
b = np.linspace(10,20,11)

a = np.array(a,dtype='f',ndmin=1)
b = np.array(b,dtype='f',ndmin=1)

import ctypes
class Supernova(ctypes.Structure):
    _fields_ = [
                ("x", ctypes.c_float),
                ("y", ctypes.c_float)]

c = Supernova()




c.x = 2.0
c.y = 10.0
print c.x,c.y


exec_call = "/home/abg6257/CFoF/fof_sne.so"
c_obj = ctypes.CDLL(exec_call)

h_out_cast=ctypes.c_float*11
H_OUT=h_out_cast()

c_obj.add_arrays(
    ctypes.c_int(11),
    a.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
    b.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
    ctypes.byref(c),
    ctypes.byref(H_OUT))

h=np.ctypeslib.as_array(H_OUT)
print h

print c.x,c.y





