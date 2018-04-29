import ctypes
import numpy as np
import time
import h5py

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


def findFoFClustering(xs,ys,zs,ids,launchTimes,coolingTimes,linkingLengths):
    NSNe = len(xs)
    return extractLinkedListValues(
        *getLinkedListHead(NSNe,xs,ys,zs,ids,launchTimes,coolingTimes,linkingLengths)
        ,write_to_filename='test')

def getLinkedListHead(NSNe,xs,ys,zs,ids,launchTimes,coolingTimes,linkingLengths):
    head = SupernovaCluster()    
    exec_call = "/home/abg6257/CFoF/fof_sne.so"
    c_obj = ctypes.CDLL(exec_call)

    h_out_cast=ctypes.c_int
    H_OUT=h_out_cast()

    print "Executing c code"
    init_time=time.time()
    numClusters = c_obj.FoFNGB(
        ctypes.c_int(NSNe),
        xs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        ys.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        zs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),

        launchTimes.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        coolingTimes.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        linkingLengths.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),

        ids.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),

        ctypes.byref(head),
        ctypes.byref(H_OUT))


    print numClusters,'many clusters found'
    print time.time()-init_time,'s elapsed'
    ## skip the empty head node
    head = head.NextCluster.contents
    return numClusters,head


def extractLinkedListValues(numClusters,head,write_to_filename=0):    
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
    if write_to_filename:
        ## overwrites
        with h5py.File(write_to_filename+".hdf5",'w') as handle:
            ## set the data group's values, lazily
            group=handle.create_group('ClusterData')
            for key,val in zip(keys[:-3],valss):
                group['flat_'+key[0]+'s']=val
                
            group["masterListIndices"]=masterListIndices
            group["numNGBs"]=numNGBs
            
            print handle.keys()
            print handle['ClusterData'].keys()
        
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
            
