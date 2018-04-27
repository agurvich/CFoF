#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


struct SupernovaCluster{
    float * xs;
    float * ys;
    float * zs;
    float * ids;
    float * linkingLengths;
    int numNGB;
    int cluster_id;
};

/*
struct SupernovaCluster * initialize_SupernovaCluster(
    float xs, float ys, float zs){
    struct LLSupernova *new_LLSN = malloc(sizeof(struct LLSupernova)); 
}
*/

typedef struct Supernova{
    float x;
    float y;
}Supernova;

struct LLSupernova{
    float x;
    float y; 
    struct LLSupernova * next_LLSN; 
};

struct LLSupernova * initialize_LLSupernova(float argx, float argy){
    struct LLSupernova *new_LLSN = malloc(sizeof(struct LLSupernova)); 

    new_LLSN->x = argx; 
    new_LLSN->y = argy; 
    new_LLSN->next_LLSN = NULL;
  //for (i = 0; i < n_part; i++) 
   // new_clump->coordinates[i] = (float *) malloc(3 * sizeof(float)); 
    return new_LLSN;
};


void * findNGBFlags(
    int Narr,
    float * dists,
    float * linkingLengths,
    float link_node,
    int * ngbflags, int * numNGB){
    for (int i=0; i<Narr;i++){
        if (dists[i] < link_node || dists[i] < linkingLengths[i]){
            ngbflags[i]=1;
            (*numNGB)++;
        }
    }
}

void calculateDists(float * point, float * xs, float * ys, float * zs,int Narr, float * dists){
    for (int i =0; i<Narr; i++){
        float dx = point[0]-xs[i];
        float dy = point[1]-ys[i];
        float dz = point[2]-zs[i];
        dists[i]=dx*dx + dy*dy + dz*dz;
    }
}

int fillFlags(
    float * point, 
    float * xs, float * ys,float * zs,
    float * linkingLengths, float link_node,
    int Narr, int * NGBFlags){

    int numNGB;
    float * dists = (float *) malloc(Narr*sizeof(float));
    // calculate the distance to each point
    calculateDists(point,xs,ys,zs,Narr,dists);

    //printArray(dists,Narr);
    findNGBFlags(
        Narr,
        dists,
        linkingLengths,link_node,
        NGBFlags,&numNGB);
}


void getIndicesFromFlags(int * NGBFlags, int Narr, int * NGBIndices){
    int filled=0;
    for (int i=0; i<Narr;i++){
        if (NGBFlags[i]){
            NGBIndices[filled]=i;
            filled++;
        }
    }
}

void extractSubarrayWithIndices(
    float * array, float * subarray, 
    int * NGBIndices, 
    int Narray, int numNGB,
    int checkNGBIndices){
    for (int i=0; i<numNGB; i++){
        int ngbindex = NGBIndices[i];
        //extract the value
        subarray[i]= array[ngbindex];

        // replace it with the value at the end
        array[ngbindex]=array[Narray-1-i];
        
        // only needs to be done once, since we'll be following the same pattern
        // of extraction/replacement for every array
        if (checkNGBIndices){
            // check if the value we just moved is actually a neighbor
            for (int j=numNGB; j>0; j--){
                // if it is, then 
                if (Narray-1-i == NGBIndices[j]){
                    NGBIndices[j]=ngbindex;
                    break;
                }
                // no need to go any further back, those indices have already been extracted
                if (i>=j)break;
            }
        }
    }
}

void printArray(float * arr,int Narr){
    for (int i=0; i< Narr; i++){
        printf("%.2f\t",arr[i]);
    }
    printf("\n");
}

void printIntArray(int * arr,int Narr){
    for (int i=0; i< Narr; i++){
        printf("%d\t",arr[i]);
    }
    printf("\n");
}

struct SupernovaCluster * findFriends(
    float * xs, float * ys, float * zs, 
    float * linkingLengths,
    float * ids,
    int Narr,
    int cluster_id){
    
    // initialize NGB variables
    int * NGBFlags;
    int * NGBIndices;
    int numNGB=0;

    // allocate ngbflags memory
    NGBFlags=(int*)malloc(Narr*sizeof(int));
    memset(NGBFlags,0,Narr*sizeof(int));

    // find all the neighbors to the first point 
    // (inclusive, so we're guaranteed to get at least one...)
    float point[3];
    point[0]=xs[0];
    point[1]=ys[0];
    point[2]=zs[0];

    numNGB=fillFlags(
        point,
        xs,ys,zs,
        linkingLengths,linkingLengths[0],
        Narr,NGBFlags);
    
    // "worst case scenario", all the particles are in this cluster
    float * buffer_xs=(float*)malloc(Narr*sizeof(float));
    float * buffer_ys=(float*)malloc(Narr*sizeof(float));
    float * buffer_zs=(float*)malloc(Narr*sizeof(float));
    float * buffer_ids=(float*)malloc(Narr*sizeof(float));
    float * buffer_linkingLengths=(float*)malloc(Narr*sizeof(float));

    // extract the sub arrays from their indices
        // only need to recalculate NGB indices on the first pass 
        // and could in principal have a separate function that does this
        // but I think this is easier to wrap one's head around
        // essentially it's just if NGBFlags[j] && NGBFlags[N-1-j] -> NGBIndices[N-1-j]=j
    extractSubarrayWithIndices(xs,buffer_xs,NGBIndices,Narr,numNGB,1);
    extractSubarrayWithIndices(ys,buffer_ys,NGBIndices,Narr,numNGB,0);
    extractSubarrayWithIndices(zs,buffer_zs,NGBIndices,Narr,numNGB,0);
    extractSubarrayWithIndices(ids,buffer_ids,NGBIndices,Narr,numNGB,0);
    extractSubarrayWithIndices(linkingLengths,buffer_linkingLengths,NGBIndices,Narr,numNGB,0);

    // new size of xs/ys/zs/...
    int Nremain=Narr-numNGB;
    int Nadded,numNewNGB;
    int cur_ngb=1; // don't need to check the first neighbor, we just did that above
    while (cur_ngb < numNGB){
        Nadded=0;
        // reset the first Nremain elements NGB flags and indices array
        // dangerous if you made a typo below, but efficient otherwise
        memset(NGBFlags,0,Nremain*sizeof(int));

        // start the loop at the first friend (of friend (of friend ...)) that we haven't checked
        //TODO pragma omp parallel accumulate Nadded accumulate private point
        for (int j=cur_ngb; j<numNGB; j++){
            point[0]=buffer_xs[j];
            point[1]=buffer_ys[j];
            point[2]=buffer_zs[j];
        
            // this is slightly inefficient in the case where one remaining point
                // is multiple neighbors' friend. Since we just set the flag to 1 
                // the race condition is not important, and we just waste a little
                // time calculating the distance a few times
            Nadded+= fillFlags(
                point,
                xs,ys,zs,
                linkingLengths,linkingLengths[j],
                Nremain,NGBFlags);

            // increment the counter that we've checked this neighbor and added its friends

        }
        
        // set cur_ngb to the end of the list
        cur_ngb=numNGB;

        if (Nadded){
            // find new NGB indices
            memset(NGBIndices,0,Nadded*sizeof(int));
            getIndicesFromFlags(NGBFlags,Nremain,NGBIndices);

            // extract all our new friends into the buffer, removing them from the main array
            // increment the pointer to point to the first open slot of the buffer
            extractSubarrayWithIndices(
                xs,buffer_xs+sizeof(float)*numNGB,
                NGBIndices,Nremain,Nadded,1);
            extractSubarrayWithIndices(
                ys,buffer_ys+sizeof(float)*numNGB,
                NGBIndices,Nremain,Nadded,0);
            extractSubarrayWithIndices(
                zs,buffer_zs+sizeof(float)*numNGB,
                NGBIndices,Nremain,Nadded,0);
            extractSubarrayWithIndices(
                ids,buffer_ids+sizeof(float)*numNGB,
                NGBIndices,Nremain,Nadded,0);
            extractSubarrayWithIndices(
                linkingLengths+sizeof(float)*numNGB,buffer_linkingLengths,
                NGBIndices,Nremain,Nadded,0);

            // update the size of the neighbors that live in the buffer
            // and the remaining particles
            numNGB+=Nadded;
            Nremain-=Nadded;
        } //if Nadded
    } //while cur_ngb < numNGB

    struct SupernovaCluster *new_cluster = malloc(sizeof(struct SupernovaCluster)); 

    // set the number of neighbors in the cluster
    new_cluster->cluster_id=cluster_id;
    new_cluster->numNGB = numNGB;

    // allocate the cluster arrays
    new_cluster->xs=(float*)malloc(numNGB*sizeof(float));
    new_cluster->ys=(float*)malloc(numNGB*sizeof(float));
    new_cluster->zs=(float*)malloc(numNGB*sizeof(float));
    new_cluster->ids=(float*)malloc(numNGB*sizeof(float));
    new_cluster->linkingLengths=(float*)malloc(numNGB*sizeof(float));
    //TODO memcopy from the buffer to the new_cluster array
 
    return new_cluster;

}// void findFriends

// Supernova * c, struct LLSupernova * d,
int FoFNGB(
    int Narr,
    float * xs, float * ys, float * zs,
    float * launchTimes, float * coolingTimes,
    float * linkingLengths,
    float * ids,
    float * arr, float * second, 
    float * H_OUT ){


    int returnVal;

    int cluster_id=0;

    while (Narr > 0){
        // take our point to be the first SNe in the list


        printf("------------------------------------------------------\n");
        printf("Working on cluster %d\n",cluster_id);
        printf("------------------------------------------------------\n");
        struct SupernovaCluster *new_cluster = findFriends(
            xs,ys,zs,
            linkingLengths,
            ids,
            Narr,
            cluster_id);

        Narr-=new_cluster->numNGB;
        printf("%d members found\n",new_cluster->numNGB);
        printArray(new_cluster->ids,new_cluster->numNGB);
        printf("%d elements remain\n",Narr);

        cluster_id++;

    } // while (Narr > 0)



    /*
    c->x=5.0;
    c->y=20.0;

    c++;

    c->x=15;
    c->y=45;


    d->x = 16;
    d->y = 18;
    //printf("Trying to set next LLSupernova... \n");
    d->next_LLSN = initialize_LLSupernova(3,7); 
    // set 3rd llsn
    d->next_LLSN->next_LLSN = initialize_LLSupernova(29,13); 

    //printf("Successfully set next LLSupernova! \n");
    */

    return 0;
}


/*
int stellargasdensity(int N_gas, int N_star,
                      float* x_gas, float* y_gas, float* z_gas, float* m_gas,
                      float* x_star, float* y_star, float* z_star,
                      int DesNgb, float Hmax, long long* id_gas, float* H_OUT)
{
    printf("Undefined");
}



struct Supernova * initializeSupernova(
    double x, double y, double z,
    float launch_time, int nsne){
    struct Supernova *new_supernova = malloc(sizeof(struct Supernova));

    // save supernova properties
    new_supernova->x = x;
    new_supernova->y = y;
    new_supernova->z = z;

    new_supernova->launch_time = launch_time;
    new_supernova->nsne = nsne;

}

void findClusters(smart_list,limit=None):
    Top level loop for finding freinds of friends clusters. Loops through each object 
     * in smart_list and finds the entire FoF cluster it belongs to recursively. 
     * Input:
     *      smart_list - a list of smart object based supernovae
     *      limit - limit to the linking length and tcool specified as a density
     *      floor in physical units, uses obj.get_cool(max(obj.density,limit)) 

    N_remain = 
    popped=0
    clusters=[]
    init_length=len(smart_list)
    while (N_remain > 0){
        N_remain--;
        obj=smart_list.pop(0)
        cluster=[obj]

        cluster+=findCluster(obj,smart_list,cnum=len(clusters),limit=limit)
        popped+=len(cluster)
        clusters+=[cluster]
    }
    return clusters


void checkMembership(Supernova obj,Supernova other)
{
     * Function to determine if an object and an other are in the same 
     * cluster, based on some membership criterion. Here we check if either 
     * are contained in the other's rcool within a time period tcool
        
    dist=distance(obj,other)
    eps2=(3e-4/5)**2

    if ((obj.rcool > dist
        and 0 < other.launch_time-obj.launch_time < obj.tcool) 
        or 
        (other.rcool > dist
        and 0 < obj.launch_time-other.launch_time < other.tcool)):
        return True

    return False
}

void findCluster(obj,smart_list,cnum,fof=1,**kwargs):

     * Recursive level call to find friends of friends cluster for object. This
        is done by finding the first object in smart_list that is a member of object's
        cluster, here called "other," then finding other's cluster (in the same way) 
        and joining the two.
        Input: 
            obj - the object whose cluster you're trying to build 
            smart_list - the list of objects you are searching for other in. 
            cnum - the number of the cluster, useful for identifying it on the
                STDOUT but not really important
            fof - a flag for if you want to use fof to find the cluster. Alternatively
                you could just loop through and get some biased result based on the order 
                of your cluster. This is not reccomended. 

    cluster=[]
    popped=0
    popped_from_the_future=0
    init_length=len(smart_list)
    for i in xrange(len(smart_list)){
        other=smart_list[i-popped]

        if checkMembership(obj,other,**kwargs){
            member=smart_list.pop(i-popped)
            cluster+=[member]
            popped+=1
            checkListPop(init_length,smart_list,popped)

            print 'found an independent member of the cluster #',cnum


            if fof{
                #friends of friends part
                new_cluster=findCluster(other,smart_list,cnum=cnum)
                popped_from_new_cluster=len(new_cluster)
                popped_from_the_future+=popped_from_new_cluster
                popped+=popped_from_new_cluster
                print len(new_cluster),'friends found'
                checkListPop(init_length,smart_list,popped)
                cluster+=new_cluster
                }
            }
        }
       
        if i < popped_from_the_future:
            #my DumbObject test tells me that this break statement 
            #stops the recursion from letting the new_cluster variable
            #be larger than smart_list
            break

    return cluster

*/

