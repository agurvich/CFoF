#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/*
struct SupernovaCluster{

}

struct SupernovaCluster * initialize_SupernovaCluster(
    float xs, float ys, float zs){
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

void getIndicesFromFlags(int * NGBFlags, int Narr, int * NGBIndices){
    int filled=0;
    for (int i=0; i<Narr;i++){
        if (NGBFlags[i]){
            NGBIndices[filled]=i;
            filled++;
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

int findFriends(
    float * point,
    float * xs, float * ys, float * zs, 
    float * linkingLengths,
    float * ids,
    int Narr,
    int * pnumNGB){

    // initialize NGB variables
    int * NGBFlags;
    int * NGBIndices;

    float * dists;

    
    float * sub_arr;
    float * subsecond;
    
    //numNGB=0;

    // allocate ngbflags memory
    NGBFlags=(int*)malloc(Narr*sizeof(int));
    memset(NGBFlags,0,Narr*sizeof(int));


    // find the neighbor indices

    dists = (float *) malloc(Narr*sizeof(float));
    // calculate the distance to each point
    calculateDists(point,xs,ys,zs,Narr,dists);
    //printArray(xs,Narr);
    //printArray(ys,Narr);
    //printArray(zs,Narr);
    //printf("------------------\n");
    printArray(dists,Narr);
    findNGBFlags(
        Narr,
        dists,
        linkingLengths,linkingLengths[0],
        NGBFlags,pnumNGB);


    printIntArray(NGBFlags,Narr);
    
    // allocate memory for, and find, NGB indices
    NGBIndices=(int*)malloc(*pnumNGB*sizeof(int));
    getIndicesFromFlags(NGBFlags,Narr,NGBIndices);
    //printIntArray(NGBIndices,numNGB);

    // allocate the buffer array
    sub_arr=(float*)malloc(*pnumNGB*sizeof(float));
    subsecond=(float*)malloc(*pnumNGB*sizeof(float));


    // extract the sub arrays from their indices
        // only need to recalculate NGB indices on the first pass 
        // and could in principal have a separate function that does this
        // but I think this is easier to wrap one's head around
        // essentially it's just if NGBFlags[j] && NGBFlags[N-1-j] -> NGBIndices[N-1-j]=j
    extractSubarrayWithIndices(dists,sub_arr,NGBIndices,Narr,*pnumNGB,1);
    extractSubarrayWithIndices(xs,subsecond,NGBIndices,Narr,*pnumNGB,0);
    extractSubarrayWithIndices(ys,subsecond,NGBIndices,Narr,*pnumNGB,0);
    extractSubarrayWithIndices(zs,subsecond,NGBIndices,Narr,*pnumNGB,0);
    extractSubarrayWithIndices(ids,subsecond,NGBIndices,Narr,*pnumNGB,0);
    
    printIntArray(NGBFlags,Narr);
    printIntArray(NGBIndices,*pnumNGB);
    printf("dists \t");
    printArray(sub_arr,*pnumNGB);

    printf("ids \t");
    printArray(subsecond,*pnumNGB);

    printArray(ids,Narr-*pnumNGB);
    return *pnumNGB;
}

// Supernova * c, struct LLSupernova * d,
int FoFNGB(
    int Narr,
    float * xs, float * ys, float * zs,
    float * launchTimes, float * coolingTimes,
    float * linkingLengths,
    float * ids,
    float * arr, float * second, 
    float * H_OUT ){
    /*
    for (int i = 0; i<N; i++){
        H_OUT[i]= a[i]+b[i];
    }


    int Narr=5;
    float arr[5],second[5];

    arr[0]=5;
    arr[1]=10;
    arr[2]=1;
    arr[3]=7;
    arr[4]=0;

    second[0]=1;
    second[1]=2;
    second[2]=3;
    second[3]=4;
    second[4]=5;
    */

    // declare all our variables
    float point[3];

    int numNGB;

    while (Narr > 0){
        numNGB=0;
        // take our point to be the first SNe in the list
        point[0]=xs[0];
        point[1]=ys[0];
        point[2]=zs[0];

        printf("------------------\n");

        numNGB = findFriends(
            point,
            xs,ys,zs,
            linkingLengths,
            ids,
            Narr,
            &numNGB);

        Narr-=numNGB;
        printf("%d many elements remain\n",Narr);

        printf("------------------\n");


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

