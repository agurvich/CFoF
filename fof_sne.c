#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


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


int popArray(int * array, int index,int *Narr){
    
    *Narr-=1;
    return array[index];

}

int twofunc(int var){
    return var > 2;
}

//
// my first function pointer, I'm so proud!!
void * findNGBFlags(
    int * array,int * Narr,
    int boolfunc(int),
    int * ngbflags, int * numNGB){
    for (int i=0; i<*Narr;i++){
        if (boolfunc(array[i])){
            ngbflags[i]=1;
            (*numNGB)++;

            // delete this particle by replacing it with the last particle
            // since the order doesn't matter, it's all good
            //array[i]=array[*Narr-*subNarr-1];
            //(*subNarr)++;
        }
    }
    // reduce the size by the number of elements popped
    //*Narr-=*subNarr;

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
    int * array, int * subarray, 
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

//int add_arrays(int N, float * a, float * b, float * H_OUT);
//
void printArray(int * arr,int Narr){
    for (int i=0; i< Narr; i++){
        printf("%d\t",arr[i]);
    }
    printf("\n");
}


int add_arrays(int N, float * a, float * b, Supernova * c, struct LLSupernova * d, float * H_OUT ){
    for (int i = 0; i<N; i++){
        H_OUT[i]= a[i]+b[i];
    }

    
    int arr[5],second[5];
    int Narr=5,subNarr=0;
    
    // make sure the array is 0'd out


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

    printArray(arr,Narr);
    printArray(second,Narr);
    printf("------------------\n");

    // initialize NGB variables
    int numNGB=0;
    int * NGBFlags;
    int * NGBIndices;
    // allocate ngbflags memory
    NGBFlags=(int*)malloc(Narr*sizeof(int));
    memset(NGBFlags,0,Narr*sizeof(int));


    // find the neighbor indices
    findNGBFlags(arr,&Narr,twofunc,NGBFlags,&numNGB);
    
    // allocate memory for, and find, NGB indices
    NGBIndices=(int*)malloc(numNGB*sizeof(int));
    getIndicesFromFlags(NGBFlags,Narr,NGBIndices);

    // allocate the buffer array
    int * sub_arr;
    int * subsecond;
    sub_arr=(int*)malloc(numNGB*sizeof(int));
    subsecond=(int*)malloc(numNGB*sizeof(int));


    printArray(NGBIndices,numNGB);
    extractSubarrayWithIndices(arr,sub_arr,NGBIndices,Narr,numNGB,1);
    extractSubarrayWithIndices(second,subsecond,NGBIndices,Narr,numNGB,0);
    
    printArray(NGBFlags,Narr);
    printArray(NGBIndices,numNGB);
    printArray(sub_arr,numNGB);
    printArray(subsecond,numNGB);

    printf("------------------\n");


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

