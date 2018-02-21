#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


int add_arrays(int num_elements, float * a, float * b, float * H_OUT);

int add_arrays(int num_elements, float * a, float * b, float * H_OUT ){
    for (int i = 0; i<num_elements; i++){
        //H_OUT[i]= a[i]+b[i];
        H_OUT[i]=5;
    }
    H_OUT[0]=5;

    printf("Hello World!");
    return 1;
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

