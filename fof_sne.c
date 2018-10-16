#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


//#define DEBUG
//#define SUPERDEBUG
//#define PRETTYPRINT

void printArray(float * arr,int Narr){
    for (int i=0; i< Narr; i++){
        printf("%.2f\t",arr[i]);
    }
    printf("\n");
}

void printArrayFixedRatio(float * arr,float divisor, int Narr){
    for (int i=0; i< Narr; i++){
        printf("%.2f\t",arr[i]/(divisor*divisor));
    }
    printf("\n");
}

void printArrayRatio(float * arr,float * arr1, int Narr){
    for (int i=0; i< Narr; i++){
        printf("%.2f\t",arr[i]/(arr1[i]*arr1[i]));
    }
    printf("\n");
}


void printIntArray(int * arr,int Narr){
    for (int i=0; i< Narr; i++){
        printf("%d\t",arr[i]);
    }
    printf("\n");
}

struct SupernovaCluster{
    float * xs;
    float * ys;
    float * zs;
    float * ids;
    float * launchTimes;
    float * coolingTimes;
    float * linkingLengths;
    int numNGB;
    int cluster_id;
    struct SupernovaCluster  * NextCluster; 
};

void * findNGBFlags(
    int Narr,
    float * dists2,
    float * launchTimes,float launch_node,
    float * coolingTimes,float cool_node,
    float * linkingLengths,float link_node,
    int * ngbflags, int * numNGB){

    float launch_leaf,cool_leaf,link_leaf;
    for (int i=0; i<Narr;i++){
        launch_leaf = launchTimes[i];
        cool_leaf = coolingTimes[i];
        link_leaf = linkingLengths[i]*linkingLengths[i];
#ifdef SUPERDEBUG
        if (i==0){
            printf("%d\t ",(dists2[i]<=(link_node*link_node)));
            printf("%d\n",((launch_leaf >= launch_node) && launch_leaf<=(launch_node+cool_node)));
            printf("%.2f %.2f %.2f \n",launch_leaf,launch_node,cool_node);
            printf("-----\n");
            printf("%d\t",(dists2[i]<=(link_leaf*link_leaf)));
            printf("%d\n",((launch_node >= launch_leaf) && launch_node<=(launch_leaf+cool_leaf)));
            printf("%.2f %.2f %.2f \n",launch_node,launch_leaf,cool_leaf);
        }
#endif
        if (
                //check if the leaf is contained in the node's hot bubble
            (   (dists2[i]<=(link_node*link_node)) && 
                //check if the leaf is launched after the node, but before it cools
                ((launch_leaf >= launch_node) && launch_leaf<=(launch_node+cool_node))
            ) || 
            (
                //check if the node is contained in the leaf's hot bubble
                (dists2[i]<=(link_leaf*link_leaf)) && 
                //check if the node is launched after the leaf, but before it cools
                ((launch_node >= launch_leaf) && launch_node<=(launch_leaf+cool_leaf))
            )
            ){
            ngbflags[i]=1;
            (*numNGB)++;
        }
    }
}

void calculateDists(float * point, float * xs, float * ys, float * zs,int Narr, float * dists2){
    for (int i =0; i<Narr; i++){
        float dx = point[0]-xs[i];
        float dy = point[1]-ys[i];
        float dz = point[2]-zs[i];
        dists2[i]=dx*dx + dy*dy + dz*dz;
    }
}

int fillFlags(
    float * point, 
    float * xs, float * ys,float * zs,
    float * ids,
    float * launchTimes, float launch_node,
    float * coolingTimes, float cool_node,
    float * linkingLengths, float link_node,
    int Narr, int * NGBFlags){

    int numNGB=0;
    float * dists2 = (float *) malloc(Narr*sizeof(float));
    // calculate the distance to each point
    calculateDists(point,xs,ys,zs,Narr,dists2);

#ifdef SUPERDEBUG
    //printArray(ids,Narr);
    printArray(dists2,Narr);
    printArray(coolingTimes,Narr);
    printArray(launchTimes,Narr);
    //printArrayRatio(dists2,linkingLengths,Narr);
    //printArrayFixedRatio(dists2,link_node,Narr);
    //printArray(point,3);
    //printArray(xs,Narr);
    //printArray(ys,Narr);
    //printArray(zs,Narr);
#endif
    findNGBFlags(
        Narr,
        dists2,
        launchTimes,launch_node,
        coolingTimes,cool_node,
        linkingLengths,link_node,
        NGBFlags,&numNGB);
#ifdef DEBUG
    printf("Found %d friends\n",numNGB);
#endif
    free(dists2);
    return numNGB;
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
            } // for (int j=numNGB; j>0; j--)
        } // if (checkNGBIndices)
    } // for (int i=0; i<numNGB; i++)
}// void extractSubarrayWithIndices

struct SupernovaCluster * findSNeFriends(
    float * xs, float * ys, float * zs, 
    float * launchTimes, float * coolingTimes, float * linkingLengths,
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
        ids,
        launchTimes,launchTimes[0],
        coolingTimes,coolingTimes[0],
        linkingLengths,linkingLengths[0],
        Narr,NGBFlags);

    // convert flags to indices
    NGBIndices=(int*)malloc(numNGB*sizeof(int));
    getIndicesFromFlags(NGBFlags,Narr,NGBIndices);
    
    // "worst case scenario", all the particles are in this cluster
    float * buffer_xs=(float*)malloc(Narr*sizeof(float));
    memset(buffer_xs,0,Narr*sizeof(int));
    float * buffer_ys=(float*)malloc(Narr*sizeof(float));
    memset(buffer_ys,0,Narr*sizeof(int));
    float * buffer_zs=(float*)malloc(Narr*sizeof(float));
    memset(buffer_zs,0,Narr*sizeof(int));
    float * buffer_ids=(float*)malloc(Narr*sizeof(float));
    memset(buffer_ids,0,Narr*sizeof(int));
    float * buffer_launchTimes=(float*)malloc(Narr*sizeof(float));
    memset(buffer_launchTimes,0,Narr*sizeof(int));
    float * buffer_coolingTimes=(float*)malloc(Narr*sizeof(float));
    memset(buffer_coolingTimes,0,Narr*sizeof(int));
    float * buffer_linkingLengths=(float*)malloc(Narr*sizeof(float));
    memset(buffer_linkingLengths,0,Narr*sizeof(int));

    // extract the sub arrays from their indices
        // only need to recalculate NGB indices on the first pass 
        // and could in principal have a separate function that does this
        // but I think this is easier to wrap one's head around
        // essentially it's just if NGBFlags[j] && NGBFlags[N-1-j] -> NGBIndices[N-1-j]=j


    extractSubarrayWithIndices(xs,buffer_xs,NGBIndices,Narr,numNGB,1);
    extractSubarrayWithIndices(ys,buffer_ys,NGBIndices,Narr,numNGB,0);
    extractSubarrayWithIndices(zs,buffer_zs,NGBIndices,Narr,numNGB,0);
    extractSubarrayWithIndices(ids,buffer_ids,NGBIndices,Narr,numNGB,0);
    extractSubarrayWithIndices(launchTimes,buffer_launchTimes,NGBIndices,Narr,numNGB,0);
    extractSubarrayWithIndices(coolingTimes,buffer_coolingTimes,NGBIndices,Narr,numNGB,0);
    extractSubarrayWithIndices(linkingLengths,buffer_linkingLengths,NGBIndices,Narr,numNGB,0);

    // new size of xs/ys/zs/...
    int Nremain=Narr-numNGB;
    int Nadded,numNewNGB;
    int cur_ngb=1; // don't need to check the first neighbor, we just did that above
    while (cur_ngb < numNGB){
#ifdef SUPERDEBUG
        //printf("Current cluster composition:\t");
        //printArray(buffer_ids,numNGB);
#endif
        // add a dot for every loop, haha
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
            fillFlags(
                point,
                xs,ys,zs,
                ids,
                launchTimes,buffer_launchTimes[j],
                coolingTimes,buffer_coolingTimes[j],
                linkingLengths,buffer_linkingLengths[j],
                Nremain,NGBFlags);
        }
        
        // find out how many we actually added
        for (int k=0; k<Nremain; k++) Nadded+=NGBFlags[k];
        
        // set cur_ngb to the end of the list
        cur_ngb=numNGB;

        if (Nadded){
#ifdef DEBUG
            printf("--------------\n");
            printf("Adding the friends we found\n");
#endif
            //reset the NGB indices 
            free(NGBIndices);
            NGBIndices=(int*)malloc(Nadded*sizeof(int));
            memset(NGBIndices,0,Nadded*sizeof(int));

            // find new NGB indices
            getIndicesFromFlags(NGBFlags,Nremain,NGBIndices);

            // extract all our new friends into the buffer, removing them from the main array
            // increment the pointer to point to the first open slot of the buffer
            extractSubarrayWithIndices(
                xs,&buffer_xs[numNGB],
                NGBIndices,Nremain,Nadded,1);
            extractSubarrayWithIndices(
                ys,&buffer_ys[numNGB],
                NGBIndices,Nremain,Nadded,0);
            extractSubarrayWithIndices(
                zs,&buffer_zs[numNGB],
                NGBIndices,Nremain,Nadded,0);
            extractSubarrayWithIndices(
                ids,&buffer_ids[numNGB],
                NGBIndices,Nremain,Nadded,0);
            extractSubarrayWithIndices(
                &launchTimes[numNGB],buffer_launchTimes,
                NGBIndices,Nremain,Nadded,0);
            extractSubarrayWithIndices(
                &coolingTimes[numNGB],buffer_coolingTimes,
                NGBIndices,Nremain,Nadded,0);
            extractSubarrayWithIndices(
                &linkingLengths[numNGB],buffer_linkingLengths,
                NGBIndices,Nremain,Nadded,0);


            // update the size of the neighbors that live in the buffer
            // and the remaining particles
            numNGB+=Nadded;
            Nremain-=Nadded;
#ifdef DEBUG
            printf("Now looking for the new friends' neighbors!\n");
            printf("--------------\n");
#endif
        } //if Nadded
    } //while cur_ngb < numNGB
#ifdef DEBUG
    printf(" finished.\n");


    printf("Building the SN cluster...");
#endif
    struct SupernovaCluster *new_cluster = malloc(sizeof(struct SupernovaCluster)); 

    // set the number of neighbors in the cluster
    new_cluster->cluster_id=cluster_id;
    new_cluster->numNGB = numNGB;

    // allocate the cluster arrays
    new_cluster->xs=(float*)malloc(numNGB*sizeof(float));
    new_cluster->ys=(float*)malloc(numNGB*sizeof(float));
    new_cluster->zs=(float*)malloc(numNGB*sizeof(float));
    new_cluster->ids=(float*)malloc(numNGB*sizeof(float));
    new_cluster->launchTimes=(float*)malloc(numNGB*sizeof(float));
    new_cluster->coolingTimes=(float*)malloc(numNGB*sizeof(float));
    new_cluster->linkingLengths=(float*)malloc(numNGB*sizeof(float));
    
    //memcopy from the buffer to the new_cluster array
    memcpy((void *)new_cluster->xs,(void *)buffer_xs, numNGB*sizeof(float));
    memcpy((void *)new_cluster->ys,(void *)buffer_ys, numNGB*sizeof(float));
    memcpy((void *)new_cluster->zs,(void *)buffer_zs, numNGB*sizeof(float));
    memcpy((void *)new_cluster->ids,(void *)buffer_ids, numNGB*sizeof(float));
    memcpy((void *)new_cluster->launchTimes,(void *)buffer_launchTimes, numNGB*sizeof(float));
    memcpy((void *)new_cluster->coolingTimes,(void *)buffer_coolingTimes, numNGB*sizeof(float));
    memcpy((void *)new_cluster->linkingLengths,(void *)buffer_linkingLengths, numNGB*sizeof(float));

    // free all those buffers
    free(buffer_xs);
    free(buffer_ys);
    free(buffer_zs);
    free(buffer_ids);
    free(buffer_launchTimes);
    free(buffer_coolingTimes);
    free(buffer_linkingLengths);

    free(NGBFlags); 

    return new_cluster;


}// void findSNeFriends

int FoFSNeNGB(
    int Narr,
    float * xs, float * ys, float * zs,
    float * launchTimes, float * coolingTimes,
    float * linkingLengths,
    float * ids,
    struct SupernovaCluster * head,
    int H_OUT ){

#ifdef SUPERDEBUG
    printArray(ids,Narr);
    printArray(xs,Narr);
    printArray(coolingTimes,Narr);
    printArray(launchTimes,Narr);
    return 0;
#endif


    int returnVal;

    int cluster_id=0;

    while (Narr > 0){
        // take our point to be the first SNe in the list
        if (cluster_id % 25 == 0){
            printf("Working on cluster %d\n",cluster_id);
        }


#ifdef PRETTYPRINT
        printf("------------------------------------------------------\n");
        printf("Working on cluster %d\n",cluster_id);
        printf("------------------------------------------------------\n");
#endif
        struct SupernovaCluster *new_cluster = findSNeFriends(
            xs,ys,zs,
            launchTimes,coolingTimes,
            linkingLengths,
            ids,
            Narr,
            cluster_id);

        Narr-=new_cluster->numNGB;

#ifdef PRETTYPRINT
        //printArray(new_cluster->ids,new_cluster->numNGB);
        printf("%d members found\n",new_cluster->numNGB);
        printf("%d elements remain\n",Narr);
#endif

        cluster_id++;
        //step the linked list
        head->NextCluster=new_cluster;
        head = new_cluster;

    } // while (Narr > 0)
    H_OUT=cluster_id;

    return cluster_id;
}



/* -----------------------------------------------------*/
struct GMCClump{
    float * xs;
    float * ys;
    float * zs;
    float * vxs;
    float * vys;
    float * vzs;
    float * masses; 
    float * sfrs;
    float * nH;
    float * ids;
    int numNGB;
    int cluster_id;
    struct GMCClump * NextCluster; 
};

void * findGMCNGBFlags(
    int Narr,
    float * dists2,
    float link_length2,
    int * ngbflags, int * numNGB){
    for (int i=0; i<Narr;i++){
        //check if the leaf is contained in the node's linking length sphere
        if ((dists2[i]<=(link_length2))){
            ngbflags[i]=1;
            (*numNGB)++;
        }
    }
}

int fillGMCFlags(
    float * point, 
    float * xs, float * ys,float * zs,
    float * ids,
    float link_length2,
    int Narr, int * NGBFlags){

    int numNGB=0;
    float * dists2 = (float *) malloc(Narr*sizeof(float));
    // calculate the distance to each point
    calculateDists(point,xs,ys,zs,Narr,dists2);

    findGMCNGBFlags(
        Narr,
        dists2,
        link_length2,
        NGBFlags,&numNGB);

    free(dists2);

    return numNGB;
}
struct GMCClump * findGMCFriends(
    float linkingLength,
    float * xs, float * ys, float * zs,  // coordinates
    float * vxs, float * vys, float * vzs, // velocities
    float * masses, float * sfrs, float * nH, // physical quantities
    float * ids, // particle ids
    int Narr, // length of array
    int cluster_id // this cluster ID
    ){
    
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

    numNGB=fillGMCFlags(
        point,
        xs,ys,zs,
        ids,
        linkingLength*linkingLength,
        Narr,NGBFlags);

    // convert flags to indices
    NGBIndices=(int*)malloc(numNGB*sizeof(int));
    getIndicesFromFlags(NGBFlags,Narr,NGBIndices);
    
    // "worst case scenario", all the particles are in this cluster
    // --- coordinates
    float * buffer_xs=(float*)malloc(Narr*sizeof(float));
    memset(buffer_xs,0,Narr*sizeof(int));
    float * buffer_ys=(float*)malloc(Narr*sizeof(float));
    memset(buffer_ys,0,Narr*sizeof(int));
    float * buffer_zs=(float*)malloc(Narr*sizeof(float));
    memset(buffer_zs,0,Narr*sizeof(int));

    // --- velocities
    float * buffer_vxs=(float*)malloc(Narr*sizeof(float));
    memset(buffer_vxs,0,Narr*sizeof(int));
    float * buffer_vys=(float*)malloc(Narr*sizeof(float));
    memset(buffer_vys,0,Narr*sizeof(int));
    float * buffer_vzs=(float*)malloc(Narr*sizeof(float));
    memset(buffer_vzs,0,Narr*sizeof(int));

    // --- physical quantities
    float * buffer_masses=(float*)malloc(Narr*sizeof(float));
    memset(buffer_masses,0,Narr*sizeof(int));
    float * buffer_sfrs=(float*)malloc(Narr*sizeof(float));
    memset(buffer_sfrs,0,Narr*sizeof(int));
    float * buffer_nH=(float*)malloc(Narr*sizeof(float));
    memset(buffer_nH,0,Narr*sizeof(int));

    // --- etc.
    float * buffer_ids=(float*)malloc(Narr*sizeof(float));
    memset(buffer_ids,0,Narr*sizeof(int));

    // extract the sub arrays from their indices
        // only need to recalculate NGB indices on the first pass 
        // and could in principal have a separate function that does this
        // but I think this is easier to wrap one's head around
        // essentially it's just if NGBFlags[j] && NGBFlags[N-1-j] -> NGBIndices[N-1-j]=j


    // --- coordinates
    extractSubarrayWithIndices(xs,buffer_xs,NGBIndices,Narr,numNGB,1);
    extractSubarrayWithIndices(ys,buffer_ys,NGBIndices,Narr,numNGB,0);
    extractSubarrayWithIndices(zs,buffer_zs,NGBIndices,Narr,numNGB,0);

    // --- velocities
    extractSubarrayWithIndices(vxs,buffer_vxs,NGBIndices,Narr,numNGB,0);
    extractSubarrayWithIndices(vys,buffer_vys,NGBIndices,Narr,numNGB,0);
    extractSubarrayWithIndices(vzs,buffer_vzs,NGBIndices,Narr,numNGB,0);
    
    // --- physical quantities
    extractSubarrayWithIndices(masses,buffer_masses,NGBIndices,Narr,numNGB,0);
    extractSubarrayWithIndices(sfrs,buffer_sfrs,NGBIndices,Narr,numNGB,0);
    extractSubarrayWithIndices(nH,buffer_nH,NGBIndices,Narr,numNGB,0);
    
    // --- etc.
    extractSubarrayWithIndices(ids,buffer_ids,NGBIndices,Narr,numNGB,0);

    // new size of xs/ys/zs/...
    int Nremain=Narr-numNGB;
    int Nadded,numNewNGB;
    int cur_ngb=1; // don't need to check the first neighbor, we just did that above
    while (cur_ngb < numNGB){
#ifdef SUPERDEBUG
        //printf("Current cluster composition:\t");
        //printArray(buffer_ids,numNGB);
#endif
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
            fillGMCFlags(
                point,
                xs,ys,zs,
                ids,
                linkingLength*linkingLength,
                Narr,NGBFlags);
        }
        
        // find out how many we actually added
        for (int k=0; k<Nremain; k++) Nadded+=NGBFlags[k];
        
        // set cur_ngb to the end of the list
        cur_ngb=numNGB;

        if (Nadded){
#ifdef DEBUG
            printf("--------------\n");
            printf("Adding the friends we found\n");
#endif
            //reset the NGB indices 
            free(NGBIndices);
            NGBIndices=(int*)malloc(Nadded*sizeof(int));
            memset(NGBIndices,0,Nadded*sizeof(int));

            // find new NGB indices
            getIndicesFromFlags(NGBFlags,Nremain,NGBIndices);

            // extract all our new friends into the buffer, removing them from the main array
            // increment the pointer to point to the first open slot of the buffer
            // --- coordinates
            extractSubarrayWithIndices(
                xs,&buffer_xs[numNGB],
                NGBIndices,Nremain,Nadded,1);
            extractSubarrayWithIndices(
                ys,&buffer_ys[numNGB],
                NGBIndices,Nremain,Nadded,0);
            extractSubarrayWithIndices(
                zs,&buffer_zs[numNGB],
                NGBIndices,Nremain,Nadded,0);
            
            // --- velocities
            extractSubarrayWithIndices(
                vxs,&buffer_vxs[numNGB],
                NGBIndices,Nremain,Nadded,0);
            extractSubarrayWithIndices(
                vys,&buffer_vys[numNGB],
                NGBIndices,Nremain,Nadded,0);
            extractSubarrayWithIndices(
                vzs,&buffer_vzs[numNGB],
                NGBIndices,Nremain,Nadded,0);

            // --- physical quantities
            extractSubarrayWithIndices(
                masses,&buffer_masses[numNGB],
                NGBIndices,Nremain,Nadded,0);
            extractSubarrayWithIndices(
                sfrs,&buffer_sfrs[numNGB],
                NGBIndices,Nremain,Nadded,0);
            extractSubarrayWithIndices(
                nH,&buffer_nH[numNGB],
                NGBIndices,Nremain,Nadded,0);

            // --- etc.
            extractSubarrayWithIndices(
                ids,&buffer_ids[numNGB],
                NGBIndices,Nremain,Nadded,0);

            // update the size of the neighbors that live in the buffer
            // and the remaining particles
            numNGB+=Nadded;
            Nremain-=Nadded;
#ifdef DEBUG
            printf("Now looking for the new friends' neighbors!\n");
            printf("--------------\n");
#endif
        } //if Nadded
    } //while cur_ngb < numNGB
#ifdef DEBUG
    printf(" finished.\n");


    printf("Building the SN cluster...");
#endif
    struct GMCClump *new_cluster = malloc(sizeof(struct GMCClump)); 

    // set the number of neighbors in the cluster
    new_cluster->cluster_id=cluster_id;
    new_cluster->numNGB = numNGB;

    // allocate the cluster arrays
    // --- coordinates
    new_cluster->xs=(float*)malloc(numNGB*sizeof(float));
    new_cluster->ys=(float*)malloc(numNGB*sizeof(float));
    new_cluster->zs=(float*)malloc(numNGB*sizeof(float));

    // --- velocities
    new_cluster->vxs=(float*)malloc(numNGB*sizeof(float));
    new_cluster->vys=(float*)malloc(numNGB*sizeof(float));
    new_cluster->vzs=(float*)malloc(numNGB*sizeof(float));

    // --- physical quantities
    new_cluster->masses=(float*)malloc(numNGB*sizeof(float));
    new_cluster->sfrs=(float*)malloc(numNGB*sizeof(float));
    new_cluster->nH=(float*)malloc(numNGB*sizeof(float));

    // --- etc.
    new_cluster->ids=(float*)malloc(numNGB*sizeof(float));

    
    //memcopy from the buffer to the new_cluster array
    // --- coordinates
    memcpy((void *)new_cluster->xs,(void *)buffer_xs, numNGB*sizeof(float));
    memcpy((void *)new_cluster->ys,(void *)buffer_ys, numNGB*sizeof(float));
    memcpy((void *)new_cluster->zs,(void *)buffer_zs, numNGB*sizeof(float));
    // --- velocities
    memcpy((void *)new_cluster->vxs,(void *)buffer_vxs, numNGB*sizeof(float));
    memcpy((void *)new_cluster->vys,(void *)buffer_vys, numNGB*sizeof(float));
    memcpy((void *)new_cluster->vzs,(void *)buffer_vzs, numNGB*sizeof(float));
    // --- physical quantities
    memcpy((void *)new_cluster->masses,(void *)buffer_masses, numNGB*sizeof(float));
    memcpy((void *)new_cluster->sfrs,(void *)buffer_sfrs, numNGB*sizeof(float));
    memcpy((void *)new_cluster->nH,(void *)buffer_nH, numNGB*sizeof(float));
    // --- etc.
    memcpy((void *)new_cluster->ids,(void *)buffer_ids, numNGB*sizeof(float));


    // free all those buffers
    free(buffer_xs);
    free(buffer_ys);
    free(buffer_zs);

    free(buffer_vxs);
    free(buffer_vys);
    free(buffer_vzs);

    free(buffer_masses);
    free(buffer_sfrs);
    free(buffer_nH);

    free(buffer_ids);


    free(NGBFlags); 

    return new_cluster;


}// void findGMCs

int FoFGMCNGB(
    int Narr,
    float * xs, float * ys, float * zs,
    float * vxs, float * vys, float * vzs,
    float * masses, float * sfrs, float * nH,
    float linkingLength,
    float * ids,
    struct GMCClump * head,
    int H_OUT ){

#ifdef SUPERDEBUG
    printArray(ids,Narr);
    printArray(xs,Narr);
    printArray(coolingTimes,Narr);
    printArray(launchTimes,Narr);
    return 0;
#endif


    int returnVal;

    int cluster_id=0;

    while (Narr > 0){
        // take our point to be the first SNe in the list

#ifdef PRETTYPRINT
        printf("------------------------------------------------------\n");
        printf("Working on cluster %d\n",cluster_id);
        printf("------------------------------------------------------\n");
#endif
        struct GMCClump *new_cluster = findGMCFriends(
            linkingLength,
            xs,ys,zs,
            vxs,vys,vzs,
            masses,sfrs,nH,
            ids,
            Narr,
            cluster_id);

        Narr-=new_cluster->numNGB;

#ifdef PRETTYPRINT
        //printArray(new_cluster->ids,new_cluster->numNGB);
        printf("%d members found\n",new_cluster->numNGB);
        printf("%d elements remain\n",Narr);
#endif

        cluster_id++;
        //step the linked list
        head->NextCluster=new_cluster;
        head = new_cluster;

    } // while (Narr > 0)
    H_OUT=cluster_id;

    return cluster_id;
}


/* -----------------------------------------------------*/
