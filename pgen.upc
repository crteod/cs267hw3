#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>
//#include <upc_atomic.h>
#include <upc_collective.h>

#include "packingDNAseq.h"
#include "kmerHash_upc.h"
#include "commonDefaults_upc.h"


int main(int argc, char *argv[]) {
  
  /** Declarations **/
  double inputTime=0.0, constrTime=0.0, traversalTime=0.0;
  const int ROOT = 0;
  shared hash_table_t * hashtable;
  
  /** Read input **/
  upc_barrier;
  inputTime -= gettime();
  ///////////////////////////////////////////
  // Your code for input file reading here //
  ///////////////////////////////////////////
  upc_barrier;
  inputTime += gettime();
  
  /** Graph construction **/
  constrTime -= gettime();
  ///////////////////////////////////////////
  // Your code for graph construction here //
  ///////////////////////////////////////////
  upc_barrier;
  constrTime += gettime();
  
  ////////////////////////////////////////////////////////////
  // TODO: after a2a, I need startNodesGlobal to be the globally defined array
  shared start_kmer_t *startNodesGlobal;
  // TODO: after a2a, I need totalStartNodes to be the size of startNodes
  int64_t totalStartNodes;
  ////////////////////////////////////////////////////////////
  
  /** Graph traversal **/
  traversalTime -= gettime();
  
  // TODO: Save your output to "output/pgen.out"
  char localOutFilename[20];
  sprintf(localOutFilename, "output/pgen%d.out", MYTHREAD);
  FILE * myOutputFile = fopen(localOutFilename, "w");
  
  /* Pick start nodes from the startNodesGlobal */
  shared int64_t *currSNIndex;
  int64_t localSNIndex = 0;
  upc_lock_t *indexLock;
  shared kmer_t *currKmerPtr;
  shared int64_t *localBases;
  shared int64_t *localContigs;
  char unpackedKmer[KMER_LENGTH+1];
  char currContig[MAXIMUM_CONTIG_SIZE];
  upc_barrier;
  indexLock = upc_all_lock_alloc();
  localBases = upc_all_alloc(THREADS, sizeof(int64_t));
  localContigs = upc_all_alloc(THREADS, sizeof(int64_t));
  upc_barrier;
  if (MYTHREAD == ROOT)
    *currSNIndex = 0;
  localContigs[MYTHREAD] = 0;
  localBases[MYTHREAD] = 0;
  upc_barrier;
  
  // Synchronization
  //while((localSNIndex = bupc_atomicI64_fetchadd_S(currSNIndex, (int64_t) 1LL)) != totalStartNodes-1) {
  while(1) {
    upc_lock(indexLock);
    (*currSNIndex)++;
    localSNIndex = *currSNIndex;
    upc_unlock(indexLock);
    if (localSNIndex >= totalStartNodes) {
      break;
    }
    
    // TODO: confirm on Piazza @203 that this does in fact work
    shared [1] start_kmer_t *currStartNode = &startNodesGlobal[localSNIndex];
    /* Unpack first seed and initialize contig */
    currKmerPtr = currStartNode->kmerPtr; // communication
    unpackSequence((unsigned char*) currKmerPtr->kmer,  (unsigned char*) unpackedKmer, KMER_LENGTH);
    memcpy(currContig, unpackedKmer, KMER_LENGTH * sizeof(char));
    int64_t posInContig = KMER_LENGTH;
    char rightExt = currKmerPtr->r_ext; // communication
    
    /* Keep adding bases until we find a terminal node */
    while (rightExt != 'F') {
      currContig[posInContig] = rightExt;
      posInContig++;
      /* The last kmer in the current contig is at position currContig[posInContig-KMER_LENGTH] */
      currKmerPtr = lookupKmer(hashtable, (const unsigned char *) &currContig[posInContig-KMER_LENGTH]);
      rightExt = currKmerPtr->r_ext; // communication
    }
    
    /* Print the contig to our local file */
    currContig[posInContig] = '\0';
    fprintf(myOutputFile,"%s\n", currContig);
    // TODO: only for debugging
    localContigs[MYTHREAD]++;
    localBases[MYTHREAD]++;
    
  }
  
  upc_barrier;
  traversalTime += gettime();
  
  /** Print timing and output info **/
  // TODO: reduce localContigs and localBases for debugging
  shared int64_t *totalBases;
  shared int64_t *totalContigs;
  upc_barrier;
  upc_all_reduceL(totalContigs, localContigs, UPC_ADD, THREADS, 1, 
		  NULL, UPC_IN_NOSYNC | UPC_OUT_NOSYNC);
  upc_all_reduceL(totalBases, localBases, UPC_ADD, THREADS, 1, 
		  NULL, UPC_IN_NOSYNC | UPC_OUT_NOSYNC);
  upc_barrier;
  upc_all_free(localContigs);
  upc_all_free(localBases);
  upc_lock_free(indexLock);
  /***** DO NOT CHANGE THIS PART ****/
  if(MYTHREAD==ROOT){
    printf("%s: Input set: %s\n", argv[0], argv[1]);
    printf("Number of UPC threads: %d\n", THREADS);
    printf("Input reading time: %f seconds\n", inputTime);
    printf("Graph construction time: %f seconds\n", constrTime);
    printf("Graph traversal time: %f seconds\n", traversalTime);
    printf("Generated %lld contigs with %lld total bases\n", *totalContigs, *totalBases);
  }
  
  return 0;
}
