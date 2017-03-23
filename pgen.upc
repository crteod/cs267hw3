#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>
//#include <upc_collective.h>
#include <bupc_collectivev.h>

#include "packingDNAseq.h"
#include "kmerHash_upc.h"
#include "commonDefaults_upc.h"


int main(int argc, char *argv[]) {
  
  /** Declarations **/
  
  // Private variables
  double inputTime=0.0, constrTime=0.0, traversalTime=0.0;
  const int ROOT = 0;
  char leftExt, rightExt;
  start_kmer_t *startKmersList = NULL, *currStartLink;
  
  // Private variables to shared memory
  shared directory_entry_t *directory;
  directory = upc_all_alloc(THREADS, sizeof(directory_entry_t));
  upc_barrier;
  directory[MYTHREAD].size = 0;
  
  /** Read input **/
  upc_barrier;
  inputTime -= gettime();
  
  /* Read the input file name */
  char *inputUFXName = argv[1];
  
  /* Initialize lookup table that will be used for the DNA packing routines */
  initLookupTable();
  
  /* Extract the number of k-mers in the input file */
  int64_t nKmers = getNumKmersInUFX(inputUFXName);
  
  /* Read the kmers from the input file and store them in the workBuffer */
  int64_t totalCharsToRead = nKmers * LINE_SIZE;
  unsigned char * workBuffer = (unsigned char*) malloc(totalCharsToRead * sizeof(unsigned char));
  FILE * inputFile = fopen(inputUFXName, "r");
  int64_t currCharsRead = fread(workBuffer, sizeof(unsigned char),totalCharsToRead , inputFile);
  fclose(inputFile);
        
  ///////////////////////////////////////////
  upc_barrier;
  inputTime += gettime();

  /** Graph construction **/
  constrTime -= gettime();
  
  int64_t heapBlockSize = nKmers / THREADS;
  if ((nKmers % THREADS) != 0)
    heapBlockSize++;
  
  /* Create a hash table */
  memory_heap_t memory_heap;
  // TODO: shared or no?
  shared hash_table_t *hashtable = createHashTable(nKmers, &memory_heap, heapBlockSize);
  
  /* Process the workBuffer and store the k-mers in the hash table */
  /* Expected format: KMER LR ,i.e. first k characters that represent the kmer, then a tab and then two chatacers, one for the left (backward) extension and one for the right (forward) extension */
  int64_t ptr = 0;
  while (ptr < currCharsRead) {
    /* workBuffer[ptr] is the start of the current k-mer                */
    /* so current left extension is at workBuffer[ptr+KMER_LENGTH+1]    */
    /* and current right extension is at workBuffer[ptr+KMER_LENGTH+2]  */
    
    leftExt = (char) workBuffer[ptr+KMER_LENGTH+1];
    rightExt = (char) workBuffer[ptr+KMER_LENGTH+2];
    
    /* Add k-mer to hash table */
    int64_t kmerIndex = addKmer(hashtable, &memory_heap, &workBuffer[ptr], leftExt, rightExt);
    
    /* Create also a list with the "start" kmers: nodes with F as left (backward) extension */
    if (leftExt == 'F') {
      addKmerToStartList(&memory_heap, &startKmersList, kmerIndex);
      directory[MYTHREAD].size++;
    }
    
    /* Move to the next k-mer in the input workBuffer */
    ptr += LINE_SIZE;
  }
  
  int my_array_size = directory[MYTHREAD].size;
  
  // TODO: check if pointer has to be shared
  directory[MYTHREAD].localStartArray = upc_alloc(my_array_size * sizeof(shared kmer_t *shared));
  
  if (directory[MYTHREAD].localStartArray == NULL) {
    fprintf(stderr, "ERROR: Could not allocate memory for the local start array: %lu bytes for thread %d\n", my_array_size * sizeof(shared kmer_t *shared), MYTHREAD);
    exit(1);
  }
  
  int currIndex = 0;
  currStartLink = startKmersList;
  while (currStartLink != NULL) {
    directory[MYTHREAD].localStartArray[currIndex] = currStartLink->kmerPtr;
    
    currStartLink = currStartLink->next;
    currIndex++;
  }
  
  
  // Shared pointer to global list of all thread's master start-node lists
  shared kmer_t * shared * globalStartNodeArray;
  
  // Local pointer to first location (with thread affinity) of each thread's master list
  //shared [1] kmer_t * startNodesGlobal;
  shared kmer_t * shared * startNodesGlobal; // TODO: was this previously
  int64_t totalStartNodes = bupc_allv_reduce(int64_t, my_array_size, 0, UPC_ADD);
  globalStartNodeArray = upc_all_alloc(THREADS, totalStartNodes * sizeof(shared kmer_t *shared));
  
  if (globalStartNodeArray == NULL) {
    fprintf(stderr, "ERROR: Could not allocate memory for the global start array: %lu bytes\n", totalStartNodes * sizeof(shared start_kmer_t *shared));
    exit(1);
  }
  
  startNodesGlobal = &globalStartNodeArray[MYTHREAD * totalStartNodes];
  
  //
  // MEGA TODO: GATHER FULL START NODE LIST AND BROADCAST!!!
  // ALL2ALL COMMUNICATION GOES HERE
  //
  
  upc_barrier;
  constrTime += gettime();
  
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
  unpackedKmer[KMER_LENGTH] = '\0';
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
    
    /* Unpack first seed and initialize contig */
    
    // TODO: confirm on Piazza @203 that shared [1] kmer_t * does in fact work
    shared kmer_t *currKmerPtr = startNodesGlobal[localSNIndex];
    unpackSequence((unsigned char*) currKmerPtr->kmer,  (unsigned char*) unpackedKmer, KMER_LENGTH);
    memcpy(currContig, unpackedKmer, KMER_LENGTH * sizeof(char));
    int64_t posInContig = KMER_LENGTH;
    rightExt = currKmerPtr->r_ext; // communication
    
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
