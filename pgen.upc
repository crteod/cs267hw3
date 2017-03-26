#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>
#include <upc_collective.h>
#include <bupc_collectivev.h>

#include "packingDNAseq.h"
#include "kmerHash_upc.h"
#include "commonDefaults_upc.h"

int main(int argc, char *argv[]) {
  
  /** Declarations **/
  
  // Private variables
  double inputTime=0.0, constrTime=0.0, traversalTime=0.0;
  char leftExt, rightExt;
  start_kmer_t *startKmersList = NULL;
  
  printf("%d init\n", MYTHREAD);
  
  // Private variables to shared memory
  shared directory_entry_t *directory = upc_all_alloc(THREADS, sizeof(directory_entry_t));
  directory[MYTHREAD].size = 0;
  
  ///////////////////////////////////////////
  /** Read input **/
  upc_barrier;
  inputTime -= gettime();
  
  /* Read the input file name */
  char *inputUFXName = argv[1];
  
  /* Initialize lookup table that will be used for the DNA packing routines */
  initLookupTable();
  
  printf("%d starting to read input file...\n", MYTHREAD);
  
  /* Extract the number of k-mers in the input file */
  int64_t nKmers = getNumKmersInUFX(inputUFXName);
  
  /* Read the kmers from the input file and store them in workBuffer */
  int64_t kmersPerThread = nKmers / THREADS;
  int64_t kmersLeftOver = nKmers - (kmersPerThread*(THREADS-1));
  int64_t charsToRead;
  if (MYTHREAD < THREADS-1) {
    charsToRead = kmersPerThread * LINE_SIZE;
  }
  else {
    charsToRead = kmersLeftOver * LINE_SIZE;
  }
  
  unsigned char * workBuffer = (unsigned char*) malloc(charsToRead * sizeof(unsigned char));
  
  FILE * inputFile = fopen(inputUFXName, "r");
  int64_t offset = MYTHREAD*kmersPerThread*LINE_SIZE;
  fseek(inputFile, offset, SEEK_SET);
  int64_t charsRead = fread(workBuffer, sizeof(unsigned char), charsToRead, inputFile);
  fclose(inputFile);
  if (charsRead != charsToRead) {
    fprintf(stderr, "ERROR: thread %d only read %ld/%ld bytes!\n", MYTHREAD, charsRead, charsToRead);
    upc_global_exit(1);
  }
  
  ///////////////////////////////////////////
  upc_barrier;
  inputTime += gettime();
  
  /** Graph construction **/
  constrTime -= gettime();
  printf("%d finished reading input file!\n", MYTHREAD);
  
  int64_t heapBlockSize = (kmersPerThread > kmersLeftOver ? kmersPerThread : kmersLeftOver);
  
  /* Create a hash table */
  memory_heap_t memoryHeap;
  hash_table_t *hashtable = createHashTable(nKmers, &memoryHeap, heapBlockSize);
  
  printf("%d created hash table!\n", MYTHREAD);
  
  /* Process the workBuffer and store the k-mers in the hash table */
  /* Expected format: KMER LR ,i.e. first k characters that represent the kmer, 
     then a tab and then two characters (one for the left (backward) extension and one for the right (forward) extension) */
  int64_t ptr = 0;
  
  while (ptr < charsRead) {
    /* workBuffer[ptr] is the start of the current k-mer                */
    /* so current left extension is at workBuffer[ptr+KMER_LENGTH+1]    */
    /* and current right extension is at workBuffer[ptr+KMER_LENGTH+2]  */
    
    leftExt = (char) workBuffer[ptr+KMER_LENGTH+1];
    rightExt = (char) workBuffer[ptr+KMER_LENGTH+2];
    
    /* Add k-mer to hash table */
    int64_t kmerIndex = addKmer(hashtable, &memoryHeap, &workBuffer[ptr], leftExt, rightExt, heapBlockSize);
    
    /* Create also a list with the "start" kmers: nodes with F as left (backward) extension */
    if (leftExt == 'F') {
      addKmerToStartList(&memoryHeap, &startKmersList, kmerIndex, heapBlockSize);
      directory[MYTHREAD].size++;
    }
    
    /* Move to the next k-mer in the input workBuffer */
    ptr += LINE_SIZE;
  }
  
  printf("%d added kmers to hash table!\n", MYTHREAD);
  upc_barrier;
  ///////////////////////////////////////////
  
  /* Create local start node array from local linked-lists */
  int64_t localArraySize = directory[MYTHREAD].size;
  int64_t *localStartArray = malloc(localArraySize * sizeof(int64_t));
  
  if (localStartArray == NULL) {
    fprintf(stderr, "ERROR: Could not allocate memory for the local start array: %lu bytes for thread %d\n", localArraySize * sizeof(shared kmer_t *shared), MYTHREAD);
    exit(1);
  }
  
  int currIndex = 0;
  start_kmer_t *currStartLink = startKmersList;
  while (currStartLink != NULL) {
    localStartArray[currIndex] = currStartLink->kmerIndex;
    
    currStartLink = currStartLink->next;
    currIndex++;
  }
  
  printf("%d populated directory of start nodes of size %ld!\n", MYTHREAD, directory[MYTHREAD].size);
  
  int64_t totalStartNodes = bupc_allv_reduce_all(int64_t, localArraySize, UPC_ADD);
  int64_t nbytesTotalStartNodes = totalStartNodes * sizeof(int64_t);
  
  printf("%d knows total start nodes = %ld!\n", MYTHREAD, totalStartNodes);
  
  /* Gather all local arrays into a global array of start kmers on the root thread */
  shared [] int64_t * rootStartNodeArray = upc_all_alloc(1, nbytesTotalStartNodes);
  if (rootStartNodeArray == NULL) {
    fprintf(stderr, "ERROR: Could not allocate memory for the root start array: %lu bytes\n", totalStartNodes * sizeof(shared kmer_t *shared));
    upc_global_exit(1);
  }
  
  int64_t rootStartNodeArrayOffset = 0;
  for (int i = 0; i < MYTHREAD; ++i) {
    // TODO: consider also gathering/broadcasting an array of sizes to minimize remote accesses for sizes
    rootStartNodeArrayOffset += directory[i].size;
  }
  
  upc_memput(&rootStartNodeArray[rootStartNodeArrayOffset], localStartArray, localArraySize * sizeof(int64_t));
  upc_barrier;
  
  // Local pointer to first location (with thread affinity) of each thread's master list
  int64_t *localStartNodeArray = malloc(nbytesTotalStartNodes);
  /* Broadcast global array of start kmers to all threads */
  upc_memget(localStartNodeArray, rootStartNodeArray, nbytesTotalStartNodes);
  
  printf("%d finishing a2a\n", MYTHREAD);
  upc_barrier;
  constrTime += gettime();
  ///////////////////////////////////////////
  
  
  /** Graph traversal **/
  traversalTime -= gettime();
  
  // TODO: Save output to "output/pgen.out"?
  char localOutFilename[20];
  sprintf(localOutFilename, "output/pgen%d.out", MYTHREAD);
  FILE * myOutputFile = fopen(localOutFilename, "w");
  
  /* Pick start nodes from the startNodesGlobal */
  shared int64_t *currSNIndex;
  int64_t localSNIndex = 0;
  int64_t localContigs = 0;
  int64_t localBases = 0;
  char unpackedKmer[KMER_LENGTH+1];
  char currContig[MAXIMUM_CONTIG_SIZE];
  if (MYTHREAD == ROOT)
    *currSNIndex = 0;
  upc_lock_t * indexLock = upc_all_lock_alloc();
  
  upc_barrier;
  
  unpackedKmer[KMER_LENGTH] = '\0';
  upc_barrier;
  printf("%d finished initialization\n", MYTHREAD);
  
  // Synchronization
  kmer_t * currKmerPtr = malloc(sizeof(kmer_t));
  //while((localSNIndex = bupc_atomicI64_fetchadd_S(*currSNIndex, (int64_t) 1LL)) != totalStartNodes-1) {  
  while(1) {
    
    // Lock to make sure only one thread updates currSNIndex
    // TODO: try atomic operation?
    upc_lock(indexLock);
    localSNIndex = *currSNIndex;
    (*currSNIndex)++;
    upc_unlock(indexLock);
    if (localSNIndex >= totalStartNodes)
      break;
    
    /* Unpack first seed and initialize contig */
    int64_t heapIndex = localStartNodeArray[localSNIndex];
    upc_memget(currKmerPtr, &memoryHeap.heap[heapIndex], sizeof(kmer_t));
    unpackSequence((unsigned char*) currKmerPtr->kmer, (unsigned char*) unpackedKmer, KMER_LENGTH);
    memcpy(currContig, unpackedKmer, KMER_LENGTH * sizeof(char));
    
    int64_t posInContig = KMER_LENGTH;
    rightExt = currKmerPtr->rExt;
    
    /* Keep adding bases until we find a terminal node */
    while (rightExt != 'F') {
      currContig[posInContig] = rightExt;
      posInContig++;
      
      /* The last kmer in the current contig is at position currContig[posInContig-KMER_LENGTH] */
      lookupKmer(hashtable, currKmerPtr, (const unsigned char *) &currContig[posInContig-KMER_LENGTH]);
      if (currKmerPtr == NULL) {
	fprintf(stderr, "ERROR: current pointer @%ld for thread=%d is null!\n", posInContig, MYTHREAD);
	exit(1);
      }
      rightExt = currKmerPtr->rExt;
    }
    
    /* Print the contig to our local file */
    currContig[posInContig] = '\0';
    fprintf(myOutputFile,"%s\n", currContig);
    // TODO: only for debugging
    localContigs++;
    localBases++;
    
  }
  
  ///////////////////////////////////////////
  upc_barrier;
  traversalTime += gettime();
  
  printf("%d finished traversal\n", MYTHREAD);
  
  /** Print timing and output info **/
  // TODO: reduce localContigs and localBases only for debugging
  upc_barrier;
  int64_t totalBases = bupc_allv_reduce(int64_t, localContigs, ROOT, UPC_ADD);
  int64_t totalContigs = bupc_allv_reduce(int64_t, localBases, ROOT, UPC_ADD);
  upc_barrier;
  
  printf("%d finished reductions\n", MYTHREAD);
  
  /** CLEAN UP */
  upc_all_free(directory);
  free(workBuffer);
  free(localStartArray);
  
  free(currKmerPtr);
  
  upc_all_lock_free(indexLock);
  
  deallocHeap(&memoryHeap);
  deallocHashtable(hashtable);
  
  /***** DO NOT CHANGE THIS PART ****/
  if(MYTHREAD==ROOT){
    printf("%s: Input set: %s\n", argv[0], argv[1]);
    printf("Number of UPC threads: %d\n", THREADS);
    printf("Input reading time: %f seconds\n", inputTime);
    printf("Graph construction time: %f seconds\n", constrTime);
    printf("Graph traversal time: %f seconds\n", traversalTime);
    printf("Generated %ld contigs with %ld total bases\n", totalContigs, totalBases);
  }
  
  return 0;
  
}
