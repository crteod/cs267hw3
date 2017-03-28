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
  
  if (workBuffer == NULL) {
    fprintf(stderr, "ERROR: Could not allocate memory for workBuffer: %lu bytes\n", charsToRead * sizeof(unsigned char));
    upc_global_exit(1);
  }
  
  FILE * inputFile = fopen(inputUFXName, "r");
  int64_t offset = MYTHREAD * kmersPerThread * LINE_SIZE;
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
  shared [1] int64_t *localPartialArraySizes = upc_all_alloc(THREADS, sizeof(int64_t));
  shared [] int64_t *rootArraySizes = upc_all_alloc(1, THREADS * sizeof(int64_t));
  int64_t *localArraySizes = malloc(THREADS * sizeof(int64_t));
  
  if ((localPartialArraySizes == NULL) || (rootArraySizes == NULL) || (localArraySizes == NULL)) {
    fprintf(stderr, "ERROR: Could not allocate memory for localPartialArraySizes or rootArraySizes or localArraySizes \n");
    upc_global_exit(1);
  }
  
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
    int64_t kmerIndex = addKmer(hashtable, &memoryHeap, &workBuffer[ptr], leftExt, rightExt);
    
    /* Create also a list with the "start" kmers: nodes with F as left (backward) extension */
    if (leftExt == 'F') {
      addKmerToStartList(&memoryHeap, &startKmersList, kmerIndex);
      localPartialArraySizes[MYTHREAD]++;
    }
    
    /* Move to the next k-mer in the input workBuffer */
    ptr += LINE_SIZE;
  }
  
  printf("%d added kmers to hash table!\n", MYTHREAD);
  upc_barrier;
  ///////////////////////////////////////////
  
  /* Create local partial start node array from local linked-lists */
  int64_t localArraySize = localPartialArraySizes[MYTHREAD];
  int64_t *localPartialSNArray = malloc(localArraySize * sizeof(int64_t));
  
  if (localPartialSNArray == NULL) {
    fprintf(stderr, "ERROR: Could not allocate memory for the local start array: %lu bytes for thread %d\n", localArraySize * sizeof(int64_t), MYTHREAD);
    upc_global_exit(1);
  }
  
  int currIndex = 0;
  start_kmer_t *currStartLink = startKmersList;
  while (currStartLink != NULL) {
    localPartialSNArray[currIndex] = currStartLink->kmerIndex;
    
    currStartLink = currStartLink->next;
    currIndex++;
  }
  
  printf("%d populated local array of start nodes of size %ld!\n", MYTHREAD, localPartialArraySizes[MYTHREAD]);
  
  // Gather all partial array sizes in one array on root thread, then broadcast a copy to all threads for local access
  upc_all_gather(rootArraySizes, localPartialArraySizes, sizeof(int64_t), UPC_IN_ALLSYNC | UPC_OUT_ALLSYNC);
  upc_memget(localArraySizes, rootArraySizes, THREADS * sizeof(int64_t));
  
  // Calculate total number of start nodes and offsets for each thread into a concatenated array of start nodes
  int totalStartNodes = 0;
  int rootStartNodeArrayOffset = 0;  
  for (int i = 0; i < THREADS; ++i)
  {
    totalStartNodes += localArraySizes[i];
    if (i < MYTHREAD)
      rootStartNodeArrayOffset += localArraySizes[i];
  }
  
  // TODO: Remove after testing gather/broadcast of size arrays vs direct remote access
  /*for (int i = 0; i < MYTHREAD; ++i) {
    rootStartNodeArrayOffset += localPartialArraySizes[i];
  }*/
  //totalStartNodes = bupc_allv_reduce_all(int64_t, localArraySize, UPC_ADD);
  
  int64_t nbytesTotalStartNodes = totalStartNodes * sizeof(int64_t);
  
  printf("%d knows total start nodes = %ld!\n", MYTHREAD, totalStartNodes);
  
  /* Gather all local partial arrays of start nodes into a global array of start kmers on the root thread */
  shared [] int64_t * rootStartNodeArray = upc_all_alloc(1, nbytesTotalStartNodes);
  if (rootStartNodeArray == NULL) {
    fprintf(stderr, "ERROR: Could not allocate memory for the root start array: %lu bytes\n", nbytesTotalStartNodes);
    upc_global_exit(1);
  }
  
  upc_memput(&rootStartNodeArray[rootStartNodeArrayOffset], localPartialSNArray, localArraySize * sizeof(int64_t));
  
  // Local array that will contain a copy of the entire set of start nodes gathered on the root thread
  int64_t *localStartNodeArray = malloc(nbytesTotalStartNodes);
  
  if (localStartNodeArray == NULL) {
    fprintf(stderr, "ERROR: Could not allocate memory for localStartNodeArray: %lu bytes\n", nbytesTotalStartNodes);
    upc_global_exit(1);
  }
  
  /* Broadcast global array of start kmers to all threads */
  upc_barrier;
  upc_memget(localStartNodeArray, rootStartNodeArray, nbytesTotalStartNodes);
  
  printf("%d finishing a2a\n", MYTHREAD);
  upc_barrier;
  constrTime += gettime();
  ///////////////////////////////////////////
  
  /** Graph traversal **/
  traversalTime -= gettime();
  
  char localOutFilename[20];
  sprintf(localOutFilename, "output/pgen-%d.out", MYTHREAD);
  FILE * myOutputFile = fopen(localOutFilename, "w");
  
  /* Pick start nodes from the startNodesGlobal */
  shared int64_t *currSNIndex = upc_all_alloc(1, sizeof(int64_t));
    
  if (currSNIndex == NULL) {
    fprintf(stderr, "ERROR: Could not allocate memory for currSNIndex\n");
    upc_global_exit(1);
  }
  
  int64_t localSNIndex = 0;
  int64_t localContigs = 0;
  int64_t localBases = 0;
  char unpackedKmer[KMER_LENGTH+1];
  char currContig[MAXIMUM_CONTIG_SIZE];
  if (MYTHREAD == ROOT)
    *currSNIndex = 0;
  
  upc_barrier;
  
  unpackedKmer[KMER_LENGTH] = '\0';
  
  printf("%d finished initialization\n", MYTHREAD);
  
  // Synchronization
  kmer_t currKmerPtr;

  while((localSNIndex = bupc_atomicI64_fetchadd_strict((shared void*)currSNIndex, (int64_t) 1)) < totalStartNodes) {  
    
    /* Unpack first seed and initialize contig */
    int64_t heapIndex = localStartNodeArray[localSNIndex];
    upc_memget(&currKmerPtr, &memoryHeap.heap[heapIndex], sizeof(kmer_t));
    unpackSequence((unsigned char*) currKmerPtr.kmer, (unsigned char*) unpackedKmer, KMER_LENGTH);
    memcpy(currContig, unpackedKmer, KMER_LENGTH * sizeof(char));
    
    int64_t posInContig = KMER_LENGTH;
    rightExt = currKmerPtr.rExt;
    
    /* Keep adding bases until we find a terminal node */
    while (rightExt != 'F') {
      currContig[posInContig] = rightExt;
      posInContig++;
      
      /* The last kmer in the current contig is at position currContig[posInContig-KMER_LENGTH] */
      int lookupFailed = lookupKmer(hashtable, &memoryHeap, &currKmerPtr, (const unsigned char *) &currContig[posInContig-KMER_LENGTH]);
      if (lookupFailed) {
	fprintf(stderr, "ERROR: Lookup failed on thread=%d!\n", MYTHREAD);
	upc_global_exit(1);
      }
      rightExt = currKmerPtr.rExt;
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
  int64_t totalBases = bupc_allv_reduce(int64_t, localContigs, ROOT, UPC_ADD);
  int64_t totalContigs = bupc_allv_reduce(int64_t, localBases, ROOT, UPC_ADD);
  
  printf("%d finished reductions\n", MYTHREAD);
  
  /** CLEAN UP */
  free(workBuffer);
  free(localPartialSNArray);
  
  deallocHeap(&memoryHeap);
  deallocHashtable(hashtable);
  
  /***** DO NOT CHANGE THIS PART ****/
  if(MYTHREAD == ROOT){
    printf("%s: Input set: %s\n", argv[0], argv[1]);
    printf("Number of UPC threads: %d\n", THREADS);
    printf("Input reading time: %f seconds\n", inputTime);
    printf("Graph construction time: %f seconds\n", constrTime);
    printf("Graph traversal time: %f seconds\n", traversalTime);
    printf("Total time: %f seconds\n", inputTime + constrTime + traversalTime);
    printf("Generated %ld contigs with %ld total bases\n", totalContigs, totalBases);
  }
  
  return 0;
  
}
