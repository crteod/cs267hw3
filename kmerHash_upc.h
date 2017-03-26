#ifndef KMER_HASH_UPC_H
#define KMER_HASH_UPC_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <upc_relaxed.h>

#include "commonDefaults_upc.h"

int64_t local2Global(int64_t local) {
  return (local*THREADS+MYTHREAD);
}

/* Creates a hash table and (pre)allocates memory for the memory heap */
hash_table_t* createHashTable(int64_t nEntries, memory_heap_t *memoryHeap, int64_t heapBlockSize) {
  hash_table_t *result;
  int64_t nBuckets = nEntries * LOAD_FACTOR;
  
  result = malloc(sizeof(hash_table_t));
  result->size = nBuckets;
  // TODO: check that sizeof shared bucket_t doesn't change things
  result->table = upc_all_alloc(nBuckets, sizeof(bucket_t));
  
  if (result->table == NULL) {
    fprintf(stderr, "ERROR: Could not allocate memory for the hash table: %ld buckets of %lu bytes\n", nBuckets, sizeof(bucket_t));
    fprintf(stderr, "ERROR: Are you sure that your input is of the correct KMER_LENGTH in Makefile?\n");
    upc_global_exit(1);
  }
  
  upc_lock_t *lock;
  upc_forall(int i = 0; i < nBuckets; ++i; i) {
    if ((lock = upc_global_lock_alloc()) == NULL) {
      fprintf(stderr, "ERROR: Failed to alloc lock\n");
      exit(1);
    }
    result->table[i].bucketLock = lock;
    result->table[i].head = NULL; // Must be set to 0 before addKmer() called
    result->table[i].size = 0;
  }
  
  memoryHeap->heap = upc_all_alloc(THREADS, heapBlockSize * sizeof(kmer_t));
  //memoryHeap->heap = upc_alloc(heapBlockSize * sizeof(kmer_t));
  
  if (memoryHeap->heap == NULL) {
    fprintf(stderr, "ERROR: Could not allocate memory for the heap!\n");
    upc_global_exit(1);
  }
  //memoryHeap->posInHeap = MYTHREAD * heapBlockSize;
  memoryHeap->posInHeap = 0;
  
  return result;
}

/* Auxiliary function for computing hash values */
int64_t hashSeq(int64_t  hashSize, char *seq, int size) {
  unsigned long hashval;
  hashval = 5381;
  for(int i = 0; i < size; i++) {
    hashval = seq[i] +  (hashval << 5) + hashval;
  }
  
  return (hashval % hashSize);
}

/* Returns the hash value of a kmer */
int64_t hashKmer(int64_t  hashtable_size, char *seq) {
  return hashSeq(hashtable_size, seq, KMER_PACKED_LENGTH);
}

/* Looks up a kmer in the hash table and returns a pointer to that entry */
void lookupKmer(hash_table_t *hashtable, kmer_t * result, const unsigned char *kmer) {
  
  char packedKmer[KMER_PACKED_LENGTH];
  packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
  int64_t hashval = hashKmer(hashtable->size, (char*) packedKmer);
  
  // TODO: I'm pretty sure that this shouldn't work?
  shared bucket_t * currBucket = &(hashtable->table[hashval]);
  upc_memget(result, currBucket->head, sizeof(kmer_t));
  
  while(result != NULL) {
    if (memcmp(packedKmer, (char *)result->kmer, KMER_PACKED_LENGTH * sizeof(char)) == 0) {
      return;
    }
    
    upc_memget(result, result->next, sizeof(kmer_t));
  }
  
  fprintf(stderr, "ERROR: %d failed lookup on contig, returning NULL \n", MYTHREAD);
  result = NULL;
  return;
  
}

/* Adds a kmer and its extensions in the hash table (note that memory heap must be preallocated!) */
int64_t addKmer(hash_table_t *hashtable, memory_heap_t *memoryHeap, const unsigned char *kmer, char leftExt, char rightExt, int64_t heapBlockSize) {
  
  /* Pack a k-mer sequence appropriately */
  char packedKmer[KMER_PACKED_LENGTH];
  packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
  int64_t hashval = hashKmer(hashtable->size, (char*) packedKmer);
  int64_t pos = local2Global(memoryHeap->posInHeap);
  
  // TODO: does this match our definitely correct syntax throughout pgen.upc?
  shared kmer_t *indexedKmer = (shared kmer_t *)(&(memoryHeap->heap[pos]));
  
  /* Add the contents to the appropriate kmer struct in the heap */
  //memcpy((memoryHeap->heap[pos]).kmer, packedKmer, KMER_PACKED_LENGTH * sizeof(char));
  //upc_memput((memoryHeap->heap[pos]).kmer, packedKmer, KMER_PACKED_LENGTH * sizeof(char));
  //(memoryHeap->heap[pos]).lExt = leftExt;
  //(memoryHeap->heap[pos]).rExt = rightExt;
  upc_memput(indexedKmer->kmer, packedKmer, KMER_PACKED_LENGTH * sizeof(char));
  indexedKmer->lExt = leftExt;
  indexedKmer->rExt = rightExt;
  
  /* Enter critical section serializing updates to bucket list */
  upc_lock((hashtable->table[hashval]).bucketLock);
  // Fix the next pointer to point to the appropriate kmer struct
  //(memoryHeap->heap[pos]).next = hashtable->table[hashval].head;
  indexedKmer->next = hashtable->table[hashval].head;
  // Fix the head pointer of the appropriate bucket to point to the current kmer
  //hashtable->table[hashval].head = (shared kmer_t *) &(memoryHeap->heap[pos]);
  hashtable->table[hashval].head = indexedKmer;
  hashtable->table[hashval].size++;
  // Exit critical section
  upc_unlock((hashtable->table[hashval]).bucketLock);
  
  // Increase the heap pointer
  memoryHeap->posInHeap++;
  
  return pos;
  
}

/* Adds a k-mer in the start list by using the memory heap */
void addKmerToStartList(memory_heap_t *memoryHeap, start_kmer_t **startKmersList, int64_t kmerIndex, int64_t heapBlockSize) {
  
  start_kmer_t *newEntry = (start_kmer_t*) malloc(sizeof(start_kmer_t));
  newEntry->next = (*startKmersList);
  newEntry->kmerIndex = kmerIndex;
  //newEntry->kmerPtr = ptrToKmer;
  (*startKmersList) = newEntry;
}

/* Deallocate heap. Call before calling deallocHashtable */
int deallocHeap(memory_heap_t *memoryHeap) {
  upc_all_free(memoryHeap->heap);
  return 0;
}

/** Deallocate hashtable */
int deallocHashtable(hash_table_t *hashtable) {
  upc_forall(int i = 0; i < hashtable->size; ++i; i) {
    upc_lock_free(hashtable->table[i].bucketLock);
  }
  upc_all_free(hashtable->table);
  free(hashtable);
  return 0;
}

#endif // KMER_HASH_H

