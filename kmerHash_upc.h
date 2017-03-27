#ifndef KMER_HASH_UPC_H
#define KMER_HASH_UPC_H
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <upc_relaxed.h>
#include "commonDefaults_upc.h"

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
  
  upc_forall(int i = 0; i < nBuckets; ++i; i) {
    result->table[i].head = -1; // Must be set to -1 before addKmer() called
  }
  
  memoryHeap->heap = upc_all_alloc(THREADS, heapBlockSize * sizeof(kmer_t));
  
  if (memoryHeap->heap == NULL) {
    fprintf(stderr, "ERROR: Could not allocate memory for the heap!\n");
    upc_global_exit(1);
  }
  
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
void lookupKmer(hash_table_t *hashtable, memory_heap_t *memoryHeap, kmer_t * result, const unsigned char *kmer) {
  
  char packedKmer[KMER_PACKED_LENGTH];
  packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
  int64_t hashval = hashKmer(hashtable->size, (char*) packedKmer);
  
  // TODO: I'm pretty sure that this shouldn't work?
  shared bucket_t * currBucket = &(hashtable->table[hashval]);
  int64_t currIndex = currBucket->head;
  
  while (currIndex != -1)
  {
    upc_memget(result, &memoryHeap->heap[currIndex], sizeof(kmer_t));
    if (memcmp(packedKmer, (char *)result->kmer, KMER_PACKED_LENGTH * sizeof(char)) == 0) {   
      return;
    }
    currIndex = result->next;
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
  int64_t pos =  memoryHeap->posInHeap * THREADS + MYTHREAD;
  
  // Atomically add kmer to bucket
  int64_t oldHead = hashtable->table[hashval].head;
  int64_t realOldHead = bupc_atomicI64_cswap_strict(&hashtable->table[hashval].head, oldHead, pos);
  while (oldHead != realOldHead)
  {
    oldHead = realOldHead;
    realOldHead = bupc_atomicI64_cswap_strict(&hashtable->table[hashval].head, oldHead, pos);
  }
  
  // TODO: does this match our definitely correct syntax throughout pgen.upc?
  shared kmer_t *indexedKmer = (shared kmer_t *)(&(memoryHeap->heap[pos]));
  kmer_t localKmer;
  
  /* Add the contents to the appropriate kmer struct in the heap */
  memcpy(localKmer.kmer, packedKmer, KMER_PACKED_LENGTH * sizeof(char));
  localKmer.lExt = leftExt;
  localKmer.rExt = rightExt;
  localKmer.next = realOldHead;
  upc_memput(indexedKmer, &localKmer, sizeof(kmer_t));
  
  // Increase the heap pointer
  memoryHeap->posInHeap++;
  
  return pos;
  
}

/* Adds a k-mer in the start list by using the memory heap */
void addKmerToStartList(memory_heap_t *memoryHeap, start_kmer_t **startKmersList, int64_t kmerIndex, int64_t heapBlockSize) {
  
  start_kmer_t *newEntry = (start_kmer_t*) malloc(sizeof(start_kmer_t));
  newEntry->next = (*startKmersList);
  newEntry->kmerIndex = kmerIndex;
  (*startKmersList) = newEntry;
}

/* Deallocate heap. Call before calling deallocHashtable */
int deallocHeap(memory_heap_t *memoryHeap) {
  upc_all_free(memoryHeap->heap);
  return 0;
}

/** Deallocate hashtable */
int deallocHashtable(hash_table_t *hashtable) {
  upc_all_free(hashtable->table);
  free(hashtable);
  return 0;
}

#endif // KMER_HASH_H

