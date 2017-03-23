#ifndef KMER_HASH_UPC_H
#define KMER_HASH_UPC_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <upc.h>

#include "commonDefaults_upc.h"

/* Creates a hash table and (pre)allocates memory for the memory heap *
hash_table_t* createHashTable(int64_t nEntries, memory_heap_t *memory_heap) {
  hash_table_t *result;
  int64_t n_buckets = nEntries * LOAD_FACTOR;
  
  result = (hash_table_t*) malloc(sizeof(hash_table_t));
  result->size = n_buckets;
  result->table = (bucket_t*) malloc(n_buckets , sizeof(bucket_t)); //need to initialize to 0
  memset(result->table, 0, n_buckets * sizeof(bucket_t));
  
  if (result->table == NULL) {
    fprintf(stderr, "ERROR: Could not allocate memory for the hash table: %lld buckets of %lu bytes\n", n_buckets, sizeof(bucket_t));
    fprintf(stderr, "ERROR: Are you sure that your input is of the correct KMER_LENGTH in Makefile?\n");
    exit(1);
  }
  
  memory_heap->heap = (kmer_t *) malloc(nEntries * sizeof(kmer_t));
  if (memory_heap->heap == NULL) {
    fprintf(stderr, "ERROR: Could not allocate memory for the heap!\n");
    exit(1);
  }
  memory_heap->posInHeap = 0;
  
  return result;
}

/* Auxiliary function for computing hash values */
int64_t hashSeq(int64_t  hashtable_size, char *seq, int size) {
  unsigned long hashval;
  hashval = 5381;
  for(int i = 0; i < size; i++) {
    hashval = seq[i] +  (hashval << 5) + hashval;
  }
  
  return (hashval % hashtable_size);
}

/* Returns the hash value of a kmer */
int64_t hashKmer(int64_t  hashtable_size, char *seq) {
  return hashSeq(hashtable_size, seq, KMER_PACKED_LENGTH);
}

/* Looks up a kmer in the hash table and returns a pointer to that entry */
shared kmer_t* lookupKmer(shared hash_table_t *hashtable, const unsigned char *kmer) {
  
  char packedKmer[KMER_PACKED_LENGTH];
  packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
  int64_t hashval = hashKmer(hashtable->size, (char*) packedKmer);
  bucket_t currBucket;
  shared kmer_t *result;
  
  currBucket = hashtable->table[hashval];
  result = currBucket.head;
  
  while(result != NULL) {
    if (memcmp(packedKmer, result->kmer, KMER_PACKED_LENGTH * sizeof(char)) == 0) {
      return result;
    }
    result = result->next;
  }
  return NULL;
  
}

/* Adds a kmer and its extensions in the hash table (note that memory heap must be preallocated!) *
int addKmer(hash_table_t *hashtable, memory_heap_t *memory_heap, const unsigned char *kmer, char left_ext, char right_ext) {
  
  // Pack a k-mer sequence appropriately 
  char packedKmer[KMER_PACKED_LENGTH];
  packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
  int64_t hashval = hashKmer(hashtable->size, (char*) packedKmer);
  int64_t pos = memory_heap->posInHeap;
  
  // Add the contents to the appropriate kmer struct in the heap
  upc_memcpy((memory_heap->heap[pos]).kmer, packedKmer, KMER_PACKED_LENGTH * sizeof(char));
  (memory_heap->heap[pos]).l_ext = left_ext;
  (memory_heap->heap[pos]).r_ext = right_ext;
  
  // Fix the next pointer to point to the appropriate kmer struct
  (memory_heap->heap[pos]).next = hashtable->table[hashval].head;
  // Fix the head pointer of the appropriate bucket to point to the current kmer
  hashtable->table[hashval].head = &(memory_heap->heap[pos]);
  
  // Increase the heap pointer
  memory_heap->posInHeap++;
  
  return 0;
  
}

/* Adds a k-mer in the start list by using the memory heap (note that the k-mer was "just added" in the memory heap at position posInHeap - 1) *
void addKmerToStartList(memory_heap_t *memory_heap, start_kmer_t **startKmersList) {
  start_kmer_t *new_entry;
  kmer_t *ptrToKmer;
  
  int64_t prevPosInHeap = memory_heap->posInHeap - 1;
  ptrToKmer = &(memory_heap->heap[prevPosInHeap]);
  new_entry = (start_kmer_t*) malloc(sizeof(start_kmer_t));
  new_entry->next = (*startKmersList);
  new_entry->kmerPtr = ptrToKmer;
  (*startKmersList) = new_entry;
}

/* Deallocate heap. Call before calling deallocHashtable */
int deallocHeap(memory_heap_t *memory_heap) {
  free(memory_heap->heap);
  return 0;
}

/** Deallocate hashtable */
int deallocHashtable(hash_table_t *hashtable) {
  free(hashtable->table);
  return 0;
}

#endif // KMER_HASH_UPC_H
