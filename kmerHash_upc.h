#ifndef KMER_HASH_H
#define KMER_HASH_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <upc.h>

#include "commonDefaults_upc.h"

/* Creates a hash table and (pre)allocates memory for the memory heap */
hash_table_t* createHashTable(int64_t nEntries, memory_heap_t *memory_heap) {
  hash_table_t *result;
  int64_t n_buckets = nEntries * LOAD_FACTOR;
  
  result = (hash_table_t*) upc_alloc(sizeof(hash_table_t));
  result->size = n_buckets;
  result->table = (bucket_t*) upc_all_alloc(n_buckets , sizeof(bucket_t)); //need to initialize to 0
  upc_memset(result->table, 0, n_buckets * sizeof(bucket_t));
  
  if (result->table == NULL) {
    fprintf(stderr, "ERROR: Could not allocate memory for the hash table: %lld buckets of %lu bytes\n", n_buckets, sizeof(bucket_t));
    fprintf(stderr, "ERROR: Are you sure that your input is of the correct KMER_LENGTH in Makefile?\n");
    exit(1);
  }
  
  memory_heap->heap = (kmer_t *) upc_alloc(nEntries * sizeof(kmer_t));
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
