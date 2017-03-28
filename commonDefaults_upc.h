#ifndef COMMON_DEFAULTS_UPC_H
#define COMMON_DEFAULTS_UPC_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <math.h>
#include <string.h>
#include <upc.h>

/** Initializes many defaults (and structs) that are used throughout the program */
#ifndef MAXIMUM_CONTIG_SIZE
#define MAXIMUM_CONTIG_SIZE 100000
#endif

#ifndef KMER_LENGTH
#define KMER_LENGTH 51
#endif

#ifndef LOAD_FACTOR
#define LOAD_FACTOR 1
#endif

#ifndef LINE_SIZE
#define LINE_SIZE (KMER_LENGTH+4)
#endif

#ifndef ROOT
#define ROOT 0
#endif

/* K-mer data structure */
typedef struct kmer_t kmer_t;
struct kmer_t{
  char kmer[KMER_PACKED_LENGTH];
  char lExt;
  char rExt;
  int64_t next;                 // Index to kmer_t in heap
};

/* Start k-mer data structure */
typedef struct start_kmer_t start_kmer_t;
struct start_kmer_t{
  int64_t kmerIndex;            // Index to kmer_t in heap
  start_kmer_t *next;
};

/* Memory heap data structure */
typedef struct memory_heap_t memory_heap_t;
struct memory_heap_t {
  shared [1] kmer_t *heap;      // Cycled heap array
  int64_t posInHeap;            // Logical thread offset/phase
};

/* Bucket data structure */
typedef struct bucket_t bucket_t;
struct bucket_t{
  int64_t head;                 // Heap index to the first entry of that bucket
};

/* Hash table data structure */
typedef struct hash_table_t hash_table_t;
struct hash_table_t {
  int64_t size;                 // Size of the hash table
  shared bucket_t *table;	// Entries of the hash table buckets
};













/* Returns the number of UFX kmers in a file, and performs some error checking on file */
int64_t getNumKmersInUFX(const char *filename) {
  FILE *f = fopen(filename, "r");
  if (f == NULL) {
    fprintf(stderr, "Could not open %s for reading!\n", filename);
    return -1;
  }
  char firstLine[ LINE_SIZE+1 ];
  firstLine[LINE_SIZE] = '\0';
  if (fread(firstLine, sizeof(char), LINE_SIZE, f) != LINE_SIZE) {
    fprintf(stderr, "Could not read %d bytes!\n", LINE_SIZE);
    return -2;
  }
  // check structure and size of kmer is correct!
  if (firstLine[LINE_SIZE] != '\0') {
    fprintf(stderr, "UFX text file is an unexpected line length for kmer length %d\n", KMER_LENGTH);
    return -3;
  }
  if (firstLine[KMER_LENGTH] != ' ' && firstLine[KMER_LENGTH] != '\t') {
    fprintf(stderr, "Unexpected format for firstLine '%s'\n", firstLine);
    return -4;
  }
  
  struct stat buf;
  int fd = fileno(f);
  if (fstat(fd, &buf) != 0) {
    fprintf(stderr, "Could not fstat %s\n", filename);
    return -5;
  }
  int64_t totalSize = buf.st_size;
  if (totalSize % LINE_SIZE != 0) {
    fprintf(stderr, "UFX file is not a multiple of %d bytes for kmer length %d\n", LINE_SIZE, KMER_LENGTH);
    return -6;
  }
  fclose(f);
  int64_t numKmers = totalSize / LINE_SIZE;
  printf("Detected %ld kmers in text UFX file: %s\n", numKmers, filename);
  return numKmers;
}

/** Utility function to get the current time */
static double gettime(void) {
  struct timeval tv;
  if (gettimeofday(&tv, NULL)) {
    perror("gettimeofday");
    abort();
  }
  return ((double)tv.tv_sec) + tv.tv_usec/1000000.0;
}

#endif // COMMON_DEFAULTS_UPC_H
