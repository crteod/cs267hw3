#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>
#include <bupc_collectivev.h>

#include "packingDNAseq.h"
#include "kmerHash_upc.h"
#include "commonDefaults_upc.h"


int main(int argc, char *argv[]) {

	/** Declarations **/
	double inputTime=0.0, constrTime=0.0, traversalTime=0.0;
        
        // Private vars
        char cur_contig[MAXIMUM_CONTIG_SIZE], unpackedKmer[KMER_LENGTH+1], left_ext, right_ext, *inputUFXName;
        int64_t posInContig, ptr = 0, nKmers, cur_chars_read, total_chars_to_read, global_start_array_size = 0;
        unpackedKmer[KMER_LENGTH] = '\0';
        start_kmer_t *startKmersList = NULL, *curStartLink;
        unsigned char *working_buffer;
        FILE *inputFile, *serialOutputFile;
        int64_t heap_block_size;
        
        
        // Private vars to shared memory
        shared int64_t contigID[THREADS];
        shared [1] int64_t totBases[THREADS];
        contigID[MYTHREAD] = 0;
        totBases[MYTHREAD] = 0;
        shared kmer_t *cur_kmer_ptr;
        shared kmer_t *shared * global_start_kmer_array;
        shared kmer_t *shared * my_kmer_array_ptr;
        shared directory_entry_t directory[THREADS];
        directory[MYTHREAD].size = 0;
        
	/** Read input **/
	upc_barrier;
	inputTime -= gettime();
	///////////////////////////////////////////
	// Your code for input file reading here //
        
        /* Read the input file name */
        inputUFXName = argv[1];
        
        /* Initialize lookup table that will be used for the DNA packing routines */
        initLookupTable();
        
        /* Extract the number of k-mers in the input file */
        nKmers = getNumKmersInUFX(inputUFXName);
        
        /* Read the kmers from the input file and store them in the working_buffer */
        total_chars_to_read = nKmers * LINE_SIZE;
        working_buffer = (unsigned char*) malloc(total_chars_to_read * sizeof(unsigned char));
        inputFile = fopen(inputUFXName, "r");
        cur_chars_read = fread(working_buffer, sizeof(unsigned char),total_chars_to_read , inputFile);
        fclose(inputFile);
        
	///////////////////////////////////////////
	upc_barrier;
	inputTime += gettime();

	/** Graph construction **/
	constrTime -= gettime();
	///////////////////////////////////////////
	// Your code for graph construction here //
        
        hash_table_t *hashtable;
        memory_heap_t memory_heap;
        
        heap_block_size = nKmers / THREADS;
        if ((nKmers % THREADS) != 0)
          heap_block_size++;
        
        /* Create a hash table */
        hashtable = createHashTable(nKmers, &memory_heap, heap_block_size);
        
        /* Process the working_buffer and store the k-mers in the hash table */
        /* Expected format: KMER LR ,i.e. first k characters that represent the kmer, then a tab and then two chatacers, one for the left (backward) extension and one for the right (forward) extension */
        
        while (ptr < cur_chars_read) {
          /* working_buffer[ptr] is the start of the current k-mer                */
          /* so current left extension is at working_buffer[ptr+KMER_LENGTH+1]    */
          /* and current right extension is at working_buffer[ptr+KMER_LENGTH+2]  */
          
          left_ext = (char) working_buffer[ptr+KMER_LENGTH+1];
          right_ext = (char) working_buffer[ptr+KMER_LENGTH+2];
          
          /* Add k-mer to hash table */
          int64_t kmer_index = addKmer(hashtable, &memory_heap, &working_buffer[ptr], left_ext, right_ext);
          
          /* Create also a list with the "start" kmers: nodes with F as left (backward) extension */
          if (left_ext == 'F') {
            addKmerToStartList(&memory_heap, &startKmersList, kmer_index);
            directory[MYTHREAD].size++;
          }
          
          /* Move to the next k-mer in the input working_buffer */
          ptr += LINE_SIZE;
        }
        
        int my_array_size = directory[MYTHREAD].size;
        
        // TODO: check if need the shared on the pointer
        directory[MYTHREAD].localStartArray = upc_alloc(my_array_size * sizeof(shared kmer_t *shared));
        
        if (directory[MYTHREAD].localStartArray == NULL) {
          fprintf(stderr, "ERROR: Could not allocate memory for the local start array: %lu bytes for thread %d\n", my_array_size * sizeof(shared kmer_t *shared), MYTHREAD);
          exit(1);
        }
        
        int curIndex = 0;
        curStartLink = startKmersList;
        while (curStartLink != NULL)
        {
          directory[MYTHREAD].localStartArray[curIndex] = curStartLink->kmerPtr;
          
          curStartLink = curStartLink->next;
          curIndex++;
        }
        
        global_start_array_size = bupc_allv_reduce(int64_t, my_array_size, 0, UPC_ADD);
        global_start_kmer_array = upc_all_alloc(THREADS, global_start_array_size * sizeof(shared kmer_t *shared));
        
        if (global_start_kmer_array == NULL) {
          fprintf(stderr, "ERROR: Could not allocate memory for the global start array: %lu bytes\n", global_start_array_size * sizeof(shared kmer_t *shared));
          exit(1);
        }
        
        my_kmer_array_ptr = &global_start_kmer_array[MYTHREAD * global_start_array_size];
        
        //
        // TODO: GATHER FULL START NODE LIST AND BROADCAST!!!
        //
        
        
        ///////////////////////////////////////////
	upc_barrier;
	constrTime += gettime();

	/** Graph traversal **/
	traversalTime -= gettime();
	////////////////////////////////////////////////////////////
	// Your code for graph traversal and output printing here //
	// Save your output to "output/pgen.out"                  //
	////////////////////////////////////////////////////////////
	upc_barrier;
	traversalTime += gettime();

	/** Print timing and output info **/
	/***** DO NOT CHANGE THIS PART ****/
	if(MYTHREAD==0){
		printf("%s: Input set: %s\n", argv[0], argv[1]);
		printf("Number of UPC threads: %d\n", THREADS);
		printf("Input reading time: %f seconds\n", inputTime);
		printf("Graph construction time: %f seconds\n", constrTime);
		printf("Graph traversal time: %f seconds\n", traversalTime);
	}
	return 0;
}
