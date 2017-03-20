#ifndef PACKING_DNA_SEQ_H
#define PACKING_DNA_SEQ_H

#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>

// Useful utility function
#define pow4(a) (1<<((a)<<1))

// Lookup table to get from packed code (0 < 255) to int value of ACGT string
unsigned int packedCodeToFourMer[256];

/** Initializes packedCodeToFourMer. Currently uses 4-mers so as to keep lookup tables small */
void initLookupTable() {
  
  int merLen = 4, i, slot, valInSlot;
  unsigned char mer[4];
  
  for ( i = 0; i < 256; i++ ) {
    // converts a packedcode to a 4-mer
    int remainder = i;
    int pos = 0;
    for( slot = merLen-1; slot >= 0; slot-- ) {
      valInSlot = remainder / pow4(slot);
      char base;
      
      if( valInSlot == 0 ) {
	base = 'A';
      }
      else if( valInSlot == 1 ) {
	base = 'C';
      }
      else if( valInSlot == 2 ) {
	base = 'G';
      }
      else if( valInSlot == 3 ) {
	base = 'T';
      }
      else{
	fprintf(stderr, "ERROR: Not A,C,G, or T!\n");
	assert(0);
      }
      
      mer[pos] = base;
      pos++;
      remainder -= valInSlot * pow4(slot);
    }
    // update lookup table to reflect this string's int value
    unsigned int *merAsUInt = (unsigned int*) mer;
    packedCodeToFourMer[i] = (unsigned int) (*merAsUInt);
  }
}

/** Returns the character code corresponding to input 4-mer */
unsigned char convertFourMerToPackedCode(const unsigned char *fourMer) {
  int retval = 0;
  int code, i;
  int pow = 64;
  
  for ( i=0; i < 4; i++) {
    char base = fourMer[i];
    switch ( base ) {
    case 'A':
      code = 0;	
      break;
    case 'C':
      code = 1;
      break;
    case 'G':
      code = 2;
      break;
    case 'T':
      code = 3;
      break;
    }
    retval += code * pow;
    pow /= 4;
  }
  return ((unsigned char) retval);
}

/* Updates the pointer to m_data to point to the result of the packing */
void packSequence(const unsigned char *seq_to_pack, unsigned char *m_data, const int m_len) {
  
  
  int ind, j = 0;     // coordinate along unpacked string ( matches with m_len )
  int i = 0;          // coordinate along packed string
  
  // Pack the leading sequence in blocks of 4
  for ( ; j <= m_len - 4; i++, j+=4 ) {
    m_data[i] = convertFourMerToPackedCode( ( unsigned char * ) ( seq_to_pack + j )) ;
  }
  
  // Last 4-block is a special case (if m_len % 4 != 0). Appends "A"s as filler
  int remainder = m_len % 4;
  unsigned char blockSeq[5] = "AAAA";
  for(ind = 0; ind < remainder; ind++) {
    blockSeq[ind] = seq_to_pack[j + ind];
  }
  m_data[i] = convertFourMerToPackedCode(blockSeq);
}

/* Unpacks an input sequence. Updates the pointer unpacked_seq to contain output sequence */
void unpackSequence(const unsigned char *seq_to_unpack, unsigned char *unpacked_seq, const int kmer_len) {
  
  int i = 0, j = 0;
  int packed_len = (kmer_len+3)/4;
  for( ; i < packed_len ; i++, j += 4 ) {
    *( ( unsigned int * )( unpacked_seq + j ) ) = packedCodeToFourMer[ seq_to_unpack[i] ];
  }
  *(unpacked_seq + kmer_len) = '\0';
}

/** Compares two packed sequences */
int comparePackedSeq(const unsigned char *seq1, const unsigned char *seq2, const int seq_len) {
  return memcmp(seq1, seq2, seq_len);
}

#endif // PACKING_DNA_SEQ_H
