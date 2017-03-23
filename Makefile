CC = CC
UPCC = upcc

KMER_LENGTH 		= 19
KMER_PACKED_LENGTH 	= $(shell echo $$((($(KMER_LENGTH)+3)/4)))

# Add -std=gnu99 to CFLAGS if use gnu compiler
CFLAGS 	= -O3 
CFLAGSUPC = -O3 -std=gnu99
DEFINE 	= -DKMER_LENGTH=$(KMER_LENGTH) -DKMER_PACKED_LENGTH=$(KMER_PACKED_LENGTH)
HEADERS	= commonDefaults.h kmerHash.h packingDNAseq.h
HEADERSUPC = commonDefaults_upc.h kmerHash_upc.h packingDNAseq.h
LIBS	=

TARGETS	= serial pgen sort

all: 	$(TARGETS)

serial: serial.c $(HEADERS)
		$(CC) $(CFLAGS) -o $@ $< -DKMER_LENGTH=$(KMER_LENGTH) -DKMER_PACKED_LENGTH=$(KMER_PACKED_LENGTH) $(LIBS)

pgen:	pgen.upc $(HEADERSUPC)
		$(UPCC) $(UPCFLAGS) -Wc,"$(CFLAGSUPC)" -o $@ $< $(DEFINE) $(LIBS)

sort:	sort.cpp
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

clean :
	rm -f *.o
	rm -rf $(TARGETS)
