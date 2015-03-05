#include <search.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "dataio.h"
#include "sequences.h"

#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?(x):(y))

int main(int argc, char ** argv) {
	int       i, j, k ;
	char*     re ;
	char*     outfile ;
	char*     fastafile ;
	int       mynseqmax = 40000 ;
    
	seqI      si ;
	char*     seq = 0;
	char*     name;
	int       size;

	int   singlestrand = 0;
	int   lenoffset = 0;

	int   l ;
	int   np;
	int   no;
	int   flank = 10 ;

	int simple = 0 ; //output binary vector


	re        	= "A[GC]GTGA" ; //get_parameter(argc, argv, "-motiffile") ;
	outfile     = "/Users/hani/Desktop/test.txt" ; //get_parameter(argc, argv, "-outfile") ;
	fastafile   = "/Users/hani/Desktop/test.fa" ; //get_parameter(argc, argv, "-fastafile") ;

	if (seqI_open(&si, fastafile) == 0) {
		printf("Error opening fastafile: %s\n", fastafile) ;
		exit(0) ;
	}
	//printf("\nLoading sequences...") ;

	int*  positions    = (int*)malloc(5000 * sizeof(char));
	int*  orientations = (int*)malloc(5000 * sizeof(char));
	int*  lengths_m    = (int*)malloc(5000 * sizeof(char));

	while ((seq = seqI_nextSequence(&si, &name, &size))) {
		if (seq && (strlen(seq) > 0)) {
			np = 0;
			no = 0;

			l = strlen(seq);
			findSites(re, name, seq, positions, &np, orientations, &no, lengths_m, singlestrand, lenoffset);
			if (simple == 1) {
				printf("%s\t%d\n", name, (np > 0 ? 1 : 0));
			} else {
				for (i = 0; i < np; i++) {
					printf("%s\t", name);
					printf("%d\t%d", positions[i], orientations[i]);
					fflush(stdout);
					int sp = max(0, positions[i] - flank);
					int sl = positions[i] - sp;
					char *ss = substr(seq, sp, sl);
					for (j = 0; j < flank; j++)
						ss[j] = tolower(ss[j]);
					printf("\t%s", ss);
					free(ss);

					fflush(stdout);

					ss = substr(seq, positions[i], lengths_m[i]);
					printf("%s", ss);
					free(ss);
					fflush(stdout);
					sp = positions[i] + lengths_m[i];

					sl = min(l, sp + flank);
					//printf("\nlen=%d, sl=%d, flank_d=%d\n", l, sl, sp+flank);

					ss = substr(seq, sp, sl - sp);
					for (j = 0; j < flank; j++)
						ss[j] = tolower(ss[j]);
					printf("%s", ss);
					fflush(stdout);
					free(ss);
					printf("\n");
				}
			}
		}
		free(seq);
    }

	return 0 ;
}
