#include <search.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <search.h>
#include <limits.h>
#include <time.h>
#include <ctype.h>

#include "dataio.h"
#include "sequences.h"

int main(int argc, char ** argv) {
	int       i, j, k ;
	char*     motiffile ;
	char*     outfile ;
	char*     fastafile ;
	int       mynseqmax = 40000 ;
    
	seqI      si ;
	char*     seq ;
	char**    seqs ;
	char*     name;
	int       size;
	int       nseqs = 0 ;
  	char**    seqnames ;
  	int       seqlen = 0 ;
    
	float     gc ;
    
	motiffile   = "/Users/hani/Life/Projects/Ongoing_Projects/tiRNAs/CIMS/3xFLAG-YBX1-LM2tr/YBX1-meme-wm.txt" ; //get_parameter(argc, argv, "-motiffile") ;
	outfile     = "/Users/hani/Desktop/test.txt" ; //get_parameter(argc, argv, "-outfile") ;
	fastafile   = "/Users/hani/Desktop/test.fa" ; //get_parameter(argc, argv, "-fastafile") ;
	gc          = 0.63636364 ; //atof(get_parameter(argc, argv, "-gc")) ;
    
	initialize_nt() ;
    
	float ** wmatrix  ;
	char *motifname ;
	int width ;
    
	readMotifProbMatrixFile (motiffile, &wmatrix, &motifname, &width, gc) ;
	
	printf("Done\n") ;
    
	if (seqI_open(&si, fastafile) == 0) {
		printf("Error opening fastafile: %s\n", fastafile) ;
		exit(0) ;
	}
	printf("\nLoading sequences...") ;
    
	seqs      = (char**) malloc (mynseqmax*sizeof(char*)) ;
	seqnames  = (char**) malloc (mynseqmax*sizeof(char*)) ;
	while ((seq = seqI_nextSequence(&si, &name, &size))) {
	    if (seqlen < strlen(seq))
	    	seqlen = strlen(seq) ;
        
	    int l = strlen(seq) ;
	    for (j=0; j<l; j++) {
            seq[j] = toupper(seq[j]) ;
    	}
    	seqs[nseqs] = strdup(seq) ;
    	seqnames[nseqs] = strdup(name) ;
    	nseqs++ ;
    	//free(seq) ;
    	free(name) ;
    }
    printf("%d sequence...Done\n", nseqs) ;
    
    FILE* f = fopen (outfile, "w") ;
    if (!f){
		printf("Couldn't open file:%s\n", outfile) ;
		exit(0) ;
	}
    
    fprintf(f, "%s\n", motifname) ;
    
    int *stars = (int*) malloc (width*sizeof(int)) ;
    int *matches_pos ;
    int num_matches ;
    char *matches_ori ;
    float *scores ;
    
    for (k=0 ; k<width ; k++)
        stars[k] = 1 ;
    
    float SX = 0 ;
    float SX2 = 0 ;
    for (j=0 ; j<nseqs ; j++){
        findAllWeightMatrixMatches(wmatrix, width, stars, 0, seqs[j], &matches_pos, &num_matches, &matches_ori, &scores, 100, 0) ;
        for(k=0 ; k<num_matches ; k++){
            fprintf(f, "%d\t%f\n", matches_pos[k], scores[k]) ;
        }
    }
    
    fprintf(f, "\n") ;
    fclose(f) ;
    
    printf("DONE\n") ;
	return 0 ;
}
