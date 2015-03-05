#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "sam.h"
#include "faidx.h"

void uctransform(char* seq) ;
char* revcmp(char* s) ;
char* my_fai_fetch(faidx_t* idx, char* c, int i, int j, char o);

int main(int argc, char ** argv) {
    char fastafile[] = "/Users/hani/Life/Projects/Applications/bowtie-0.12.7/indexes/hg19.fa" ;

    faidx_t* fidx = 0;

    // load index
    fidx = fai_load(fastafile);
    if (fidx == 0) {
        printf("Cannot load fasta index\n");
        return -1 ;
    }

    char c[] = "chr10" ;
    int be = 132317295 ;
    int en = 132317552 ;
    char st = '-' ;

    char* ss =  my_fai_fetch(fidx, c, be, en, st);
    if (ss == 0) {
        printf("Problem: could not extract %s:%d-%d (%c)\n", c, be, en, st) ;
        exit(-1);
    }
    printf("%s\n", ss) ;

    return 0;
}

char* my_fai_fetch(faidx_t* idx, char* c, int i, int j, char o) {
    int len;
    char id[1000];
    sprintf(id, "%s:%d-%d", c, i, j);
    if(o == '+')
        return fai_fetch(idx, id, &len) ;
    else
        return revcmp(fai_fetch(idx, id, &len)) ;
}

void uctransform(char* seq) {
    int l;
    int i;
    l = strlen(seq);
    for (i=0; i<l; i++) {
        seq[i] = toupper(seq[i]);
    }
}

char* revcmp(char* s) {
    int l = strlen(s);
    char* c;
    char  d;
    int i;

    c = (char*)calloc(l+1, sizeof(char));

    for (i=l-1; i>=0; i--) {
        if (s[i] == 'A')
            d = 'T';
        else if (s[i] == 'T')
            d = 'A';
        else if (s[i] == 'G')
            d = 'C';
        else if (s[i] == 'C')
            d = 'G';
        else if (s[i] == '[')
            d = ']';
        else if (s[i] == ']')
            d = '[';
        else if (s[i] == '.')
            d = '.';
        else
            d = 'N';
        strncat(c, &d, 1);
    }

    return c;

}
