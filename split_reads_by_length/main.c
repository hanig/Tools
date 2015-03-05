#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include "mystrio.h"

KSEQ_INIT(gzFile, gzread) // STEP 1: declare the type of file handler and the read() function

int main(int argc, char *argv[]){
    int i, j ;

    gzFile fp;
    gzFile out1, out2;
    kseq_t *seq;

    int t1=15, t2 = 25 ; //thresholding

    int l;
    char *fastqfile=strdup(argv[1]) ;
    printf("%s\n", fastqfile) ;

    //build outfile names
    char outfile1[1000], outfile2[1000] ;
    char *root ;
    root = replace_str(fastqfile,".fastq.gz", "") ;
    sprintf(outfile1, "%s_st%i.fastq.gz", root, t2) ;
    sprintf(outfile2, "%s_lt%i.fastq.gz", root, t2) ;

    fp = gzopen(fastqfile, "rb"); // STEP 2: open the file handler
    out1 = gzopen(outfile1, "wb");
    out2 = gzopen(outfile2, "wb");

    seq = kseq_init(fp); // STEP 3: initialize seq

    int total=0 ;
    int cnt=0 ;
    while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
        total++ ;
        if (total % 1000000 == 0){
            printf(".") ;
            fflush(stdout) ;
        }
        if (strlen(seq->seq.s)<=t2 && strlen(seq->seq.s)>t1) {
            cnt++ ;
            char *tmpseq = (char*) calloc(10000,sizeof(char)) ;
            sprintf(tmpseq, "@%s\n%s\n+\n%s\n", seq->name.s, seq->seq.s, seq->qual.s) ;
            gzwrite(out1, (char*)tmpseq, (unsigned)strlen(tmpseq)) ;
            free(tmpseq) ;
        }else if(strlen(seq->seq.s)>t2){
            char *tmpseq = (char*) calloc(10000,sizeof(char)) ;
            sprintf(tmpseq, "@%s\n%s\n+\n%s\n", seq->name.s, seq->seq.s, seq->qual.s) ;
            gzwrite(out2, (char*)tmpseq, (unsigned)strlen(tmpseq)) ;
            free(tmpseq) ;
        }
    }
    printf("\n%d/%d\n", cnt, total) ;
    kseq_destroy(seq); // STEP 5: destroy seq
    gzclose(fp) ; // STEP 6: close the file handler
    gzclose(out1) ;
    gzclose(out2) ;
    return 0;
}