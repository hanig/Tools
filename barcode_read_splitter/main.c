#include <zlib.h>
#include <stdio.h>

#include "kseq.h"
#include "mystrio.h"
#include "hashtable.h"

KSEQ_INIT(gzFile, gzread) // STEP 1: declare the type of file handler and the read() function

int main(int argc, char *argv[]){
    int i ;
    int bclen=8 ;

    gzFile fp;
    kseq_t *seq;
    int l;

    char*** barcodes ;
    int nbarcodes ;
    int ncols ;

    char *barcodefile=strdup(argv[1]) ;
    char *fastqfile=strdup(argv[2]) ;
    char *outdir=strdup(argv[3]) ;

    readStringTable(barcodefile, &barcodes, &nbarcodes, &ncols) ;

    HASH    h_barcodes ;
    HASH_init(&h_barcodes, 1000) ;

    gzFile *outfiles = (gzFile*) malloc(nbarcodes*sizeof(gzFile)) ;
    char  **outfilenames = (char**) malloc (nbarcodes*sizeof(char*)) ;

    for (i=0 ; i<nbarcodes ; i++){
        char *str = strdup(barcodes[i][0]) ;
        outfilenames[i] = (char*) malloc (1000*sizeof(char)) ;
        //sprintf(outfilenames[i], "%s/%s-%s.fastq.gz", outdir, barcodes[i][1], str) ;
        sprintf(outfilenames[i], "%s/%s.fastq.gz", outdir, barcodes[i][1]) ;
        outfiles[i] = gzopen(outfilenames[i], "wb") ;
        HASH_enter(&h_barcodes, str, i);
    }

    fp = gzopen(fastqfile, "rb"); // STEP 2: open the file handler
    seq = kseq_init(fp); // STEP 3: initialize seq
    int total=0 ;
    int *cnt=(int*) calloc(nbarcodes,sizeof(int)) ;
    while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
        total++ ;
        char *name = strdup(seq->name.s) ;
        char *read = strdup(seq->seq.s) ;
        char *bc = strndup (read, bclen) ;

        int index ;
        if (HASH_find(&h_barcodes, bc, &index)){
            cnt[index]++ ;
            char *nread = strndup (read+bclen, strlen(read)-bclen) ;
            char *nqual = strndup (seq->qual.s+bclen, strlen(seq->qual.s)-bclen) ;

            char *tmpseq = (char*) calloc(10000,sizeof(char)) ;
            sprintf(tmpseq, "@%s#1\n%s\n+\n%s\n", name, nread, nqual) ;
            gzwrite(outfiles[index], (char*)tmpseq, (unsigned)strlen(tmpseq)) ;
            free(tmpseq) ;
            if (total%1000000==0){
                printf(".") ;
                fflush(stdout) ;
            }

            free(nread) ;
            free(nqual) ;
        }
        free(name) ;
        free(read) ;
    }

    for (i=0 ; i<nbarcodes ; i++){
        gzclose(outfiles[i]);
        printf("%s\t%d\n", barcodes[i][0], cnt[i]) ;
    }

    kseq_destroy(seq); // STEP 5: destroy seq
    gzclose(fp); // STEP 6: close the file handler
    HASH_destroy(&h_barcodes) ;
    return 0;
}