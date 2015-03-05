#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <string>
#include <vector>

extern "C" {
#include "dataio.h"
#include "statistics.h"
#include "sequences.h"
#include "readio.h"
}
#include "interval_tree.h"

#define min(a,b) (a<b?a:b)
#define max(a,b) (a>b?a:b)

#define NO 0
#define TOTAL 1
#define NBED 2

using namespace std;

class GenInt : public Interval
{
public:
    int i;			// start position
    int j;			// end position
    string c;		// chromosome
    string id;		// id
    float score;	// score
    int strand;		// strand
    string str;
    int matches;
    int matchesplus;
    int num;
    int GetLowPoint() const
    {
        return i;
    }
    int GetHighPoint() const
    {
        return j;
    }
};

vector<GenInt> a_int;

int main(int argc, char** argv) {
    char bamfile[] = "/Volumes/Hani_MyBook_Th/seqfiles/09242013_small-RNA_MDA_LM2s_normal_hypoxic/S1_CGATGT_sorted.bam" ;
    char outfile[] = "/Volumes/Hani_MyBook_Th/seqfiles/09242013_small-RNA_MDA_LM2s_normal_hypoxic/S1_CGATGT_sorted.txt" ;
    char bedfile[] = "/Users/hani/Life/Projects/Databases/hg19_beds/hg19_miR.bed" ;
    char chrdata[] = "/Users/hani/Life/Projects/Applications/ChIPseeqer2.0/data/hg19.chrdata" ;

    // interval tree
    // store intervals
    IntervalTree* chr_tree;
    vector<GenInt> a_int2;
    GenInt tmpInt;

    // read iterators
    int   numreads = 0;
    MappedRead     r;
    ReadI          ri;
    char* formattext;
    int format = BAMUNQ;

    // chrdata
    char** chrnames;
    int* chrlens;
    int numchroms;
    HASH hc; //chromosome index

    readChrData(chrdata, &chrnames, &chrlens, &numchroms);
    HASH_create(&hc, numchroms);
    for (int i=0; i<numchroms; i++) {
        HASH_enter(&hc, chrnames[i], i);
    }

    int nbeds = nbLinesInFile(bedfile);

    HASH hm; // hash to store miRNAs
    HASH_create(&hm, nbeds);

    int midx = -1;
    FILE* fp = fopen(bedfile, "r");

    int  mynmax = 10000;
    char* buff = 0;
    buff  = (char*)malloc(mynmax * sizeof(char));
    int m;
    char** a;
    char chr[1000];
    int  cidx = -1;
    int numint = 0;

    // norm
    int normalize       = NBED;
    char* normalize_txt = 0;
    int normto = 1000000;

    chr_tree = new IntervalTree[numchroms];

    while (fgets(buff, mynmax, fp) != 0) {
        chomp(buff);
        if (buff[0] == '#') {
            continue;
        }
        split_line_delim(buff, (char*)"\t", &a, &m);
        // chr
        sprintf(chr, "%s", a[0]);

        // chr idx
        cidx = -1;
        HASH_find(&hc, chr, &cidx);
        if (cidx == -1) {
            printf("Cannot find %s\n", chr);
            exit(1);
        }

        //printf("%s\t%d\n", chr, cidx);

        if (HASH_find(&hm, a[3], &midx)) { //check if the name already exists
            continue;
        }

        HASH_enter(&hm, a[3], numint);
        //printf ("%s\t%s\t%s\t%s\n", a[0], a[1], a[2], a[3]) ;
        a_int.push_back(tmpInt);

        a_int[numint].id = string(a[3]);

        a_int[numint].c = string(a[0]);
        a_int[numint].i = atoi(a[1]);
        a_int[numint].j = atoi(a[2]);
        a_int[numint].matches = 0;
        a_int[numint].matchesplus = 0;
        a_int[numint].num = numint;

        chr_tree[cidx].Insert(new GenInt(a_int[numint]));
        numint++;
        //printf ("%s\t%d\t%d\t%s\n", (a_int[numint-1].c).c_str(), a_int[numint-1].i, a_int[numint-1].j, (a_int[numint-1].id).c_str()) ;
        free(a) ;
    }
    fclose(fp);

    if (!ReadIopen(&ri, bamfile, format)) {
        printf("Cannot open file %s\n", bamfile);
    }

    ri.verbose = 0;
    int nbedmatches = 0;
    numreads  = 0;

    while (nextRead(&ri, &r)) {
        int x1 = r.pos ;
        int x2 = r.pos + r.lenread;

        // chr idx
        cidx = -1;
        HASH_find(&hc, r.seqname, &cidx);
        if (cidx == -1) {
            printf("Cannot find %s\n", chr);
            exit(1);
        }

        TemplateStack<GenInt*> *res = (TemplateStack<GenInt*>*)chr_tree[cidx].Enumerate(x1, x2);

        /* Number of overlapping intervals */
        int mycnt		= res->Size();
        if (mycnt > 0) {
            //printf("Found %d overlapping miRNAs: %s\n", mycnt, (*res)[0]->id.c_str());

            int ov = sequencesOverlap(x1,x2,a_int[ (*res)[0]->num ].i, a_int[ (*res)[0]->num ].j);
            if (ov < 15) {
                continue;
            }

            a_int[ (*res)[0]->num ].matches ++;

            if (r.st == 1) {
                a_int[ (*res)[0]->num ].matchesplus ++;
            }

            nbedmatches ++;
        }
        numreads ++;

        if (numreads % 10000 == 0) {
            fprintf(stderr, "# %d reads read       \r", numreads);
        }

    }

    ReadIclose(&ri);

    FILE* fpo = fopen(outfile, "w");
    if (!fpo) {
        die("cabbot open outfile\n");
    }

    for (int i=0; i<numint; i++) {
        fprintf(fpo, "%s",  a_int[i].id.c_str());
        if (normalize == NO) {
            fprintf(fpo, "\t%d",  a_int[i].matches);
        } else if (normalize == TOTAL) {
            fprintf(fpo, "\t%3.2f",  normto * (double)(a_int[i].matches)/numreads);
        } else if (normalize == NBED) {
            fprintf(fpo, "\t%3.2f",  normto * (double)(a_int[i].matches)/nbedmatches);
        }
        fprintf(fpo, "\n");
    }

    return 0;
}