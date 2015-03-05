#define min(x,y) (x<y?x:y)
#define max(x,y) (x>y?x:y)

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sequences.h"
#include "dataio.h"
#include "pcre.h"


void initialize_nt()
{

  N = (int*)calloc(256, sizeof(int));

  N['A'] = 0;
  N['C'] = 1;
  N['T'] = 2;
  N['G'] = 3;
  N['N'] = 4;
  N['X'] = 5;

  C = (int*)calloc(256, sizeof(int));

  C['A'] = 2;
  C['C'] = 3;
  C['T'] = 0;
  C['G'] = 1;
  C['N'] = 4;
  C['X'] = 4;


}


void printWM(float** wm, int w)
{
  int i,j;
  for (j=0; j<w; j++) {

    for (i=0; i<4; i++) {
      printf("\t%3.2f", wm[i][j]);
    }
    printf("\n");
  }


}

//
//  read a single WM from a file
//
//  transfo = 0 -> integer, 1 -> float, 2 -> log
//
void readWM(char* wmfile, float*** wm, int* w, char*** sites, int* n, int** stars, int transfo, float* bkg)
{

  FILE*   fp3;
  char*   buff;
  int     len = 100000;
  float** mywm;
  int     cnt_sites;
  int     l;
  int     col;
  int     myw;
  int     i;
  float   epsilon = 0.38;
  int     star = 0;
  char**  mysites;
  int*    mystars;

  buff = (char*)malloc(1000*sizeof(char));

  //
  //  read in WM
  //
  fp3 = fopen(wmfile, "r");
  //getline(&buff, &len, fp3);
  fgets(buff, len, fp3);

  //getline(&buff, &len, fp3);
  fgets(buff, len, fp3);

  myw = strlen(buff)-1;
  cnt_sites = 1;
  while (!feof(fp3)) {
    //getline(&buff, &len, fp3);
    fgets(buff, len, fp3);
    if (feof(fp3)) {
      break;
    }
    cnt_sites ++;
  }
  rewind(fp3);


  //printf("mw=%d, cnt_sites=%d\n", myw, cnt_sites);

  //
  // input wm, log weights
  //
  mywm = (float**)malloc(myw * sizeof(float*));
  for (col= 0; col < myw; col++) {
    mywm[col] = (float*)malloc(5 * sizeof(float));
    for (i=0; i<5; i++) {
      mywm[col][i] = 0.0;
    }
  }

  mystars = (int*)malloc(myw * sizeof(int));

  //
  //  read in each site
  //
  mysites = (char**)malloc(cnt_sites * sizeof(char*));

  cnt_sites = 0;
  // getline(&buff, &len, fp3);
  fgets(buff, len, fp3);

  while (!feof(fp3)) {
    fgets(buff, len, fp3);
    //getline(&buff, &len, fp3);
    if (feof(fp3)) {
      break;
    }
    chomp(buff);
    l = strlen(buff);

    //printf("read %s, l=%d\n", buff,l);

    star = 0;
    for (i=0; i<l; i++) {
      if ((buff[i] == ' ') || (buff[i] == '*')) {
	star = 1;  break;
      }
    }

    if (star == 1) {
      for (i=0; i<l; i++) {
	if (buff[i] == '*')
	  mystars[i] = 1;
	else
	  mystars[i] = 0;

      }

      break;
    }

    for (i=0; i<l; i++) {
      //printf("buff[i] = %c\n", buff[i]);
      mywm[i][ N[ (int)(buff[i]) ] ] += 1.0;
    }

    mysites[ cnt_sites ] = strdup(buff);

    cnt_sites++;

    //printf("loop end\n");
  }


  /*
  for (col= 0; col < myw; col++) {
     for (i=0; i<5; i++) {
       printf("%5.2f\t", mywm[col][i]);

     }
     printf("\n");
  }
  */

  //
  //  normalize by the number of sites
  //
  for (col= 0; col < myw; col++) {
      for (i=0; i<4; i++) {
	//+= bkg[i]; //epsilon;
	int ini = mywm[col][i];
	mywm[col][i]  = log2( ((mywm[col][i] + bkg[i]) / (float)(cnt_sites + 1)) / bkg[i] );

	//if (mystars[col] == 1)
	  //printf("col = %d, i = %d, ini = %d, si = %d, gc=%f, w = %f\n", col, i, ini, cnt_sites, bkg[i], mywm[col][i]);
	//mywm[col][i]  = log2(mywm[col][i]);
      }
  }



  /*
  printf("\n");
  for (col= 0; col < myw; col++) {
     for (i=0; i<5; i++) {
       printf("%5.2f\t", mywm[col][i]);

     }
     printf("\n");
  }
  */

  //
  //  output
  //
  *w     = myw;
  *wm    = mywm;
  *sites = mysites;
  *n     = cnt_sites;
  *stars = mystars;

}

void readMotifProbMatrixFile (char* pfile, float*** wm, char** name, int* w, float gc) {
  FILE*   fp ;
  char*   buff ;
  int     len = 10000 ;
  int     i,j = 0 ;
  char*   s ;


  buff = (char*)malloc(1000*sizeof(char)) ;

  //open file
  fp = fopen(pfile, "r") ;
  if (!fp){
  	die ("readMotifProbMatrixFile: wrong file name\n") ;
  }
  fgets(buff, len, fp) ;
  chomp(buff) ;
  *name = strdup(buff) ;

  fgets(buff, len, fp) ;
  s = mystrtok(buff, '\t') ; free(s) ;
  *w = 0 ;
  while ((s = mystrtok(0, '\t'))) {
	(*w)++;
	free(s);
  }
  rewind (fp) ;
  fgets(buff, len, fp) ;
  printf("motif name: %s (width = %d)\n", *name, *w) ;

  *wm = (float**) malloc (4*sizeof(float*)) ;
  for (i=0 ; i<4 ; i++){
  	(*wm)[i] = (float*) malloc ((*w)*sizeof(float)) ;
  }

  while (!feof(fp)){
  	fgets(buff, len, fp) ;
  	if (feof(fp))
     break;

  	s = mystrtok(buff, '\t') ;
   	int id = N [(int)(s[0])] ;
   	free(s) ;

   	int pos = 0 ;
   	while ((s = mystrtok(0, '\t'))) {
		(*wm)[id][pos] = atof(s) ;
		free(s);
		pos++ ;
  	}
  }

  float bkg[4] ;
  bkg[0] = (1.0 - gc) / 2.0 ;
  bkg[1] =        gc  / 2.0 ;
  bkg[2] = (1.0 - gc) / 2.0 ;
  bkg[3] =        gc  / 2.0 ;
  for (i=0 ; i<4 ; i++){
  	for (j=0 ; j<*w ; j++){
  		(*wm)[i][j] = log2( (*wm)[i][j] / bkg[i] );
  	}
  }
  free (buff) ;
  fclose (fp) ;
}

//
//  log score
//
float getScoreOnSeq(float** wml, int w, int* stars, char* seq, int i, int* strand, int* hasn, int complement)
{
  float score1 = 0.0;
  float score2 = 0.0;

  int  j, k;

  *hasn = 0;
  for (j=0; j<w; j++) {
    if (stars[j] == 1) {
      if ((seq[i + j] != 'A') && (seq[i + j] != 'T') && (seq[i + j] != 'C') && (seq[i + j] != 'G')) {
		*hasn = 1;
		return 0.0;
      }
        float test = wml[ N[ (int)(seq[i + j]) ] ][j] ;
      score1 += wml[ N[ (int)(seq[i + j]) ] ][j] ;
    }
  }
  if (complement == 1){
  	*strand = 1 ;
    return score1 ;
  }

  for (j=w-1,k=0; j>=0; j--,k++) {
    if (stars[k] == 1) {
      score2 += wml[ C[ (int)(seq[i + j]) ] ][k] ;
    }
  }

  if (complement == -1){
  	*strand = -1 ;
    return score2 ;
  }

  if (score1 >= score2) {
    *strand = 1;
    return score1;
  } else {
    *strand = -1;
    return score2;
  }
}

float findWeightMatrixMaxScore(float** wml, int w, int* stars, char* seq, int* match_pos, int complement, int low, int high){
  int i;
  float s;
  int strand;
  int  hasn;

  float max = -1000 ;
  for (i=low; i<high; i++) {
    s = getScoreOnSeq(wml, w, stars, seq, i, &strand, &hasn, complement) ;
    if (!hasn && s>max){
    	*match_pos = i ;
    	max = s ;
    }
  }
  return max ;
}

void findAllWeightMatrixMatches(float** wml, int w, int* stars, float t, char* seq, int** matches_pos, int* num_matches, char** matches_ori, float **scores, int max_num_matches, int complement)
{

  int i;
  int l = strlen(seq);
  float s;
  int strand;
  int my_num_matches = 0;
  int my_max_num_matches = max_num_matches;
  void* ptr;
  int  hasn;

  *matches_pos = (int*) malloc(my_max_num_matches * sizeof(int));
  *scores = (float*) malloc(my_max_num_matches * sizeof(float));
  if (!*matches_pos) {
    die("findAllWeightMatrixMatches: not enough memory for matches_pos\n");
  }
  //*matches_ori = (char*)malloc(my_max_num_matches * sizeof(char));
  //if (!*matches_ori) {
  //  die("findAllWeightMatrixMatches: not enough memory for matches_ori\n");
  //}



  for (i=0; i<l; i++) {
    s = getScoreOnSeq(wml, w, stars, seq, i, &strand, &hasn, complement);

    if (!hasn && (s >= t)) {

      /*
      if (strand > 0)
	printf("%s => %4.3f (%d)\n", substr(seq, i, w), s, strand);
      else
	printf("%s => %4.3f (%d)\n", complement(substr(seq, i, w)), s, strand);
	//printf("%s => %4.3f (%d)\n", , s, strand);
	*/

      (*matches_pos)[ my_num_matches ] = i;
      (*scores)[ my_num_matches ] = s;
      //(*matches_ori)[ my_num_matches ] = (char)strand;
      my_num_matches ++;

      if (my_num_matches == my_max_num_matches) {
	my_max_num_matches += max_num_matches;
	ptr = realloc( *matches_pos, my_max_num_matches * sizeof(int));
	if (ptr == 0)
	  die("findAllWeightMatrixMatches: not enough memory for matches_pos (realloc)\n");
	//ptr = realloc( *matches_ori, my_max_num_matches * sizeof(char));
	//if (ptr == 0)
	//  die("findAllWeightMatrixMatches: not enough memory for matches_ori (realloc)\n");

      }

    }
  }

  *num_matches = my_num_matches;

}



//
//  get the score threshold out of a WM (log or not)
//
float getWMThreshold(float** wm, int w, int* stars, float* bkg, int transfo, char** sites, int n) {

  float score;
  int   i, j;
  //char* T    = "ACTG";
  float SX   = 0.0;
  float SX2  = 0.0;
  float AVG;
  float STD  = 0.0;
  float t;
  float* scores;
  int   cnt_above;

  scores = (float*)malloc(n * sizeof(float));

  for (i=0; i<n; i++) {

    //printf("site = %s\n", sites[i]);

    if (transfo == 0) {
      score = 1.0;
    } else {
      score = 0.0;
    }

    for (j=0; j<w; j++) {
      if (stars[j] == 0)
      	continue;
      if (transfo == 0) {
	//printf(" %c mult pby %4.3f / %4.3f\n",sites[i][j], wm[j][ N[sites[i][j]] ],
	//     bkg[ N[ sites[i][j] ] ]);

	score *= wm[j][ N[  (int)(sites[i][j]) ] ] / bkg[ N[ (int)(sites[i][j]) ] ];
	//printf(" score = %4.3f\n", score);
      } else {
	score += wm[j][ N[ (int)(sites[i][j]) ] ]; // - bkg[ N[ (int)(sites[i][j]) ] ];
	//printf("j=%d %c mult pby %4.3f (%c) - %4.3f\n", j, sites[i][j], wm[j][ N[sites[i][j]] ], sites[i][j], bkg[ N[ sites[i][j]] ]);

      }
    }

    ///printf("%s\t%4.3f\n", sites[i], score);

    SX  += score;
    SX2 += score * score;

    scores[i] = score;

  }

  AVG =  SX/n;
  STD =  (float)sqrt(( SX2 - SX*SX / n )  / (n - 1));
  //printf("AVG = %3.2f, S = %3.2f\n", AVG, STD);

  t = AVG - 2.0 * STD;

  cnt_above = 0;
  for (i=0; i<n; i++) {
    if (scores[i] >= t) {
      cnt_above ++;
    }
  }

  //printf("%d/%d sites score better than t\n", cnt_above, n);

  return t;

}




char* getGappedKmer(char* kmer, int gap)
{

  char*    stmp;
  int      kmersize;
  int      j;
  int      half1, half2;


  kmersize = strlen(kmer);
  half1 = kmersize/2;
  half2 = kmersize-half1;

  stmp = (char*)calloc((kmersize+gap+1), sizeof(char));

  strncpy(stmp, kmer, half1);
  for (j=0; j<gap; j++) {
    stmp[half1+j] = '.';
  }
  strncpy(stmp+half1+gap, kmer + half1, half2);
  stmp[kmersize+gap] = '\0';


  return stmp;
}


char* getGappedMotif(char* motif, int motifsize, int gap)
{

  char*    newmotif;
  int      half1;
  int      l = strlen(motif);
  int      i = 0;  // cnt in ungapped motif
  int      j = 0;  // cnt in   gapped motif
  int      k = 0;  // cnt characters
  int      h = 0;
  int      p_b;

  newmotif = (char*)calloc((l+gap+1), sizeof(char));

  half1 = motifsize/2;

  while (i < l) {

    if (motif[i] == '[') {

      p_b = i;
      while (motif[p_b] != ']') {
	p_b++;
      }

      while (i <= p_b) {
	newmotif[j] = motif[i];
	i++; j++;
      }

      k++;

    } else {

      newmotif[j] = motif[i];
      j++; i++; k++;

    }

    if (k == half1) {
      for (h=0; h<gap; h++) {
	newmotif[j] = '.';
	j++;
      }
    }

  }

  newmotif[l+gap] = '\0';

  return newmotif;
}



int getRegexpMotifLength(char* motif, int gap)
{


  int      l = strlen(motif);
  int      i = 0;  // cnt in ungapped motif
  int      j = 0;  // cnt in   gapped motif
  int      k = 0;  // cnt characters
  //int      h = 0;
  int      p_b;


  while (i < l) {

    if (motif[i] == '[') {

      p_b = i;
      while (motif[p_b] != ']') {
	p_b++;
      }

      while (i <= p_b) {
	i++; j++;
      }

      k++;
    } else {
      j++; i++; k++;
    }
  }

  return k + gap;
}

/** REGEXP, returns 1 if a site was found, 0 otherwise **/
int findSites(char* r, char* name, char* seq, int* positions, int* np, int* orientations, int* no, int* lengths_m, int singlestrand, int lenoffset)
{

  int i, j; //, j, n;
  int   score1, score2;
  pcre* re;
  const char* error;
  int   erroffset;
  int   rc;
  int   ovector[30];
  char* cr = 0;
  char*  substring_start;
  int  substring_length;
  int startoffset = 0;

  char* seq_pos;
  int pos;

  //
  //  store the position where a RE has been found, to avoid counting twice the palindromes
  //
  seq_pos = (char*)calloc(strlen(seq), sizeof(char));

  *np = 0;
  *no = 0;



  //printf("searching for %s in %s\n", r, seq);


  //
  // ONE STRAND
  //


  re = pcre_compile(r,
          PCRE_CASELESS,
          &error,
          &erroffset,
          NULL);

  while ( (rc = pcre_exec(re,
          NULL,
          seq,
          strlen(seq),
          startoffset,
          0,
          ovector,
          30) ) > 0) {


    substring_start    = seq         + ovector[0];
    substring_length   = ovector[1]  - ovector[0];

    // save the position
    pos = strlen(seq) - ovector[0] - 1;



    //if (reldist == 1) {
    //  positions[ *np ] = pos;
    //} else {
    positions[ *np ] = ovector[0];
    //}

    if (lengths_m != 0)
      lengths_m[ *np ] = substring_length;


    //printf("%s\n", substr(seq, pos, substring_length));
    //printf("%.*s\n", substring_length, substring_start);

    // save the orientation
    orientations[ *no ] = 1;

    if (lenoffset == 0)
      startoffset        = ovector[0] + substring_length;
    else
      startoffset        = ovector[0] + lenoffset;

    seq_pos[ pos ] = 1;



    (*np)++;
    (*no)++;

  }


  //return 0;

  if (singlestrand == 0) {

    // OTHER STRAND

    cr = complement(r);

    //printf("Complement of %s is %s\n", r, cr);

    startoffset = 0;

    re = pcre_compile(cr,
            PCRE_CASELESS,
            &error,
            &erroffset,
            NULL);

    while ( (rc = pcre_exec(re,
            NULL,
            seq,
            strlen(seq),
            startoffset,
            0,
            ovector,
            30)) > 0) {
      substring_start = seq + ovector[0];
      substring_length = ovector[1] - ovector[0];
      //printf("%.*s\n", substring_length, substring_start);
      //  if not a palindromic version
      pos = strlen(seq) - ovector[0] - 1;
      if ((seq_pos[ pos ] != 1) || ((seq_pos[ pos ] == 1) )) {
        // save the position
        //if (reldist == 1) {
        //  positions[ *np ] = pos;
        // } else {
        positions[ *np ] = ovector[0];
        // }

        if (lengths_m != 0)
          lengths_m[ *np ] = substring_length;

        //positions[ *np ] = pos;
        (*np)++;

      }

      // save the orientation in all cases
      orientations[ *no ] = -1;
      (*no)++;

      if (lenoffset == 0)
        startoffset        = ovector[0] + substring_length;
      else
        startoffset        = ovector[0] + lenoffset;

    }

  }


  if (cr != 0)
    free(cr);
  free(seq_pos);
  pcre_free(re);
  return 0;
}

