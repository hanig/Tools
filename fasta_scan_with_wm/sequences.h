#ifndef SEQUENCES_H_
#define SEQUENCES_H_

int* C;
int* N;

char* getGappedKmer(char* kmer, int gap);
char* getGappedMotif(char* motif, int motifsize, int gap);
int getRegexpMotifLength(char* motif, int gap);
void readMotifProbMatrixFile (char* pfile, float*** wm, char** name, int* w, float gc) ;
float findWeightMatrixMaxScore(float** wml, int w, int* stars, char* seq, int* matches_pos, int complement, int low, int high) ;

void  initialize_nt();
float getWMThreshold(float** wm, int w, int* stars, float* bkg, int log, char** sites, int n);
float getScoreOnSeq(float** wml, int w, int* stars, char* seq, int i, int* strand, int* hasn, int complement);
void readWM(char* wmfile, float*** wm, int* w, char*** sites, int* n, int** stars, int transfo, float* bkgraw);
void findAllWeightMatrixMatches(float** wml, int w, int* stars, float t, char* seq, int** matches_pos, int* num_matches, char** matches_ori, float **scores, int max_num_matches, int complement);
void printWM(float** wm, int w);

#endif /* SEQUENCES_H_ */
