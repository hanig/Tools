#include <stdio.h>
#include "mystrio.h"
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

char *replace_str(const char *str, const char *old, const char *new)
{
    char *ret, *r;
    const char *p, *q;
    size_t oldlen = strlen(old);
    size_t count, retlen, newlen = strlen(new);
    int samesize = (oldlen == newlen);

    if (!samesize) {
        for (count = 0, p = str; (q = strstr(p, old)) != NULL; p = q + oldlen)
            count++;
        /* This is undefined if p - str > PTRDIFF_MAX */
        retlen = p - str + strlen(p) + count * (newlen - oldlen);
    } else
        retlen = strlen(str);

    if ((ret = malloc(retlen + 1)) == NULL)
        return NULL;

    r = ret, p = str;
    while (1) {
        /* If the old and new strings are different lengths - in other
         * words we have already iterated through with strstr above,
         * and thus we know how many times we need to call it - then we
         * can avoid the final (potentially lengthy) call to strstr,
         * which we already know is going to return NULL, by
         * decrementing and checking count.
         */
        if (!samesize && !count--)
            break;
        /* Otherwise i.e. when the old and new strings are the same
         * length, and we don't know how many times to call strstr,
         * we must check for a NULL return here (we check it in any
         * event, to avoid further conditions, and because there's
         * no harm done with the check even when the old and new
         * strings are different lengths).
         */
        if ((q = strstr(p, old)) == NULL)
            break;
        /* This is undefined if q - p > PTRDIFF_MAX */
        ptrdiff_t l = q - p;
        memcpy(r, p, l);
        r += l;
        memcpy(r, new, newlen);
        r += newlen;
        p = q + oldlen;
    }
    strcpy(r, p);

    return ret;
}

void die(char* s)
{
    printf("%s", s);
    exit(0);
}

void chomp(char* s)
{
    
    int j;
    
    j = strlen(s);
    while (j >= 0) {
        if (s[j] == '\n') {
            s[j] = '\0';
            return;
        }
        j--;
    }
    
    
}

char* mystrtok(char* s, char delim)
{
    
    int posc;
    int l;
    char* myretstr;
    
    static char* sp;
    static int   posi;
    
    
    
    if (s != NULL) {
        sp   = s;
        posi = 0;
    }
    
    l = strlen(sp);
    
    if (posi > l)
        return NULL;
    
    
    // move to the next delim
    posc = posi;
    while ((posc < l) && (sp[posc] != delim)) {
        ///printf("posc++(becomes %d,c=%c)\n", posc+1, sp[posc+1]);
        posc++;
    }
    
    
    // alloc enough memory
    myretstr = (char*)calloc(posc - posi + 1, sizeof(char));
    
    // cpy
    strncat(myretstr, sp + posi, posc - posi);
    myretstr[posc - posi] = '\0';
    
    // update posi for next call  (posc points to the current delimiter)
    posi = posc+1;
    
    
    //printf("posi updated to %d\n", posi);
    
    return myretstr;
    
}

void readStringTable(char* filename, char**** data, int* n, int* m) {
    char* s;
    int   mym;
    int   myn;
    char* buff;
    FILE* fp;
    int   i;
    int   mynmax = 100000;
    
    buff  = (char*)malloc(100000 * sizeof(char));
    fp    = fopen(filename, "r");
    if (!fp) {
        die("readFloatTable: please enter a valid filename\n");
    }
    
    //
    // get the first line and estimate the number of columns
    //
    
    fgets(buff, mynmax, fp);
    chomp(buff);
    
    s = mystrtok(buff, '\t'); free(s);
    mym = 1;
    while ((s = mystrtok(0, '\t'))) {
        mym ++; free(s);
    }
    
    //
    // allocate the right number of columns
    //
    *data      = (char***)malloc(mynmax * sizeof(char**));
    (*data)[0] = (char**)malloc(mym * sizeof(char*));
    chomp(buff);
    s   = mystrtok(buff, '\t');
    (*data)[0][0] = strdup(s); free(s);
    i   = 1;
    while ((s = mystrtok(0, '\t'))) {
        (*data)[0][i] = strdup(s); free(s);
        i++;
    }
    
    myn = 1;
    while (1) {
        fgets(buff, mynmax, fp);
        
        if (feof(fp)) {
            break;
        }
        
        chomp(buff);
        
        
        (*data)[myn] = (char**)malloc(mym * sizeof(char*));
        
        s = mystrtok(buff, '\t');
        (*data)[myn][0] = strdup(s);  free(s);
        i = 1;
        while ((s = mystrtok(0, '\t'))) {
            
            (*data)[myn][i] = strdup(s);
            free(s);
            i ++;
            
        }
        
        myn ++;
        
        if (myn == mynmax) {
            die("readFloatTable: running out of memory, please recompile ..\n");
        }
    }
    
    *n = myn;
    *m = mym;
    
    free(buff); fclose(fp);
}
