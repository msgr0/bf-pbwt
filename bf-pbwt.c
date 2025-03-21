#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
// #include <qsort.h>

#define N_ROWS 5
#define M_COLUMNS 6
#define WINDOW 1

typedef struct pbwt pbwt;


struct pbwt{
    size_t k;
    size_t dim; // n& of rows
    size_t* pref;
    size_t* div;
};


void compute_pbwt(pbwt* p, uint8_t mat[][N_ROWS], size_t w) {
    uint8_t allele;
    size_t* zeros;
    size_t* ones; //1
    zeros = malloc(sizeof(*zeros) * (p->dim));
    ones = malloc(sizeof(*ones) * (p->dim));
    size_t r = 0;
    size_t q = 0;

    for (size_t i = 0; i < p->dim; i++) {
        size_t row = p->pref[i];
        allele = mat[p->k][row];

        if (allele) {
            ones[q++] = row;
        }
        else {
            zeros[r++] = row;
        }

    }

    assert((r + q)== p->dim);
    for (size_t i = 0; i < r; i ++) {
        p->pref[i] = zeros[i];
    }

    for (size_t i = 0; i < q; i ++) {
        p->pref[r+i] = ones[i];
    }
    p->k++;
    free(ones);
    free(zeros);

}

int main(int argc, char* argv[]) {

    pbwt *p;
    p = malloc(sizeof(p));

    p->k = 0;
    p->dim = N_ROWS;
    p->pref = (size_t[N_ROWS]){0, 1, 2, 3, 4};
    p->div = (size_t[N_ROWS]){0, 0, 0, 0, 0};

    uint8_t mat[M_COLUMNS][ N_ROWS ] = {
        {0, 1, 0, 1, 0},
        {1, 0, 0, 1, 1},
        {0, 0, 1, 0, 0},
        {1, 1, 0, 1, 1},
        {1, 0, 1, 0, 0},
        {1, 0, 1, 1, 0}
    };
    
    for (int j = 0; j < M_COLUMNS ; j=j+WINDOW) {
        compute_pbwt(p, mat, WINDOW);
        printf("pref at %d:", j);
        for (int i = 0; i < p->dim; i++) {
            printf(" %zu", p->pref[i]);
        }
        printf("\n");
        
    }





    return EXIT_SUCCESS;
}