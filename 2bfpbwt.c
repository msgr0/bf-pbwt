#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define FREE(x)                                                                \
  do {                                                                         \
    free((x));                                                                 \
    (x) = NULL;                                                                \
  } while (0)

#define parr(n, a, fmt)                                                        \
  do {                                                                         \
    for (int i = 0; i < (n); i++) {                                            \
      printf((fmt), (a)[i]);                                                   \
    }                                                                          \
    puts("");                                                                  \
  } while (0)

// get row and column count from file
void fgetrc(FILE *fd, size_t *restrict nr, size_t *restrict nc) {
  fseek(fd, 0, SEEK_SET);

  *nc = *nr = 0;
  int c, prev;
  while ((c = fgetc(fd)) != 0xA)
    (*nc)++;

  (*nr)++;
  prev = c;

  while ((c = fgetc(fd)) != EOF) {
    *nr += (c == 0xA) & (prev != 0xA);
    prev = c;
  }
}

// get column i from file
// c[n] is a pointer to store the column,
// nc is total number of columns, needed for fseek
static inline void fgetcoli(FILE *fd, size_t i, size_t n, uint8_t *restrict c,
                            size_t nc) {
  // NOTE: this assumes ASCII text file, offset are computed assuming
  // 1-byte size for each character
  int x;
  fseek(fd, i, SEEK_SET);
  for (size_t r = 0; r < n; r++) {
    x = fgetc(fd) - 48;
    fseek(fd, nc, SEEK_CUR);
    /*printf("%d\n", x);*/
    c[r] = x;
  }
}

typedef struct pbwtad pbwtad;
struct pbwtad {
  size_t *a;
  size_t *d;
};

// msgr0's version
// c[n] is a pointer to the column
// p is the `pbwtad` of A and D arrays of the previous column
// return: `pbwtad` of the current column
pbwtad *cpbwt(size_t n, uint8_t *restrict c, pbwtad *restrict p) {
  // NOTE: these two arrays do not need be allocated and free'd each time.
  // it would be possible to allocate it once (per process)
  // and re-use them each time.
  size_t *z = malloc(n * sizeof *z);
  size_t *o = malloc(n * sizeof *o);

  pbwtad *ret = malloc(sizeof *ret);
  size_t *a = malloc(n * sizeof *a);

  if (!z || !o || !a) {
    // FIXME: error
    return NULL;
  }
  size_t r = 0, q = 0;

  size_t i;
  for (i = 0; i < n; i++) {
    if (c[i] == 1) {
      o[q++] = p->a[i];
    } else {
      z[r++] = p->a[i];
    }
  }

  assert(r + q == n);
  for (i = 0; i < r; i++) {
    a[i] = z[i];
  }
  for (i = 0; i < q; i++) {
    a[r + i] = o[i];
  }

  FREE(o);
  FREE(z);

  ret->a = a;
  return ret;
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s FILE\n", argv[0]);
    return EXIT_FAILURE;
  }

  FILE *fin = fopen(argv[1], "r");
  if (!fin) {
    perror("[main]");
    return EXIT_FAILURE;
  }

  size_t nrow, ncol;
  printf("row: %zu, col:%zu\n", nrow, ncol);
  fgetrc(fin, &nrow, &ncol);
  printf("row: %zu, col:%zu\n", nrow, ncol);

  uint8_t *tcol = malloc(nrow * sizeof *tcol);

  pbwtad *p = malloc(sizeof *p);
  p->a = malloc(nrow * sizeof *(p->a));
  for (int j = 0; j < nrow; j++) {
    p->a[j] = j;
  }

  parr(nrow, p->a, "%zu ");
  for (int j = 0; j < ncol; j++) {
    fgetcoli(fin, j, nrow, tcol, ncol);
    p = cpbwt(nrow, tcol, p);
    parr(nrow, p->a, "%zu ");
  }

  fclose(fin);

  return EXIT_SUCCESS;
}
