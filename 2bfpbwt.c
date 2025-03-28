#include "lib/quadsort/quadsort.h"
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

#define FGETCOLIW_DECLARE(W)                                                   \
  static inline void fgetcoliw##W(FILE *fd, size_t i, size_t n,                \
                                  uint64_t *restrict c, size_t nc) {           \
    int x;                                                                     \
    fseek(fd, i *W, SEEK_SET);                                                 \
    for (size_t r = 0; r < n; r++) {                                           \
      c[r] = 0;                                                                \
      /* NOTE: should be tried to do unrolling */                              \
      /*       check commented version below for how */                        \
      for (size_t s = 0; s < W; s++) {                                         \
        x = fgetc(fd) - 48;                                                    \
        c[r] = (c[r] << 1) | x;                                                \
      }                                                                        \
      fseek(fd, nc - W + 1, SEEK_CUR);                                         \
    }                                                                          \
  }                                                                            \
  static inline void w##W##mrgsi(size_t n, uint64_t const *wc,                 \
                                 uint64_t const *wp, uint64_t *restrict c,     \
                                 size_t i) {                                   \
    uint64_t c1;                                                               \
    for (size_t r = 0; r < n; r++) {                                           \
      c1 = wp[r] & ((1 << i) - 1);                                             \
      c[r] = (c1 << (W - i)) | (wc[r] >> i);                                   \
    }                                                                          \
  }                                                                            \
  static inline void fgetcoliw##W##r(FILE *fd, size_t i, size_t n,             \
                                     uint64_t *restrict c, size_t nc) {        \
    int x;                                                                     \
    fseek(fd, i *W, SEEK_SET);                                                 \
    for (size_t r = 0; r < n; r++) {                                           \
      c[r] = 0;                                                                \
      for (size_t s = 0; s < W; s++) {                                         \
        x = fgetc(fd) - 48;                                                    \
        c[r] = (x << (s)) | c[r];                                              \
      }                                                                        \
      fseek(fd, nc - W + 1, SEEK_CUR);                                         \
    }                                                                          \
  }                                                                            \
  static inline void wr##W##mrgsi(size_t n, uint64_t const *wc,                \
                                  uint64_t const *wp, uint64_t *restrict c,    \
                                  size_t i) {                                  \
    uint64_t c1;                                                               \
    static const uint64_t mask = (UINT64_MAX >> (64 - W));                     \
    for (size_t r = 0; r < n; r++) {                                           \
      c[r] = (wp[r] >> (W - i)) | ((wc[r] << i) & mask);                       \
    }                                                                          \
  }
/*static inline void fgetcoliw8(FILE *fd, size_t i, size_t n,*/
/*                              uint64_t *restrict c, size_t nc) {*/
/*  int x;*/
/*  fseek(fd, i * 8, SEEK_SET);*/
/*  for (size_t r = 0; r < n; r++) {*/
/*    c[r] = 0;*/
/*#ifdef __GNUC__*/
/*#ifdef __clang__*/
/*#pragma unroll*/
/*#else*/
/*#pragma GCC unroll 8*/
/*#endif*/
/*#endif*/
/*    for (size_t s = 0; s < 8; s++) {*/
/*      x = fgetc(fd) - 48;*/
/*      c[r] = (c[r] << 1) | x;*/
/*    }*/
/*    fseek(fd, nc - 8 + 1, SEEK_CUR);*/
/*  }*/
/*}*/

FGETCOLIW_DECLARE(8)
FGETCOLIW_DECLARE(16)
FGETCOLIW_DECLARE(32)
FGETCOLIW_DECLARE(64)

// Window MeRGe Shift I:
// merge windows `wc[n]` and `wp[n]`, storing in `c`
// with `n` number of rows
static inline void wmrgsi(size_t n, uint64_t const *wc, uint64_t const *wp,
                          uint64_t *restrict c, size_t i, uint8_t w) {
  // How it works
  // |    wp    |    wc    |
  // | '. . . .'| '. . . .'|
  // |    w1    |    w2    |
  // i = 2
  // |  . .'. . |  . .'. . |
  // c1:
  //   w1 & ((1<<i) -1)
  //            c2:
  //              w2 >> i
  // c = (c1 << (w-i)) | c2
  uint64_t c1;
  for (size_t r = 0; r < n; r++) {
    c1 = wp[r] & ((1 << i) - 1);
    c[r] = (c1 << (w - i)) | (wc[r] >> i);
  }
}

// Window-Reverse MeRGe Shift I:
// merge windows `wc[n]` and `wp[n]`, storing in `c`
// with `n` number of rows
static inline void wrmrgsi(size_t n, uint64_t const *wc, uint64_t const *wp,
                           uint64_t *restrict c, size_t i, uint8_t w) {
  // How it works: similary to `wmrgsi`, but now windows are reversed
  //  |   wp  |   wc  |
  //  | abcde | 12345 |
  //  |  w1   |  w2   |
  //
  //  Represented as
  //    edcba   54321
  //  Suppose w=5 i=2, result would be
  //    de123
  //  represented as
  //    321ed
  // c1:
  //    w1 >> (w-i)
  //            c2:
  //              w2 << i
  // c = c1 | c2
  uint64_t c1;
  uint64_t mask = (UINT64_MAX >> (64 - w));
  for (size_t r = 0; r < n; r++) {
    c[r] = (wp[r] >> (w - i)) | ((wc[r] << i) & mask);
  }
}

static inline void _sort_window(size_t n, uint64_t *restrict c) {
  /*quadsort_uint64(c, n, CMPFUNC *cmp);*/
  quadsort_prim(c, n, sizeof(long long) + 1);
}

/*
 * Sort `c[n]`, sorted permutations will be in saved in `s[n]`,
 * using externally allocated `aux[n]` auxiliary array.
 * This version assumes the `s` array to be already initialized.
 */
void rrsortx(size_t n, uint64_t *c, uint64_t *s, uint64_t *aux) {
  uint64_t *tmp;
  size_t j;
  uint64_t *pre = s;
  uint64_t *post = aux;
  uint8_t b;

  for (size_t i = 0; i < 8; i++) {
    size_t cnt[256] = {0};

    // frequencies
    for (j = 0; j < n; j++) {
      /*cnt[mask2(c[j], i)]++;*/
      b = (c[j] >> (8 * i)) & 0xFFULL;
      cnt[b]++;
    }
    // prefix sum
    for (size_t j = 1; j < 256; j++)
      cnt[j] += cnt[j - 1];
    // sorting
    for (ssize_t j = n - 1; j >= 0; --j) {
      /*cnt[mask2(c[pre[j]], i)]--;*/
      /*post[cnt[mask2(c[pre[j]], i)]] = pre[j];*/

      b = (c[pre[j]] >> (8 * i)) & 0xFFULL;
      cnt[b]--;
      post[cnt[b]] = pre[j];
    }
    // swap s and aux
    tmp = pre;
    pre = post;
    post = tmp;
  }
}

/*
 * Sort `c[n]`, sorted permutations will be in saved in `s[n]`,
 * using externally allocated `aux[n]` auxiliary array.
 * This version initialize the sorting from position 0,
 * meaning that there will be a pass of setting the initial
 * positions array to 1..n
 */
void rrsort0(size_t n, uint64_t *c, uint64_t *s, uint64_t *aux) {
  uint64_t *tmp;
  size_t j;
  uint64_t *pre = s;
  uint64_t *post = aux;
  uint8_t b;

  // this is needed if:
  // 1. we want to sort numbers
  // 2. this is the first iteration
  //
  // In normal BWT cases we assume to have
  // `s` equal to the sorting of the previous "column"
  for (size_t i = 0; i < n; ++i)
    pre[i] = i;

  for (size_t i = 0; i < 8; i++) {
    size_t cnt[256] = {0};

    // frequencies
    for (j = 0; j < n; j++) {
      /*cnt[mask2(c[j], i)]++;*/
      b = (c[j] >> (8 * i)) & 0xFFULL;
      cnt[b]++;
    }
    // prefix sum
    for (size_t j = 1; j < 256; j++)
      cnt[j] += cnt[j - 1];
    // sorting
    for (ssize_t j = n - 1; j >= 0; --j) {
      /*cnt[mask2(c[pre[j]], i)]--;*/
      /*post[cnt[mask2(c[pre[j]], i)]] = pre[j];*/

      b = (c[pre[j]] >> (8 * i)) & 0xFFULL;
      cnt[b]--;
      post[cnt[b]] = pre[j];
    }
    // swap s and aux
    tmp = pre;
    pre = post;
    post = tmp;
  }
}

static inline void _sort_window_brian(size_t n, uint64_t *restrict c) {
  size_t cnt[256] = {0};
  for (size_t i = 0; i < 8; i++) {
    uint8_t b = (c[n] >> (64 - (8 * (i + 1)))) & 0xFF;
    printf("%zu ", 64 - (8 * (i + 1)));
  }
  puts("");
  for (ssize_t i = 8; i > 0; i--) {
    uint8_t b = (c[n] >> (8 * (i - 1))) & 0xFF;
    printf("%zu ", 8 * (i - 1));
  }
  puts("");
}

static inline void fgetcoliwg(FILE *fd, size_t i, size_t n,
                              uint64_t *restrict c, size_t nc, uint8_t w) {
  // NOTE: this assumes ASCII text file, offset are computed assuming
  // 1-byte size for each character
  int x;
  fseek(fd, i * w, SEEK_SET);
  for (size_t r = 0; r < n; r++) {
    c[r] = 0;
    for (size_t s = 0; s < w; s++) {
      x = fgetc(fd) - 48;
      /*printf("%d", x);*/
      c[r] = (c[r] << 1) | x;
    }
    /*printf("=%llu", c[r]);*/
    /*puts("");*/
    fseek(fd, nc - w + 1, SEEK_CUR);
  }
}
static inline void fgetcoliwgr(FILE *fd, size_t i, size_t n,
                               uint64_t *restrict c, size_t nc, uint8_t w) {
  // NOTE: this assumes ASCII text file, offset are computed assuming
  // 1-byte size for each character
  int x;
  fseek(fd, i * w, SEEK_SET);
  for (size_t r = 0; r < n; r++) {
    c[r] = 0;
    /*printf("r(");*/
    for (size_t s = 0; s < w; s++) {
      x = fgetc(fd) - 48;
      /*printf("%d", x);*/
      c[r] = (x << (s)) | c[r];
    }
    /*printf(")=%llu\n", c[r]);*/
    fseek(fd, nc - w + 1, SEEK_CUR);
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
  puts("\n---------- Column reading ----------------------------------\n");

  uint64_t *tmcol = malloc(nrow * sizeof *tmcol);

  fgetcoliwg(fin, 0, nrow, tmcol, ncol, 8);
  parr(nrow, tmcol, "%llu ");
  fgetcoliw8(fin, 0, nrow, tmcol, ncol);
  parr(nrow, tmcol, "%llu ");

  puts("   .... reverse ....");
  fgetcoliwgr(fin, 0, nrow, tmcol, ncol, 8);
  parr(nrow, tmcol, "%llu ");
  fgetcoliw8r(fin, 0, nrow, tmcol, ncol);
  parr(nrow, tmcol, "%llu ");
  puts("\n---------- Shifting of windows -----------------------------\n");

  uint64_t *twc = malloc(nrow * sizeof *twc);
  uint64_t *twp = malloc(nrow * sizeof *twp);
  uint64_t *tc = malloc(nrow * sizeof *tc);

  fgetcoliw8(fin, 0, nrow, twp, ncol);
  fgetcoliw8(fin, 1, nrow, twc, ncol);
  parr(nrow, twp, "%3llu ");
  parr(nrow, twc, "%3llu ");

  w8mrgsi(nrow, twc, twp, tc, 3);
  parr(nrow, tc, "%3llu ");

  puts("   .... reverse ....");
  fgetcoliw8r(fin, 0, nrow, twp, ncol);
  fgetcoliw8r(fin, 1, nrow, twc, ncol);
  parr(nrow, twp, "%3llu ");
  parr(nrow, twc, "%3llu ");

  wrmrgsi(nrow, twc, twp, tc, 2, 8);
  parr(nrow, tc, "%3llu ");

  puts("\n------------------------------------------------------------\n");
  _sort_window(nrow, tc);
  parr(nrow, tc, "%3llu ");

  _sort_window_brian(nrow, tc);

  fclose(fin);
  return EXIT_SUCCESS;
}
