#include "lib/quadsort/quadsort.h"
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

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
    uint64_t x;                                                                \
    fseek(fd, i *W, SEEK_SET);                                                 \
    for (size_t r = 0; r < n; r++) {                                           \
      c[r] = 0;                                                                \
      for (size_t s = 0; s < W; s++) {                                         \
        x = fgetc(fd) - 48;                                                    \
        c[r] = (x << s) | c[r];                                                \
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
void rrsortx(size_t n, uint64_t *c, size_t *s, size_t *aux) {
  size_t *tmp;
  size_t j;
  size_t *pre = s;
  size_t *post = aux;
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
void rrsort0(size_t n, uint64_t *c, size_t *s, size_t *aux) {
  size_t *tmp;
  size_t j;
  size_t *pre = s;
  size_t *post = aux;
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
  uint64_t x;
  fseek(fd, i * w, SEEK_SET);
  for (size_t r = 0; r < n; r++) {
    c[r] = 0;
    /*printf("r(");*/
    for (size_t s = 0; s < w; s++) {
      x = fgetc(fd) - 48;
      /*printf("%llu", x);*/
      c[r] = (x << s) | c[r];
    }
    /*printf(")=%llu\n", c[r]);*/
    fseek(fd, nc - w + 1, SEEK_CUR);
  }
}

// This version does not use window size to move in the file.
// It starts reading from column `i` as-is without window-offset computation
static inline void fgetcoliwgri(FILE *fd, size_t i, size_t n,
                                uint64_t *restrict c, size_t nc, uint8_t w) {
  // NOTE: this assumes ASCII text file, offset are computed assuming
  // 1-byte size for each character
  uint64_t x;
  fseek(fd, i, SEEK_SET);
  for (size_t r = 0; r < n; r++) {
    c[r] = 0;
    /*printf("r(");*/
    for (size_t s = 0; s < w; s++) {
      x = fgetc(fd) - 48;
      /*printf("%llu", x);*/
      c[r] = (x << s) | c[r];
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
// TODO: try with using external o and z arrays
pbwtad *cpbwt_0(size_t n, uint8_t *restrict c, pbwtad *restrict p) {
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
    if (c[p->a[i]] == 1) {
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
pbwtad *cpbwt(size_t n, uint8_t *restrict c, pbwtad *restrict p,
              size_t *restrict z, size_t *restrict o) {
  // NOTE: these two arrays do not need be allocated and free'd each time.
  // it would be possible to allocate it once (per process)
  // and re-use them each time.

  pbwtad *ret = malloc(sizeof *ret);
  size_t *a = malloc(n * sizeof *a);

  size_t r = 0, q = 0;

  size_t i;
  for (i = 0; i < n; i++) {
    if (c[p->a[i]] == 1) {
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

  ret->a = a;
  return ret;
}

#define TEST_LOG

#ifdef TEST_LOG
#define DPRINT(format, args...)                                                \
  do {                                                                         \
    fprintf(stderr, format, ##args);                                           \
  } while (0)

#define DPARR(n, a, fmt)                                                       \
  do {                                                                         \
    for (int i = 0; i < (n); i++) {                                            \
      fprintf(stderr, (fmt), (a)[i]);                                          \
    }                                                                          \
    fputc(0xA, stderr);                                                        \
  } while (0)
#else
#define DPRINT(args...)
#define DPARR(args...)
#endif

pbwtad **linc(FILE *fin, size_t nrow, size_t ncol) {
  uint8_t *c0 = malloc(nrow * sizeof *c0);
  size_t *o = malloc(nrow * sizeof *o);
  size_t *z = malloc(nrow * sizeof *z);
  // NOTE: right now I don't know what I need, so I'm keeping
  // everything in memory, we'll see later
  pbwtad **pl = malloc(ncol * sizeof(pbwtad *));
  if (!pl)
    return NULL;

  pbwtad *p0 = malloc(sizeof *p0);
  p0->a = malloc(nrow * sizeof *(p0->a));
  for (int j = 0; j < nrow; j++) {
    p0->a[j] = j;
  }
  fgetcoli(fin, 0, nrow, c0, ncol);
  pl[0] = cpbwt(nrow, c0, p0, z, o);
  FREE(p0->a);
  FREE(p0);

  for (size_t j = 1; j < ncol; j++) {
    fgetcoli(fin, j, nrow, c0, ncol);
    pl[j] = cpbwt(nrow, c0, pl[j - 1], z, o);
  }

#if 0
  for (size_t j = 0; j < ncol; j++) {
    DPRINT("lin %3zu: ", j);
    DPARR(nrow, pl[j]->a, "%zu ");
  }
#endif
  FREE(c0);
  FREE(o);
  FREE(z);
  return pl;
}

typedef enum { APPROX_MODE_LAST_WINDOW, APPROX_MODE_LAST_LIN } APPROX_MODE;
pbwtad **wapproxc_rrs(FILE *fin, size_t nrow, size_t ncol,
                      APPROX_MODE lastmode) {
  // NOTE: right now I don't know what I need, so I'm keeping
  // everything in memory, we'll see later
  pbwtad **pb = malloc(ncol * sizeof(pbwtad *));
#define W 64

  // Compute the bit-packed windows
  uint64_t *pw = malloc(nrow * sizeof *pw);
  size_t *aux = malloc(nrow * sizeof *aux);

  pbwtad *ps = malloc(nrow * sizeof *ps);
  ps->a = malloc(nrow * sizeof *(ps->a));
  pb[W - 1] = ps;
  fgetcoliw64r(fin, 0, nrow, pw, ncol);
  rrsort0(nrow, pw, ps->a, aux);

  size_t j;
  for (j = 1; j * W <= ncol - W; j++) {
    pbwtad *ps = malloc(nrow * sizeof *ps);
    ps->a = malloc(nrow * sizeof *(ps->a));
    pb[W * (j + 1) - 1] = ps;
    memcpy(ps->a, pb[W * j - 1]->a, nrow * sizeof *(ps->a));
    fgetcoliw64r(fin, j, nrow, pw, ncol);
    rrsortx(nrow, pw, ps->a, aux);
  }

  uint8_t *c0 = NULL;
  switch (lastmode) {
  case APPROX_MODE_LAST_LIN:
    c0 = malloc(nrow * sizeof *c0);
    for (j = j * W; j < ncol; j++) {
      fgetcoli(fin, j, nrow, c0, ncol);
      pb[j] = cpbwt_0(nrow, c0, pb[j - 1]);
    }
    FREE(c0);
    break;
  case APPROX_MODE_LAST_WINDOW:
    j *= W;
    fgetcoliwgri(fin, j, nrow, pw, ncol, ncol - j);
    pbwtad *ps = malloc(nrow * sizeof *ps);
    ps->a = malloc(nrow * sizeof *(ps->a));
    pb[ncol - 1] = ps;
    memcpy(ps->a, pb[j - 1]->a, nrow * sizeof *(ps->a));
    rrsortx(nrow, pw, ps->a, aux);
    break;
  }

#if 0
  for (size_t j = 0; j < ncol; j++) {
    DPRINT("ars %3zu: ", j);
    if (pb[j])
      DPARR(nrow, pb[j]->a, "%zu ");
    else
      DPRINT(" ---\n");
  }
#endif
  FREE(c0);
  return pb;
}

// WARN: I believe that this might be wrong.
// I am not sure it works, as to use the order based on the previous column
// I tried with pw and pw0, and now it works.
// The next idea was to use aux to save the value read in order of the previous
// ordering but it doesn't work. Now I have fever and I am not in the best
// condition to understand what is happening.
// However this might not be needed at all.
pbwtad **wapproxc_qs(FILE *fin, size_t nrow, size_t ncol,
                     APPROX_MODE lastmode) {
  // NOTE: right now I don't know what I need, so I'm keeping
  // everything in memory, we'll see later
  pbwtad **pb = malloc(ncol * sizeof(pbwtad *));
#define W 64

  // Compute the bit-packed windows
  uint64_t *pw = malloc(nrow * sizeof *pw);
  uint64_t *pw0;
  /*uint64_t *aux = malloc(nrow * sizeof *aux);*/

  pbwtad *ps = malloc(nrow * sizeof *ps);
  ps->a = malloc(nrow * sizeof *(ps->a));
  pb[W - 1] = ps;
  fgetcoliw64r(fin, 0, nrow, pw, ncol);
  uint64_t **qsp = quadsort_u64_ix(pw, nrow, NULL);
  pw0 = pw;

  size_t j;
  for (j = 1; j * W <= ncol - W; j++) {
    fgetcoliw64r(fin, j, nrow, pw, ncol);
    for (size_t i = 0; i < nrow; i++) {
      ps->a[i] = qsp[i] - pw0;
      /*aux[i] = pw[ps->a[i]];*/
    }
    ps = malloc(nrow * sizeof *ps);
    ps->a = malloc(nrow * sizeof *(ps->a));
    pb[W * (j + 1) - 1] = ps;
    qsp = quadsort_u64_ix(pw, nrow, NULL);
    pw0 = pw;
  }
  for (size_t i = 0; i < nrow; i++) {
    ps->a[i] = qsp[i] - pw0;
  }

  uint8_t *c0 = NULL;
  switch (lastmode) {
  case APPROX_MODE_LAST_LIN:
    c0 = malloc(nrow * sizeof *c0);
    for (j = j * W; j < ncol; j++) {
      fgetcoli(fin, j, nrow, c0, ncol);
      pb[j] = cpbwt_0(nrow, c0, pb[j - 1]);
    }
    FREE(c0);
    break;
  case APPROX_MODE_LAST_WINDOW:
    j *= W;
    fgetcoliwgri(fin, j, nrow, pw, ncol, ncol - j);
    ps = malloc(nrow * sizeof *ps);
    ps->a = malloc(nrow * sizeof *(ps->a));
    pb[ncol - 1] = ps;
    qsp = quadsort_u64_ix(pw, nrow, NULL);
    for (size_t i = 0; i < nrow; i++) {
      ps->a[i] = qsp[i] - pw0;
    }
    break;
  }

#if 0
  for (size_t j = 0; j < ncol; j++) {
    DPRINT("lin %3zu: ", j);
    if (pb[j])
      DPARR(nrow, pb[j]->a, "%zu ");
    else
      DPRINT(" ---\n");
  }
#endif
  FREE(c0);
  /*FREE(pw);*/
  FREE(pw0);
  return pb;
}

#define SWAP(x, y)                                                             \
  do {                                                                         \
    typeof((x)) tmp = (x);                                                     \
    (x) = (y);                                                                 \
    (y) = tmp;                                                                 \
  } while (0)

pbwtad **wparc_rrs(FILE *fin, size_t nrow, size_t ncol) {
  // NOTE: right now I don't know what I need, so I'm keeping
  // everything in memory, we'll see later
  pbwtad **pb = malloc(ncol * sizeof(pbwtad *));
#define W 64

  // first W=(64 for now), must be computed linearly
  // TODO: ask if true

  uint8_t *c0 = malloc(nrow * sizeof *c0);
  pbwtad *p0 = malloc(sizeof *p0);
  p0->a = malloc(nrow * sizeof *(p0->a));
  for (int j = 0; j < nrow; j++) {
    p0->a[j] = j;
  }
  fgetcoli(fin, 0, nrow, c0, ncol);
  pb[0] = cpbwt_0(nrow, c0, p0);
  FREE(p0->a);
  FREE(p0);

  for (int j = 1; j < W; j++) {
    fgetcoli(fin, j, nrow, c0, ncol);
    pb[j] = cpbwt_0(nrow, c0, pb[j - 1]);
  }

  uint64_t *pw0 = malloc(nrow * sizeof *pw0); // the first one is not freed
  uint64_t *pw1 = malloc(nrow * sizeof *pw1);
  size_t *aux = malloc(nrow * sizeof *aux);
  fgetcoliw64r(fin, 0, nrow, pw0, ncol);

  size_t j;

  /*for (j = 1; j * W <= W * 2; j++) {*/
  for (j = 1; j * W <= ncol - W; j++) {
    pbwtad *ps = malloc(nrow * sizeof *ps);
    ps->a = malloc(nrow * sizeof *(ps->a));
    pb[W * (j + 1) - 1] = ps;
    memcpy(ps->a, pb[W * j - 1]->a, nrow * sizeof *(ps->a));
    fgetcoliw64r(fin, j, nrow, pw1, ncol);
    rrsortx(nrow, pw1, ps->a, aux);

// okay now we should the inner loops, this is the part that
// can be parallel
#pragma omp parallel for num_threads(8) shared(pw1, pw0, pb, j)
    for (size_t x = 1; x < W; x++) {
      uint64_t *w = malloc(nrow * sizeof *w);
      size_t J = W * (j + 1) - 1;

      wr64mrgsi(nrow, pw1, pw0, w, x);
      pbwtad *ps = malloc(nrow * sizeof *ps);
      ps->a = malloc(nrow * sizeof *(ps->a));
      pb[J - x] = ps;
      memcpy(ps->a, pb[J - W - x]->a, nrow * sizeof *(ps->a));
      rrsortx(nrow, w, ps->a, aux);

      FREE(w);
    }

    SWAP(pw0, pw1);
  }

  c0 = malloc(nrow * sizeof *c0);
  for (j = j * W; j < ncol; j++) {
    fgetcoli(fin, j, nrow, c0, ncol);
    pb[j] = cpbwt_0(nrow, c0, pb[j - 1]);
  }

#if 0
  for (size_t j = 0; j < ncol; j++) {
    DPRINT("prs %3zu: ", j);
    if (pb[j])
      DPARR(nrow, pb[j]->a, "%zu ");
    else
      DPRINT(" ---\n");
  }
#endif
  FREE(pw0);
  FREE(pw1);
  FREE(aux);
  FREE(c0);
  return pb;
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s [lin|ars|aqs|prs] FILE\n", argv[0]);
    return EXIT_FAILURE;
  }

  FILE *fin = fopen(argv[2], "r");
  if (!fin) {
    perror("[main]");
    return EXIT_FAILURE;
  }

  size_t nrow, ncol;
  fgetrc(fin, &nrow, &ncol);
  DPRINT("[%s] row: %5zu, col: %5zu\n", __func__, nrow, ncol);

  if (strcmp(argv[1], "lin") == 0) {
    linc(fin, nrow, ncol);
  } else if (strcmp(argv[1], "ars") == 0) {
    wapproxc_rrs(fin, nrow, ncol, APPROX_MODE_LAST_WINDOW);
  } else if (strcmp(argv[1], "aqs") == 0) {
    wapproxc_qs(fin, nrow, ncol, APPROX_MODE_LAST_WINDOW);
  } else if (strcmp(argv[1], "prs") == 0) {
    pbwtad **r = wparc_rrs(fin, nrow, ncol);
    for (size_t i = 0; i < ncol; i++) {
      if (r[i]) {
        FREE(r[i]->a);
        FREE(r[i]);
      }
    }
  }

  return EXIT_FAILURE;
}

int main2(int argc, char *argv[]) {
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

  DPRINT("[%s] row: %5zu, col: %5zu\n", __func__, nrow, ncol);

  uint8_t *c0 = malloc(nrow * sizeof *c0);
  // NOTE: right now I don't know what I need, so I'm keeping
  // everything in memory, we'll see later
  pbwtad **pl = malloc(ncol * sizeof(pbwtad *));
  pbwtad **pb = malloc(ncol * sizeof(pbwtad *));
#define W 64
  // first W=(64 for now), must be computed linearly
  // TODO: ask if true

  pbwtad *p0 = malloc(sizeof *p0);
  p0->a = malloc(nrow * sizeof *(p0->a));
  for (int j = 0; j < nrow; j++) {
    p0->a[j] = j;
  }
  fgetcoli(fin, 0, nrow, c0, ncol);
  pl[0] = cpbwt_0(nrow, c0, p0);
  FREE(p0->a);
  FREE(p0);

  for (int j = 1; j < W; j++) {
    fgetcoli(fin, j, nrow, c0, ncol);
    pl[j] = cpbwt_0(nrow, c0, pl[j - 1]);
    pb[j] = pl[j];
  }
  DPRINT("\tlin %3d: ", W - 1);
  DPARR(nrow, pl[W - 1]->a, "%zu ");

  for (int j = W; j < 2 * W; j++) {
    fgetcoli(fin, j, nrow, c0, ncol);
    pl[j] = cpbwt_0(nrow, c0, pl[j - 1]);
  }
  DPRINT("\tlin %3d: ", 2 * W - 1);
  DPARR(nrow, pl[2 * W - 1]->a, "%zu ");

  // now let's start with the others.
  // First I need to get the bitpack of the first (computed) W cols
  uint64_t *pw0 = malloc(nrow * sizeof *pw0);
  uint64_t *pw1 = malloc(nrow * sizeof *pw1);
  size_t *aux = malloc(nrow * sizeof *aux);
  fgetcoliw64r(fin, 0, nrow, pw0, ncol);
  /*fgetcoliw64r(fin, 1, nrow, pw1, ncol);*/
  // FIXME: for non W-divisible input

  size_t j;
  /*for (j = 1; j < 2; j++) {*/
  /*  pbwtad *ps = malloc(nrow * sizeof *ps);*/
  /*  ps->a = malloc(nrow * sizeof *(ps->a));*/
  /*  pb[W * (j + 1) - 1] = ps;*/
  /*  memcpy(ps->a, pb[W * j - 1]->a, nrow * sizeof *(ps->a));*/
  /*  fgetcoliw64r(fin, j, nrow, pw1, ncol);*/
  /*  rrsortx(nrow, pw1, ps->a, aux);*/
  /**/
  /*  pw0 = pw1;*/
  /*}*/
  /*DPRINT("\trrs %3d: ", 2 * W - 1);*/
  /*DPARR(nrow, pb[2 * W - 1]->a, "%zu ");*/

  for (j = 1; j * W <= ncol - W; j++) {
    pbwtad *ps = malloc(nrow * sizeof *ps);
    ps->a = malloc(nrow * sizeof *(ps->a));
    pb[W * (j + 1) - 1] = ps;
    memcpy(ps->a, pb[W * j - 1]->a, nrow * sizeof *(ps->a));
    fgetcoliw64r(fin, j, nrow, pw1, ncol);
    rrsortx(nrow, pw1, ps->a, aux);

    pw0 = pw1;
    DPRINT("\trrs %3zu: ", (j + 1) * W - 1);
    DPARR(nrow, pb[(j + 1) * W - 1]->a, "%zu ");
  }
  /*for (j = j * W; j < ncol; j++) {*/
  /*  fgetcoli(fin, j, nrow, c0, ncol);*/
  /*  pb[j] = cpbwt(nrow, c0, pb[j - 1]);*/
  /*  DPRINT("\tlin %3zu: ", j);*/
  /*  DPARR(nrow, pb[j]->a, "%zu ");*/
  /*}*/

  j *= W;
  fgetcoliwgri(fin, j, nrow, pw1, ncol, ncol - j);
  pbwtad *ps = malloc(nrow * sizeof *ps);
  ps->a = malloc(nrow * sizeof *(ps->a));
  pb[ncol - 1] = ps;
  memcpy(ps->a, pb[j - 1]->a, nrow * sizeof *(ps->a));
  rrsortx(nrow, pw1, ps->a, aux);
  DPRINT("\tlin %3zu: ", ncol - 1);
  DPARR(nrow, pb[ncol - 1]->a, "%zu ");

#if 0 // NOTE: just a test that window work as expected
  uint64_t *sorted = malloc(nrow * sizeof *sorted);
  uint64_t *aux = malloc(nrow * sizeof *aux);
  rrsort0(nrow, pw0, sorted, aux);
  DPRINT("\trrs: ");
  DPARR(nrow, sorted, "%llu ");

  uint64_t **qsp = quadsort_u64_ix(pw0, nrow, NULL);
  DPRINT("\tqds: ");
  for (size_t i = 0; i < nrow; i++) {
    DPRINT("%ld ", qsp[i] - pw0);
  }
  puts("");
#endif

  return EXIT_SUCCESS;
}

int maintest(int argc, char *argv[]) {
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
    p = cpbwt_0(nrow, tcol, p);
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
  fgetcoliwgri(fin, 1, nrow, twc, ncol, 6);
  parr(nrow, twp, "%3llu ");

  _sort_window_brian(nrow, tc);

  fclose(fin);
  return EXIT_SUCCESS;
}
