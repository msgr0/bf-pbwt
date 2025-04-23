#include "io.h"
#include "lib/quadsort/quadsort.h"
#include <assert.h>
#include <fcntl.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef BF2IOMODE_BCF
#include "htslib/synced_bcf_reader.h"
#endif

#define W 64

#define DBDUMP
uint8_t DO_DUMP = 0;
#ifdef DBDUMP
#define PDUMP(i, p)                                                            \
  do {                                                                         \
    if (DO_DUMP) {                                                             \
      printf("%zu:", (size_t)(i));                                             \
      size_t pdump_j__;                                                        \
      for (pdump_j__ = 0; pdump_j__ < nrow - 1; pdump_j__++)                   \
        printf("%zu ", (p)->a[pdump_j__]);                                     \
      printf("%zu", (p)->a[pdump_j__]);                                        \
      fputc('|', stdout);                                                      \
      for (pdump_j__ = 0; pdump_j__ < nrow - 1; pdump_j__++)                   \
        printf("%zu ", (p)->d[pdump_j__]);                                     \
      printf("%zu", (p)->d[pdump_j__]);                                        \
      fputc(0xA, stdout);                                                      \
    }                                                                          \
  } while (0)

#define PDUMP_SEQ(s, e, p)                                                     \
  do {                                                                         \
    for (size_t pdump_ix__ = (s); pdump_ix__ < (e); pdump_ix__++) {            \
      if (DO_DUMP) {                                                           \
        printf("%zu:", (size_t)(pdump_ix__));                                  \
        size_t pdump_j__;                                                      \
        for (pdump_j__ = 0; pdump_j__ < nrow - 1; pdump_j__++)                 \
          printf("%zu ", (p)[pdump_ix__]->a[pdump_j__]);                       \
        printf("%zu", (p)[pdump_ix__]->a[pdump_j__]);                          \
        fputc('|', stdout);                                                    \
        for (pdump_j__ = 0; pdump_j__ < nrow - 1; pdump_j__++)                 \
          printf("%zu ", (p)[pdump_ix__]->d[pdump_j__]);                       \
        printf("%zu", (p)[pdump_ix__]->d[pdump_j__]);                          \
        fputc(0xA, stdout);                                                    \
      }                                                                        \
    }                                                                          \
  } while (0)
#else
#define PDUMP(p)
#define PDUMP_SEQ(s, e, p)
#endif

#define FREE(x)                                                                \
  do {                                                                         \
    free((x));                                                                 \
    (x) = NULL;                                                                \
  } while (0)

#define SWAP(x, y)                                                             \
  do {                                                                         \
    typeof((x)) tmp = (x);                                                     \
    (x) = (y);                                                                 \
    (y) = tmp;                                                                 \
  } while (0)

#define parr(n, a, fmt)                                                        \
  do {                                                                         \
    for (int parr_i__ = 0; parr_i__ < (n); parr_i__++) {                       \
      printf((fmt), (a)[parr_i__]);                                            \
    }                                                                          \
    puts("");                                                                  \
  } while (0)

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
 * This version assumes the `s` array to be already initialized.
 */
void rrsortx_noaux(size_t n, uint64_t *c, size_t *s) {
  size_t *tmp;
  size_t j;
  size_t *pre = s;
  size_t *post = malloc(n * sizeof *post);
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
  FREE(post);
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

  // TODO: this should also consider for D[-1]

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

typedef struct pbwtad pbwtad;
struct pbwtad {
  size_t *a;
  size_t *d;
};

static inline pbwtad *pbwtad_new(size_t n) {
  pbwtad *p = malloc(sizeof *p);
  p->a = malloc(n * sizeof *(p->a));
  p->d = malloc(n * sizeof *(p->d));
  // NOTE: maybe we want to have a look at continuos array for both a and d and
  // allocate as follows, in which case the struct might be changed a bit:
  // pbwtad *p = malloc(sizeof *p + 2*n * sizeof *(p->data));
  return p;
}

#define PBWTAD_FREE(p)                                                         \
  do {                                                                         \
    FREE((p)->a);                                                              \
    FREE((p)->d);                                                              \
    FREE(p);                                                                   \
  } while (0)

// msgr0's version
// c[n] is a pointer to the column
// p is the `pbwtad` of A and D arrays of the previous column
// return: `pbwtad` of the current column
pbwtad *cpbwt_0(size_t n, uint8_t *restrict c, pbwtad *restrict p) {
  // NOTE: these two arrays do not need be allocated and free'd each time.
  // it would be possible to allocate it once (per process)
  // and re-use them each time.
  size_t *z = malloc(n * sizeof *z);
  size_t *o = malloc(n * sizeof *o);
  size_t *zd = malloc(n * sizeof *zd);
  size_t *od = malloc(n * sizeof *od);

  pbwtad *ret = malloc(sizeof *ret);
  size_t *a = malloc(n * sizeof *a);
  size_t *d = malloc(n * sizeof *d);

  if (!z || !o || !a) {
    // FIXME: error
    return NULL;
  }
  size_t r = 0, q = 0;
  size_t f = n, g = n;

  size_t i;
  for (i = 0; i < n; i++) {

    if (c[p->d[i]] > f) {
      f = c[p->d[i]];
    }
    if (c[p->d[i]] > g) {
      g = c[p->d[i]];
    }

    if (c[p->a[i]] == 1) {
      o[q] = p->a[i];
      od[q++] = g = 0;
      g = 0;
    } else {
      z[r] = p->a[i];
      zd[r++] = f;
      f = 0;
    }
  }

  assert(r + q == n);
  for (i = 0; i < r; i++) {
    a[i] = z[i];
    d[i] = zd[i];
  }
  for (i = 0; i < q; i++) {
    a[r + i] = o[i];
    d[r + i] = od[i];
  }

  FREE(o);
  FREE(z);
  FREE(zd);
  FREE(od);

  ret->a = a;
  ret->d = d;
  return ret;
}

static pbwtad *cpbwt(size_t n, uint8_t *restrict c, pbwtad *restrict p) {
  static size_t *o = NULL;
  static size_t *h = NULL;
  static size_t k = 1;

  if (!o)
    o = malloc(n * sizeof *o);
  if (!h)
    h = malloc(n * sizeof *h);

  pbwtad *ret = malloc(sizeof *ret);
  ret->a = malloc(n * sizeof *(ret->a));
  ret->d = malloc(n * sizeof *(ret->d));

  size_t r = 0, q = 0;
  size_t f = k, g = k;

  size_t i;
#if 0
  for (i = 0; i < n; i++) {
    /*printf("i: %6zu - p->a[i]: %zu\n", i, p->a[i]);*/
    if (c[p->a[i]] == 1) {
      o[q++] = p->a[i];
    } else {
      ret->a[r++] = p->a[i];
      /*z[r++] = p->a[i];*/
    }
  }
#else
  for (i = 0; i < n; i++) {
    size_t idx = p->a[i];
    size_t ddx = p->d[i];

    f = (ddx > f) ? ddx : f;
    g = (ddx > g) ? ddx : g;

    size_t mask = c[idx];
    o[q] = idx;
    ret->a[r] = idx;
#if 0
    if (mask) {
      h[q] = g;
      g = 0;
    } else {
      ret->d[r] = f;
      f = 0;
    }
#else
    h[q] = g;
    ret->d[r] = f;

    f &= -mask;       // f = 0 if mask == 0, unchanged if mask == 1
    g &= -(1 - mask); // g = 0 if mask == 1, unchanged if mask == 0
#endif
    q += mask;     // Increment q if mask is 1
    r += mask ^ 1; // Increment r if mask is 0
  }
#endif

  memcpy(ret->a + r, o, q * sizeof(size_t));
  memcpy(ret->d + r, h, q * sizeof(size_t));

  k++;
  return ret;
}

static pbwtad *cpbwtk(size_t n, uint8_t *restrict c, pbwtad *restrict p,
                      size_t k) {
  static size_t *o = NULL;
  static size_t *h = NULL;

  if (!o)
    o = malloc(n * sizeof *o);
  if (!h)
    h = malloc(n * sizeof *h);

  pbwtad *ret = malloc(sizeof *ret);
  ret->a = malloc(n * sizeof *(ret->a));
  ret->d = malloc(n * sizeof *(ret->d));

  size_t r = 0, q = 0;
  size_t f = k, g = k;

  size_t i;
#if 0
  for (i = 0; i < n; i++) {
    /*printf("i: %6zu - p->a[i]: %zu\n", i, p->a[i]);*/
    if (c[p->a[i]] == 1) {
      o[q++] = p->a[i];
    } else {
      ret->a[r++] = p->a[i];
      /*z[r++] = p->a[i];*/
    }
  }
#else
  for (i = 0; i < n; i++) {
    size_t idx = p->a[i];
    size_t ddx = p->d[i];

    f = (ddx > f) ? ddx : f;
    g = (ddx > g) ? ddx : g;

    size_t mask = c[idx];
    o[q] = idx;
    ret->a[r] = idx;

    if (mask) {
      h[q] = g;
      g = 0;
    } else {
      ret->d[r] = f;
      f = 0;
    }
    q += mask;     // Increment q if mask is 1
    r += mask ^ 1; // Increment r if mask is 0
  }
#endif

  memcpy(ret->a + r, o, q * sizeof(size_t));
  memcpy(ret->d + r, h, q * sizeof(size_t));

  return ret;
}

#define TEST_LOG

#ifdef TEST_LOG
#define DPRINT(format, args...)                                                \
  do {                                                                         \
    fprintf(stderr, format, ##args);                                           \
  } while (0)

#define DPARR(n, a, fmt, ...)                                                  \
  do {                                                                         \
    __VA_OPT__(fprintf(stderr, __VA_ARGS__);)                                  \
    for (int parr_i__ = 0; parr_i__ < (n); parr_i__++) {                       \
      fprintf(stderr, (fmt), (a)[parr_i__]);                                   \
    }                                                                          \
    fputc(0xA, stderr);                                                        \
  } while (0)
#else
#define DPRINT(args...)
#define DPARR(args...)
#endif

pbwtad **linc(void *fin, size_t nrow, size_t ncol) {
  uint8_t *c0 = malloc(nrow * sizeof *c0);
  size_t *o = malloc(nrow * sizeof *o);
  size_t *z = malloc(nrow * sizeof *z);

  // NOTE: right now I don't know what I need, so I'm keeping
  // everything in memory, we'll see later
  pbwtad **pl = malloc(ncol * sizeof(pbwtad *));
  if (!pl)
    return NULL;

  pbwtad *p0 = pbwtad_new(nrow);
  for (int j = 0; j < nrow; j++) {
    p0->a[j] = j;
    p0->d[j] = 0;
  }
  fgetcoli(fin, 0, nrow, c0, ncol);
  // pl[0] = cpbwtk(nrow, c0, p0, 1);
  pl[0] = cpbwt(nrow, c0, p0);
  PDUMP(0, p0);
  PBWTAD_FREE(p0);

  for (size_t j = 1; j < ncol; j++) {
    fgetcoli(fin, j, nrow, c0, ncol);
    // pl[j] = cpbwtk(nrow, c0, pl[j - 1], j + 1);
    pl[j] = cpbwt(nrow, c0, pl[j - 1]);
    PDUMP(j, pl[j]);
  }

  FREE(c0);
  FREE(o);
  FREE(z);
  return pl;
}
pbwtad **blinc(void *fin, size_t nrow, size_t ncol) {
  uint8_t *c0 = malloc(nrow * sizeof *c0);
  // NOTE: right now I don't know what I need, so I'm keeping
  // everything in memory, we'll see later
  pbwtad **pl = malloc(ncol * sizeof(pbwtad *));
  if (!pl)
    return NULL;

  pbwtad *p0 = pbwtad_new(nrow);
  for (int j = 0; j < nrow; j++) {
    p0->a[j] = j;
    p0->d[j] = 0;
  }
  bfgetcoln(fin, nrow, c0, ncol);
  pl[0] = cpbwt(nrow, c0, p0);
  PDUMP(0, p0);
  PBWTAD_FREE(p0);

  for (size_t j = 1; j < ncol; j++) {
    bfgetcoln(fin, nrow, c0, ncol);
    pl[j] = cpbwt(nrow, c0, pl[j - 1]);
    PDUMP(j, pl[j]);
  }

  FREE(c0);
  return pl;
}
pbwtad **sblinc(int fin, size_t nrow, size_t ncol) {
  uint8_t *c0 = malloc(nrow * sizeof *c0);
  // NOTE: right now I don't know what I need, so I'm keeping
  // everything in memory, we'll see later
  pbwtad **pl = malloc(ncol * sizeof(pbwtad *));
  if (!pl)
    return NULL;

  pbwtad *p0 = pbwtad_new(nrow);
  for (int j = 0; j < nrow; j++) {
    p0->a[j] = j;
    p0->d[j] = 0;
  }
  sbfgetcoln(fin, nrow, c0, ncol);
  pl[0] = cpbwt(nrow, c0, p0);
  PDUMP(0, p0);
  PBWTAD_FREE(p0);

  for (size_t j = 1; j < ncol; j++) {
    sbfgetcoln(fin, nrow, c0, ncol);
    pl[j] = cpbwt(nrow, c0, pl[j - 1]);
    PDUMP(j, pl[j]);
  }

  FREE(c0);
  return pl;
}

pbwtad **mblinc(int fin, size_t nrow, size_t ncol) {
  uint8_t *c0 = malloc(nrow * sizeof *c0);
  // NOTE: right now I don't know what I need, so I'm keeping
  // everything in memory, we'll see later
  pbwtad **pl = malloc(ncol * sizeof(pbwtad *));
  if (!pl)
    return NULL;

  pbwtad *p0 = pbwtad_new(nrow);
  for (int j = 0; j < nrow; j++) {
    p0->a[j] = j;
    p0->d[j] = 0;
  }
  mbfgetcoln(fin, nrow, c0, ncol);
  pl[0] = cpbwt(nrow, c0, p0);
  PDUMP(0, p0);
  PBWTAD_FREE(p0);

  for (size_t j = 1; j < ncol; j++) {
    mbfgetcoln(fin, nrow, c0, ncol);
    pl[j] = cpbwt(nrow, c0, pl[j - 1]);
    PDUMP(j, pl[j]);
  }

  FREE(c0);
  return pl;
}

pbwtad **wapproxc_rrs(void *fin, size_t nrow, size_t ncol) {
  // NOTE: right now I don't know what I need, so I'm keeping
  // everything in memory, we'll see later
  pbwtad **pb = malloc(ncol * sizeof(pbwtad *));
  memset(pb, 0, ncol * sizeof(pbwtad *));

  // Compute the bit-packed windows
  uint64_t *pw = malloc(nrow * sizeof *pw);
  size_t *aux = malloc(nrow * sizeof *aux);

  pbwtad *ps = pbwtad_new(nrow);
  pb[W - 1] = ps;
  fgetcoliw64r(fin, 0, nrow, pw, ncol);
  rrsort0(nrow, pw, ps->a, aux);
  PDUMP(W - 1, ps);

  size_t j;
  for (j = 1; j * W <= ncol - W; j++) {
    pbwtad *ps = pbwtad_new(nrow);
    pb[W * (j + 1) - 1] = ps;
    memcpy(ps->a, pb[W * j - 1]->a, nrow * sizeof *(ps->a));
    fgetcoliw64r(fin, j, nrow, pw, ncol);
    rrsortx(nrow, pw, ps->a, aux);
    PDUMP(W * (j + 1) - 1, ps);
  }

  uint8_t *c0 = NULL;

#if defined(BF2IOMODE_BM)
  j *= W;
  fgetcolwgri(fin, j, nrow, pw, ncol, ncol - j);
#elif defined(BF2IOMODE_BCF)
  fgetcoliw64r(fin, j, nrow, pw, ncol);
  j *= W;
#else
#error UNDEFINED BEHAVIOUR
#endif

  // last column needs special handling, since it is < W
  ps = pbwtad_new(nrow);
  pb[ncol - 1] = ps;
  memcpy(ps->a, pb[j - 1]->a, nrow * sizeof *(ps->a));
  rrsortx(nrow, pw, ps->a, aux);
  PDUMP(ncol - 1, ps);

  FREE(c0);
  return pb;
}
pbwtad **wbapproxc_rrs(void *fin, size_t nrow, size_t ncol) {
  // NOTE: right now I don't know what I need, so I'm keeping
  // everything in memory, we'll see later
  pbwtad **pb = malloc(ncol * sizeof(pbwtad *));

  // Compute the bit-packed windows
  uint64_t *pw = malloc(nrow * sizeof *pw);
  size_t *aux = malloc(nrow * sizeof *aux);

  pbwtad *ps = pbwtad_new(nrow);
  pb[W - 1] = ps;
  bfgetcolw64rn(fin, nrow, pw, ncol);
  rrsort0(nrow, pw, ps->a, aux);
  PDUMP(W - 1, ps);

  size_t j;
  for (j = 1; j * W <= ncol - W; j++) {
    pbwtad *ps = pbwtad_new(nrow);
    pb[W * (j + 1) - 1] = ps;
    memcpy(ps->a, pb[W * j - 1]->a, nrow * sizeof *(ps->a));
    bfgetcolw64rn(fin, nrow, pw, ncol);
    rrsortx(nrow, pw, ps->a, aux);
    PDUMP(W * (j + 1) - 1, ps);
  }

  uint8_t *c0 = NULL;
  j *= W;
  fgetcolwgri(fin, j, nrow, pw, ncol, ncol - j);
  ps = pbwtad_new(nrow);
  pb[ncol - 1] = ps;
  memcpy(ps->a, pb[j - 1]->a, nrow * sizeof *(ps->a));
  rrsortx(nrow, pw, ps->a, aux);
  PDUMP(ncol - 1, ps);

  FREE(c0);
  return pb;
}
pbwtad **swbapproxc_rrs(int fin, size_t nrow, size_t ncol) {
  // NOTE: right now I don't know what I need, so I'm keeping
  // everything in memory, we'll see later
  pbwtad **pb = malloc(ncol * sizeof(pbwtad *));

  // Compute the bit-packed windows
  uint64_t *pw = malloc(nrow * sizeof *pw);
  size_t *aux = malloc(nrow * sizeof *aux);

  pbwtad *ps = pbwtad_new(nrow);
  pb[W - 1] = ps;
  sbfgetcolw64rn(fin, nrow, pw, ncol);
  rrsort0(nrow, pw, ps->a, aux);
  PDUMP(W - 1, ps);

  size_t j;
  for (j = 1; j * W <= ncol - W; j++) {
    pbwtad *ps = pbwtad_new(nrow);
    pb[W * (j + 1) - 1] = ps;
    memcpy(ps->a, pb[W * j - 1]->a, nrow * sizeof *(ps->a));
    sbfgetcolw64rn(fin, nrow, pw, ncol);
    rrsortx(nrow, pw, ps->a, aux);
    PDUMP(W * (j + 1) - 1, ps);
  }

  uint8_t *c0 = NULL;
  j *= W;
  sfgetcolwgri(fin, j, nrow, pw, ncol, ncol - j);
  ps = pbwtad_new(nrow);
  pb[ncol - 1] = ps;
  memcpy(ps->a, pb[j - 1]->a, nrow * sizeof *(ps->a));
  rrsortx(nrow, pw, ps->a, aux);
  PDUMP(ncol - 1, ps);

  FREE(c0);
  return pb;
}
pbwtad **mwbapproxc_rrs(int fin, size_t nrow, size_t ncol) {
  // NOTE: right now I don't know what I need, so I'm keeping
  // everything in memory, we'll see later
  pbwtad **pb = malloc(ncol * sizeof(pbwtad *));

  // Compute the bit-packed windows
  uint64_t *pw = malloc(nrow * sizeof *pw);
  size_t *aux = malloc(nrow * sizeof *aux);

  pbwtad *ps = pbwtad_new(nrow);
  pb[W - 1] = ps;
  sbfgetcolw64rn_mmap(fin, nrow, pw, ncol);
  rrsort0(nrow, pw, ps->a, aux);
  PDUMP(W - 1, ps);

  size_t j;
  for (j = 1; j * W <= ncol - W; j++) {
    pbwtad *ps = pbwtad_new(nrow);
    pb[W * (j + 1) - 1] = ps;
    memcpy(ps->a, pb[W * j - 1]->a, nrow * sizeof *(ps->a));
    sbfgetcolw64rn_mmap(fin, nrow, pw, ncol);
    rrsortx(nrow, pw, ps->a, aux);
    PDUMP(W * (j + 1) - 1, ps);
  }

  uint8_t *c0 = NULL;
  j *= W;
  sfgetcolwgri(fin, j, nrow, pw, ncol, ncol - j);
  ps = pbwtad_new(nrow);
  pb[ncol - 1] = ps;
  memcpy(ps->a, pb[j - 1]->a, nrow * sizeof *(ps->a));
  rrsortx(nrow, pw, ps->a, aux);
  PDUMP(ncol - 1, ps);
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
pbwtad **wapproxc_qs(FILE *fin, size_t nrow, size_t ncol) {
  // NOTE: right now I don't know what I need, so I'm keeping
  // everything in memory, we'll see later
  pbwtad **pb = malloc(ncol * sizeof(pbwtad *));

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
  j *= W;
  fgetcolwgri(fin, j, nrow, pw, ncol, ncol - j);
  ps = malloc(nrow * sizeof *ps);
  ps->a = malloc(nrow * sizeof *(ps->a));
  pb[ncol - 1] = ps;
  qsp = quadsort_u64_ix(pw, nrow, NULL);
  for (size_t i = 0; i < nrow; i++) {
    ps->a[i] = qsp[i] - pw0;
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

pbwtad **wparc_rrs(void *fin, size_t nrow, size_t ncol) {
  // NOTE: right now I don't know what I need, so I'm keeping
  // everything in memory, we'll see later

  pbwtad **pb = malloc(ncol * sizeof(pbwtad *));
  uint8_t *c0 = malloc(nrow * sizeof *c0);
  pbwtad *p0 = pbwtad_new(nrow);
  for (int j = 0; j < nrow; j++) {
    p0->a[j] = j;
    p0->d[j] = 0;
  }
#ifdef BF2IOMODE_BCF
  void *_tfin = fin;
  bcf_srs_t *_sr = bcf_sr_init();
  bcf_sr_add_reader(_sr, ((bcf_srs_t *)fin)->readers[0].fname);
  fin = _sr;
#endif

  fgetcoli(fin, 0, nrow, c0, ncol);
  pb[0] = cpbwt_0(nrow, c0, p0);
  PDUMP(0, pb[0]);
  PBWTAD_FREE(p0);

  for (int j = 1; j < W; j++) {
    fgetcoli(fin, j, nrow, c0, ncol);
    pb[j] = cpbwt_0(nrow, c0, pb[j - 1]);
    PDUMP(j, pb[j]);
  }

#ifdef BF2IOMODE_BCF
  bcf_sr_destroy(_sr);
  fin = _tfin;
#endif
  uint64_t *pw0 = malloc(nrow * sizeof *pw0);
  uint64_t *pw1 = malloc(nrow * sizeof *pw1);
  size_t *aux = malloc(nrow * sizeof *aux);
  fgetcoliw64r(fin, 0, nrow, pw0, ncol);

  size_t j;

  for (j = 1; j * W <= ncol - W; j++) {
    pbwtad *ps = pbwtad_new(nrow);
    pb[W * (j + 1) - 1] = ps;
    memcpy(ps->a, pb[W * j - 1]->a, nrow * sizeof *(ps->a));
    fgetcoliw64r(fin, j, nrow, pw1, ncol);
    rrsortx(nrow, pw1, ps->a, aux);

#pragma omp parallel for shared(pw1, pw0, pb, j)
    for (size_t x = 1; x < W; x++) {
      uint64_t *w = malloc(nrow * sizeof *w);
      size_t J = W * (j + 1) - 1;

      wr64mrgsi(nrow, pw1, pw0, w, x);
      pbwtad *ps = pbwtad_new(nrow);
      pb[J - x] = ps;
      memcpy(ps->a, pb[J - W - x]->a, nrow * sizeof *(ps->a));
      rrsortx_noaux(nrow, w, ps->a);
      FREE(w);
    }
    PDUMP_SEQ(W * j - 1, W * (j + 1) - 1, pb);
    SWAP(pw0, pw1);
  }

  for (j = j * W; j < ncol; j++) {
#if defined(BF2IOMODE_BM)
    fgetcoli(fin, j, nrow, c0, ncol);
#elif defined(BF2IOMODE_BCF)
    fgetcoli(fin, j, nrow, c0, 0);
#else
#error UNDEFINED BEHAVIOUR
#endif
    pb[j] = cpbwt_0(nrow, c0, pb[j - 1]);
    PDUMP(j, pb[j]);
  }

  FREE(pw0);
  FREE(pw1);
  FREE(aux);
  FREE(c0);
  return pb;
}
pbwtad **bwparc_rrs(void *fin, size_t nrow, size_t ncol) {
  // NOTE: right now I don't know what I need, so I'm keeping
  // everything in memory, we'll see later
  pbwtad **pb = malloc(ncol * sizeof(pbwtad *));

  uint8_t *c0 = malloc(nrow * sizeof *c0);
  pbwtad *p0 = pbwtad_new(nrow);
  for (int j = 0; j < nrow; j++) {
    p0->a[j] = j;
    p0->d[j] = 0;
  }
  fgetcoli(fin, 0, nrow, c0, ncol);
  pb[0] = cpbwt_0(nrow, c0, p0);
  PDUMP(0, pb[0]);
  PBWTAD_FREE(p0);

  for (int j = 1; j < W; j++) {
    fgetcoli(fin, j, nrow, c0, ncol);
    pb[j] = cpbwt_0(nrow, c0, pb[j - 1]);
    PDUMP(j, pb[j]);
  }

  uint64_t *pw0 = malloc(nrow * sizeof *pw0);
  uint64_t *pw1 = malloc(nrow * sizeof *pw1);
  size_t *aux = malloc(nrow * sizeof *aux);
  bfgetcolw64rn(fin, nrow, pw0, ncol);

  size_t j;

  for (j = 1; j * W <= ncol - W; j++) {
    pbwtad *ps = pbwtad_new(nrow);
    pb[W * (j + 1) - 1] = ps;
    memcpy(ps->a, pb[W * j - 1]->a, nrow * sizeof *(ps->a));
    bfgetcolw64rn(fin, nrow, pw1, ncol);
    rrsortx(nrow, pw1, ps->a, aux);

#pragma omp parallel for shared(pw1, pw0, pb, j)
    for (size_t x = 1; x < W; x++) {
      uint64_t *w = malloc(nrow * sizeof *w);
      size_t J = W * (j + 1) - 1;

      wr64mrgsi(nrow, pw1, pw0, w, x);
      pbwtad *ps = pbwtad_new(nrow);
      pb[J - x] = ps;
      memcpy(ps->a, pb[J - W - x]->a, nrow * sizeof *(ps->a));
      rrsortx_noaux(nrow, w, ps->a);

      FREE(w);
    }

    PDUMP_SEQ(W * j - 1, W * (j + 1) - 1, pb);
    SWAP(pw0, pw1);
  }

  for (j = j * W; j < ncol; j++) {
    fgetcoli(fin, j, nrow, c0, ncol);
    pb[j] = cpbwt_0(nrow, c0, pb[j - 1]);
    PDUMP(j, pb[j]);
  }

  FREE(pw0);
  FREE(pw1);
  FREE(aux);
  FREE(c0);
  return pb;
}

pbwtad **wstagparc_rrs(char *fpath, size_t nrow, size_t ncol) {
  // NOTE: right now I don't know what I need, so I'm keeping
  // everything in memory, we'll see later
  pbwtad **pb = malloc(ncol * sizeof(pbwtad *));

  // first W must be computed linearly
  FILE *fin = fopen(fpath, "r");
  if (!fin) {
    perror("[spr]");
    exit(32);
  }
  uint8_t *c0 = malloc(nrow * sizeof *c0);
  pbwtad *p0 = pbwtad_new(nrow);
  for (int j = 0; j < nrow; j++) {
    p0->a[j] = j;
    p0->d[j] = 0;
  }
  fgetcoli(fin, 0, nrow, c0, ncol);
  pb[0] = cpbwt_0(nrow, c0, p0);
  PDUMP(0, pb[0]);
  PBWTAD_FREE(p0);
  /*printf("[man] c:0 (<1..N)\n");*/

  for (size_t j = 1; j < W; j++) {
    fgetcoli(fin, j, nrow, c0, ncol);
    pb[j] = cpbwt_0(nrow, c0, pb[j - 1]);
    /*printf("[lin] c:%zu (<%zu)\n", j, j - 1);*/
  }

// This works differently from previous versions of windows.
// here I am computing (j+w) using value in j
#pragma omp parallel
  {
    FILE *fin = fopen(fpath, "r");
    if (!fin) {
      perror("[spr]");
      exit(33);
    }
    int tid = omp_get_thread_num();
    int nthreads = omp_get_num_threads();

    int base = W / nthreads;
    int rem = W % nthreads;

    int start = tid * base + (tid < rem ? tid : rem);
    int count = base + (tid < rem ? 1 : 0);

    for (int offset = 0; offset < count; offset++) {
      int lane = start + offset;

      for (size_t j = lane; j + W < ncol; j += W) {

        uint64_t *pw = malloc(nrow * sizeof *pw);
        pbwtad *ps = pbwtad_new(nrow);
        pb[j + W] = ps;
        memcpy(ps->a, pb[j]->a, nrow * sizeof *(ps->a));

        fgetcolwgri(fin, j + 1, nrow, pw, ncol, W);
        // maybe change to aux for each thread
        rrsortx_noaux(nrow, pw, ps->a);

#ifdef DBDUMP
#pragma omp critical
        { PDUMP(j + W, ps); }
#endif

        FREE(pw);
      }
    }
  }

  FREE(c0);
  return pb;
}

pbwtad **wseq_rrs(FILE *fin, size_t nrow, size_t ncol) {
  // NOTE: right now I don't know what I need, so I'm keeping
  // everything in memory, we'll see later
  pbwtad **pb = malloc(ncol * sizeof(pbwtad *));

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

  uint64_t *pw0 = malloc(nrow * sizeof *pw0);
  uint64_t *pw1 = malloc(nrow * sizeof *pw1);
  size_t *aux = malloc(nrow * sizeof *aux);
  fgetcoliw64r(fin, 0, nrow, pw0, ncol);

  size_t j;

  /*for (j = 1; j * W <= W * 2; j++) {*/
  for (j = 1; j * W <= ncol - W; j++) {
    fprintf(stderr, "\r%10zu/%zu", (j * W) + 1, ncol);
    pbwtad *ps = malloc(nrow * sizeof *ps);
    ps->a = malloc(nrow * sizeof *(ps->a));
    pb[W * (j + 1) - 1] = ps;
    memcpy(ps->a, pb[W * j - 1]->a, nrow * sizeof *(ps->a));
    fgetcoliw64r(fin, j, nrow, pw1, ncol);
    rrsortx(nrow, pw1, ps->a, aux);

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

  FREE(pw0);
  FREE(pw1);
  FREE(aux);
  FREE(c0);
  return pb;
}

int main(int argc, char *argv[]) {
  char _usage_args_[] = "[lin|bli[s|m]|ars|aqs|bar[s|m]|prs|bpr|spr] FILE\n";
  if (argc < 2) {
    fprintf(stderr, "Usage: %s %s FILE\n", argv[0], _usage_args_);
    return EXIT_FAILURE;
  }

#if defined(BF2IOMODE_BM)
  FILE *fin = fopen(argv[2], "r");
  int fd = open(argv[2], O_RDONLY);
  if (!fin) {
    perror("[main]");
    return EXIT_FAILURE;
  }
#elif defined(BF2IOMODE_BCF)
  int fd = -1;
  void *fin = NULL;

#if 0
  htsFile *fp = hts_open(argv[2], "rb");
  bcf_hdr_t *hdr = bcf_hdr_read(fp);
  bcf1_t *rec = bcf_init();
  size_t n = 0;
  while (bcf_read(fp, hdr, rec) >= 0) {
    bcf_unpack(rec, BCF_UN_ALL); // NOTE: ?
    int32_t *gt_arr = NULL, ngt_arr = 0;
    int i, j, ngt, nsmpl = bcf_hdr_nsamples(hdr);
    ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
    int max_ploidy = ngt / nsmpl;
    // printf("\nmax_ploidy = %d\n", max_ploidy);
    // parr(ngt, gt_arr, "%d ");
    uint8_t *col = malloc(nsmpl * 2 * sizeof *col);
    size_t icol = 0;
    for (i = 0; i < nsmpl; i++) {
      int32_t *ptr = gt_arr + i * max_ploidy;
      // printf("  %d", *ptr);
      for (j = 0; j < max_ploidy; j++) {
        // if true, the sample has smaller ploidy
        // If a sample has less genotypes than max_ploidy,
        // the "vector" retains the size of max_ploidy, but
        // missing ploid are filled with `bcf_int32_vector_end` macro
        if (ptr[j] == bcf_int32_vector_end)
          break;

        // missing allele
        if (bcf_gt_is_missing(ptr[j]))
          exit(-1);

        // the VCF 0-based allele index
        col[icol] = bcf_gt_allele(ptr[j]);
        icol++;
      }
    }
    printf("%zu: %d\n", n, col[1]);
    n++;
    // putchar(0xA);
    // parr(nsmpl * 2, col, "%d ");
  }
#else
  bcf_srs_t *sr = bcf_sr_init();
  bcf_sr_add_reader(sr, argv[2]);

#endif

  fin = sr;

#else
#error BF2IOMODE is not specified
#endif

  size_t nrow, ncol;
#ifdef BF2IOMODE_BCF
  bcf_srs_t *_sr = bcf_sr_init();
  bcf_sr_add_reader(_sr, argv[2]);
  fgetrc(_sr, &nrow, &ncol);
  // WARN: seeking does not work for unknown reasons.
  // Therefore here I need to duplicate the input.
  // bcf_sr_seek(sr, NULL, 0);
  bcf_sr_destroy(_sr);
#else
  fgetrc(fin, &nrow, &ncol);
#endif
  DPRINT("[%s] row: %5zu, col: %5zu\n", __func__, nrow, ncol);
  // uint8_t *cc = malloc(nrow * sizeof *cc);
  // fgetcoli(fin, 0, nrow, cc, 0);
  // parr(nrow, cc, "%d ");
  // fgetcoli(fin, 1, nrow, cc, 0);
  // parr(nrow, cc, "%d ");
  // fgetcoli(fin, 2, nrow, cc, 0);
  // parr(nrow, cc, "%d ");
  // fgetcoli(fin, 3, nrow, cc, 0);
  // parr(nrow, cc, "%d ");
  // fgetcoli(fin, 4, nrow, cc, 0);
  // parr(nrow, cc, "%d ");
  // fgetcoli(fin, 5, nrow, cc, 0);
  // parr(nrow, cc, "%d ");
  // fgetcoli(fin, 7, nrow, cc, 0);

  // bcf_hdr_t *hdr = sr->readers[0].header;
  // int nsmpl = bcf_hdr_nsamples(hdr);
  // printf("nsmpl: %d\n", nsmpl);

  // uint8_t *col = malloc(nsmpl * 2 * sizeof *col); // NOTE: assume diploid
  // size_t n = 0;
  // while (bcf_sr_next_line(sr)) {
  //   bcf1_t *line = bcf_sr_get_line(sr, 0);
  //   int32_t *gt_arr = NULL, ngt_arr = 0;
  //   int ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);
  //
  //   size_t icol = 0;
  //   for (size_t i = 0; i < nsmpl; i++) {
  //     int32_t *ptr = gt_arr + i * 2;
  //     // hap 1
  //     // if (ptr[0] == bcf_int32_vector_end)
  //     //   exit(-2);
  //     // if (bcf_gt_is_missing(ptr[0]))
  //     //   exit(-1);
  //     col[2 * i] = bcf_gt_allele(ptr[0]);
  //
  //     // hap 2
  //     // if (ptr[1] == bcf_int32_vector_end)
  //     //   exit(-2);
  //     // if (bcf_gt_is_missing(ptr[1]))
  //     //   exit(-1);
  //     col[2 * i + 1] = bcf_gt_allele(ptr[1]);
  //   }
  //   // parr(nsmpl * 2, col, "%d ");
  //   printf("%zu: %d\n", n, col[1]);
  //   n++;
  // }
  pbwtad **r;

  if (argc > 3) {
    if (strcmp(argv[3], "DUMP") == 0) {
      DO_DUMP = 1;
    }
  }
  if (strcmp(argv[1], "lin") == 0) {
    r = linc(fin, nrow, ncol);
  } else if (strcmp(argv[1], "bli") == 0) {
    r = blinc(fin, nrow, ncol);
  } else if (strcmp(argv[1], "blis") == 0) {
    r = sblinc(fd, nrow, ncol);
  } else if (strcmp(argv[1], "blim") == 0) {
    r = mblinc(fd, nrow, ncol);
  } else if (strcmp(argv[1], "ars") == 0) {
    r = wapproxc_rrs(fin, nrow, ncol);
  } else if (strcmp(argv[1], "aqs") == 0) {
    fprintf(stderr, "\e[0;33mWARNING: this is not been properly tested and "
                    "might not work.\e[0m\n");
    r = wapproxc_qs(fin, nrow, ncol);
  } else if (strcmp(argv[1], "bar") == 0) {
    r = wbapproxc_rrs(fin, nrow, ncol);
  } else if (strcmp(argv[1], "bars") == 0) {
    r = swbapproxc_rrs(fd, nrow, ncol);
  } else if (strcmp(argv[1], "barm") == 0) {
    r = mwbapproxc_rrs(fd, nrow, ncol);
  } else if (strcmp(argv[1], "prs") == 0) {
    r = wparc_rrs(fin, nrow, ncol);
  } else if (strcmp(argv[1], "bpr") == 0) {
    r = bwparc_rrs(fin, nrow, ncol);
  } else if (strcmp(argv[1], "spr") == 0) {
    r = wstagparc_rrs(argv[2], nrow, ncol);
  } else if (strcmp(argv[1], "srs") == 0) {
    r = wseq_rrs(fin, nrow, ncol);
  } else {
    fprintf(stderr, "Usage: %s %s FILE\n", argv[0], _usage_args_);
    return EXIT_FAILURE;
  }

  /*fclose(fin);*/

  // if (r != NULL) {
  //   for (size_t i = 0; i < ncol; i++) {
  //     if (r[i]) {
  //       PBWTAD_FREE(r[i]);
  //     }
  //   }
  //   FREE(r);
  // }

  return EXIT_SUCCESS;
}

int maint(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s [lin|ars|aqs|prs] FILE\n", argv[0]);
    return EXIT_FAILURE;
  }

  FILE *fin = fopen(argv[2], "r");
  if (!fin) {
    perror("[main]");
    return EXIT_FAILURE;
  }
  int ifin = open(argv[2], O_RDONLY);

  size_t nrow, ncol;
  fgetrc(fin, &nrow, &ncol);
  DPRINT("[%s] row: %5zu, col: %5zu\n", __func__, nrow, ncol);

  /*uint8_t *c0 = malloc(nrow * sizeof *c0);*/
  /*sbfgetcoln(ifin, nrow, c0, ncol);*/
  /*parr(nrow, c0, "%d ");*/
  /*fgetcoli(fin, 0, nrow, c0, ncol);*/
  /*parr(nrow, c0, "%d ");*/
  /*sbfgetcoln(ifin, nrow, c0, ncol);*/
  /*parr(nrow, c0, "%d ");*/
  /*fgetcoli(fin, 1, nrow, c0, ncol);*/
  /*parr(nrow, c0, "%d ");*/
  /*sbfgetcoln(ifin, nrow, c0, ncol);*/
  /*parr(nrow, c0, "%d ");*/
  /*fgetcoli(fin, 2, nrow, c0, ncol);*/
  /*parr(nrow, c0, "%d ");*/
  /*sbfgetcoln(ifin, nrow, c0, ncol);*/
  /*parr(nrow, c0, "%d ");*/
  /*fgetcoli(fin, 3, nrow, c0, ncol);*/
  /*parr(nrow, c0, "%d ");*/
  /*sbfgetcoln(ifin, nrow, c0, ncol);*/
  /*parr(nrow, c0, "%d ");*/
  /*fgetcoli(fin, 4, nrow, c0, ncol);*/
  /*parr(nrow, c0, "%d ");*/
  /*sbfgetcoln(ifin, nrow, c0, ncol);*/
  /*parr(nrow, c0, "%d ");*/
  /*fgetcoli(fin, 5, nrow, c0, ncol);*/
  /*parr(nrow, c0, "%d ");*/

  uint64_t *w0 = malloc(nrow * sizeof *w0);
  sbfgetcolw64rn(ifin, nrow, w0, ncol);
  parr(nrow, w0, "%llu ");
  fgetcoliw64r(fin, 0, nrow, w0, ncol);
  parr(nrow, w0, "%llu ");
  sbfgetcolw64rn(ifin, nrow, w0, ncol);
  parr(nrow, w0, "%llu ");
  fgetcoliw64r(fin, 1, nrow, w0, ncol);
  parr(nrow, w0, "%llu ");
  sbfgetcolw64rn(ifin, nrow, w0, ncol);
  parr(nrow, w0, "%llu ");
  fgetcoliw64r(fin, 2, nrow, w0, ncol);
  parr(nrow, w0, "%llu ");
  return EXIT_SUCCESS;
}
