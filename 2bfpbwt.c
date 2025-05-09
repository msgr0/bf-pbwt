#include "io.h"
#include "tracing.h"
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
#define PDUMPR(i, p)                                                           \
  do {                                                                         \
    if (DO_DUMP) {                                                             \
      printf("%zu:", (size_t)(i));                                             \
      size_t pdump_j__;                                                        \
      for (pdump_j__ = 0; pdump_j__ < nrow - 1; pdump_j__++)                   \
        printf("%zu ", (p)->a[pdump_j__]);                                     \
      printf("%zu", (p)->a[pdump_j__]);                                        \
      fputc('|', stdout);                                                      \
      for (pdump_j__ = 0; pdump_j__ < nrow - 1; pdump_j__++)                   \
        printf("%zu ", 1 + (i) - (p)->d[pdump_j__]);                           \
      printf("%zu", 1 + (i) - (p)->d[pdump_j__]);                              \
      fputc(0xA, stdout);                                                      \
    }                                                                          \
  } while (0)
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
#define PDUMP_SEQ_OFFSET(s, e, p, offset)                                      \
  do {                                                                         \
    for (size_t pdump_ix__ = (s); pdump_ix__ < (e); pdump_ix__++) {            \
      if (DO_DUMP) {                                                           \
        printf("%zu:", (size_t)(offset) + (size_t)(pdump_ix__));               \
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

typedef struct pbwtad pbwtad;
struct pbwtad {
  size_t *a;
  size_t *d;
};

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
void reversec(pbwtad *p, pbwtad *rev, size_t n) {
  for (size_t i = 0; i < n; i++) {
    rev->a[p->a[i]] = i;
    rev->d[p->a[i]] = p->d[i]; // == p->d[rev->a[i]]
    assert(rev->d[p->a[i]] == p->d[rev->a[p->a[i]]]);
  }
}

size_t recover_div(size_t n, size_t w, size_t i, size_t i0, uint64_t *c,
                   pbwtad *p, pbwtad *ppr, pbwtad *prev, pbwtad *pprrev) {
  size_t d;
  // printf("recovering...");
  //   printf("pprA:");
  // for (size_t _n = 0; _n < n; _n++) {
  //   printf("%zu,", ppr->a[_n]);
  // }
  // printf(" -- pprD:");
  // for (size_t _n = 0; _n < n; _n++) {
  //   printf("%zu,", ppr->d[_n]);
  // }
  //   printf("||  ppprevA:");
  // for (size_t _n = 0; _n < n; _n++) {
  //   printf("%zu,", pprrev->a[_n]);
  // }
  // printf(" -- pprrevD:");
  // for (size_t _n = 0; _n < n; _n++) {
  //   printf("%zu,", pprrev->d[_n]);
  // }
  // printf(" -- i-1: %zu, i: %zu", i0, i);
  //
  // printf("\n");

  // if (pprrev->a[p->a[i]] == pprrev->a[p->a[i-1]] +1 ) { // if the current
  // a[i] and a[i-1] values are one after the other in the previous a (ppr), do
  // ... ;  ( so pprrev->a[p->a[i]] )
  //   d = w + prev->d[p->a[i]];
  // }
  // else {
  size_t min = ppr->d[pprrev->a[i]];
  // printf("ai0: %zu, ai1: %zu", pprrev->a[p->a[i-1]], pprrev->a[p->a[i]]);
  // devo scorrere
  for (size_t j = (pprrev->a[i0]) + 1; j < (pprrev->a[i]); j++) {
    // printf("j: %zu\n", j);
    // printf("iterating j -- ");
    if (ppr->d[j] < min) {
      // printf("updating min... :%zu \n", min);
      min = ppr->d[j];
      // printf("min at j(%zu): %zu", j, min);
    }
    // printf("\n");
  }
  d = w + min;
  // }
  // find min in range
  // printf("d: %zu\n", d);
  return d;
}

void divc0(size_t n, uint64_t *c, pbwtad *p) {
  uint64_t x = 0;
  size_t div = 0;
  p->d[0] = 0;
  for (size_t i = 1; i < n; i++) {
    x = c[p->a[i]] ^ c[p->a[i - 1]];
    p->d[i] = x ? __builtin_clzll(x) : 64;
  }
}

void divc_last(size_t n, uint64_t *c, pbwtad *p, pbwtad *ppr, pbwtad *prev,
               pbwtad *pprrev, uint64_t *_x, size_t w) {
  // c contains 64bit-encoded ints
  // xor of each c[s[i]] and its preceeding;
  // x[0] contains no information, previous x information is discarded;
  // here 64 is the size of the window
  static int8_t kk = 0;
  uint64_t x = 0;
  // size_t w = 64;
  size_t div;
  p->d[0] = 0;

#if 1
  for (size_t i = 1; i < n; i++) {
    x = c[p->a[i]] ^ c[p->a[i - 1]];
    div = x ? __builtin_clzll(x) - w : w;
    // if (div > w) div = 64 - div;
    p->d[i] = (div >= w) ? recover_div(n, w, p->a[i], p->a[i - 1], c, p, ppr,
                                       prev, pprrev)
                         : div;
  }

#else
#if 0 
  for (size_t i = 1; i < n; i++) {
    x = c[p->a[i]] ^ c[p->a[i - 1]];
    div = x ? __builtin_clzll(x) : 64; // if x == 0, it means that are completly
                                       // equal --> div = 64. otherwise ctzll
    p->d[i] = (div == 64) ? (prev->d[p->a[i]] + div)
                          : div; // BUG: when div == 64, it is not always true
                                 // that (( div = (p_rev->d[p->a[i]] + div) ))
                                 // but div could be less than above value
  }
#else
  // reading pbwtPrev
  //
  size_t bi = 0;        // Backward Index;
  size_t pbi = 0;       // Previous Backward Index;
  size_t mbdiv = 65535; // Minimum Backward DIVergence;
                        // Stack di minimi;
  // GOAL: mbdiv[#<64]
  // mbdiv
  //
  for (size_t i = 0; i < n; i++) {

    bi = prev->a[ppr->a[i]];

    if (bi == 0) {
      // p->d[0] = 0;
      div = 0;
    } else {
      x = (c[p->a[bi]] ^ c[p->a[bi - 1]]);
      div = (x ? __builtin_clzll(x)
               : 64); // quando div < 64 --> ho una nuova finestra; # < 64
    }
    // printf("%zu, %zu, %zu, +++++++\n", i, bi, div);

    if (div < 64) { // mismatch, nuovo carattere
      p->d[bi] = div;
      // mbdiv = i ?  ((ppr->d[i] < mbdiv) ? ppr->d[i] : mbdiv) : mbdiv;

    } else { // match, stesso carattere.
      assert(div == 64);
      if (bi - 1 == pbi) {
        // mbdiv = i ?  ((ppr->d[i] < mbdiv) ? ppr->d[i] : mbdiv) : mbdiv;
        p->d[bi] = div + ppr->d[i];
      } else {
        mbdiv = i ? ((ppr->d[i] < mbdiv) ? ppr->d[i] : mbdiv) : mbdiv;
        p->d[bi] = div + mbdiv;
        mbdiv = 65535;
      }
    }
    pbi = bi;
  }
#endif
#endif
  kk += W;
}

void divc(size_t n, uint64_t *c, pbwtad *p, pbwtad *ppr, pbwtad *prev,
          pbwtad *pprrev, uint64_t *_x) {
  // c contains 64bit-encoded ints
  // xor of each c[s[i]] and its preceeding;
  // x[0] contains no information, previous x information is discarded;
  // here 64 is the size of the window
  static int8_t kk = 0;
  uint64_t x = 0;
  size_t w = 64;
  size_t div;
  p->d[0] = 0;

#if 1
  for (size_t i = 1; i < n; i++) {
    x = c[p->a[i]] ^ c[p->a[i - 1]];
    div = x ? __builtin_clzll(x) : w;
    p->d[i] = (div == w) ? recover_div(n, w, p->a[i], p->a[i - 1], c, p, ppr,
                                       prev, pprrev)
                         : div;
  }

#else
#if 0 
  for (size_t i = 1; i < n; i++) {
    x = c[p->a[i]] ^ c[p->a[i - 1]];
    div = x ? __builtin_clzll(x) : 64; // if x == 0, it means that are completly
                                       // equal --> div = 64. otherwise ctzll
    p->d[i] = (div == 64) ? (prev->d[p->a[i]] + div)
                          : div; // BUG: when div == 64, it is not always true
                                 // that (( div = (p_rev->d[p->a[i]] + div) ))
                                 // but div could be less than above value
  }
#else
  // reading pbwtPrev
  //
  size_t bi = 0;        // Backward Index;
  size_t pbi = 0;       // Previous Backward Index;
  size_t mbdiv = 65535; // Minimum Backward DIVergence;
                        // Stack di minimi;
  // GOAL: mbdiv[#<64]
  // mbdiv
  //
  for (size_t i = 0; i < n; i++) {

    bi = prev->a[ppr->a[i]];

    if (bi == 0) {
      // p->d[0] = 0;
      div = 0;
    } else {
      x = (c[p->a[bi]] ^ c[p->a[bi - 1]]);
      div = (x ? __builtin_clzll(x)
               : 64); // quando div < 64 --> ho una nuova finestra; # < 64
    }
    // printf("%zu, %zu, %zu, +++++++\n", i, bi, div);

    if (div < 64) { // mismatch, nuovo carattere
      p->d[bi] = div;
      // mbdiv = i ?  ((ppr->d[i] < mbdiv) ? ppr->d[i] : mbdiv) : mbdiv;

    } else { // match, stesso carattere.
      assert(div == 64);
      if (bi - 1 == pbi) {
        // mbdiv = i ?  ((ppr->d[i] < mbdiv) ? ppr->d[i] : mbdiv) : mbdiv;
        p->d[bi] = div + ppr->d[i];
      } else {
        mbdiv = i ? ((ppr->d[i] < mbdiv) ? ppr->d[i] : mbdiv) : mbdiv;
        p->d[bi] = div + mbdiv;
        mbdiv = 65535;
      }
    }
    pbi = bi;
  }
#endif
#endif
  kk += W;
}

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

static int cpbwti(size_t n, uint8_t *restrict c, pbwtad *restrict pp,
                  pbwtad *restrict pc) {
  static size_t *o = NULL;
  static size_t *h = NULL;
  static size_t k = 1;

  if (!o)
    o = malloc(n * sizeof *o);
  if (!h)
    h = malloc(n * sizeof *h);

  size_t r = 0, q = 0;
  size_t f = k, g = k;

  size_t i;
  for (i = 0; i < n; i++) {
    size_t idx = pp->a[i];
    size_t ddx = pp->d[i];

    f = (ddx > f) ? ddx : f;
    g = (ddx > g) ? ddx : g;

    size_t mask = c[idx];
    o[q] = idx;
    pc->a[r] = idx;
    h[q] = g;
    pc->d[r] = f;

    f &= -mask;       // f = 0 if mask == 0, unchanged if mask == 1
    g &= -(1 - mask); // g = 0 if mask == 1, unchanged if mask == 0
    q += mask;        // Increment q if mask is 1
    r += mask ^ 1;    // Increment r if mask is 0
  }

  memcpy(pc->a + r, o, q * sizeof(size_t));
  memcpy(pc->d + r, h, q * sizeof(size_t));

  k++;
  return 1;
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

  pbwtad *p0 = pbwtad_new(nrow);
  pbwtad *p1 = pbwtad_new(nrow);
  for (int j = 0; j < nrow; j++) {
    p0->a[j] = j;
    p0->d[j] = 0;
  }
  fgetcoli(fin, 0, nrow, c0, ncol);
  // pl[0] = cpbwtk(nrow, c0, p0, 1);
  cpbwti(nrow, c0, p0, p1);
  PDUMP(0, p1);
  SWAP(p0, p1);

#if defined(BF2IOMODE_BM) || defined(BF2IOMODE_ENC)
  for (size_t j = 1; j < ncol;) {
    fgetcoli(fin, j, nrow, c0, ncol);
#elif defined(BF2IOMODE_BCF)
  size_t j = 1;
  while (fgetcoli(fin, j, nrow, c0, 1)) {
#else
#error UNDEFINED BEHAVIOUR
#endif
    // for (size_t j = 1; j < ncol; j++) {
    // fgetcoli(fin, j, nrow, c0, ncol);
    cpbwti(nrow, c0, p0, p1);
    PDUMP(j, p1);
    SWAP(p0, p1);
    j++;
  }

  PBWTAD_FREE(p0);
  PBWTAD_FREE(p1);
  FREE(c0);
  return NULL;
}
pbwtad **blinc(void *fin, size_t nrow, size_t ncol) {
  uint8_t *c0 = malloc(nrow * sizeof *c0);

  pbwtad *p0 = pbwtad_new(nrow);
  pbwtad *p1 = pbwtad_new(nrow);
  for (int j = 0; j < nrow; j++) {
    p0->a[j] = j;
    p0->d[j] = 0;
  }
  bfgetcoln(fin, nrow, c0, ncol);
  cpbwti(nrow, c0, p0, p1);
  PDUMP(0, p1);
  SWAP(p0, p1);

  for (size_t j = 1; j < ncol; j++) {
    bfgetcoln(fin, nrow, c0, ncol);
    cpbwti(nrow, c0, p0, p1);
    PDUMP(j, p1);
    SWAP(p0, p1);
  }

  PBWTAD_FREE(p0);
  PBWTAD_FREE(p1);
  FREE(c0);
  return NULL;
}
pbwtad **sblinc(int fin, size_t nrow, size_t ncol) {
  uint8_t *c0 = malloc(nrow * sizeof *c0);

  pbwtad *p0 = pbwtad_new(nrow);
  pbwtad *p1 = pbwtad_new(nrow);
  for (int j = 0; j < nrow; j++) {
    p0->a[j] = j;
    p0->d[j] = 0;
  }
  sbfgetcoln(fin, nrow, c0, ncol);
  cpbwti(nrow, c0, p0, p1);
  PDUMP(0, p1);
  SWAP(p0, p1);

  for (size_t j = 1; j < ncol; j++) {
    sbfgetcoln(fin, nrow, c0, ncol);
    cpbwti(nrow, c0, p0, p1);
    PDUMP(j, p1);
    SWAP(p0, p1);
  }

  PBWTAD_FREE(p0);
  PBWTAD_FREE(p1);
  FREE(c0);
  return NULL;
}

pbwtad **mblinc(int fin, size_t nrow, size_t ncol) {
  uint8_t *c0 = malloc(nrow * sizeof *c0);

  pbwtad *p0 = pbwtad_new(nrow);
  pbwtad *p1 = pbwtad_new(nrow);
  for (int j = 0; j < nrow; j++) {
    p0->a[j] = j;
    p0->d[j] = 0;
  }
  mbfgetcoln(fin, nrow, c0, ncol);
  cpbwti(nrow, c0, p0, p1);
  PDUMP(0, p1);
  SWAP(p0, p1);

  for (size_t j = 1; j < ncol; j++) {
    mbfgetcoln(fin, nrow, c0, ncol);
    cpbwti(nrow, c0, p0, p1);
    PDUMP(j, p1);
    SWAP(p0, p1);
  }

  PBWTAD_FREE(p0);
  PBWTAD_FREE(p1);
  FREE(c0);
  return NULL;
}

pbwtad **wapproxc_rrs(void *fin, size_t nrow, size_t ncol) {
  // Compute the bit-packed windows
  uint64_t *w64 = malloc(nrow * sizeof *w64);
  uint64_t *xor = malloc(nrow * sizeof *xor);
  size_t *aux = malloc(nrow * sizeof *aux);

  pbwtad *pbwt = pbwtad_new(nrow);
  pbwtad *pbwtPr = pbwtad_new(nrow);
  pbwtad *pbwtRev = pbwtad_new(nrow);
  pbwtad *pbwtPrRev = pbwtad_new(nrow);
  // pw is the actual windows data

  fgetcoliw64r(fin, 0, nrow, w64, ncol);
  rrsort0(nrow, w64, pbwt->a, aux);
  memcpy(pbwtPrRev->a, pbwtRev->a, nrow * sizeof *(pbwtRev->a));
  memcpy(pbwtPrRev->d, pbwtRev->d, nrow * sizeof *(pbwtRev->d));
  reversec(pbwt, pbwtRev, nrow);
  divc0(nrow, w64, pbwt);
  // divc(nrow, w64, pbwtPr, pbwt, pbwtRev, xor);

  PDUMPR(W - 1, pbwt);
  // memcpy(pbwtPr->a , pbwt->a, nrow * sizeof *(pbwt->a));
  // memcpy(pbwtPr->d , pbwt->d, nrow * sizeof *(pbwt->a));
  // printf("--------fine prima it\n");
  // SWAP(p0, p1);
  size_t j;
  size_t k = 1;
#if defined(BF2IOMODE_BM) || defined(BF2IOMODE_ENC)
  for (j = 1; j * W <= ncol - W;) {
    // memcpy(p1->a, p0->a, nrow * sizeof *(p1->a));
    fgetcoliw64r(fin, j, nrow, w64, ncol);
#elif defined(BF2IOMODE_BCF)
  j = 1;
  ncol = W;
  size_t _ncol = 0;
  while ((_ncol = fgetcoliw64r(fin, j, nrow, w64, 0)) == W) {
    ncol += _ncol;
#else
#error UNDEFINED BEHAVIOUR
#endif
    memcpy(pbwtPr->a, pbwt->a, nrow * sizeof *(pbwt->a));
    memcpy(pbwtPr->d, pbwt->d, nrow * sizeof *(pbwt->d));
    rrsortx(nrow, w64, pbwt->a, aux); // BUG: pbwtPr->d doesnt contain anything.
    // FIXME: moved reversec after rrsortx, could now be integrated in rrsortx
    // as #1 pull request;
    memcpy(pbwtPrRev->a, pbwtRev->a, nrow * sizeof *(pbwtRev->a));
    memcpy(pbwtPrRev->d, pbwtRev->d, nrow * sizeof *(pbwtRev->d));
    reversec(pbwt, pbwtRev, nrow);
    // PDUMPR(W * (j + 1) - 1, pbwt);
    divc(nrow, w64, pbwt, pbwtPr, pbwtRev, pbwtPrRev,
         xor); // FIXME: xor here can be eliminated, thus also the allocation.

    // reversec(p0, prev, nrow);
    PDUMPR(W * (j + 1) - 1, pbwt);
    // SWAP(p0, p1);
    k++;
    j++;
  }

  uint8_t *c0 = NULL;

#if defined(BF2IOMODE_BM) || defined(BF2IOMODE_ENC)
  j *= W;
  fgetcolwgri(fin, j, nrow, w64, ncol, ncol - j);
#elif defined(BF2IOMODE_BCF)
  // no need to read here as it is already updated in failed condition
  // of the reading while
  ncol += _ncol;
  j *= W;
#else
#error UNDEFINED BEHAVIOUR
#endif
  // last column needs special handling, since it is < W
  // memcpy(p1->a, p0->a, nrow * sizeof *(p1->a));
  memcpy(pbwtPr->a, pbwt->a, nrow * sizeof *(pbwt->a));
  memcpy(pbwtPr->d, pbwt->d, nrow * sizeof *(pbwt->d));
  rrsortx(nrow, w64, pbwt->a, aux);
  memcpy(pbwtPrRev->a, pbwtRev->a, nrow * sizeof *(pbwtRev->a));
  memcpy(pbwtPrRev->d, pbwtRev->d, nrow * sizeof *(pbwtRev->d));
  reversec(pbwt, pbwtRev, nrow);
  // PDUMPR(W * (j + 1) - 1, pbwt);
  printf("last w has width:%zu\n", ncol - j);
  divc_last(nrow, w64, pbwt, pbwtPr, pbwtRev, pbwtPrRev, xor, ncol - j);
  PDUMPR(ncol - 1, pbwt);

  PBWTAD_FREE(pbwt);
  FREE(c0);
  FREE(pbwt);
  FREE(pbwtRev);
  FREE(pbwtPr);
  FREE(pbwtPrRev);
  FREE(w64);
  FREE(xor);
  return NULL;
}

pbwtad **wbapproxc_rrs(void *fin, size_t nrow, size_t ncol) {
  // Compute the bit-packed windows
  uint64_t *pw = malloc(nrow * sizeof *pw);
  size_t *aux = malloc(nrow * sizeof *aux);

  pbwtad *p0 = pbwtad_new(nrow);
  pbwtad *p1 = pbwtad_new(nrow);
  bfgetcolw64rn(fin, nrow, pw, ncol);
  rrsort0(nrow, pw, p1->a, aux);
  PDUMP(W - 1, p1);
  SWAP(p0, p1);

  size_t j;
  for (j = 1; j * W <= ncol - W; j++) {
    memcpy(p1->a, p0->a, nrow * sizeof *(p1->a));
    bfgetcolw64rn(fin, nrow, pw, ncol);
    rrsortx(nrow, pw, p1->a, aux);
    PDUMP(W * (j + 1) - 1, p1);
    SWAP(p0, p1);
  }

  uint8_t *c0 = NULL;
  j *= W;
  fgetcolwgri(fin, j, nrow, pw, ncol, ncol - j);
  memcpy(p1->a, p0->a, nrow * sizeof *(p1->a));
  rrsortx(nrow, pw, p1->a, aux);
  PDUMP(ncol - 1, p1);

  PBWTAD_FREE(p0);
  PBWTAD_FREE(p1);
  FREE(c0);
  FREE(pw);
  FREE(aux);
  return NULL;
}

pbwtad **swbapproxc_rrs(int fin, size_t nrow, size_t ncol) {
  // Compute the bit-packed windows
  uint64_t *w64 = malloc(nrow * sizeof *w64);
  uint64_t *xor = malloc(nrow * sizeof *xor);
  size_t *aux = malloc(nrow * sizeof *aux);

  pbwtad *pbwt = pbwtad_new(nrow);
  pbwtad *pbwtPr = pbwtad_new(nrow);
  pbwtad *pbwtRev = pbwtad_new(nrow);
  pbwtad *pbwtPrRev = pbwtad_new(nrow);

  sbfgetcolw64rn(fin, nrow, w64, ncol);
  rrsort0(nrow, w64, pbwt->a, aux);
  memcpy(pbwtPrRev->a, pbwtRev->a, nrow * sizeof *(pbwtRev->a));
  memcpy(pbwtPrRev->d, pbwtRev->d, nrow * sizeof *(pbwtRev->d));
  reversec(pbwt, pbwtRev, nrow);
  divc0(nrow, w64, pbwt);
  PDUMPR(W - 1, pbwt);

  size_t j, k = 1;
  for (j = 1; j * W <= ncol - W; j++) {
    sbfgetcolw64rn(fin, nrow, w64, ncol);
    memcpy(pbwtPr->a, pbwt->a, nrow * sizeof *(pbwt->a));
    memcpy(pbwtPr->d, pbwt->d, nrow * sizeof *(pbwt->d));
    rrsortx(nrow, w64, pbwt->a, aux); // BUG: pbwtPr->d doesnt contain anything.
    // FIXME: moved reversec after rrsortx, could now be integrated in rrsortx
    // as #1 pull request;
    memcpy(pbwtPrRev->a, pbwtRev->a, nrow * sizeof *(pbwtRev->a));
    memcpy(pbwtPrRev->d, pbwtRev->d, nrow * sizeof *(pbwtRev->d));
    reversec(pbwt, pbwtRev, nrow);
    // PDUMPR(W * (j + 1) - 1, pbwt);
    divc(nrow, w64, pbwt, pbwtPr, pbwtRev, pbwtPrRev,
         NULL); // FIXME: xor here can be eliminated, thus also the allocation.

    // reversec(p0, prev, nrow);
    PDUMPR(W * (j + 1) - 1, pbwt);
    // SWAP(p0, p1);
    k++;
  }

  uint8_t *c0 = NULL;
  j *= W;
  sfgetcolwgri(fin, j, nrow, w64, ncol, ncol - j);
  
  // last column needs special handling, since it is < W
  // memcpy(p1->a, p0->a, nrow * sizeof *(p1->a));
  memcpy(pbwtPr->a, pbwt->a, nrow * sizeof *(pbwt->a));
  memcpy(pbwtPr->d, pbwt->d, nrow * sizeof *(pbwt->d));
  rrsortx(nrow, w64, pbwt->a, aux);
  memcpy(pbwtPrRev->a, pbwtRev->a, nrow * sizeof *(pbwtRev->a));
  memcpy(pbwtPrRev->d, pbwtRev->d, nrow * sizeof *(pbwtRev->d));
  reversec(pbwt, pbwtRev, nrow);
  // PDUMPR(W * (j + 1) - 1, pbwt);
  divc_last(nrow, w64, pbwt, pbwtPr, pbwtRev, pbwtPrRev, xor, ncol - j);
  PDUMPR(ncol - 1, pbwt);

  PBWTAD_FREE(pbwt);
  FREE(c0);
  FREE(pbwt);
  FREE(pbwtRev);
  FREE(pbwtPr);
  FREE(pbwtPrRev);
  FREE(w64);
  FREE(xor);
  return NULL;
}

pbwtad **mwbapproxc_rrs(int fin, size_t nrow, size_t ncol) {
  // Compute the bit-packed windows
  uint64_t *pw = malloc(nrow * sizeof *pw);
  size_t *aux = malloc(nrow * sizeof *aux);

  pbwtad *p0 = pbwtad_new(nrow);
  pbwtad *p1 = pbwtad_new(nrow);
  sbfgetcolw64rn_mmap(fin, nrow, pw, ncol);
  rrsort0(nrow, pw, p1->a, aux);
  PDUMP(W - 1, p1);
  SWAP(p0, p1);

  size_t j;
  for (j = 1; j * W <= ncol - W; j++) {
    memcpy(p1->a, p0->a, nrow * sizeof *(p1->a));
    sbfgetcolw64rn_mmap(fin, nrow, pw, ncol);
    rrsortx(nrow, pw, p1->a, aux);
    PDUMP(W * (j + 1) - 1, p1);
  }

  uint8_t *c0 = NULL;
  j *= W;
  sfgetcolwgri(fin, j, nrow, pw, ncol, ncol - j);
  memcpy(p1->a, p0->a, nrow * sizeof *(p1->a));
  rrsortx(nrow, pw, p1->a, aux);
  PDUMP(ncol - 1, p1);

  PBWTAD_FREE(p0);
  PBWTAD_FREE(p1);
  FREE(c0);
  FREE(pw);
  FREE(aux);
  return NULL;
}

pbwtad **wparc_rrs(void *fin, size_t nrow, size_t ncol) {
  // NOTE: here it is necessary to keep in memory the entire windows
  pbwtad **pb0 = malloc(W * sizeof(pbwtad *));
  pbwtad **pb1 = malloc(W * sizeof(pbwtad *));
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
  pb0[0] = cpbwt(nrow, c0, p0);
  PDUMP(0, pb0[0]);
  PBWTAD_FREE(p0);

  for (int j = 1; j < W; j++) {
    fgetcoli(fin, j, nrow, c0, ncol);
    pb0[j] = cpbwt(nrow, c0, pb0[j - 1]);
    PDUMP(j, pb0[j]);
  }

#ifdef BF2IOMODE_BCF
  bcf_sr_destroy(_sr);
  fin = _tfin;
#endif
  uint64_t *pw0 = malloc(nrow * sizeof *pw0);
  uint64_t *pw1 = malloc(nrow * sizeof *pw1);
  size_t *aux = malloc(nrow * sizeof *aux);
  fgetcoliw64r(fin, 0, nrow, pw0, ncol);

  // pb0 is now filled with computed values,
  // to allow reusing I need to fill pb1 with empty values
  for (int j = 0; j < W; j++) {
    pb1[j] = pbwtad_new(nrow);
  }

  size_t j;
#if defined(BF2IOMODE_BM) || defined(BF2IOMODE_ENC)
  for (j = 1; j * W <= ncol - W;) {
    fgetcoliw64r(fin, j, nrow, pw1, ncol);
#elif defined(BF2IOMODE_BCF)
  j = 1;
  ncol = W;
  size_t _ncol = 0;
  while ((_ncol = fgetcoliw64r(fin, j, nrow, pw1, 0)) == W) {
    ncol += _ncol;
#else
#error UNDEFINED BEHAVIOUR
#endif
    pbwtad *ps = pb1[W - 1];
    memcpy(ps->a, pb0[W - 1]->a, nrow * sizeof *(ps->a));
    rrsortx(nrow, pw1, ps->a, aux);

#pragma omp parallel for shared(pw1, pw0, pb0, pb1, j)
    for (size_t x = 1; x < W; x++) {
      uint64_t *w = malloc(nrow * sizeof *w);
      // size_t J = W * (j + 1) - 1;
      size_t J = W - 1;

      wr64mrgsi(nrow, pw1, pw0, w, x);
      // pbwtad *ps = pbwtad_new(nrow);
      pbwtad *ps = pb1[J - x];
      // pb0[J - x] = ps;
      memcpy(ps->a, pb0[J - x]->a, nrow * sizeof *(ps->a));
      rrsortx_noaux(nrow, w, ps->a);
      FREE(w);
    }
    PDUMP_SEQ_OFFSET(0, W, pb1, W * j - 1);
    SWAP(pw0, pw1);
    SWAP(pb0, pb1);
    j++;
  }

  pbwtad *pp0, *pp1;
  pp0 = pb0[W - 1];
  pp1 = pb0[W - 2];

#ifdef BF2IOMODE_BCF
  ncol += _ncol;
  size_t wix = 0;
#endif

  for (j = j * W; j < ncol; j++) {
#if defined(BF2IOMODE_BM) || defined(BF2IOMODE_ENC)
    fgetcoli(fin, j, nrow, c0, ncol);
#elif defined(BF2IOMODE_BCF)
    // no need to read here as it is already updated in failed condition
    // of the reading while
    // WARN: if something goes wrong this should be the first place to
    // investigate, it seems correct but I'm not 100% sure
    for (size_t _i = 0; _i < nrow; _i++) {
      c0[_i] = (pw1[_i] >> wix) & 0x1;
    }
    wix++;
#else
#error UNDEFINED BEHAVIOUR
#endif
    cpbwti(nrow, c0, pp0, pp1);
    PDUMP(j, pp1);
    SWAP(pp0, pp1);
  }

  return NULL;
  for (int j = 0; j < W; j++) {
    PBWTAD_FREE(pb0[j]);
    PBWTAD_FREE(pb1[j]);
  }
  FREE(pb0);
  FREE(pb1);
  FREE(pw0);
  FREE(pw1);
  FREE(aux);
  FREE(c0);
  return NULL;
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
        {
          PDUMP(j + W, ps);
        }
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
  char _usage_args_[] = "[lin|bli[s|m]|ars|bar[s|m]|prs|bpr|spr] FILE\n";
  if (argc < 2) {
    fprintf(stderr, "Usage: %s %s FILE\n", argv[0], _usage_args_);
    return EXIT_FAILURE;
  }

#if defined(BF2IOMODE_BM) || defined(BF2IOMODE_ENC)
  FILE *fin = fopen(argv[2], "r");
  int fd = open(argv[2], O_RDONLY);
  if (!fin) {
    perror("[main]");
    return EXIT_FAILURE;
  }
#elif defined(BF2IOMODE_BCF)
  bcf_srs_t *sr = bcf_sr_init();
  bcf_sr_add_reader(sr, argv[2]);

  int fd = -1;
  void *fin = sr;
#else
#error BF2IOMODE is not specified
#endif

  size_t nrow, ncol;
  TRACE(fgetrc(fin, &nrow, &ncol));
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
    // r = linc(fin, nrow, ncol);
    TRACE(linc(fin, nrow, ncol), r);
  } else if (strcmp(argv[1], "bli") == 0) {
    // r = blinc(fin, nrow, ncol);
    TRACE(blinc(fin, nrow, ncol), r);
  } else if (strcmp(argv[1], "blis") == 0) {
    // r = sblinc(fd, nrow, ncol);
    TRACE(sblinc(fd, nrow, ncol), r);
  } else if (strcmp(argv[1], "blim") == 0) {
    // r = mblinc(fd, nrow, ncol);
    TRACE(mblinc(fd, nrow, ncol), r);
  } else if (strcmp(argv[1], "ars") == 0) {
    // r = wapproxc_rrs(fin, nrow, ncol);
    TRACE(wapproxc_rrs(fin, nrow, ncol), r);
  } else if (strcmp(argv[1], "bar") == 0) {
    // r = wbapproxc_rrs(fin, nrow, ncol);
    TRACE(wbapproxc_rrs(fin, nrow, ncol), r);
  } else if (strcmp(argv[1], "bars") == 0) {
    // r = swbapproxc_rrs(fd, nrow, ncol);
    TRACE(swbapproxc_rrs(fd, nrow, ncol), r);
  } else if (strcmp(argv[1], "barm") == 0) {
    // r = mwbapproxc_rrs(fd, nrow, ncol);
    TRACE(mwbapproxc_rrs(fd, nrow, ncol), r);
  } else if (strcmp(argv[1], "prs") == 0) {
    // r = wparc_rrs(fin, nrow, ncol);
    TRACE(wparc_rrs(fin, nrow, ncol), r);
  } else if (strcmp(argv[1], "bpr") == 0) {
    // r = bwparc_rrs(fin, nrow, ncol);
    TRACE(bwparc_rrs(fin, nrow, ncol), r);
  } else if (strcmp(argv[1], "spr") == 0) {
    // r = wstagparc_rrs(argv[2], nrow, ncol);
    TRACE(wstagparc_rrs(argv[2], nrow, ncol), r);
  } else if (strcmp(argv[1], "srs") == 0) {
    // r = wseq_rrs(fin, nrow, ncol);
    TRACE(wseq_rrs(fin, nrow, ncol), r);
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

#if defined(BF2IOMODE_BM) || defined(BF2IOMODE_ENC)
  // TODO: file cleanup
#elif defined(BF2IOMODE_BCF)
  bcf_sr_destroy(sr);
#else
#endif

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
