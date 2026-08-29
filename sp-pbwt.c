/* vim: set ft=c */
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

/*
 * pbidx_t is the element type of PBWT prefix (a) and divergence (d) arrays.
 * Row counts (~2*nsamples) and column counts (up to a few million, e.g.
 * ~6.2M for the largest 1000G phase3 chromosome) both fit comfortably in
 * 32 bits, so uint32_t roughly halves the memory traffic of every memcpy,
 * radix scatter, and reversec pass versus the previous size_t width.
 * Switch back to a 64-bit width by changing this typedef/macros if inputs
 * ever exceed UINT32_MAX rows or columns.
 */
typedef uint32_t pbidx_t;
#define PBIDX_FMT "%u"
#define PBIDX_MAX UINT32_MAX

/*
 * Defensive guard: column counters that get stored into a pbidx_t
 * divergence/position value (see cpbwt/cpbwti/cpbwtiLCP) must not exceed
 * PBIDX_MAX or they would silently wrap. Not expected to trigger for the
 * target datasets (columns top out around 6.2M), but BCF input's column
 * count is only known at EOF, so this is checked incrementally.
 */
static inline void pbidx_guard(size_t v) {
  if (v > PBIDX_MAX) {
    fprintf(stderr,
            "sp-pbwt: column index %zu exceeds pbidx_t capacity (%u); "
            "rebuild with a wider pbidx_t\n",
            v, (unsigned)PBIDX_MAX);
    exit(EXIT_FAILURE);
  }
}

#define DBDUMP
uint8_t DO_DUMP = 0;
#ifdef DBDUMP
#define CDUMP8(i, c)                                                           \
  do {                                                                         \
    if (DO_DUMP) {                                                             \
      printf("%zu:", (size_t)(i));                                             \
      size_t cdump_j__;                                                        \
      for (cdump_j__ = 0; cdump_j__ < nrow - 1; cdump_j__++)                   \
        printf("%u ", c[cdump_j__]);                                           \
      printf("%u", c[cdump_j__]);                                              \
      fputc(0xA, stdout);                                                      \
    }                                                                          \
  } while (0)

#define CDUMP(i, c)                                                            \
  do {                                                                         \
    if (DO_DUMP) {                                                             \
      printf("%zu:", (size_t)(i));                                             \
      size_t cdump_j__;                                                        \
      for (cdump_j__ = 0; cdump_j__ < nrow - 1; cdump_j__++)                   \
        printf("%llu ", c[cdump_j__]);                                         \
      printf("%llu", c[cdump_j__]);                                            \
      fputc(0xA, stdout);                                                      \
    }                                                                          \
  } while (0)
#define PDUMPR(i, p)                                                           \
  do {                                                                         \
    if (DO_DUMP) {                                                             \
      printf("%zu:", (size_t)(i));                                             \
      size_t pdump_j__;                                                        \
      for (pdump_j__ = 0; pdump_j__ < nrow - 1; pdump_j__++)                   \
        printf(PBIDX_FMT " ", (p)->a[pdump_j__]);                              \
      printf(PBIDX_FMT, (p)->a[pdump_j__]);                                    \
      fputc('|', stdout);                                                      \
      for (pdump_j__ = 0; pdump_j__ < nrow - 1; pdump_j__++)                   \
        printf("%zu ", (size_t)(i) + 1 - (p)->d[pdump_j__]);                   \
      printf("%zu", (size_t)(i) + 1 - (p)->d[pdump_j__]);                      \
      fputc(0xA, stdout);                                                      \
    }                                                                          \
  } while (0)
#define PDUMP(i, p)                                                            \
  do {                                                                         \
    if (DO_DUMP) {                                                             \
      printf("%zu:", (size_t)(i));                                             \
      size_t pdump_j__;                                                        \
      for (pdump_j__ = 0; pdump_j__ < nrow - 1; pdump_j__++)                   \
        printf(PBIDX_FMT " ", (p)->a[pdump_j__]);                              \
      printf(PBIDX_FMT, (p)->a[pdump_j__]);                                    \
      fputc('|', stdout);                                                      \
      for (pdump_j__ = 0; pdump_j__ < nrow - 1; pdump_j__++)                   \
        printf(PBIDX_FMT " ", (p)->d[pdump_j__]);                              \
      printf(PBIDX_FMT, (p)->d[pdump_j__]);                                    \
      fputc(0xA, stdout);                                                      \
    }                                                                          \
  } while (0)

#define PDUMP_SEQR(s, e, p)                                                    \
  do {                                                                         \
    for (size_t pdump_ix__ = (s); pdump_ix__ < (e); pdump_ix__++) {            \
      if (DO_DUMP) {                                                           \
        printf("%zu:", (size_t)(pdump_ix__));                                  \
        size_t pdump_j__;                                                      \
        for (pdump_j__ = 0; pdump_j__ < nrow - 1; pdump_j__++)                 \
          printf(PBIDX_FMT " ", (p)[pdump_ix__]->a[pdump_j__]);                \
        printf(PBIDX_FMT, (p)[pdump_ix__]->a[pdump_j__]);                      \
        fputc('|', stdout);                                                    \
        for (pdump_j__ = 0; pdump_j__ < nrow - 1; pdump_j__++)                 \
          printf("%zu ", 1 + (pdump_ix__) - (p)[pdump_ix__]->d[pdump_j__]);    \
        printf("%zu", 1 + (pdump_ix__) - (p)[pdump_ix__]->d[pdump_j__]);       \
        fputc(0xA, stdout);                                                    \
      }                                                                        \
    }                                                                          \
  } while (0)
#define PDUMP_SEQ(s, e, p)                                                     \
  do {                                                                         \
    for (size_t pdump_ix__ = (s); pdump_ix__ < (e); pdump_ix__++) {            \
      if (DO_DUMP) {                                                           \
        printf("%zu:", (size_t)(pdump_ix__));                                  \
        size_t pdump_j__;                                                      \
        for (pdump_j__ = 0; pdump_j__ < nrow - 1; pdump_j__++)                 \
          printf(PBIDX_FMT " ", (p)[pdump_ix__]->a[pdump_j__]);                \
        printf(PBIDX_FMT, (p)[pdump_ix__]->a[pdump_j__]);                      \
        fputc('|', stdout);                                                    \
        for (pdump_j__ = 0; pdump_j__ < nrow - 1; pdump_j__++)                 \
          printf(PBIDX_FMT " ", (p)[pdump_ix__]->d[pdump_j__]);                \
        printf(PBIDX_FMT, (p)[pdump_ix__]->d[pdump_j__]);                      \
        fputc(0xA, stdout);                                                    \
      }                                                                        \
    }                                                                          \
  } while (0)
#define PDUMP_SEQ_OFFSETR(s, e, p, offset)                                     \
  do {                                                                         \
    for (size_t pdump_ix__ = (s); pdump_ix__ < (e); pdump_ix__++) {            \
      if (DO_DUMP) {                                                           \
        printf("%zu:", (size_t)(offset) + (size_t)(pdump_ix__));               \
        size_t pdump_j__;                                                      \
        for (pdump_j__ = 0; pdump_j__ < nrow - 1; pdump_j__++)                 \
          printf(PBIDX_FMT " ", (p)[pdump_ix__]->a[pdump_j__]);                \
        printf(PBIDX_FMT, (p)[pdump_ix__]->a[pdump_j__]);                      \
        fputc('|', stdout);                                                    \
        for (pdump_j__ = 0; pdump_j__ < nrow - 1; pdump_j__++)                 \
          printf("%zu ",                                                       \
                 offset + pdump_ix__ + 1 - (p)[pdump_ix__]->d[pdump_j__]);     \
        printf("%zu",                                                          \
               offset + pdump_ix__ + 1 - (p)[pdump_ix__]->d[pdump_j__]);       \
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
          printf(PBIDX_FMT " ", (p)[pdump_ix__]->a[pdump_j__]);                \
        printf(PBIDX_FMT, (p)[pdump_ix__]->a[pdump_j__]);                      \
        fputc('|', stdout);                                                    \
        for (pdump_j__ = 0; pdump_j__ < nrow - 1; pdump_j__++)                 \
          printf(PBIDX_FMT " ", (p)[pdump_ix__]->d[pdump_j__]);                \
        printf(PBIDX_FMT, (p)[pdump_ix__]->d[pdump_j__]);                      \
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
  pbidx_t *a;
  pbidx_t *d;
};

void rrsortx(size_t n, uint64_t *c, pbidx_t *s, pbidx_t *aux) {
  pbidx_t *tmp;
  size_t j;
  pbidx_t *pre = s;
  pbidx_t *post = aux;
  uint8_t b;

  size_t cnt[8][256] = {0};
  for (size_t j = 0; j < n; j++) {
    uint64_t val = c[j];
    cnt[0][(val) & 0xFFULL]++;
    cnt[1][(val >> 8) & 0xFFULL]++;
    cnt[2][(val >> 16) & 0xFFULL]++;
    cnt[3][(val >> 24) & 0xFFULL]++;
    cnt[4][(val >> 32) & 0xFFULL]++;
    cnt[5][(val >> 40) & 0xFFULL]++;
    cnt[6][(val >> 48) & 0xFFULL]++;
    cnt[7][(val >> 56) & 0xFFULL]++;
  }

  for (size_t i = 0; i < 8; i++) {
    // prefix sum
    for (size_t j = 1; j < 256; j++)
      cnt[i][j] += cnt[i][j - 1];
    // sorting
    for (ssize_t j = n - 1; j >= 0; --j) {
      b = (c[pre[j]] >> (8 * i)) & 0xFFULL;
      cnt[i][b]--;
      post[cnt[i][b]] = pre[j];
    }
    // swap s and aux
    tmp = pre;
    pre = post;
    post = tmp;
  }
}

/*
 * Same as rrsortx, but reads the starting permutation from `src` instead of
 * sorting in place: `src` is left untouched, and the final sorted result
 * lands in `dst` after all 8 passes (never in `src` or `aux`). Lets a caller
 * double-buffer by pointer-swapping whole pbwtad structs across windows
 * instead of memcpy-ing a snapshot of the previous window's `a` array before
 * every sort.
 */
void rrsortx_src(size_t n, uint64_t *c, pbidx_t *src, pbidx_t *dst,
                 pbidx_t *aux) {
  pbidx_t *tmp;
  pbidx_t *pre = src;
  pbidx_t *post = aux;
  uint8_t b;

  size_t cnt[8][256] = {0};
  for (size_t j = 0; j < n; j++) {
    uint64_t val = c[j];
    cnt[0][(val) & 0xFFULL]++;
    cnt[1][(val >> 8) & 0xFFULL]++;
    cnt[2][(val >> 16) & 0xFFULL]++;
    cnt[3][(val >> 24) & 0xFFULL]++;
    cnt[4][(val >> 32) & 0xFFULL]++;
    cnt[5][(val >> 40) & 0xFFULL]++;
    cnt[6][(val >> 48) & 0xFFULL]++;
    cnt[7][(val >> 56) & 0xFFULL]++;
  }

  for (size_t i = 0; i < 8; i++) {
    // prefix sum
    for (size_t j = 1; j < 256; j++)
      cnt[i][j] += cnt[i][j - 1];
    // sorting
    for (ssize_t j = n - 1; j >= 0; --j) {
      b = (c[pre[j]] >> (8 * i)) & 0xFFULL;
      cnt[i][b]--;
      post[cnt[i][b]] = pre[j];
    }
    if (i == 0) {
      // pass 0 read from `src` (never written); from here on ping-pong
      // between `aux` and `dst` only, so `src` stays intact.
      pre = aux;
      post = dst;
    } else {
      // swap pre and post
      tmp = pre;
      pre = post;
      post = tmp;
    }
  }
}

/*
 * Sort `c[n]`, sorted permutations will be in saved in `s[n]`,
 * without using externally allocated `aux[n]` auxiliary array.
 * This version assumes the `s` array to be already initialized.
 */
void rrsortx_noaux(size_t n, uint64_t *c, pbidx_t *s) {
  pbidx_t *tmp;
  size_t j;
  pbidx_t *pre = s;
  pbidx_t *post = malloc(n * sizeof *post);
  uint8_t b;

  for (size_t i = 0; i < 8; i++) {
    size_t cnt[256] = {0};

    // frequencies
    for (j = 0; j < n; j++) {
      b = (c[j] >> (8 * i)) & 0xFFULL;
      cnt[b]++;
    }
    // prefix sum
    for (size_t j = 1; j < 256; j++)
      cnt[j] += cnt[j - 1];
    // sorting
    for (ssize_t j = n - 1; j >= 0; --j) {
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
void rrsort0(size_t n, uint64_t *c, pbidx_t *s, pbidx_t *aux) {
  pbidx_t *tmp;
  size_t j;
  pbidx_t *pre = s;
  pbidx_t *post = aux;
  uint8_t b;

  // this is needed if:
  // 1. we want to sort numbers
  // 2. this is the first iteration
  //
  // In normal BWT cases we assume to have
  // `s` equal to the sorting of the previous "column"
  for (size_t i = 0; i < n; ++i)
    pre[i] = (pbidx_t)i;

  for (size_t i = 0; i < 8; i++) {
    size_t cnt[256] = {0};

    // frequencies
    for (j = 0; j < n; j++) {
      b = (c[j] >> (8 * i)) & 0xFFULL;
      cnt[b]++;
    }
    // prefix sum
    for (size_t j = 1; j < 256; j++)
      cnt[j] += cnt[j - 1];
    // sorting
    for (ssize_t j = n - 1; j >= 0; --j) {
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

/* Compute *p's reverse auxiliary pbwt arrays in *rev
 * rev->a[i] contains the position of row #i in in p->a[]
 * rev->a[p->a[i]] = i = p->a[rev->a[i]]
 *
 * rev->d[i] instead contains the divergence of row i in p->a[]
 * rev->d[p->a[i]] = p->d[i] = p->d[rev->a[p->a[i]]]
 *
 * run this after rrsorting and before computing div on windows
 */
void reversec(pbwtad *p, pbwtad *rev, size_t n) {
  for (size_t i = 0; i < n; i++) {
    rev->a[p->a[i]] = (pbidx_t)i;
    rev->d[p->a[i]] = p->d[i]; // == p->d[rev->a[i]]
    assert(rev->d[p->a[i]] == p->d[rev->a[p->a[i]]]);
  }
}
void reversecprev(pbwtad *p, pbwtad *pp, pbwtad *rev, size_t n) {
  for (size_t i = 0; i < n; i++) {
    rev->a[p->a[i]] = (pbidx_t)i;
    rev->d[p->a[i]] = pp->d[i]; // == p->d[rev->a[i]]
  }
}

/*
 * Recover divergence of a match possibly longer than W.
 * Iterates over the range between previous and current row in p->a
 * using the previous pbwt array. Compute correct divergence using
 * reverse arrays computed by reversec
 * WARN: prev (current pbwt reverse) is not used, consider not memcpy in
 * wapprox computation if not needed.
 */
pbidx_t recover_div(size_t n, size_t w, size_t i, size_t i0, uint64_t *c,
                   pbwtad *p, pbwtad *ppr, pbwtad *prev, pbwtad *pprrev) {

  pbidx_t d;
  size_t j = pprrev->a[i];
  pbidx_t min = ppr->d[j];

  for (size_t j = (pprrev->a[i0]) + 1; j < (pprrev->a[i]); j++) {
    if (ppr->d[j] < min) {
      min = ppr->d[j];
    }
  }
  d = (pbidx_t)w + min;
  return d;
}

/*
 * compute the divergence of the first w64 windows
 */
void divc0(size_t n, uint64_t *c, pbwtad *p) {
  uint64_t x = 0;
  p->d[0] = 0;
  for (size_t i = 1; i < n; i++) {
    x = c[p->a[i]] ^ c[p->a[i - 1]];
    p->d[i] = x ? __builtin_clzll(x) : 64;
  }
}

/*
 * Computes the divergence of a generic w64 window;
 * LCP values equal to the window size get recoverd by recover_div function
 */
void divc(size_t n, uint64_t *c, pbwtad *p, pbwtad *ppr, pbwtad *prev,
          pbwtad *pprrev, size_t wi) {
  // c contains 64bit-encoded ints
  // xor of each c[s[i]] and its preceeding;
  // x[0] contains no information, previous x information is discarded;
  // here 64 is the size of the window
  /* (a write-only `static int8_t kk` debug counter used to live here; it was
   * never read, but it was written from every thread of the blockpar/stagpar
   * teams — ThreadSanitizer flagged it as a real data race. Removed.) */
  uint64_t x = 0;
  size_t w = wi ? wi : W;
  pbidx_t div;
  p->d[0] = 0;

  for (size_t i = 1; i < n; i++) {
    x = c[p->a[i]] ^ c[p->a[i - 1]];
    div = x ? __builtin_clzll(x) : w;
    p->d[i] = (div >= w) ? recover_div(n, w, p->a[i], p->a[i - 1], c, p, ppr,
                                       prev, pprrev)
                         : div;
  }
}

static inline pbwtad *pbwtad_new(size_t n) {
  pbwtad *p = malloc(sizeof *p);
  // a and d are allocated as a single contiguous 2*n block (a is the base,
  // d is the second half); this is two mallocs instead of three, and keeps
  // a/d adjacent in memory. PBWTAD_FREE relies on `a` being the block's
  // base pointer for every pbwtad, however constructed (see cpbwt()).
  pbidx_t *block = malloc(2 * n * sizeof *block);
  p->a = block;
  p->d = block + n;
  return p;
}

#define PBWTAD_FREE(p)                                                         \
  do {                                                                         \
    /* a/d share one allocation with `a` as its base pointer (see            \
     * pbwtad_new/cpbwt); free only `a`, and just clear `d` to avoid a       \
     * dangling pointer being used after the block is gone. */               \
    FREE((p)->a);                                                              \
    (p)->d = NULL;                                                             \
    FREE(p);                                                                   \
  } while (0)

/*
 * Given the current column index, swaps divergence values
 * between LCP and starting position of a match.
 *
 * Window computation is currently written using LCP values,
 * i.e. length of the actual longest co-lexicographical match
 * between p->a[i] and p->a[i-1],
 * while linear computation relis on "classical" divergence
 * definition of "starting position of the longest match
 * between p->a[i] and p->a[i-1]
 *
 * For tesing only during development phase
 */
void swapdiv(pbwtad *p, size_t n, size_t k) {
  for (size_t t = 0; t < n; t++) {
    p->d[t] = 1 + k - p->d[t];
  }
}

/* Number of zeros in c[0..n-1].  Since the prefix array is a permutation of
 * 0..n-1, this equals the number of loop iterations below whose mask is 0,
 * i.e. the final value of `r` in the old two-staging-array formulation.
 * Sequential byte pass, auto-vectorizes. */
static inline size_t count_zeros(size_t n, const uint8_t *restrict c) {
  size_t n0 = 0;
  for (size_t t = 0; t < n; t++)
    n0 += (c[t] == 0);
  return n0;
}

static pbwtad *cpbwt(size_t n, uint8_t *restrict c, pbwtad *restrict p) {
  static pbidx_t k = 1;

  pbwtad *ret = malloc(sizeof *ret);
  // single 2*n block, `a` as base — see pbwtad_new/PBWTAD_FREE.
  pbidx_t *block = malloc(2 * n * sizeof *block);
  ret->a = block;
  ret->d = block + n;

  // zeros land at 0.., ones directly at n0.. — no staging, no trailing memcpy
  size_t r = 0, q = count_zeros(n, c);
  pbidx_t f = k, g = k;

  size_t i;
  for (i = 0; i < n; i++) {
    pbidx_t idx = p->a[i];
    pbidx_t ddx = p->d[i];

    f = (ddx > f) ? ddx : f;
    g = (ddx > g) ? ddx : g;

    size_t mask = c[idx];
    size_t pos = mask ? q : r; // cmov
    pbidx_t dv = mask ? g : f; // cmov
    ret->a[pos] = idx;
    ret->d[pos] = dv;

    f &= -mask;       // f = 0 if mask == 0, unchanged if mask == 1
    g &= -(1 - mask); // g = 0 if mask == 1, unchanged if mask == 0
    q += mask;        // Increment q if mask is 1
    r += mask ^ 1;    // Increment r if mask is 0
  }

  k++;
  pbidx_guard(k);
  return ret;
}

/* BUG: (possibily?, needs testing)
 * use cpbwt std version (withtout *LCP) and then swap divergence
 * with swapdiv if necessary.
 */
static int cpbwtiLCP(size_t n, size_t k, uint8_t *restrict c,
                     pbwtad *restrict pp, pbwtad *restrict pc) {
  pbidx_guard(k);
  swapdiv(pp, n, k - 1);

  // zeros land at 0.., ones directly at n0.. — no staging, no trailing memcpy
  size_t r = 0, q = count_zeros(n, c);
  pbidx_t f = (pbidx_t)k + 1, g = (pbidx_t)k + 1;

  size_t i;
  for (i = 0; i < n; i++) {
    pbidx_t idx = pp->a[i];
    pbidx_t ddx = pp->d[i];

    f = (ddx > f) ? ddx : f;
    g = (ddx > g) ? ddx : g;

    size_t mask = c[idx];
    size_t pos = mask ? q : r; // cmov
    pbidx_t dv = mask ? g : f; // cmov
    pc->a[pos] = idx;
    pc->d[pos] = dv;

    f &= -mask;       // f = 0 if mask == 0, unchanged if mask == 1
    g &= -(1 - mask); // g = 0 if mask == 1, unchanged if mask == 0
    q += mask;        // Increment q if mask is 1
    r += mask ^ 1;    // Increment r if mask is 0
  }

  swapdiv(pc, n, k);
  return 1;
}

static int cpbwti(size_t n, uint8_t *restrict c, pbwtad *restrict pp,
                  pbwtad *restrict pc) {
  static pbidx_t k = 1;

  // zeros land at 0.., ones directly at n0.. — no staging, no trailing memcpy
  size_t r = 0, q = count_zeros(n, c);
  pbidx_t f = k, g = k;

  size_t i;
  for (i = 0; i < n; i++) {
    pbidx_t idx = pp->a[i];
    pbidx_t ddx = pp->d[i];

    f = (ddx > f) ? ddx : f;
    g = (ddx > g) ? ddx : g;

    size_t mask = c[idx];
    size_t pos = mask ? q : r; // cmov
    pbidx_t dv = mask ? g : f; // cmov
    pc->a[pos] = idx;
    pc->d[pos] = dv;

    f &= -mask;       // f = 0 if mask == 0, unchanged if mask == 1
    g &= -(1 - mask); // g = 0 if mask == 1, unchanged if mask == 0
    q += mask;        // Increment q if mask is 1
    r += mask ^ 1;    // Increment r if mask is 0
  }

  k++;
  pbidx_guard(k);
  return 1;
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

pbwtad **wapproxc_rrs(void *fin, size_t nrow, size_t ncol) { // ARS
  // Compute the bit-packed windows
  uint64_t *w64 =
      malloc(nrow * sizeof *w64); // window data collected by fgetcoliw64r
  pbidx_t *aux = malloc(nrow * sizeof *aux);

  pbwtad *pbwt = pbwtad_new(nrow);      // curr pbwt
  pbwtad *pbwtPr = pbwtad_new(nrow);    // prev pbwt
  pbwtad *pbwtRev = pbwtad_new(nrow);   // curr pbwt REVERSE
  pbwtad *pbwtPrRev = pbwtad_new(nrow); // prev pbwt REVERSE

#if defined(BF2IOMODE_BM) || defined(BF2IOMODE_BCF)
  fgetcoliw64r(fin, 0, nrow, w64, ncol);
  // CDUMP(0, w64);
#elif defined(BF2IOMODE_ENC)
  fgetcoliwg(fin, 0, nrow, w64, ncol, W);
  // parr(nrow, w64, "%llu,");
#else
#error UNDEFINED BEHAVIOUR
#endif

  rrsort0(nrow, w64, pbwt->a, aux);
  // WARN: following 2 memcpy(s) dump current pbwtRev (empty??) into pbwtPrRev
  // probably useless, both arrays should be already initialized.
  // memcpy(pbwtPrRev->a, pbwtRev->a, nrow * sizeof *(pbwtRev->a));
  // memcpy(pbwtPrRev->d, pbwtRev->d, nrow * sizeof *(pbwtRev->d));
  reversec(pbwt, pbwtRev, nrow);
  divc0(nrow, w64, pbwt);

  PDUMPR(W - 1, pbwt);
  size_t j;
  size_t k = 1;
#if defined(BF2IOMODE_BM)
  for (j = 1; j * W <= ncol - W;) {
    // memcpy(p1->a, p0->a, nrow * sizeof *(p1->a));
    fgetcoliw64r(fin, j, nrow, w64, ncol);

#elif defined(BF2IOMODE_ENC)
  for (j = 1; j * W <= ncol - W;) {
    fgetcoliwg(fin, j, nrow, w64, ncol, W);

#elif defined(BF2IOMODE_BCF)
  j = 1;
  ncol = W;
  size_t _ncol = 0;
  while ((_ncol = fgetcoliw64r(fin, j, nrow, w64, 0)) == W) {
    ncol += _ncol;
#else
#error UNDEFINED BEHAVIOUR
#endif
    // Double-buffer instead of memcpy-snapshotting: pbwtPr->a and
    // pbwtPrRev->d are never read (divc only reads ppr->d and pprrev->a, see
    // divc/recover_div), so those two directions need no data at all, live
    // or otherwise. For the two that are read, swapping the whole pbwtad*
    // (never individual a/d — they share one allocation, see pbwtad_new)
    // makes last iteration's already-computed values show up as "prev" for
    // free, instead of copying them.
    SWAP(pbwt, pbwtPr);
    SWAP(pbwtRev, pbwtPrRev);
    rrsortx_src(nrow, w64, pbwtPr->a, pbwt->a,
                aux); // sort pbwtPr's (prev window's) permutation into pbwt->a
    reversec(pbwt, pbwtRev,
             nrow); // computing reversec after sorting the new array
    divc(nrow, w64, pbwt, pbwtPr, pbwtRev, pbwtPrRev, W);
    PDUMPR(W * (j + 1) - 1, pbwt);
    k++;
    j++;
  }

  uint8_t *c0 = NULL;

  // LAST WINDOW
#if defined(BF2IOMODE_BM)
  j *= W;
  fgetcolwgri(fin, j, nrow, w64, ncol, ncol - j);
#elif defined(BF2IOMODE_ENC)
  j *= W;
  fgetcoliwg(fin, j, nrow, w64, ncol, W);
#elif defined(BF2IOMODE_BCF)
  // no need to read here as it is already updated in failed condition
  // of the reading while
  ncol += _ncol;
  j *= W;
#else
#error UNDEFINED BEHAVIOUR
#endif
  // last column needs special handling, since it is < W
  SWAP(pbwt, pbwtPr);
  SWAP(pbwtRev, pbwtPrRev);
  rrsortx_src(nrow, w64, pbwtPr->a, pbwt->a, aux);
  reversec(pbwt, pbwtRev, nrow);
  divc(nrow, w64, pbwt, pbwtPr, pbwtRev, pbwtPrRev, ncol - j);
  PDUMPR(ncol - 1, pbwt);

  PBWTAD_FREE(pbwt);
  PBWTAD_FREE(pbwtRev);
  PBWTAD_FREE(pbwtPr);
  PBWTAD_FREE(pbwtPrRev);
  FREE(w64);
  FREE(aux);
  return NULL;
}

pbwtad **swbapproxc_rrs(int fin, size_t nrow, size_t ncol) { // BARS
  // Compute the bit-packed windows
  uint64_t *w64 = malloc(nrow * sizeof *w64);
  pbidx_t *aux = malloc(nrow * sizeof *aux);

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
  PDUMPR(W - 1, pbwtRev);

  PDUMPR(W - 1, pbwtPr);
  PDUMPR(W - 1, pbwtPrRev);

  size_t j, k = 1;
  for (j = 1; j * W <= ncol - W; j++) {
    sbfgetcolw64rn(fin, nrow, w64, ncol);
    memcpy(pbwtPr->a, pbwt->a, nrow * sizeof *(pbwt->a));
    memcpy(pbwtPr->d, pbwt->d, nrow * sizeof *(pbwt->d));
    rrsortx(nrow, w64, pbwt->a, aux);
    memcpy(pbwtPrRev->a, pbwtRev->a, nrow * sizeof *(pbwtRev->a));
    memcpy(pbwtPrRev->d, pbwtRev->d, nrow * sizeof *(pbwtRev->d));
    reversec(pbwt, pbwtRev, nrow);
    divc(nrow, w64, pbwt, pbwtPr, pbwtRev, pbwtPrRev, W);
    PDUMPR(W * (j + 1) - 1, pbwt);
    k++;
  }

  uint8_t *c0 = NULL;
  j *= W;
  sfgetcolwgri(fin, j, nrow, w64, ncol, ncol - j);

  // last column needs special handling, since it is < W
  memcpy(pbwtPr->a, pbwt->a, nrow * sizeof *(pbwt->a));
  memcpy(pbwtPr->d, pbwt->d, nrow * sizeof *(pbwt->d));
  rrsortx(nrow, w64, pbwt->a, aux);
  memcpy(pbwtPrRev->a, pbwtRev->a, nrow * sizeof *(pbwtRev->a));
  memcpy(pbwtPrRev->d, pbwtRev->d, nrow * sizeof *(pbwtRev->d));
  reversec(pbwt, pbwtRev, nrow);
  divc(nrow, w64, pbwt, pbwtPr, pbwtRev, pbwtPrRev, ncol - j);
  PDUMPR(ncol - 1, pbwt);

  PBWTAD_FREE(pbwt);
  FREE(c0);
  FREE(pbwt);
  FREE(pbwtRev);
  FREE(pbwtPr);
  FREE(pbwtPrRev);
  FREE(w64);
  return NULL;
}

pbwtad **mwbapproxc_rrs(int fin, size_t nrow, size_t ncol) { // BARM
  // Compute the bit-packed windows
  uint64_t *pw = malloc(nrow * sizeof *pw);
  pbidx_t *aux = malloc(nrow * sizeof *aux);

  pbwtad *p0 = pbwtad_new(nrow);
  pbwtad *p1 = pbwtad_new(nrow);
  pbwtad *p0rev = pbwtad_new(nrow);
  pbwtad *p1rev = pbwtad_new(nrow);
  sbfgetcolw64rn_mmap(fin, nrow, pw, ncol);
  rrsort0(nrow, pw, p0->a, aux);

  reversec(p0, p0rev, nrow);
  divc0(nrow, pw, p0);
  PDUMPR(W - 1, p0);

  size_t j;
  for (j = 1; j * W <= ncol - W; j++) {
    sbfgetcolw64rn_mmap(fin, nrow, pw, ncol); // read next window
    memcpy(p1->a, p0->a,
           nrow * sizeof *(p1->a)); // copy previous pbwt in P1 (a)
    memcpy(p1->d, p0->d,
           nrow * sizeof *(p1->d)); // copy previous pbwt in P1 (d)

    rrsortx(nrow, pw, p0->a, aux); // radix sorting p0->a with auxiliary array
    memcpy(p1rev->a, p0rev->a,
           nrow * sizeof *(p0rev->a)); // copy previous pbwtReverse in P1rev (a)
    memcpy(p1rev->d, p0rev->d,
           nrow * sizeof *(p0rev->d)); // copy previous pbwtReverse in P1rev (d)
    reversec(p0, p0rev, nrow); // computing reversec after sorting the new array
                               // in P0 where P0[i] = P0rev[ P0->a[i] ]
    divc(nrow, pw, p0, p1, p0rev, p1rev,
         W); // compute divergence from p1 to p0, using also reverse arrays
    PDUMPR(W * (j + 1) - 1, p1);
  }
  // same stuff for last window of size < W
  uint8_t *c0 = NULL;
  j *= W;
  sfgetcolwgri(fin, j, nrow, pw, ncol, ncol - j);
  memcpy(p1->a, p0->a, nrow * sizeof *(p1->a));
  memcpy(p1->d, p0->d, nrow * sizeof *(p1->d));
  rrsortx(nrow, pw, p0->a, aux);
  memcpy(p1rev->a, p0rev->a, nrow * sizeof *(p1rev->a));
  memcpy(p1rev->d, p0rev->d, nrow * sizeof *(p1rev->d));
  reversec(p0, p0rev, nrow);
  divc(nrow, pw, p0, p1, p0rev, p1rev, ncol - j);

  PDUMPR(ncol - 1, p0);

  PBWTAD_FREE(p0);
  PBWTAD_FREE(p1);
  PBWTAD_FREE(p0rev);
  PBWTAD_FREE(p1rev);
  FREE(c0);
  FREE(pw);
  FREE(aux);
  return NULL;
}

/* parallel mixed-windows pbwt
 *
 */
pbwtad **wparc_rrs(void *fin, size_t nrow, size_t ncol) { // PRS
  // NOTE: here it is necessary to keep in memory the entire windows
  pbwtad **pb0 = malloc(W * sizeof(pbwtad *));
  pbwtad **pb0rev = malloc(W * sizeof(pbwtad *));
  pbwtad **pb1 = malloc(W * sizeof(pbwtad *));
  pbwtad **pb1rev = malloc(W * sizeof(pbwtad *));
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
  // printf("first col guard\n");
  /* The invariant maintained for the whole window array is
   *   pb0rev[j]->a[row] == position of `row` in pb0[j]->a
   * (see the reversecprev() call below for j>=1, and reversec() in the
   * parallel body). pb0rev[0] used to be aliased to p0, i.e. the identity
   * permutation, which is the reverse of the state *before* column 0, not
   * of pb0[0]. recover_div() indexes ppr->d through pprrev->a, and the
   * (ppr, pprrev) pair for the first column of every block is
   * (pb0[0], pb0rev[0]) — so that off-by-one column produced wrong
   * divergences at exactly the block-boundary columns j*W. */
  pb0rev[0] = pbwtad_new(nrow);
  reversecprev(pb0[0], p0, pb0rev[0], nrow);
  PBWTAD_FREE(p0);
  PDUMP(0, pb0[0]);
  for (int j = 1; j < W; j++) {
    fgetcoli(fin, j, nrow, c0, ncol);
    pb0[j] = cpbwt(nrow, c0, pb0[j - 1]);
    pb0rev[j] = pbwtad_new(nrow);
    reversecprev(pb0[j], pb0[j - 1], pb0rev[j], nrow);
    PDUMP(j, pb0[j]);
  }

#ifdef BF2IOMODE_BCF
  bcf_sr_destroy(_sr);
  fin = _tfin;
#endif
  /* swapping div values with LCP for window computation */
  for (int j = 0; j < W; j++) {
    swapdiv(pb0[j], nrow, j);
    swapdiv(pb0rev[j], nrow, j);
    // PDUMPR(j, pb0[j]);
  }
  /* Three window buffers: pwPrev (window j-1), pwCur (window j, being
   * computed on) and pwNext (window j+1, being read *concurrently* by the
   * master thread while the rest of the team works on pwPrev/pwCur). They
   * are rotated once per round. `pw0`/`pw1` below alias pwPrev/pwCur after
   * the loop, for the sequential tail. */
  uint64_t *pwPrev = malloc(nrow * sizeof *pwPrev);
  uint64_t *pwCur = malloc(nrow * sizeof *pwCur);
  uint64_t *pwNext = malloc(nrow * sizeof *pwNext);
  uint64_t *pw0, *pw1;
  fgetcoliw64r(fin, 0, nrow, pwPrev, ncol);

  // pb0 is now filled with computed values,
  // to allow reusing I need to fill pb1 with empty values
  for (int j = 0; j < W; j++) {
    pb1[j] = pbwtad_new(nrow);
    pb1rev[j] = pbwtad_new(nrow);
  }

  size_t j = 1;
  size_t _ncol = 0;
  /* WPARC_FETCH(idx, buf): read window `idx` into `buf`; evaluates to 1 if a
   * full window was obtained (so another round can run), 0 otherwise. This
   * reproduces exactly the loop conditions of the two original backend
   * variants (BM/ENC: bounded by the known `ncol`; BCF: EOF-driven). */
#if defined(BF2IOMODE_BM) || defined(BF2IOMODE_ENC)
#define WPARC_FETCH(idx, buf)                                                  \
  (((idx) * W <= ncol - W)                                                     \
       ? (fgetcoliw64r(fin, (idx), nrow, (buf), ncol), 1)                      \
       : 0)
#elif defined(BF2IOMODE_BCF)
  ncol = W;
#define WPARC_FETCH(idx, buf)                                                  \
  ((_ncol = fgetcoliw64r(fin, (idx), nrow, (buf), 0)) == W ? (ncol += _ncol, 1)\
                                                          : 0)
#else
#error UNDEFINED BEHAVIOUR
#endif

  int active = WPARC_FETCH(j, pwCur);
  int nxt_ok = 0;

  /* One persistent team for the whole column loop: previously a fresh
   * `omp parallel for` was forked per 64-column window (~100k team
   * fork/joins on a whole chromosome) and every iteration malloc'd its own
   * merge buffer while rrsortx_noaux malloc'd a radix scratch internally.
   * Now the team lives across all windows, each thread owns its merge
   * buffer and its radix `aux` for its whole lifetime, and the sort is
   * rrsortx (one fused histogram pass) rather than rrsortx_noaux (eight). */
#pragma omp parallel
  {
    uint64_t *w = malloc(nrow * sizeof *w);
    pbidx_t *taux = malloc(nrow * sizeof *taux);

    while (active) {
      /* I/O overlap: the master thread issues the next window's read and
       * then joins the (dynamically scheduled) work-sharing loop below, so
       * the serial BCF decode is hidden behind the team's compute on the
       * current window. The reader writes only pwNext; the team reads only
       * pwPrev/pwCur. The `omp for`'s implicit barrier orders the read
       * against the buffer rotation in the `single` that follows. Pinning
       * the read to the master thread (rather than `single nowait`) keeps
       * iobcf.c's per-thread decode staging buffers on one single thread. */
#pragma omp master
      { nxt_ok = WPARC_FETCH(j + 1, pwNext); }

      /* x == 0 is the unshifted window (pwCur as-is); x >= 1 needs the
       * merge of pwCur with the tail of pwPrev. Folding the old serial
       * "head" in as x == 0 balances all W columns across the team.
       * Only ps->a is copied in: rrsortx needs it as the starting
       * permutation, whereas ps->d is fully overwritten by divc() and both
       * of psrev's arrays are fully overwritten by reversec(), so the three
       * other memcpys the old code did per column were dead. */
#pragma omp for schedule(dynamic, 1)
      for (size_t x = 0; x < W; x++) {
        size_t J = W - 1;
        uint64_t *cw;
        if (x == 0) {
          cw = pwCur;
        } else {
          wr64mrgsi(nrow, pwCur, pwPrev, w, x);
          cw = w;
        }
        pbwtad *ps = pb1[J - x];
        pbwtad *psrev = pb1rev[J - x];
        memcpy(ps->a, pb0[J - x]->a, nrow * sizeof *(ps->a));
        rrsortx(nrow, cw, ps->a, taux);
        reversec(ps, psrev, nrow);
        divc(nrow, cw, ps, pb0[J - x], psrev, pb0rev[J - x], W);
      }

#pragma omp single
      {
        PDUMP_SEQ_OFFSETR(0, W, pb1, W * j);
        SWAP(pb0, pb1);
        SWAP(pb0rev, pb1rev);
        /* rotate: prev <- cur <- next <- (recycled prev) */
        uint64_t *_spare = pwPrev;
        pwPrev = pwCur;
        pwCur = pwNext;
        pwNext = _spare;
        j++;
        active = nxt_ok;
      }
      /* implicit barrier at the end of `single`: every thread sees the new
       * `active`, `j`, buffer pointers and pb0/pb1 before the next round. */
    }

    FREE(w);
    FREE(taux);
  }
#undef WPARC_FETCH
  pw0 = pwPrev;
  pw1 = pwCur;
#if 1
  pbwtad *pp0, *pp1;
  pp0 = pb0[W - 1];
  pp1 = pb0[W - 2];

#ifdef BF2IOMODE_BCF
  ncol += _ncol;
  size_t wix = 0;
#endif
  for (j = j * W; j < ncol; j++) {
    // printf("entering last w at j %zu\n", j);
#if defined(BF2IOMODE_BM) || defined(BF2IOMODE_ENC)
    fgetcoli(fin, j, nrow, c0, ncol);
#elif defined(BF2IOMODE_BCF)
    for (size_t _i = 0; _i < nrow; _i++) {
      c0[_i] = (pw1[_i] >> wix) & 0x1;
    }
    wix++;
#else
#error UNDEFINED BEHAVIOUR
#endif
    cpbwtiLCP(nrow, j, c0, pp0, pp1);
    PDUMPR(j, pp1);
    SWAP(pp0, pp1);
  }

#endif
  /* The `return NULL;` used to sit *before* this block, so every call
   * leaked all 4*W pbwtads and every window buffer. */
  for (int j = 0; j < W; j++) {
    PBWTAD_FREE(pb0[j]);
    PBWTAD_FREE(pb1[j]);
    PBWTAD_FREE(pb0rev[j]);
    PBWTAD_FREE(pb1rev[j]);
  }
  FREE(pb0);
  FREE(pb1);
  FREE(pb0rev);
  FREE(pb1rev);
  FREE(pwPrev);
  FREE(pwCur);
  FREE(pwNext);
  pw0 = pw1 = NULL;
  (void)pw0;
  (void)pw1;
  FREE(c0);
  return NULL;
}

pbwtad **sbwparc_rrs(int fin, size_t nrow, size_t ncol) {
  // NOTE: right now I don't know what I need, so I'm keeping
  // everything in memory, we'll see later
  pbwtad **pb0 = malloc(W * sizeof(pbwtad *));
  pbwtad **pb0rev = malloc(W * sizeof(pbwtad *));
  pbwtad **pb1 = malloc(W * sizeof(pbwtad *));
  pbwtad **pb1rev = malloc(W * sizeof(pbwtad *));

  FILE *ffin = fdopen(fin, "r");
  uint8_t *c0 = malloc(nrow * sizeof *c0);
  pbwtad *p0 = pbwtad_new(nrow);
  for (int j = 0; j < nrow; j++) {
    p0->a[j] = j;
    p0->d[j] = 0;
  }

  // sbfgetcoln(fin, nrow, c0, ncol);
  fgetcoli(ffin, 0, nrow, c0, ncol);
  pb0[0] = cpbwt(nrow, c0, p0);
  pb0rev[0] = p0;
  PDUMP(0, pb0[0]);
  // PBWTAD_FREE(p0);

  for (int j = 1; j < W; j++) {
    fgetcoli(ffin, j, nrow, c0, ncol);
    pb0[j] = cpbwt(nrow, c0, pb0[j - 1]);
    pb0rev[j] = pbwtad_new(nrow);
    reversecprev(pb0[j], pb0[j - 1], pb0rev[j], nrow);
    PDUMP(j, pb0[j]);
  }

  /* swapping div values with LCP for window computation */
  for (int j = 0; j < W; j++) {
    swapdiv(pb0[j], nrow, j);
    swapdiv(pb0rev[j], nrow, j);
  }

  uint64_t *pw0 = malloc(nrow * sizeof *pw0);
  uint64_t *pw1 = malloc(nrow * sizeof *pw1);
  pbidx_t *aux = malloc(nrow * sizeof *aux);
  sbfgetcolw64rn(fin, nrow, pw0, ncol);

  // pb0 is now filled with computed values,
  // to allow reusing I need to fill pb1 with empty values
  for (int j = 0; j < W; j++) {
    pb1[j] = pbwtad_new(nrow);
    pb1rev[j] = pbwtad_new(nrow);
  }

  size_t j;

  for (j = 1; j * W <= ncol - W; j++) {
    sbfgetcolw64rn(fin, nrow, pw1, ncol);
    pbwtad *ps = pb1[W - 1];
    pbwtad *psrev = pb1rev[W - 1];
    memcpy(ps->a, pb0[W - 1]->a, nrow * sizeof *(ps->a));
    memcpy(ps->d, pb0[W - 1]->d, nrow * sizeof *(ps->d));
    memcpy(psrev->a, pb0rev[W - 1]->a, nrow * sizeof *(ps->a));
    memcpy(psrev->d, pb0rev[W - 1]->d, nrow * sizeof *(ps->d));
    rrsortx(nrow, pw1, ps->a, aux);
    reversec(ps, psrev, nrow);
    divc(nrow, pw1, ps, pb0[W - 1], psrev, pb0rev[W - 1], W);

#pragma omp parallel for shared(pw1, pw0, pb0, pb1, j)
    for (size_t x = 1; x < W; x++) {
      uint64_t *w = malloc(nrow * sizeof *w);
      size_t J = W - 1;

      wr64mrgsi(nrow, pw1, pw0, w, x);
      pbwtad *ps = pb1[J - x];
      pbwtad *psrev = pb1rev[J - x];
      memcpy(ps->a, pb0[J - x]->a, nrow * sizeof *(ps->a));
      memcpy(ps->d, pb0[J - x]->d, nrow * sizeof *(ps->d));
      memcpy(psrev->a, pb0[J - x]->a, nrow * sizeof *(psrev->a));
      memcpy(psrev->d, pb0rev[J - x]->a, nrow * sizeof *(psrev->d));
      rrsortx_noaux(nrow, w, ps->a);
      reversec(ps, psrev, nrow);
      divc(nrow, w, ps, pb0[J - x], psrev, pb0rev[J - x], W);
      FREE(w);
    }
    PDUMP_SEQ_OFFSETR(0, W, pb1, W * j);
    SWAP(pw0, pw1);
    SWAP(pb0, pb1);
    SWAP(pb0rev, pb1rev);
  }

  pbwtad *pp0, *pp1;
  pp0 = pb0[W - 1];
  pp1 = pb0[W - 2];

  for (j = j * W; j < ncol; j++) {
    fgetcoli(ffin, j, nrow, c0, ncol);
    cpbwtiLCP(nrow, j, c0, pp0, pp1);
    PDUMPR(j, pp1);
    SWAP(pp0, pp1);
  }
  for (int j = 0; j < W; j++) {
    PBWTAD_FREE(pb0[j]);
    PBWTAD_FREE(pb1[j]);
  }
  FREE(pb0);
  FREE(pb1);
  FREE(pw0);
  FREE(aux);
  FREE(c0);
  return NULL;
}

pbwtad **mbwparc_rrs(int fin, size_t nrow, size_t ncol) { // BPRM
  pbwtad **pb0 = malloc(W * sizeof(pbwtad *));
  pbwtad **pb0rev = malloc(W * sizeof(pbwtad *));
  pbwtad **pb1 = malloc(W * sizeof(pbwtad *));
  pbwtad **pb1rev = malloc(W * sizeof(pbwtad *));

  uint8_t *c0 = malloc(nrow * sizeof *c0);
  pbwtad *p0 = pbwtad_new(nrow);
  for (int j = 0; j < nrow; j++) {
    p0->a[j] = j;
    p0->d[j] = 0;
  }

  FILE *ffin = fdopen(fin, "r");
  fgetcoli(ffin, 0, nrow, c0, ncol);
  pb0[0] = cpbwt(nrow, c0, p0);
  pb0rev[0] = p0;
  PDUMP(0, pb0[0]);

  for (int j = 1; j < W; j++) {
    fgetcoli(ffin, j, nrow, c0, ncol);
    pb0[j] = cpbwt(nrow, c0, pb0[j - 1]);
    pb0rev[j] = pbwtad_new(nrow);
    reversecprev(pb0[j], pb0[j - 1], pb0rev[j], nrow);
    PDUMP(j, pb0[j]);
  }

  /* swapping div values with LCP for window computation */
  for (int j = 0; j < W; j++) {
    swapdiv(pb0[j], nrow, j);
    swapdiv(pb0rev[j], nrow, j);
  }

  uint64_t *pw0 = malloc(nrow * sizeof *pw0);
  uint64_t *pw1 = malloc(nrow * sizeof *pw1);
  pbidx_t *aux = malloc(nrow * sizeof *aux);
  // bfgetcolw64rn(fin, nrow, pw0, ncol);
  sbfgetcolw64rn_mmap(fin, nrow, pw0, ncol);

  // pb0 is now filled with computed values,
  // to allow reusing I need to fill pb1 with empty values
  for (int j = 0; j < W; j++) {
    pb1[j] = pbwtad_new(nrow);
    pb1rev[j] = pbwtad_new(nrow);
  }

  size_t j;

  for (j = 1; j * W <= ncol - W; j++) {
    sbfgetcolw64rn_mmap(fin, nrow, pw1, ncol);
    pbwtad *ps = pb1[W - 1];
    pbwtad *psrev = pb1rev[W - 1];
    memcpy(ps->a, pb0[W - 1]->a, nrow * sizeof *(ps->a));
    memcpy(ps->d, pb0[W - 1]->d, nrow * sizeof *(ps->d));
    memcpy(psrev->a, pb0rev[W - 1]->a, nrow * sizeof *(ps->a));
    memcpy(psrev->d, pb0rev[W - 1]->d, nrow * sizeof *(ps->d));
    rrsortx(nrow, pw1, ps->a, aux);
    reversec(ps, psrev, nrow);
    divc(nrow, pw1, ps, pb0[W - 1], psrev, pb0rev[W - 1], W);

#pragma omp parallel for shared(pw1, pw0, pb0, pb1, j)
    for (size_t x = 1; x < W; x++) {
      uint64_t *w = malloc(nrow * sizeof *w);
      size_t J = W - 1;

      wr64mrgsi(nrow, pw1, pw0, w, x);
      pbwtad *ps = pb1[J - x];
      pbwtad *psrev = pb1rev[J - x];
      memcpy(ps->a, pb0[J - x]->a, nrow * sizeof *(ps->a));
      memcpy(ps->d, pb0[J - x]->d, nrow * sizeof *(ps->d));
      memcpy(psrev->a, pb0[J - x]->a, nrow * sizeof *(psrev->a));
      memcpy(psrev->d, pb0rev[J - x]->a, nrow * sizeof *(psrev->d));
      rrsortx_noaux(nrow, w, ps->a);
      reversec(ps, psrev, nrow);
      divc(nrow, w, ps, pb0[J - x], psrev, pb0rev[J - x], W);
      FREE(w);
    }
    PDUMP_SEQ_OFFSETR(0, W, pb1, W * j);
    SWAP(pw0, pw1);
    SWAP(pb0, pb1);
    SWAP(pb0rev, pb1rev);
  }

  pbwtad *pp0, *pp1;
  pp0 = pb0[W - 1];
  pp1 = pb0[W - 2];

  for (j = j * W; j < ncol; j++) {
    fgetcoli(ffin, j, nrow, c0, ncol);
    cpbwtiLCP(nrow, j, c0, pp0, pp1);
    PDUMPR(j, pp1);
    SWAP(pp0, pp1);
  }
  for (int j = 0; j < W; j++) {
    PBWTAD_FREE(pb0[j]);
    PBWTAD_FREE(pb1[j]);
  }
  FREE(pb0);
  FREE(pb1);
  FREE(pw0);
  FREE(aux);
  FREE(c0);
  return NULL;
}

/* Staggered-parallel PBWT (SPR).
 *
 * Restructured in Step 5 of the optimization plan. Previously every one of
 * the W=64 "lanes" opened its *own* bcf_srs_t and streamed the whole file
 * from the beginning, so a BCF was BGZF-decompressed 64 times per run, and
 * the loop nest was lane-outer / round-inner.
 *
 * The observation that removes all of that: lane `l` at round `r` consumes
 * the W columns [l + rW, l + rW + W). Across all W lanes at a *fixed* round
 * `r` that is exactly the columns of blocks `r` and `r+1` — precisely the
 * two-window working set `wparc_rrs` (blockpar) already keeps. So the loop
 * nest is inverted (round outer, lanes as the inner `omp for`) and the 64
 * private readers are replaced by ONE sequential reader driving the same
 * three-buffer ring `wparc_rrs` uses, with `fgetcoliw64r` + `wr64mrgsi`.
 *
 * Window for lane `l` at round `r` starts at column rW + l:
 *   l == 0 : block r as-is                          (pwCur)
 *   l >= 1 : wr64mrgsi(pwNext=block r+1, pwPrev=block r, shift W-l)
 * which is bit-for-bit the same packing the old per-lane fgetcolwgri produced.
 *
 * Because a lane is no longer pinned to a thread, its carried state (the old
 * per-thread `pt0`/`pt0rev`) is now per LANE: pb[l] and pbrev[l]. pt1/pt1rev
 * stay per-thread — they are only a one-iteration snapshot of the previous
 * state consumed by divc().
 */
pbwtad **wstagparc_rrs(char *fpath, size_t nrow, size_t ncol) { // SPR
#if defined(BF2IOMODE_BCF)
  ncol = W;
#endif

  /* pb[0..W-1] are the lane states (pb[l] seeded with the PBWT after columns
   * 0..l-1); pb[W] is only the tail of the prologue's cpbwti chain and is not
   * a lane. pbrev[l] is lane l's reverse array, carried across rounds. */
  pbwtad **pb = malloc((W + 1) * sizeof(pbwtad *));
  pbwtad **pbrev = malloc(W * sizeof(pbwtad *));

  for (int j = 0; j < W + 1; j++) {
    pb[j] = pbwtad_new(nrow);
  }
  for (int j = 0; j < W; j++) {
    pbrev[j] = pbwtad_new(nrow);
  }
  for (size_t i = 0; i < nrow; i++) {
    pb[0]->a[i] = i;
    pb[0]->d[i] = 0;
  }

#ifdef BF2IOMODE_BCF
  bcf_srs_t *_pr = bcf_sr_init();
  bcf_sr_add_reader(_pr, fpath);
  void *fin = _pr;
#elif defined(BF2IOMODE_BM) || defined(BF2IOMODE_ENC)
  FILE *fin = fopen(fpath, "r");
  if (!fin) {
    perror("[spr]");
    exit(32);
  }
#else
#error UNDEFINED BEHAVIOUR
#endif

  uint8_t *c0 = malloc(nrow * sizeof *c0);
  fgetcoli(fin, 0, nrow, c0, ncol);
  cpbwti(nrow, c0, pb[0], pb[1]);

  for (size_t j = 1; j < W; j++) {
    fgetcoli(fin, j, nrow, c0, ncol);

    cpbwti(nrow, c0, pb[j], pb[j + 1]);

    PDUMP(j, pb[j + 1]);
  }
  for (size_t j = 0; j < W; j++) {
    swapdiv(pb[j], nrow, j - 1);
  }
  FREE(c0);

#ifdef BF2IOMODE_BCF
  /* The prologue consumed columns 0..W-1 from _pr, but the window scan has to
   * start again at block 0, so it gets a fresh reader — still exactly ONE for
   * the whole run, versus W of them before. */
  bcf_sr_destroy(_pr);
  bcf_srs_t *_sr = bcf_sr_init();
  bcf_sr_add_reader(_sr, fpath);
  fin = _sr;
#endif

  /* STAG_FETCH(idx, buf): read block `idx` (columns [idx*W, idx*W+W)) into
   * `buf`, returning how many columns were actually available (W for a full
   * block, less at EOF). Mirrors wparc_rrs's WPARC_FETCH but has to return the
   * *count*, not a boolean: at the last round lanes 1..k still have a complete
   * window when the trailing block holds only k columns. */
#if defined(BF2IOMODE_BM) || defined(BF2IOMODE_ENC)
  /* NOTE: unreachable in the BM binary (BM `stagpar` dispatches to
   * swstagparc_rrs/mwstagparc_rrs); kept compiling only. Full blocks only —
   * the BM readers have no partial-window/EOF signal to key off. */
#define STAG_FETCH(idx, buf)                                                   \
  ((((idx) + 1) * W <= ncol) ? (fgetcoliw64r(fin, (idx), nrow, (buf), ncol), W)\
                             : 0)
#elif defined(BF2IOMODE_BCF)
#define STAG_FETCH(idx, buf) ((int)fgetcoliw64r(fin, (idx), nrow, (buf), 0))
#else
#error UNDEFINED BEHAVIOUR
#endif

  /* Three window buffers rotated once per round: pwCur = block r,
   * pwNext = block r+1 (both consumed by this round's lanes), pwSpare = the
   * block r+2 the master thread reads *while* the team computes round r. */
  uint64_t *pwCur = malloc(nrow * sizeof *pwCur);
  uint64_t *pwNext = malloc(nrow * sizeof *pwNext);
  uint64_t *pwSpare = malloc(nrow * sizeof *pwSpare);

  int nCur = STAG_FETCH(0, pwCur);
  int nNext = (nCur == W) ? STAG_FETCH(1, pwNext) : 0;
  int nSpare = 0;
  size_t r = 0;
  int active = (nCur == W);

#pragma omp parallel
  {
    uint64_t *w = malloc(nrow * sizeof *w);
    pbidx_t *aux = malloc(nrow * sizeof *aux);
    pbwtad *pt1 = pbwtad_new(nrow);
    pbwtad *pt1rev = pbwtad_new(nrow);

    /* Seed every lane's reverse array once (the old code did this per lane at
     * the top of its private stream). Implicit barrier publishes them. */
#pragma omp for schedule(static)
    for (size_t lane = 0; lane < W; lane++) {
      reversec(pb[lane], pbrev[lane], nrow);
    }

    while (active) {
      /* I/O overlap, same shape as wparc_rrs: the master issues the read of
       * block r+2 into pwSpare and then joins the dynamically scheduled lane
       * loop. Readers write only pwSpare, lanes read only pwCur/pwNext, and
       * the `omp for`'s implicit barrier orders that write against the buffer
       * rotation in the `single` below. Pinned to master (not `single nowait`)
       * so iobcf.c's per-thread decode staging buffer stays on one thread. */
#pragma omp master
      { nSpare = (nNext == W) ? STAG_FETCH(r + 2, pwSpare) : 0; }

#pragma omp for schedule(dynamic, 1)
      for (size_t lane = 0; lane < W; lane++) {
        /* lane l needs columns up to rW + l + W - 1, i.e. the first `l`
         * columns of block r+1; skip the lanes the trailing block can't
         * complete. (Identical to the old per-lane loop's EOF condition.) */
        if (lane > (size_t)nNext)
          continue;
        uint64_t *cw;
        if (lane == 0) {
          cw = pwCur; /* shift 0 would be a shift-by-64 in wr64mrgsi (UB) */
        } else {
          wr64mrgsi(nrow, pwNext, pwCur, w, W - lane);
          cw = w;
        }
        pbwtad *p = pb[lane];
        pbwtad *prv = pbrev[lane];
        /* divc() reads only ppr->d and pprrev->a of the previous state; the
         * other two memcpys the old body did (pt1->a, pt1rev->d) were dead. */
        memcpy(pt1->d, p->d, nrow * sizeof *(p->d));
        memcpy(pt1rev->a, prv->a, nrow * sizeof *(prv->a));

        rrsortx(nrow, cw, p->a, aux);
        reversec(p, prv, nrow);
        divc(nrow, cw, p, pt1, prv, pt1rev, W);
      }

#pragma omp single
      {
#ifdef DBDUMP
        /* Emitted here, serially and in lane order, instead of from an
         * `omp critical` inside the parallel body: the dump is now
         * deterministic (columns strictly increasing across rounds) rather
         * than interleaved differently on every run. */
        if (DO_DUMP) {
          size_t lmax = (nNext < W - 1) ? (size_t)nNext : (size_t)(W - 1);
          for (size_t lane = 0; lane <= lmax; lane++) {
            PDUMPR(r * W + lane + W - 1, pb[lane]);
          }
        }
#endif
        r++;
        nCur = nNext;
        nNext = nSpare;
        uint64_t *_sp = pwCur;
        pwCur = pwNext;
        pwNext = pwSpare;
        pwSpare = _sp;
        active = (nCur == W);
      }
      /* implicit barrier of `single`: everyone sees the rotated buffers,
       * the new counts and `active` before the next round. */
    }

    PBWTAD_FREE(pt1);
    PBWTAD_FREE(pt1rev);
    FREE(aux);
    FREE(w);
  }
#undef STAG_FETCH

#ifdef BF2IOMODE_BCF
  bcf_sr_destroy(_sr);
#else
  fclose(fin);
#endif

  for (int j = 0; j < W + 1; j++) {
    PBWTAD_FREE(pb[j]);
  }
  for (int j = 0; j < W; j++) {
    PBWTAD_FREE(pbrev[j]);
  }
  FREE(pb);
  FREE(pbrev);
  FREE(pwCur);
  FREE(pwNext);
  FREE(pwSpare);
  /* Consistent with linc/wparc_rrs: main() only ever assigns the return value
   * and never reads it, so nothing is returned and nothing is leaked. The old
   * code returned `pb` holding just the 64 *seed* states while the real
   * per-lane results died with the thread-locals. */
  return NULL;
}

pbwtad **swstagparc_rrs(char *fpath, size_t nrow, size_t ncol) { // SPRS
#if defined(BF2IOMODE_BCF)
  fprintf(stderr, "Not valid version for BCF/VCF\n");
  exit(3);
#endif

  pbwtad **pb = malloc((W + 1) * sizeof(pbwtad *));
  pbwtad *pbprev = pbwtad_new(nrow);

  for (int j = 0; j < W + 1; j++) {
    pb[j] = pbwtad_new(nrow);
  }
  for (size_t i = 0; i < nrow; i++) {
    pb[0]->a[i] = i;
    pb[0]->d[i] = 0;
  }

  FILE *fin = fopen(fpath, "r");
  if (!fin) {
    perror("[spr]");
    exit(32);
  }

  uint8_t *c0 = malloc(nrow * sizeof *c0);
  fgetcoli(fin, 0, nrow, c0, ncol);
  cpbwti(nrow, c0, pb[0], pb[1]);

  for (size_t j = 1; j < W; j++) {
    fgetcoli(fin, j, nrow, c0, ncol);
    cpbwti(nrow, c0, pb[j], pb[j + 1]);
    PDUMP(j, pb[j + 1]);
  }
  for (size_t j = 0; j < W; j++) {
    swapdiv(pb[j], nrow, j - 1);
  }

#pragma omp parallel
  {

    int fd = open(fpath, O_RDONLY);

    size_t tid = omp_get_thread_num();
    size_t nthreads = omp_get_num_threads();

    size_t base = W / nthreads;
    size_t rem = W % nthreads;

    size_t start = tid * base + (tid < rem ? tid : rem);
    size_t count = base + (tid < rem ? 1 : 0);
    pbwtad *pt0 = pbwtad_new(nrow);
    pbwtad *pt0rev = pbwtad_new(nrow);
    pbwtad *pt1 = pbwtad_new(nrow);
    pbwtad *pt1rev = pbwtad_new(nrow);
    pbidx_t *aux = malloc(nrow * sizeof *aux);
    uint64_t *pw = malloc(nrow * sizeof *pw);

    for (size_t offset = 0; offset < count; offset++) {
      size_t lane = start + offset;
      memcpy(pt0->a, pb[lane]->a, nrow * sizeof *(pb[lane]->a));
      memcpy(pt0->d, pb[lane]->d, nrow * sizeof *(pb[lane]->d));
      reversec(pt0, pt0rev, nrow);

      for (size_t j = lane; j + W <= ncol; j += W) {
        spfgetcolwgri(fd, j, nrow, pw, ncol, W);
        memcpy(pt1->a, pt0->a, nrow * sizeof *(pt0->a));
        memcpy(pt1->d, pt0->d, nrow * sizeof *(pt0->d));
        memcpy(pt1rev->a, pt0rev->a, nrow * sizeof *(pt0rev->a));
        memcpy(pt1rev->d, pt0rev->d, nrow * sizeof *(pt0rev->d));

        rrsortx(nrow, pw, pt0->a, aux);
        reversec(pt0, pt0rev, nrow);
        divc(nrow, pw, pt0, pt1, pt0rev, pt1rev, W);

#ifdef DBDUMP
#pragma omp critical
        {
          PDUMPR(j + W - 1, pt0);
        }
#endif
      }
    }

    PBWTAD_FREE(pt0);
    PBWTAD_FREE(pt1);
    PBWTAD_FREE(pt0rev);
    PBWTAD_FREE(pt1rev);
    FREE(pt0);
    FREE(pt1);
    FREE(pt0rev);
    FREE(pt1rev);
    FREE(aux);
    FREE(pw);
  }
  FREE(c0);
  return pb;
}

pbwtad **mwstagparc_rrs(char *fpath, size_t nrow, size_t ncol) { // SPRM
#if defined(BF2IOMODE_BCF)
  fprintf(stderr, "Not valid version for BCF/VCF\n");
  exit(3);
#endif

  pbwtad **pb = malloc((W + 1) * sizeof(pbwtad *));
  pbwtad *pbprev = pbwtad_new(nrow);

  for (int j = 0; j < W + 1; j++) {
    pb[j] = pbwtad_new(nrow);
  }
  for (size_t i = 0; i < nrow; i++) {
    pb[0]->a[i] = i;
    pb[0]->d[i] = 0;
  }

  FILE *fin = fopen(fpath, "r");
  if (!fin) {
    perror("[spr]");
    exit(32);
  }

  uint8_t *c0 = malloc(nrow * sizeof *c0);
  fgetcoli(fin, 0, nrow, c0, ncol);
  cpbwti(nrow, c0, pb[0], pb[1]);
  // PDUMP(0, pb[1]);

  for (size_t j = 1; j < W; j++) {
    fgetcoli(fin, j, nrow, c0, ncol);
    cpbwti(nrow, c0, pb[j], pb[j + 1]);
    PDUMP(j, pb[j + 1]);
  }
  for (size_t j = 0; j < W; j++) {
    swapdiv(pb[j], nrow, j - 1);
  }

  int fd = open(fpath, O_RDONLY);
  if (fd < 0) {
    perror("open");
    exit(EXIT_FAILURE);
  }
  struct stat st;
  if (fstat(fd, &st) < 0) {
    perror("fstat");
    exit(EXIT_FAILURE);
  }
  if (st.st_size == 0) {
    fprintf(stderr, "Error: File is empty\n");
    close(fd);
    exit(EXIT_FAILURE);
  }

  static uint8_t *fdmm = NULL;
  fdmm = mmap(NULL, st.st_size, PROT_READ, __MMAP_FLAGS, fd, 0);
  if (fdmm == MAP_FAILED) {
    perror("mmap");
    exit(EXIT_FAILURE);
  }
  close(fd);

#pragma omp parallel
  {
    size_t tid = omp_get_thread_num();
    size_t nthreads = omp_get_num_threads();

    size_t base = W / nthreads;
    size_t rem = W % nthreads;

    size_t start = tid * base + (tid < rem ? tid : rem);
    size_t count = base + (tid < rem ? 1 : 0);
    pbwtad *pt0 = pbwtad_new(nrow);
    pbwtad *pt0rev = pbwtad_new(nrow);
    pbwtad *pt1 = pbwtad_new(nrow);
    pbwtad *pt1rev = pbwtad_new(nrow);
    pbidx_t *aux = malloc(nrow * sizeof *aux);
    uint64_t *pw = malloc(nrow * sizeof *pw);

    for (size_t offset = 0; offset < count; offset++) {
      size_t lane = start + offset;
      memcpy(pt0->a, pb[lane]->a, nrow * sizeof *(pb[lane]->a));
      memcpy(pt0->d, pb[lane]->d, nrow * sizeof *(pb[lane]->d));
      reversec(pt0, pt0rev, nrow);

      for (size_t j = lane; j + W <= ncol; j += W) {

        fgetcolwgri_mmap(fdmm, j, nrow, pw, ncol, W);
        memcpy(pt1->a, pt0->a, nrow * sizeof *(pt0->a));
        memcpy(pt1->d, pt0->d, nrow * sizeof *(pt0->d));
        memcpy(pt1rev->a, pt0rev->a, nrow * sizeof *(pt0rev->a));
        memcpy(pt1rev->d, pt0rev->d, nrow * sizeof *(pt0rev->d));

        rrsortx(nrow, pw, pt0->a, aux);
        reversec(pt0, pt0rev, nrow);
        divc(nrow, pw, pt0, pt1, pt0rev, pt1rev, W);

#ifdef DBDUMP
#pragma omp critical
        {
          PDUMPR(j + W - 1, pt0);
        }
#endif
      }
    }

    PBWTAD_FREE(pt0);
    PBWTAD_FREE(pt1);
    PBWTAD_FREE(pt0rev);
    PBWTAD_FREE(pt1rev);
    FREE(pt0);
    FREE(pt1);
    FREE(pt0rev);
    FREE(pt1rev);
    FREE(aux);
    FREE(pw);
  }
  FREE(c0);
  return pb;
}

int main(int argc, char *argv[]) {
#if defined(BF2IOMODE_BM)
  char _usage_args_[] =
      "[sampled|linear|blockpar|stagpar]-[syscall|mmap] FILE\n";
  if (strcmp(argv[1], "linear-syscall") == 0) {
  } else if (strcmp(argv[1], "linear-mmap") == 0) {
  } else if (strcmp(argv[1], "sample-syscall") == 0) {
  } else if (strcmp(argv[1], "sample-mmap") == 0) {
  } else if (strcmp(argv[1], "blockpar-syscall") == 0) {
  } else if (strcmp(argv[1], "blockpar-mmap") == 0) {
  } else if (strcmp(argv[1], "stagpar-syscall") == 0) {
  } else if (strcmp(argv[1], "stagpar-mmap") == 0) {
#elif defined(BF2IOMODE_BCF)
  char _usage_args_[] = "[sampled|linear|blockpar|stagpar] FILE\n";
  if (strcmp(argv[1], "linear") == 0) {
  } else if (strcmp(argv[1], "sampled") == 0) {
  } else if (strcmp(argv[1], "blockpar") == 0) {
  } else if (strcmp(argv[1], "stagpar") == 0) {
#endif
  } else {
    fprintf(stderr, "Wrong mode \"%s\"\nUsage: %s %s", argv[1], argv[0],
            _usage_args_);
    return EXIT_FAILURE;
  }
  if (argc < 3) {
    fprintf(stderr, "Missing input file\nUsage: %s %s", argv[0], _usage_args_);
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
  if (nrow > PBIDX_MAX) {
    fprintf(stderr,
            "sp-pbwt: nrow (%zu) exceeds pbidx_t capacity (%u); rebuild "
            "with a wider pbidx_t\n",
            nrow, (unsigned)PBIDX_MAX);
    return EXIT_FAILURE;
  }
  pbwtad **r;

  if (argc > 3 && strcmp(argv[3], "DUMP") == 0) {
    DO_DUMP = 1;
  }

#if defined(BF2IOMODE_BM)
  if (strcmp(argv[1], "linear-syscall") == 0) {
    TRACE(sblinc(fd, nrow, ncol), r);
  } else if (strcmp(argv[1], "linear-mmap") == 0) {
    TRACE(mblinc(fd, nrow, ncol), r);
  } else if (strcmp(argv[1], "sample-syscall") == 0) {
    TRACE(swbapproxc_rrs(fd, nrow, ncol), r);
  } else if (strcmp(argv[1], "sample-mmap") == 0) {
    TRACE(mwbapproxc_rrs(fd, nrow, ncol), r);
  } else if (strcmp(argv[1], "blockpar-syscall") == 0) {
    TRACE(sbwparc_rrs(fd, nrow, ncol), r);
  } else if (strcmp(argv[1], "blockpar-mmap") == 0) {
    TRACE(mbwparc_rrs(fd, nrow, ncol), r);
  } else if (strcmp(argv[1], "stagpar-syscall") == 0) {
    TRACE(swstagparc_rrs(argv[2], nrow, ncol), r);
  } else if (strcmp(argv[1], "stagpar-mmap") == 0) {
    TRACE(mwstagparc_rrs(argv[2], nrow, ncol), r);
#elif defined(BF2IOMODE_BCF)
  if (strcmp(argv[1], "linear") == 0) {
    TRACE(linc(fin, nrow, ncol), r);
  } else if (strcmp(argv[1], "sampled") == 0) {
    TRACE(wapproxc_rrs(fin, nrow, ncol), r);
  } else if (strcmp(argv[1], "blockpar") == 0) {
    TRACE(wparc_rrs(fin, nrow, ncol), r);
  } else if (strcmp(argv[1], "stagpar") == 0) {
    TRACE(wstagparc_rrs(argv[2], nrow, ncol), r);
#endif
  }
#if defined(BF2IOMODE_BCF)
  bcf_sr_destroy(sr);
#endif
  return EXIT_SUCCESS;
}
