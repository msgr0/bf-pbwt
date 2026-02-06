/* vim: set ft=c */
#include "io.h"
#include "tracing.h"
#include <assert.h>
#include <fcntl.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#ifdef BF2IOMODE_BCF
#include "htslib/synced_bcf_reader.h"
#endif

#define W 64

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

#define PDUMP_SEQR(s, e, p)                                                    \
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
#define PDUMP_SEQ_OFFSETR(s, e, p, offset)                                     \
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

typedef struct {
    size_t *tree;   // The array holding the tree nodes
    size_t n;       // Original array size
    size_t size;    // Power of 2 size (capacity)
} SegTree;

#define BLOCK_SIZE 32

typedef struct {
    size_t *summary_table; // Sparse Table for the blocks
    size_t *block_mins;    // Min value of each block
    size_t *original_data; // Pointer to input data (not owned)
    size_t n_blocks;       // Number of blocks
    size_t n;              // Original size
    size_t *logs;          // Precomputed logs
    size_t st_cols;        // Width of sparse table
} LinearRMQ;

// 1. ALLOCATE (Call ONCE at startup)
LinearRMQ* rmq_alloc(size_t max_n) {
    LinearRMQ *rmq = malloc(sizeof(LinearRMQ));
    
    // Max blocks needed
    size_t max_blocks = (max_n + BLOCK_SIZE - 1) / BLOCK_SIZE;
    
    rmq->block_mins = malloc(max_blocks * sizeof(size_t));
    
    // Precompute logs for the summary table
    rmq->logs = malloc((max_blocks + 1) * sizeof(size_t));
    rmq->logs[1] = 0;
    for (size_t i = 2; i <= max_blocks; i++)
        rmq->logs[i] = rmq->logs[i / 2] + 1;
    
    // Max log height
    size_t K = rmq->logs[max_blocks];
    rmq->st_cols = K + 1;
    
    // Sparse table is flattened: [col][row] -> index = col * max_blocks + row
    // We allocate worst-case size.
    rmq->summary_table = malloc(rmq->st_cols * max_blocks * sizeof(size_t));
    
    return rmq;
}

// 2. BUILD (Call every column) - O(N)
void rmq_build(LinearRMQ *rmq, size_t *data, size_t n) {
    rmq->original_data = data;
    rmq->n = n;
    rmq->n_blocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;

    // A. Compute Block Minimums - Linear Scan O(N)
    for (size_t b = 0; b < rmq->n_blocks; b++) {
        size_t start = b * BLOCK_SIZE;
        size_t end = start + BLOCK_SIZE;
        if (end > n) end = n;
        
        size_t min = SIZE_MAX;
        for (size_t i = start; i < end; i++) {
            if (data[i] < min) min = data[i];
        }
        rmq->block_mins[b] = min;
    }

    // B. Build Sparse Table on Block Minimums - O(N/B log(N/B)) -> Very Fast
    // Fill level 0 (original block mins)
    for (size_t i = 0; i < rmq->n_blocks; i++) {
        rmq->summary_table[i] = rmq->block_mins[i]; // row 0
    }

    size_t K = rmq->logs[rmq->n_blocks];
    size_t stride = rmq->n_blocks; // stride for accessing next level in flattened array

    for (size_t j = 1; j <= K; j++) {
        size_t prev_offset = (j - 1) * stride;
        size_t curr_offset = j * stride;
        size_t range = 1 << (j - 1);
        
        for (size_t i = 0; i + (1 << j) <= rmq->n_blocks; i++) {
            size_t val1 = rmq->summary_table[prev_offset + i];
            size_t val2 = rmq->summary_table[prev_offset + i + range];
            rmq->summary_table[curr_offset + i] = (val1 < val2) ? val1 : val2;
        }
    }
}

// Helper: Standard Sparse Table Query on the blocks
static inline size_t query_summary(LinearRMQ *rmq, size_t L_block, size_t R_block) {
    size_t j = rmq->logs[R_block - L_block + 1];
    size_t stride = rmq->n_blocks;
    size_t range = 1 << j;
    
    size_t val1 = rmq->summary_table[j * stride + L_block];
    size_t val2 = rmq->summary_table[j * stride + R_block - range + 1];
    
    return (val1 < val2) ? val1 : val2;
}

// 3. QUERY (Call inside recover_div) - O(1) amortized
size_t rmq_query(LinearRMQ *rmq, size_t L, size_t R) {
    size_t b_L = L / BLOCK_SIZE;
    size_t b_R = R / BLOCK_SIZE;
    size_t min_val = SIZE_MAX;

    if (b_L == b_R) {
        // Case 1: Range inside single block -> Linear Scan (tiny, max 32 ops)
        for (size_t i = L; i <= R; i++) {
            if (rmq->original_data[i] < min_val) min_val = rmq->original_data[i];
        }
    } else {
        // Case 2: Range spans blocks
        
        // 1. Scan Suffix of first block
        size_t end_L = (b_L + 1) * BLOCK_SIZE;
        for (size_t i = L; i < end_L; i++) {
            if (rmq->original_data[i] < min_val) min_val = rmq->original_data[i];
        }

        // 2. Scan Prefix of last block
        size_t start_R = b_R * BLOCK_SIZE;
        for (size_t i = start_R; i <= R; i++) {
            if (rmq->original_data[i] < min_val) min_val = rmq->original_data[i];
        }

        // 3. Query Summary Table for blocks strictly in between
        if (b_L + 1 <= b_R - 1) {
            size_t block_min = query_summary(rmq, b_L + 1, b_R - 1);
            if (block_min < min_val) min_val = block_min;
        }
    }
    
    return min_val;
}

// 4. CLEANUP
void rmq_free(LinearRMQ *rmq) {
    free(rmq->summary_table);
    free(rmq->block_mins);
    free(rmq->logs);
    free(rmq);
}

// 1. ALLOCATE (Call this ONCE at program start)
SegTree* st_alloc(size_t max_n) {
    SegTree *st = malloc(sizeof(SegTree));
    
    // Calculate next power of 2 >= max_n
    st->size = 1;
    while (st->size < max_n) st->size <<= 1;
    
    // Allocate 2*size. 
    // We use calloc to ensure safety, or just malloc since we overwrite.
    st->tree = malloc(2 * st->size * sizeof(size_t));
    return st;
}

void st_build(SegTree *st, size_t *data, size_t n) {
    st->n = n;
    
    // A. Fill the leaves (starting at index 'size')
    for (size_t i = 0; i < n; i++) {
        st->tree[st->size + i] = data[i];
    }
    // Fill the rest of the padding leaves with MAX value (infinity)
    for (size_t i = n; i < st->size; i++) {
        st->tree[st->size + i] = SIZE_MAX;
    }

    // B. Build internal nodes (working backwards from parents of leaves)
    for (size_t i = st->size - 1; i > 0; i--) {
        size_t left_child  = st->tree[2 * i];
        size_t right_child = st->tree[2 * i + 1];
        st->tree[i] = (left_child < right_child) ? left_child : right_child;
    }
}
// 3. QUERY (Call this inside 'recover_div')
// Returns min value in range [L, R] inclusive.
size_t st_query(SegTree *st, size_t L, size_t R) {
    size_t min_val = SIZE_MAX;
    // Shift indices to leaf level
    L += st->size;
    R += st->size;

    while (L <= R) {
        // If L is a right child, it means the parent includes values to the left 
        // that we DON'T want. So take L and move right.
        if (L % 2 == 1) {
            if (st->tree[L] < min_val) min_val = st->tree[L];
            L++;
        }
        
        // If R is a left child, take R and move left.
        if (R % 2 == 0) {
            if (st->tree[R] < min_val) min_val = st->tree[R];
            R--;
        }
        
        // Move up to parents
        L /= 2;
        R /= 2;
    }
    return min_val;
}
// 4. CLEANUP (Call ONCE at program end)
void st_free(SegTree *st) {
    free(st->tree);
    free(st);
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
 * without using externally allocated `aux[n]` auxiliary array.
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
    rev->a[p->a[i]] = i;
    rev->d[p->a[i]] = p->d[i]; // == p->d[rev->a[i]]
    assert(rev->d[p->a[i]] == p->d[rev->a[p->a[i]]]);
  }
}
void reversecprev(pbwtad *p, pbwtad *pp, pbwtad *rev, size_t n) {
  for (size_t i = 0; i < n; i++) {
    rev->a[p->a[i]] = i;
    rev->d[p->a[i]] = pp->d[i]; // == p->d[rev->a[i]]
    // assert(rev->d[p->a[i]] == pp->d[rev->a[p->a[i]]]);
  }
}

size_t recover_div_st(size_t w, size_t i, size_t i0, 
                   pbwtad *pprrev, LinearRMQ *st) {
    
    // Get ranks in previous sort
    size_t rank1 = pprrev->a[i];
    size_t rank2 = pprrev->a[i0];

    // Identify range [L, R]
    // We want the min strictly BETWEEN these ranks.
    // Logic: min(rank1, rank2) + 1  TO  max(rank1, rank2)
    
    size_t L, R;
    if (rank1 < rank2) {
        L = rank1 + 1;
        R = rank2;
    } else {
        L = rank2 + 1;
        R = rank1;
    }
    // L = L + 1; 

    // Safety check (should not trigger if lo// Safety check (should not trigger if logic is correct)
    if (L > R) return w; // No gap? Identical items

    size_t min_d = rmq_query(st, L, R);

    // If min_d is still SIZE_MAX, it means we queried out of bounds 
    // or the array wasn't filled. Fallback to w.
    if (min_d == SIZE_MAX) return w;

    return w + min_d;}
/*
 * Recover divergence of a match possibly longer than W.
 * Iterates over the range between previous and current row in p->a
 * using the previous pbwt array. Compute correct divergence using
 * reverse arrays computed by reversec
 * WARN: prev (current pbwt reverse) is not used, consider not memcpy in
 * wapprox computation if not needed.
 */
size_t recover_div(size_t n, size_t w, size_t i, size_t i0, uint64_t *c,
                   pbwtad *p, pbwtad *ppr, pbwtad *prev, pbwtad *pprrev) {

  size_t d;
  size_t j = pprrev->a[i];
  size_t min = ppr->d[j];

  for (size_t j = (pprrev->a[i0]) + 1; j < (pprrev->a[i]); j++) {
    if (ppr->d[j] < min) {
      min = ppr->d[j];
    }
  }
  d = w + min;
  return d;
}

/*
 * compute the divergence of the first w64 windows
 */
void divc0(size_t n, uint64_t *c, pbwtad *p) {
  uint64_t x = 0;
  size_t div = 0;
  p->d[0] = 0;
  for (size_t i = 1; i < n; i++) {
    x = c[p->a[i]] ^ c[p->a[i - 1]];
    p->d[i] = x ? __builtin_clzll(x) : 64;
  }
}

/*
 * computes the divergence of the last <=64 window (hence it could be padded
 */
// void divc_last(size_t n, uint64_t *c, pbwtad *p, pbwtad *ppr, pbwtad *prev,
//                pbwtad *pprrev,
//                size_t w) { // WARN: afteer some correction code seems to be
//                            // similar to divc(), check and merge
//   static int8_t kk = 0;
//   uint64_t x = 0;
//   // size_t w = 64;
//   size_t div;
//   p->d[0] = 0;
//   for (size_t i = 1; i < n; i++) {
//     x = c[p->a[i]] ^ c[p->a[i - 1]];
//     div = x ? __builtin_clzll(x) - (W - w) : w;
//     // if (div > w) div = 64 - div;
//     p->d[i] = (div >= w) ? recover_div(n, w, p->a[i], p->a[i - 1], c, p, ppr,
//                                        prev, pprrev)
//                          : div;
//   }
//   kk += W;
// }
//

void divc_st(size_t n, uint64_t *c, pbwtad *p, pbwtad *ppr, pbwtad *prev,
          pbwtad *pprrev, size_t wi, LinearRMQ *st) {
  // c contains 64bit-encoded ints
  // xor of each c[s[i]] and its preceeding;
  // x[0] contains no information, previous x information is discarded;
  // here 64 is the size of the window
  static int8_t kk = 0;
  uint64_t x = 0;
  size_t w = wi ? wi : W;
  size_t div;
  p->d[0] = 0;

  for (size_t i = 1; i < n; i++) {
    x = c[p->a[i]] ^ c[p->a[i - 1]];
    div = x ? __builtin_clzll(x) : w;
    p->d[i] = (div >= w) ? recover_div_st(w, p->a[i], p->a[i - 1], pprrev, st): div;
  }
  kk += W;
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
  static int8_t kk = 0;
  uint64_t x = 0;
  size_t w = wi ? wi : W;
  size_t div;
  p->d[0] = 0;

  for (size_t i = 1; i < n; i++) {
    x = c[p->a[i]] ^ c[p->a[i - 1]];
    div = x ? __builtin_clzll(x) : w;
    p->d[i] = (div >= w) ? recover_div(n, w, p->a[i], p->a[i - 1], c, p, ppr,
                                       prev, pprrev)
                         : div;
  }

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

// msgr0's version (Durbin's linear pbwt)
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

/* BUG: (possibily?, needs testing)
 * use cpbwt std version (withtout *LCP) and then swap divergence
 * with swapdiv if necessary.
 */
static pbwtad *cpbwtLCP(size_t n, uint8_t *restrict c, pbwtad *restrict p) {
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
  for (i = 0; i < n; i++) {
    size_t idx = p->a[i];
    size_t ddx = p->d[i];

    f = (ddx > f) ? ddx : f;
    g = (ddx > g) ? ddx : g;

    size_t mask = c[idx];
    o[q] = idx;
    ret->a[r] = idx;
    h[q] = g;
    ret->d[r] = f;

    f &= -mask;       // f = 0 if mask == 0, unchanged if mask == 1
    g &= -(1 - mask); // g = 0 if mask == 1, unchanged if mask == 0
    q += mask;        // Increment q if mask is 1
    r += mask ^ 1;    // Increment r if mask is 0
  }
  memcpy(ret->a + r, o, q * sizeof(size_t));
  memcpy(ret->d + r, h, q * sizeof(size_t));

  // correct div computation for DIV as LCP
  //
  for (size_t t = 0; t < n; t++) {
    ret->d[t] = 1 + k - ret->d[t];
  }
  k++;

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

/* BUG: (possibily?, needs testing)
 * use cpbwt std version (withtout *LCP) and then swap divergence
 * with swapdiv if necessary.
 */
static int cpbwtiLCP(size_t n, size_t k, uint8_t *restrict c,
                     pbwtad *restrict pp, pbwtad *restrict pc) {
  static size_t *o = NULL;
  static size_t *h = NULL;
  // static size_t k = 1;

  swapdiv(pp, n, k - 1);
  size_t nrow = n;
  // PDUMP(k, pp);
  if (!o)
    o = malloc(n * sizeof *o);
  if (!h)
    h = malloc(n * sizeof *h);

  size_t r = 0, q = 0;
  size_t f = k + 1, g = k + 1;

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

  // PDUMP(k, pc);
  swapdiv(pc, n, k);
  return 1;
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
pbwtad **wapproxcst_rrs(void *fin, size_t nrow, size_t ncol) { // ARS
  LinearRMQ *st = rmq_alloc(nrow); // Allocate ONCE

  // Compute the bit-packed windows
  uint64_t *w64 =
      malloc(nrow * sizeof *w64); // window data collected by fgetcoliw64r
  size_t *aux = malloc(nrow * sizeof *aux);

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
    memcpy(pbwtPr->a, pbwt->a, nrow * sizeof *(pbwt->a));
    memcpy(pbwtPr->d, pbwt->d, nrow * sizeof *(pbwt->d));
    memcpy(pbwtPrRev->a, pbwtRev->a, nrow * sizeof *(pbwtRev->a));
    memcpy(pbwtPrRev->d, pbwtRev->d, nrow * sizeof *(pbwtRev->d));
    rmq_build(st, pbwtPr->d, nrow);
    rrsortx(nrow, w64, pbwt->a,
            aux); // radix sorting pbwt->a with auxiliary array
    reversec(pbwt, pbwtRev,
             nrow); // computing reversec after sorting the new array
    divc_st(nrow, w64, pbwt, pbwtPr, pbwtRev, pbwtPrRev, W, st);
    // divc(nrow, w64, pbwt, pbwtPr, pbwtRev, pbwtPrRev,
         // W); // FIXME: check divc comment, pbwtRev could be removed.
    PDUMPR(W * (j + 1) - 1, pbwt);
    // CDUMP(W * (j+1) - 1, w64);
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
  // memcpy(p1->a, p0->a, nrow * sizeof *(p1->a));
  memcpy(pbwtPr->a, pbwt->a, nrow * sizeof *(pbwt->a));
  memcpy(pbwtPr->d, pbwt->d, nrow * sizeof *(pbwt->d));

  rmq_build(st, pbwtPr->d, nrow);
  rrsortx(nrow, w64, pbwt->a, aux);
  memcpy(pbwtPrRev->a, pbwtRev->a, nrow * sizeof *(pbwtRev->a));
  memcpy(pbwtPrRev->d, pbwtRev->d, nrow * sizeof *(pbwtRev->d));
  reversec(pbwt, pbwtRev, nrow);
  // PDUMPR(W * (j + 1) - 1, pbwt);
  // printf("last w has width:%zu\n", ncol - j);
  divc_st(nrow, w64, pbwt, pbwtPr, pbwtRev, pbwtPrRev, ncol-j, st);
  // divc(nrow, w64, pbwt, pbwtPr, pbwtRev, pbwtPrRev, ncol - j);
  PDUMPR(ncol - 1, pbwt);

  PBWTAD_FREE(pbwt);
  // FREE(c0);
  FREE(pbwt);
  FREE(pbwtRev);
  FREE(pbwtPr);
  FREE(pbwtPrRev);
  FREE(w64);
  FREE(aux);
  free(st);
  return NULL;
}

pbwtad **wapproxc_rrs(void *fin, size_t nrow, size_t ncol) { // ARS

  // Compute the bit-packed windows
  uint64_t *w64 =
      malloc(nrow * sizeof *w64); // window data collected by fgetcoliw64r
  size_t *aux = malloc(nrow * sizeof *aux);

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
    memcpy(pbwtPr->a, pbwt->a, nrow * sizeof *(pbwt->a));
    memcpy(pbwtPr->d, pbwt->d, nrow * sizeof *(pbwt->d));
    memcpy(pbwtPrRev->a, pbwtRev->a, nrow * sizeof *(pbwtRev->a));
    memcpy(pbwtPrRev->d, pbwtRev->d, nrow * sizeof *(pbwtRev->d));
    rrsortx(nrow, w64, pbwt->a,
            aux); // radix sorting pbwt->a with auxiliary array
    reversec(pbwt, pbwtRev,
             nrow); // computing reversec after sorting the new array
    divc(nrow, w64, pbwt, pbwtPr, pbwtRev, pbwtPrRev,
         W); // FIXME: check divc comment, pbwtRev could be removed.
    PDUMPR(W * (j + 1) - 1, pbwt);
    // CDUMP(W * (j+1) - 1, w64);
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
  // memcpy(p1->a, p0->a, nrow * sizeof *(p1->a));
  memcpy(pbwtPr->a, pbwt->a, nrow * sizeof *(pbwt->a));
  memcpy(pbwtPr->d, pbwt->d, nrow * sizeof *(pbwt->d));

  rrsortx(nrow, w64, pbwt->a, aux);
  memcpy(pbwtPrRev->a, pbwtRev->a, nrow * sizeof *(pbwtRev->a));
  memcpy(pbwtPrRev->d, pbwtRev->d, nrow * sizeof *(pbwtRev->d));
  reversec(pbwt, pbwtRev, nrow);
  // PDUMPR(W * (j + 1) - 1, pbwt);
  // printf("last w has width:%zu\n", ncol - j);
  divc(nrow, w64, pbwt, pbwtPr, pbwtRev, pbwtPrRev, ncol - j);
  PDUMPR(ncol - 1, pbwt);

  PBWTAD_FREE(pbwt);
  // FREE(c0);
  FREE(pbwt);
  FREE(pbwtRev);
  FREE(pbwtPr);
  FREE(pbwtPrRev);
  FREE(w64);
  FREE(aux);
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

pbwtad **swbapproxc_rrs(int fin, size_t nrow, size_t ncol) { // BARS
  // Compute the bit-packed windows
  uint64_t *w64 = malloc(nrow * sizeof *w64);
  // uint64_t *xor = malloc(nrow * sizeof *xor);
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
  PDUMPR(W - 1, pbwtRev);

  PDUMPR(W - 1, pbwtPr);
  PDUMPR(W - 1, pbwtPrRev);

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
    divc(nrow, w64, pbwt, pbwtPr, pbwtRev, pbwtPrRev, W);
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
  size_t *aux = malloc(nrow * sizeof *aux);

  pbwtad *p0 = pbwtad_new(nrow);
  pbwtad *p1 = pbwtad_new(nrow);
  pbwtad *p0rev = pbwtad_new(nrow);
  pbwtad *p1rev = pbwtad_new(nrow);
  sbfgetcolw64rn_mmap(fin, nrow, pw, ncol);
  rrsort0(nrow, pw, p0->a, aux);

  // memcpy(p1rev->a, p0rev->a, nrow * sizeof *(p0rev->a));
  // memcpy(p1rev->d, p0rev->d, nrow * sizeof *(p0rev->d));
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
  pb0rev[0] = p0;
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
  uint64_t *pw0 = malloc(nrow * sizeof *pw0);
  uint64_t *pw1 = malloc(nrow * sizeof *pw1);
  size_t *aux = malloc(nrow * sizeof *aux);
  fgetcoliw64r(fin, 0, nrow, pw0, ncol);

  // pb0 is now filled with computed values,
  // to allow reusing I need to fill pb1 with empty values
  for (int j = 0; j < W; j++) {
    pb1[j] = pbwtad_new(nrow);
    pb1rev[j] = pbwtad_new(nrow);
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
      // size_t J = W * (j + 1) - 1;
      size_t J = W - 1;

      wr64mrgsi(nrow, pw1, pw0, w, x);
      // pbwtad *ps = pbwtad_new(nrow);
      pbwtad *ps = pb1[J - x];
      pbwtad *psrev = pb1rev[J - x];
      // pb0[J - x] = ps;
      memcpy(ps->a, pb0[J - x]->a, nrow * sizeof *(ps->a));
      memcpy(ps->d, pb0[J - x]->d, nrow * sizeof *(ps->d));
      memcpy(psrev->a, pb0[J - x]->a, nrow * sizeof *(psrev->a));
      memcpy(psrev->d, pb0rev[J - x]->a, nrow * sizeof *(psrev->d));
      rrsortx_noaux(nrow, w, ps->a);
      reversec(ps, psrev, nrow);
      divc(nrow, w, ps, pb0[J - x], psrev, pb0rev[J - x], W);
      FREE(w);
    }
    PDUMP_SEQ_OFFSETR(0, W, pb1,
                      W * j); // FIXME: print here should be rewritten
                              // for LCP display as divergence
    SWAP(pw0, pw1);
    SWAP(pb0, pb1);
    SWAP(pb0rev, pb1rev);

    j++;
  }
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
    // if ((j/W) % W == 3) swapdiv(pp1, nrow, j);
    // if ((j/W) % W == 2) swapdiv(pp0, nrow, j);
    cpbwtiLCP(nrow, j, c0, pp0, pp1);
    PDUMPR(j, pp1);
    SWAP(pp0, pp1);
  }

#endif
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
pbwtad **bwparc_rrs(void *fin, size_t nrow, size_t ncol) { // BPR
  // NOTE: right now I don't know what I need, so I'm keeping
  // everything in memory, we'll see later
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

  fgetcoli(fin, 0, nrow, c0, ncol);
  pb0[0] = cpbwt(nrow, c0, p0);
  pb0rev[0] = p0;
  PDUMP(0, pb0[0]);
  // PBWTAD_FREE(p0);

  for (int j = 1; j < W; j++) {
    fgetcoli(fin, j, nrow, c0, ncol);
    pb0[j] = cpbwt(nrow, c0, pb0[j - 1]);
    pb0rev[j] = pbwtad_new(nrow);
    reversecprev(pb0[j], pb0[j - 1], pb0rev[j], nrow);
    PDUMP(j, pb0[j]);
  }

  /* swapping div values with LCP for window computation */
  for (int j = 0; j < W; j++) {
    swapdiv(pb0[j], nrow, j);
    swapdiv(pb0rev[j], nrow, j);
    // PDUMP(j, pb0[j]);
  }

  uint64_t *pw0 = malloc(nrow * sizeof *pw0);
  uint64_t *pw1 = malloc(nrow * sizeof *pw1);
  size_t *aux = malloc(nrow * sizeof *aux);
  bfgetcolw64rn(fin, nrow, pw0, ncol);

  // pb0 is now filled with computed values,
  // to allow reusing I need to fill pb1 with empty values
  for (int j = 0; j < W; j++) {
    pb1[j] = pbwtad_new(nrow);
    pb1rev[j] = pbwtad_new(nrow);
  }

  size_t j;

  for (j = 1; j * W <= ncol - W; j++) {
    bfgetcolw64rn(fin, nrow, pw1, ncol);
    pbwtad *ps = pb1[W - 1];
    pbwtad *psrev = pb1rev[W - 1];
    memcpy(ps->a, pb0[W - 1]->a, nrow * sizeof *(ps->a));
    memcpy(ps->d, pb0[W - 1]->d, nrow * sizeof *(ps->d));
    memcpy(psrev->a, pb0rev[W - 1]->a, nrow * sizeof *(ps->a));
    memcpy(psrev->d, pb0rev[W - 1]->d, nrow * sizeof *(ps->d));
    rrsortx(nrow, pw1, ps->a, aux);
    reversec(ps, psrev, nrow);
    divc(nrow, pw1, ps, pb0[W - 1], psrev, pb0rev[W - 1], W);

    // pbwtad *ps = pbwtad_new(nrow);
    // pb[W * (j + 1) - 1] = ps;
    // memcpy(ps->a, pb[W * j - 1]->a, nrow * sizeof *(ps->a));
    // rrsortx(nrow, pw1, ps->a, aux);

#pragma omp parallel for shared(pw1, pw0, pb0, pb1, j)
    for (size_t x = 1; x < W; x++) {
      uint64_t *w = malloc(nrow * sizeof *w);
      size_t J = W - 1;

      wr64mrgsi(nrow, pw1, pw0, w, x);
      pbwtad *ps = pb1[J - x];
      pbwtad *psrev = pb1rev[J - x];
      // pbwtad *ps = pbwtad_new(nrow);
      // pb[J - x] = ps;
      memcpy(ps->a, pb0[J - x]->a, nrow * sizeof *(ps->a));
      memcpy(ps->d, pb0[J - x]->d, nrow * sizeof *(ps->d));
      memcpy(psrev->a, pb0[J - x]->a, nrow * sizeof *(psrev->a));
      memcpy(psrev->d, pb0rev[J - x]->a, nrow * sizeof *(psrev->d));
      // memcpy(ps->a, pb[J - W - x]->a, nrow * sizeof *(ps->a));
      rrsortx_noaux(nrow, w, ps->a);
      reversec(ps, psrev, nrow);
      divc(nrow, w, ps, pb0[J - x], psrev, pb0rev[J - x], W);
      FREE(w);
    }
    PDUMP_SEQ_OFFSETR(0, W, pb1, W * j);
    // PDUMP_SEQ(W * j - 1, W * (j + 1) - 1, pb1);
    SWAP(pw0, pw1);
    SWAP(pb0, pb1);
    SWAP(pb0rev, pb1rev);
    // j++;
  }

  pbwtad *pp0, *pp1;
  pp0 = pb0[W - 1];
  pp1 = pb0[W - 2];

  for (j = j * W; j < ncol; j++) {
    fgetcoli(fin, j, nrow, c0, ncol);
    cpbwtiLCP(nrow, j, c0, pp0, pp1);
    PDUMPR(j, pp1);
    SWAP(pp0, pp1);
    // pb[j] = cpbwt_0(nrow, c0, pb[j - 1]);
    // PDUMP(j, pb[j]);
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

pbwtad **wstagparc_rrs(char *fpath, size_t nrow, size_t ncol) { // SPR
#if defined(BF2IOMODE_BCF)
  ncol = 0;
  htsFile *fp = hts_open(fpath, "r");
  if (fp) {
    // Read the header to advance the file pointer to the data
    bcf_hdr_t *h = bcf_hdr_read(fp);
    if (h) {
      bcf1_t *line = bcf_init();

      // bcf_read is faster than bcf_sr_next_line as it skips
      // synchronization logic and deep unpacking
      while (bcf_read(fp, h, line) == 0) {
        ncol++;
      }

      bcf_destroy(line);
      bcf_hdr_destroy(h);
    }
    hts_close(fp);
  }
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

#ifdef BF2IOMODE_BCF
  bcf_srs_t *_sr = bcf_sr_init();
  bcf_sr_add_reader(_sr, fpath);
  void *fin = _sr;
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
  memcpy(pbprev->a, pb[0]->a, nrow * sizeof(*pb[0]->a));
  memcpy(pbprev->d, pb[0]->d, nrow * sizeof(*pb[0]->d));

  cpbwtiLCP(nrow, 1, c0, pbprev, pb[1]);

  PDUMPR(0, pb[1]);

  for (size_t j = 1; j < W;) {
    fgetcoli(fin, j, nrow, c0, ncol);
    memcpy(pbprev->a, pb[j]->a, nrow * sizeof(*pb[j]->a));
    memcpy(pbprev->d, pb[j]->d, nrow * sizeof(*pb[j]->d));

    cpbwtiLCP(nrow, j, c0, pbprev, pb[j + 1]);
    PDUMPR(j, pb[j + 1]);
    j++;
  }

#pragma omp parallel
  {

#ifdef BF2IOMODE_BCF
    bcf_srs_t *_sr = bcf_sr_init();
    bcf_sr_add_reader(_sr, fpath);
    void *fin = _sr;
    size_t lastrowread = 0;
#elif defined(BF2IOMODE_BM) || defined(BF2IOMODE_ENC)
    FILE *fin = fopen(fpath, "r");

    if (!fin) {
      perror("[spr]");
      exit(32);
    }

#else
#error UNDEFINED BEHAVIOUR
#endif

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
    size_t *aux = malloc(nrow * sizeof *aux);
    uint64_t *pw = malloc(nrow * sizeof *pw);

    for (size_t offset = 0; offset < count; offset++) {
      size_t lane = start + offset;
      memcpy(pt0->a, pb[lane]->a, nrow * sizeof *(pb[lane]->a));
      memcpy(pt0->d, pb[lane]->d, nrow * sizeof *(pb[lane]->d));

      reversec(pt0, pt0rev, nrow);

      for (size_t j = lane; j + W <= ncol; j += W) {

#ifdef BF2IOMODE_BCF
        lastrowread = fgetcolwgri(fin, j + 1, nrow, pw, lastrowread, W);
#else
        fgetcolwgri(fin, j, nrow, pw, ncol, W);
#endif

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
  // for (int j = 0; j < W + 1; j++) {
  //   PBWTAD_FREE(pb[j]);
  //   FREE(pb[j]);
  // }
  // FREE(pb);
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
    TRACE(sblinc(fd, nrow, ncol), r); // BLIS
  } else if (strcmp(argv[1], "blim") == 0) {
    // r = mblinc(fd, nrow, ncol);
    TRACE(mblinc(fd, nrow, ncol), r); // BLIM
  } else if (strcmp(argv[1], "arst") == 0) {
    // r = wapproxc_rrs(fin, nrow, ncol);
    TRACE(wapproxcst_rrs(fin, nrow, ncol), r);
  } else if (strcmp(argv[1], "ars") == 0) {
    // r = wapproxc_rrs(fin, nrow, ncol);
    TRACE(wapproxc_rrs(fin, nrow, ncol), r);
  } else if (strcmp(argv[1], "bar") == 0) {
    // r = wbapproxc_rrs(fin, nrow, ncol);
    TRACE(wbapproxc_rrs(fin, nrow, ncol), r);
  } else if (strcmp(argv[1], "bars") == 0) {
    // r = swbapproxc_rrs(fd, nrow, ncol);
    TRACE(swbapproxc_rrs(fd, nrow, ncol), r);
  } else if (strcmp(argv[1], "barm") == 0) { // fixed
    // r = mwbapproxc_rrs(fd, nrow, ncol);
    TRACE(mwbapproxc_rrs(fd, nrow, ncol), r);
  } else if (strcmp(argv[1], "prs") == 0) {
    // r = wparc_rrs(fin, nrow, ncol);
    TRACE(wparc_rrs(fin, nrow, ncol), r);
  } else if (strcmp(argv[1], "bpr") == 0) { // fixed
    // r = bwparc_rrs(fin, nrow, ncol);
    TRACE(bwparc_rrs(fin, nrow, ncol), r);
  } else if (strcmp(argv[1], "spr") == 0) { // broken
    // r = wstagparc_rrs(argv[2], nrow, ncol);
    TRACE(wstagparc_rrs(argv[2], nrow, ncol), r);
  } else if (strcmp(argv[1], "srs") == 0) { // hidle
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
