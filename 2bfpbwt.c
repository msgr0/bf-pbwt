#include "lib/quadsort/quadsort.h"
#include <assert.h>
#include <fcntl.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

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

// get column i from file
// c[n] is a pointer to store the column,
// nc is total number of columns, needed for fseek
static inline void bfgetcoln(FILE *fd, size_t n, uint8_t *restrict c,
                             size_t nc) {
  // NOTE: this assumes ASCII text file, offset are computed assuming
  // 1-byte size for each character

#define BFGETCOLI_BUF_SIZE 128
  static char rbuf[BFGETCOLI_BUF_SIZE];
  static size_t i = 0;
  static size_t bufn = BFGETCOLI_BUF_SIZE;
  static uint8_t *buf = NULL;
  if (!buf)
    buf = malloc(BFGETCOLI_BUF_SIZE * n * sizeof *buf);

#if 0
  if (bufn == BFGETCOLI_BUF_SIZE) {
    int x;
    fseek(fd, i, SEEK_SET);
    for (size_t r = 0; r < n; r++) {

      fread(rbuf, 1, BFGETCOLI_BUF_SIZE, fd);
      for (size_t s = 0; s < BFGETCOLI_BUF_SIZE; s++) {
        x = rbuf[s] - 48;
        buf[(n * s) + r] = x;
      }
      c[r] = buf[r];
      fseek(fd, nc - BFGETCOLI_BUF_SIZE + 1, SEEK_CUR);
    }
    bufn = 1;
  } else {
    for (size_t r = 0; r < n; r++) {
      c[r] = buf[(n * bufn) + r];
    }
    bufn++;
  }
#else
  if (bufn == BFGETCOLI_BUF_SIZE) {
    int x;
    fseek(fd, i, SEEK_SET);
    for (size_t r = 0; r < n; r++) {
      fread(&buf[r * BFGETCOLI_BUF_SIZE], 1, BFGETCOLI_BUF_SIZE, fd);
      c[r] = buf[r * BFGETCOLI_BUF_SIZE] - 48;
      fseek(fd, nc - BFGETCOLI_BUF_SIZE + 1, SEEK_CUR);
    }
    bufn = 1;
  } else {
    for (size_t r = 0; r < n; r++) {
      c[r] = buf[r * BFGETCOLI_BUF_SIZE + bufn] - 48;
    }
    bufn++;
  }
#endif
  i++;
}

static inline void sbfgetcoln(int fd, size_t n, uint8_t *restrict c,
                              size_t nc) {
  // NOTE: this assumes ASCII text file, offset are computed assuming
  // 1-byte size for each character

#define BFGETCOLI_BUF_SIZE 128
  static char rbuf[BFGETCOLI_BUF_SIZE];
  static size_t i = 0;
  static size_t bufn = BFGETCOLI_BUF_SIZE;
  static uint8_t *buf = NULL;
  if (!buf)
    buf = malloc(BFGETCOLI_BUF_SIZE * n * sizeof *buf);

  size_t offset = i;
  if (bufn == BFGETCOLI_BUF_SIZE) {
    for (size_t r = 0; r < n; r++) {
      pread(fd, &buf[r * BFGETCOLI_BUF_SIZE], BFGETCOLI_BUF_SIZE, offset);
      c[r] = buf[r * BFGETCOLI_BUF_SIZE] - 48;
      offset += nc + 1;
    }
    bufn = 1;
  } else {
    for (size_t r = 0; r < n; r++) {
      c[r] = buf[r * BFGETCOLI_BUF_SIZE + bufn] - 48;
    }
    bufn++;
  }
  i++;
}

#define FGETCOLIW_DECLARE(W)                                                   \
  static inline void fgetcoliw##W(FILE *fd, size_t i, size_t n,                \
                                  uint64_t *restrict c, size_t nc) {           \
    int x;                                                                     \
    fseek(fd, i *W, SEEK_SET);                                                 \
    for (size_t r = 0; r < n; r++) {                                           \
      c[r] = 0;                                                                \
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
  }                                                                            \
  static inline void bfgetcolw##W##rn(FILE *fd, size_t n,                      \
                                      uint64_t *restrict c, size_t nc) {       \
    static size_t i = 0;                                                       \
    static size_t bufn = BFGETCOLWR_BUF_SIZE;                                  \
    static uint64_t *buf = NULL;                                               \
    if (!buf)                                                                  \
      buf = malloc(BFGETCOLWR_BUF_SIZE * n * sizeof *buf);                     \
    if (bufn == BFGETCOLWR_BUF_SIZE) {                                         \
      uint64_t x;                                                              \
      fseek(fd, i *W, SEEK_SET);                                               \
      for (size_t r = 0; r < n; r++) {                                         \
        for (size_t s = 0; s < BFGETCOLWR_BUF_SIZE; s++) {                     \
          buf[(n * s) + r] = 0;                                                \
          for (size_t j = 0; j < W; j++) {                                     \
            x = fgetc(fd) - 48;                                                \
            buf[(n * s) + r] = (x << j) | buf[(n * s) + r];                    \
          }                                                                    \
        }                                                                      \
        c[r] = buf[r];                                                         \
        fseek(fd, nc - (BFGETCOLWR_BUF_SIZE * W) + 1, SEEK_CUR);               \
      }                                                                        \
      bufn = 1;                                                                \
    } else {                                                                   \
      for (size_t r = 0; r < n; r++) {                                         \
        c[r] = buf[(n * bufn) + r];                                            \
      }                                                                        \
      bufn++;                                                                  \
    }                                                                          \
    i++;                                                                       \
  }

#define BFGETCOLWR_BUF_SIZE 2
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
// get column i from file
// c[n] is a pointer to store the column,
// nc is total number of columns, needed for fseek
static void bfgetcolwgrn(FILE *fd, size_t n, uint64_t *restrict c, size_t nc,
                         uint8_t w) {
  // NOTE: this assumes ASCII text file, offset are computed assuming
  // 1-byte size for each character

  static size_t i = 0;
  static size_t bufn = BFGETCOLWR_BUF_SIZE;
  static uint64_t *buf = NULL;
  if (!buf)
    buf = malloc(BFGETCOLWR_BUF_SIZE * n * sizeof *buf);

  if (bufn == BFGETCOLWR_BUF_SIZE) {
    uint64_t x;
    fseek(fd, i * w, SEEK_SET);
    for (size_t r = 0; r < n; r++) {
      for (size_t s = 0; s < BFGETCOLWR_BUF_SIZE; s++) {
        buf[(n * s) + r] = 0;
        for (size_t j = 0; j < w; j++) {
          x = fgetc(fd) - 48;
          buf[(n * s) + r] = (x << j) | buf[(n * s) + r];
        }
      }
      c[r] = buf[r];
      fseek(fd, nc - (BFGETCOLWR_BUF_SIZE * w) + 1, SEEK_CUR);
    }
    bufn = 1;
  } else {
    for (size_t r = 0; r < n; r++) {
      c[r] = buf[(n * bufn) + r];
    }
    bufn++;
  }
  i++;
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

static pbwtad *cpbwt(size_t n, uint8_t *restrict c, pbwtad *restrict p,
                     size_t *_o, size_t *_z) {
  static size_t *o = NULL;
  if (!o)
    o = malloc(n * sizeof *o);

  pbwtad *ret = malloc(sizeof *ret);
  ret->a = malloc(n * sizeof *(ret->a));

  size_t r = 0, q = 0;
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
    size_t mask = c[idx]; // 1 if true, 0 if false
    o[q] = idx;
    ret->a[r] = idx;
    q += mask;     // Increment q if mask is 1
    r += mask ^ 1; // Increment r if mask is 0
  }
#endif

#if 0
  for (i = 0; i < r; i++) {
    ret->a[i] = z[i];
  }
  for (i = 0; i < q; i++) {
    ret->a[r + i] = o[i];
  }
#else
  memcpy(ret->a + r, o, q * sizeof(size_t));
#endif

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
    /*fprintf(stderr, "\r%10zu/%zu", j + 1, ncol);*/
    fgetcoli(fin, j, nrow, c0, ncol);
    pl[j] = cpbwt(nrow, c0, pl[j - 1], z, o);
  }
  fputc(0xA, stderr);

#if 0
  for (size_t j = 0; j < ncol; j++) {
    /*DPRINT("lin %3zu: ", j);*/
    DPRINT("%3zu: ", j);
    DPARR(nrow, pl[j]->a, "%zu ");
  }
#endif
  FREE(c0);
  FREE(o);
  FREE(z);
  return pl;
}
pbwtad **blinc(FILE *fin, size_t nrow, size_t ncol) {
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
  bfgetcoln(fin, nrow, c0, ncol);
  /*parr(nrow, c0, "%d ");*/
  pl[0] = cpbwt(nrow, c0, p0, z, o);
  FREE(p0->a);
  FREE(p0);
  /*for (size_t i = 0; i < nrow; i++) {*/
  /*  if (c0[i] != 0 && c0[i] != 1) {*/
  /*    printf("c[%zu]=%hhu\n", i, c0[i]);*/
  /*    exit(122);*/
  /*  }*/
  /*}*/

  for (size_t j = 1; j < ncol; j++) {
    /*fprintf(stderr, "\r%10zu/%zu", j + 1, ncol);*/
    /*fprintf(stderr, "\r%10zu/%zu\n", j + 1, ncol);*/
    bfgetcoln(fin, nrow, c0, ncol);
    /*for (size_t i = 0; i < nrow; i++) {*/
    /*  if (c0[i] != 0 && c0[i] != 1) {*/
    /*    printf("c0[%zu]=%hhu\n", i, c0[i]);*/
    /*    exit(122);*/
    /*  }*/
    /*}*/
    pl[j] = cpbwt(nrow, c0, pl[j - 1], z, o);
  }
  fputc(0xA, stderr);

#if 0
  for (size_t j = 0; j < ncol; j++) {
    /*DPRINT("bli %3zu: ", j);*/
    DPRINT("%3zu: ", j);
    DPARR(nrow, pl[j]->a, "%zu ");
  }
#endif
  FREE(c0);
  FREE(o);
  FREE(z);
  return pl;
}
pbwtad **sblinc(int fin, size_t nrow, size_t ncol) {
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
  sbfgetcoln(fin, nrow, c0, ncol);
  pl[0] = cpbwt(nrow, c0, p0, z, o);
  FREE(p0->a);
  FREE(p0);

  for (size_t j = 1; j < ncol; j++) {
    sbfgetcoln(fin, nrow, c0, ncol);
    pl[j] = cpbwt(nrow, c0, pl[j - 1], z, o);
  }
  fputc(0xA, stderr);

#if 0
  for (size_t j = 0; j < ncol; j++) {
    /*DPRINT("bli %3zu: ", j);*/
    DPRINT("%3zu: ", j);
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
    /*fprintf(stderr, "\r%10zu/%zu", (j * W) + 1, ncol);*/
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
  fputc(0xA, stderr);

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
pbwtad **wbapproxc_rrs(FILE *fin, size_t nrow, size_t ncol,
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
  bfgetcolw64rn(fin, nrow, pw, ncol);
  rrsort0(nrow, pw, ps->a, aux);

  size_t j;
  for (j = 1; j * W <= ncol - W; j++) {
    pbwtad *ps = malloc(nrow * sizeof *ps);
    ps->a = malloc(nrow * sizeof *(ps->a));
    pb[W * (j + 1) - 1] = ps;
    memcpy(ps->a, pb[W * j - 1]->a, nrow * sizeof *(ps->a));
    bfgetcolw64rn(fin, nrow, pw, ncol);
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
  fputc(0xA, stderr);

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

  uint64_t *pw0 = malloc(nrow * sizeof *pw0);
  uint64_t *pw1 = malloc(nrow * sizeof *pw1);
  size_t *aux = malloc(nrow * sizeof *aux);
  fgetcoliw64r(fin, 0, nrow, pw0, ncol);

  size_t j;

  for (j = 1; j * W <= ncol - W; j++) {
    /*fprintf(stderr, "\r%10zu/%zu", (j * W) + 1, ncol);*/
    pbwtad *ps = malloc(nrow * sizeof *ps);
    ps->a = malloc(nrow * sizeof *(ps->a));
    pb[W * (j + 1) - 1] = ps;
    memcpy(ps->a, pb[W * j - 1]->a, nrow * sizeof *(ps->a));
    fgetcoliw64r(fin, j, nrow, pw1, ncol);
    rrsortx(nrow, pw1, ps->a, aux);

#pragma omp parallel for shared(pw1, pw0, pb, j)
    for (size_t x = 1; x < W; x++) {
      uint64_t *w = malloc(nrow * sizeof *w);
      size_t J = W * (j + 1) - 1;

      wr64mrgsi(nrow, pw1, pw0, w, x);
      pbwtad *ps = malloc(nrow * sizeof *ps);
      ps->a = malloc(nrow * sizeof *(ps->a));
      pb[J - x] = ps;
      memcpy(ps->a, pb[J - W - x]->a, nrow * sizeof *(ps->a));
      rrsortx_noaux(nrow, w, ps->a);

      FREE(w);
    }

    SWAP(pw0, pw1);
  }

  c0 = malloc(nrow * sizeof *c0);
  for (j = j * W; j < ncol; j++) {
    fgetcoli(fin, j, nrow, c0, ncol);
    pb[j] = cpbwt_0(nrow, c0, pb[j - 1]);
  }
  fputc(0xA, stderr);

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

pbwtad **wseq_rrs(FILE *fin, size_t nrow, size_t ncol) {
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
  fputc(0xA, stderr);

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
  char _usage_args_[] = "[lin|bli|ars|aqs|bar|prs] FILE\n";
  if (argc < 2) {
    fprintf(stderr, "Usage: %s %s FILE\n", argv[0], _usage_args_);
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
  pbwtad **r;

  if (strcmp(argv[1], "lin") == 0) {
    r = linc(fin, nrow, ncol);
  } else if (strcmp(argv[1], "bli") == 0) {
    r = blinc(fin, nrow, ncol);
  } else if (strcmp(argv[1], "blis") == 0) {
    /*fclose(fin);*/
    /*fin = NULL;*/
    int fd = open(argv[2], O_RDONLY);
    r = sblinc(fd, nrow, ncol);
  } else if (strcmp(argv[1], "ars") == 0) {
    r = wapproxc_rrs(fin, nrow, ncol, APPROX_MODE_LAST_WINDOW);
  } else if (strcmp(argv[1], "aqs") == 0) {
    fprintf(stderr, "\e[0;33mWARNING: this is not been properly tested and "
                    "might not work.\e[0m\n");
    r = wapproxc_qs(fin, nrow, ncol, APPROX_MODE_LAST_WINDOW);
  } else if (strcmp(argv[1], "bar") == 0) {
    r = wbapproxc_rrs(fin, nrow, ncol, APPROX_MODE_LAST_WINDOW);
  } else if (strcmp(argv[1], "prs") == 0) {
    r = wparc_rrs(fin, nrow, ncol);
  } else if (strcmp(argv[1], "srs") == 0) {
    r = wseq_rrs(fin, nrow, ncol);
  } else {
    fprintf(stderr, "Usage: %s %s FILE\n", argv[0], _usage_args_);
    return EXIT_FAILURE;
  }

  /*fclose(fin);*/

  if (r)
    for (size_t i = 0; i < ncol; i++) {
      if (r[i]) {
        FREE(r[i]->a);
        FREE(r[i]);
      }
    }

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

  uint8_t *c0 = malloc(nrow * sizeof *c0);
  sbfgetcoln(ifin, nrow, c0, ncol);
  parr(nrow, c0, "%d ");
  fgetcoli(fin, 0, nrow, c0, ncol);
  parr(nrow, c0, "%d ");
  sbfgetcoln(ifin, nrow, c0, ncol);
  parr(nrow, c0, "%d ");
  fgetcoli(fin, 1, nrow, c0, ncol);
  parr(nrow, c0, "%d ");
  sbfgetcoln(ifin, nrow, c0, ncol);
  parr(nrow, c0, "%d ");
  fgetcoli(fin, 2, nrow, c0, ncol);
  parr(nrow, c0, "%d ");
  sbfgetcoln(ifin, nrow, c0, ncol);
  parr(nrow, c0, "%d ");
  fgetcoli(fin, 3, nrow, c0, ncol);
  parr(nrow, c0, "%d ");
  sbfgetcoln(ifin, nrow, c0, ncol);
  parr(nrow, c0, "%d ");
  fgetcoli(fin, 4, nrow, c0, ncol);
  parr(nrow, c0, "%d ");
  sbfgetcoln(ifin, nrow, c0, ncol);
  parr(nrow, c0, "%d ");
  fgetcoli(fin, 5, nrow, c0, ncol);
  parr(nrow, c0, "%d ");

  /*uint64_t *w0 = malloc(nrow * sizeof *w0);*/
  /*bfgetcolw64rn(fin, nrow, w0, ncol);*/
  /*parr(nrow, w0, "%llu ");*/
  /*fgetcoliw64r(fin, 0, nrow, w0, ncol);*/
  /*parr(nrow, w0, "%llu ");*/
  /*bfgetcolw64rn(fin, nrow, w0, ncol);*/
  /*parr(nrow, w0, "%llu ");*/
  /*fgetcoliw64r(fin, 1, nrow, w0, ncol);*/
  /*parr(nrow, w0, "%llu ");*/
  /*bfgetcolw64rn(fin, nrow, w0, ncol);*/
  /*parr(nrow, w0, "%llu ");*/
  /*fgetcoliw64r(fin, 2, nrow, w0, ncol);*/
  /*parr(nrow, w0, "%llu ");*/
  return EXIT_SUCCESS;
}
