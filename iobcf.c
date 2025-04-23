// vim:ft=c
#include "htslib/synced_bcf_reader.h"
#include "htslib/vcf.h"
#include "io.h"
#include <stdio.h>
#include <string.h>

#ifndef IOBCF_UNUSED_EXITCODE
#define IOBCF_UNUSED_EXITCODE 3
#endif

#define IOBCF_ASSUME_WELLFORMED
#ifndef IOBCF_ASSUME_WELLFORMED
#define IOBCF_WELLFORMED_CHECK(ptr)                                            \
  do {                                                                         \
    if ((ptr) == bcf_int32_vector_end)                                         \
      exit(-2);                                                                \
    if (bcf_gt_is_missing((ptr)))                                              \
      exit(-1);                                                                \
  } while (0)
#else
#define IOBCF_WELLFORMED_CHECK(ptr)
#endif

void fgetrc(void *fd, size_t *nr, size_t *nc) {
  bcf_srs_t *sr = bcf_sr_init();
  bcf_sr_add_reader(sr, ((bcf_srs_t *)fd)->readers[0].fname);
  bcf_hdr_t *hdr = sr->readers[0].header;

  *nr = bcf_hdr_nsamples(hdr) * 2;
  *nc = 0;

  // NOTE: maybe some parts might be rewritten to avoid doing this,
  // however there is not much overhead. For chr10 it takes ~5 secs
  while (bcf_sr_next_line(sr)) {
    (*nc)++;
  }

  bcf_sr_destroy(sr);
}

void fgetcoli(void *fd, size_t i, size_t n, uint8_t *restrict c, size_t nc) {
  // NOTE: nc is used as a flag. If nc == 0, ignore _li and assume sequential
  // read
  bcf_srs_t *sr = fd;
  static bcf_hdr_t *hdr = NULL;
  if (!hdr)
    hdr = sr->readers[0].header;

  static ssize_t _li = -1;
  // NOTE: I am doing a bit of trickery here assuming that
  // 1. htslib reads lines with incrementing iterator;
  // 2. tipically cols (BCF-row) are read sequentially;
  // I am keeping track of the last position requested and if
  // i == _li+1, then I can just read the next position,
  // otherwise real seeking is necessary
  if (i == _li + 1 || !nc) {
    // this check should not be necessary
    if (!bcf_sr_next_line(sr))
      exit(124);
    bcf1_t *line = bcf_sr_get_line(sr, 0);
    int32_t *gt_arr = NULL, ngt_arr = 0;
    int ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);

    for (size_t i = 0; i < n / 2; i++) {
      int32_t *ptr = gt_arr + i * 2;
      // hap 1
      IOBCF_WELLFORMED_CHECK(ptr[0]);
      c[2 * i] = bcf_gt_allele(ptr[0]);

      // hap 2
      IOBCF_WELLFORMED_CHECK(ptr[1]);
      c[2 * i + 1] = bcf_gt_allele(ptr[1]);
    }
    free(gt_arr);
    // bcf_destroy(line);
  } else {
    // WARN: this does not work, at the moment.
    // If we do not fix this, it will be a bloodbath
    // bcf_sr_seek(sr, NULL, 0);
    errno = EPERM;
    perror("NOT IMPLMENTED YET");
    exit(24);
  }
  _li = i;
}

void bfgetcoln(void *fd, size_t n, uint8_t *restrict c, size_t nc) {
  // NOTE: nc is not used here.
  bcf_srs_t *sr = fd;
  bcf_hdr_t *hdr = sr->readers[0].header;

  static size_t bufn = BFGETCOLI_BUF_SIZE;
  static uint8_t *buf = NULL;
  if (!buf)
    buf = malloc(BFGETCOLI_BUF_SIZE * n * sizeof *buf);

  if (bufn == BFGETCOLI_BUF_SIZE) {
    int x;
    for (size_t r = 0; r < BFGETCOLI_BUF_SIZE; r++) {
      if (!bcf_sr_next_line(sr))
        break;
      bcf1_t *line = bcf_sr_get_line(sr, 0);
      int32_t *gt_arr = NULL, ngt_arr = 0;
      int ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);
      for (size_t i = 0; i < n / 2; i++) {
        int32_t *ptr = gt_arr + i * 2;
        // hap 1
        IOBCF_WELLFORMED_CHECK(ptr[0]);
        buf[r * n + 2 * i] = bcf_gt_allele(ptr[0]);

        // hap 2
        IOBCF_WELLFORMED_CHECK(ptr[1]);
        buf[r * n + 2 * i + 1] = bcf_gt_allele(ptr[1]);
      }
    }
    bufn = 0;
  }
  memcpy(c, buf + (bufn * n), n);
  bufn++;
}

void sbfgetcoln(int fd, size_t n, uint8_t *restrict c, size_t nc) {
  fputs("\e[0;33mMode not used for this type of file. Exiting.\e[0m\n", stderr);
  exit(IOBCF_UNUSED_EXITCODE);
}

void mbfgetcoln(int fd, size_t n, uint8_t *restrict c, size_t nc) {
  fputs("\e[0;33mMode not used for this type of file. Exiting.\e[0m\n", stderr);
  exit(IOBCF_UNUSED_EXITCODE);
}

#define FGETCOLIW_IMPL(W)                                                      \
  void fgetcoliw##W(void *fd, size_t i, size_t n, uint64_t *restrict c,        \
                    size_t nc) {                                               \
    fprintf(stderr, "\e[0;33m[%s] Not Implemented Yet.\e[0m\n", __func__);     \
  }                                                                            \
  void w##W##mrgsi(size_t n, uint64_t const *wc, uint64_t const *wp,           \
                   uint64_t *restrict c, size_t i) {                           \
    uint64_t c1;                                                               \
    for (size_t r = 0; r < n; r++) {                                           \
      c1 = wp[r] & ((1 << i) - 1);                                             \
      c[r] = (c1 << (W - i)) | (wc[r] >> i);                                   \
    }                                                                          \
  }                                                                            \
  void fgetcoliw##W##r(void *fd, size_t i, size_t n, uint64_t *restrict c,     \
                       size_t nc) {                                            \
    bcf_srs_t *sr = fd;                                                        \
    bcf_hdr_t *hdr = sr->readers[0].header;                                    \
                                                                               \
    static ssize_t _li = -1;                                                   \
    if (i == _li + 1) {                                                        \
      memset(c, 0, n * sizeof *c);                                             \
      for (size_t wix = 0; wix < W; wix++) {                                   \
        if (!bcf_sr_next_line(sr))                                             \
          return;                                                              \
        bcf1_t *line = bcf_sr_get_line(sr, 0);                                 \
        int32_t *gt_arr = NULL, ngt_arr = 0;                                   \
        int ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);             \
                                                                               \
        for (size_t i = 0; i < n / 2; i++) {                                   \
          int32_t *ptr = gt_arr + i * 2;                                       \
          IOBCF_WELLFORMED_CHECK(ptr[0]);                                      \
          c[2 * i] |= ((uint64_t)bcf_gt_allele(ptr[0]) << wix);                \
                                                                               \
          IOBCF_WELLFORMED_CHECK(ptr[1]);                                      \
          c[2 * i + 1] |= ((uint64_t)bcf_gt_allele(ptr[1]) << wix);            \
        }                                                                      \
      }                                                                        \
    } else {                                                                   \
      errno = EPERM;                                                           \
      perror("NOT IMPLMENTED YET");                                            \
      exit(26);                                                                \
    }                                                                          \
    _li = i;                                                                   \
  }                                                                            \
  void wr##W##mrgsi(size_t n, uint64_t const *wc, uint64_t const *wp,          \
                    uint64_t *restrict c, size_t i) {                          \
    uint64_t c1;                                                               \
    static const uint64_t mask = (UINT64_MAX >> (64 - W));                     \
    for (size_t r = 0; r < n; r++) {                                           \
      c[r] = (wp[r] >> (W - i)) | ((wc[r] << i) & mask);                       \
    }                                                                          \
  }                                                                            \
  void bfgetcolw##W##rn(void *fd, size_t n, uint64_t *restrict c, size_t nc) { \
    fprintf(stderr, "\e[0;33m[%s] Not Implemented Yet.\e[0m\n", __func__);     \
  }                                                                            \
  void sbfgetcolw##W##rn(int fd, size_t n, uint64_t *restrict c, size_t nc) {  \
    fputs("\e[0;33mMode not used for this type of file. Exiting.\e[0m\n",      \
          stderr);                                                             \
    exit(IOBCF_UNUSED_EXITCODE);                                               \
  }                                                                            \
  void sbfgetcolw##W##rn_mmap(int fd, size_t n, uint64_t *restrict c,          \
                              size_t nc) {                                     \
    fputs("\e[0;33mMode not used for this type of file. Exiting.\e[0m\n",      \
          stderr);                                                             \
    exit(IOBCF_UNUSED_EXITCODE);                                               \
  }

FGETCOLIW_IMPL(8)
FGETCOLIW_IMPL(16)
FGETCOLIW_IMPL(32)
FGETCOLIW_IMPL(64)

void fgetcoliwg(void *fd, size_t i, size_t n, uint64_t *restrict c, size_t nc,
                uint8_t w) {}

void fgetcoliwgr(void *fd, size_t i, size_t n, uint64_t *restrict c, size_t nc,
                 uint8_t w) {
  // NOTE: nc is not used here.
  bcf_srs_t *sr = fd;
  bcf_hdr_t *hdr = sr->readers[0].header;

  static ssize_t _li = -1;
  // NOTE: Same trickery here as in `fgetcoli`, however
  // _li is not the last index or row, but the last index of window.
  // Current BCF row is (_li * w)
  if (i == _li + 1) {
    memset(c, 0, n * sizeof *c);
    // this check should not be necessary
    for (size_t wix = 0; wix < w; wix++) {
      // printf("wix:%zu\n", wix);
      if (!bcf_sr_next_line(sr))
        exit(124);
      bcf1_t *line = bcf_sr_get_line(sr, 0);
      int32_t *gt_arr = NULL, ngt_arr = 0;
      int ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);

      for (size_t i = 0; i < n / 2; i++) {
        int32_t *ptr = gt_arr + i * 2;
        // hap 1
        IOBCF_WELLFORMED_CHECK(ptr[0]);
        c[2 * i] |= ((uint64_t)bcf_gt_allele(ptr[0]) << wix);

        // hap 2
        IOBCF_WELLFORMED_CHECK(ptr[1]);
        c[2 * i + 1] |= ((uint64_t)bcf_gt_allele(ptr[1]) << wix);
      }
    }
  } else {
    // WARN: this does not work, at the moment.
    // If we do not fix this, it will be a bloodbath
    // bcf_sr_seek(sr, NULL, 0);
    errno = EPERM;
    perror("NOT IMPLMENTED YET");
    exit(25);
  }
  _li = i;
}

void bfgetcolwgrn(void *fd, size_t n, uint64_t *restrict c, size_t nc,
                  uint8_t w) {
  fprintf(stderr, "\e[0;33m[%s] Not Implemented Yet.\e[0m\n", __func__);
}
void sbfgetcolwgrn(int fd, size_t n, uint64_t *restrict c, size_t nc,
                   uint8_t w) {
  fputs("\e[0;33mMode not used for this type of file. Exiting.\e[0m\n", stderr);
  exit(IOBCF_UNUSED_EXITCODE);
}

void fgetcolwgri(void *fd, size_t i, size_t n, uint64_t *restrict c, size_t nc,
                 uint8_t w) {
  fprintf(stderr, "\e[0;33m[%s] Not Implemented Yet.\e[0m\n", __func__);
}

void sfgetcolwgri(int fd, size_t i, size_t n, uint64_t *restrict c, size_t nc,
                  uint8_t w) {
  fputs("\e[0;33mMode not used for this type of file. Exiting.\e[0m\n", stderr);
  exit(IOBCF_UNUSED_EXITCODE);
}
