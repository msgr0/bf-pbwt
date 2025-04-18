// vim:ft=c
#include "htslib/synced_bcf_reader.h"
#include "htslib/vcf.h"
#include "io.h"
#include <string.h>

#define IOBCF_ASSUME_WELLFORMED

void fgetrc(void *fd, size_t *nr, size_t *nc) {
  bcf_srs_t *sr = fd;
  bcf_hdr_t *hdr = sr->readers[0].header;
  *nr = bcf_hdr_nsamples(hdr) * 2;
  *nc = 0;

  // NOTE: maybe some parts might be rewritten to avoid doing this,
  // however there is not much overhead. For chr10 it takes ~5 secs
  while (bcf_sr_next_line(sr)) {
    (*nc)++;
  }
  // WARN: this does not work
  // bcf_sr_seek(sr, NULL, 0);
}

void fgetcoli(void *fd, size_t i, size_t n, uint8_t *restrict c, size_t nc) {
  // NOTE: nc is not used here.
  bcf_srs_t *sr = fd;
  bcf_hdr_t *hdr = sr->readers[0].header;

  static ssize_t _li = -1;
  // NOTE: I am doing a bit of trickery here assuming that
  // 1. htslib reads lines with incrementing iterator;
  // 2. tipically cols (BCF-row) are read sequentially;
  // I am keeping track of the last position requested and if
  // i == _li+1, then I can just read the next position,
  // otherwise real seeking is necessary
  if (i == _li + 1) {
    // this check should not be necessary
    if (!bcf_sr_next_line(sr))
      exit(124);
    bcf1_t *line = bcf_sr_get_line(sr, 0);
    int32_t *gt_arr = NULL, ngt_arr = 0;
    int ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);

    size_t icol = 0;
    for (size_t i = 0; i < n / 2; i++) {
      int32_t *ptr = gt_arr + i * 2;
      // hap 1
#ifndef IOBCF_ASSUME_WELLFORMED
      if (ptr[0] == bcf_int32_vector_end)
        exit(-2);
      if (bcf_gt_is_missing(ptr[0]))
        exit(-1);
#endif
      c[2 * i] = bcf_gt_allele(ptr[0]);

      // hap 2
#ifndef IOBCF_ASSUME_WELLFORMED
      // if (ptr[1] == bcf_int32_vector_end)
      //   exit(-2);
      // if (bcf_gt_is_missing(ptr[1]))
      //   exit(-1);
#endif
      c[2 * i + 1] = bcf_gt_allele(ptr[1]);
    }
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
#ifndef IOBCF_ASSUME_WELLFORMED
        if (ptr[0] == bcf_int32_vector_end)
          exit(-2);
        if (bcf_gt_is_missing(ptr[0]))
          exit(-1);
#endif
        buf[r * n + 2 * i] = bcf_gt_allele(ptr[0]);

        // hap 2
#ifndef IOBCF_ASSUME_WELLFORMED
        // if (ptr[1] == bcf_int32_vector_end)
        //   exit(-2);
        // if (bcf_gt_is_missing(ptr[1]))
        //   exit(-1);
#endif
        buf[r * n + 2 * i + 1] = bcf_gt_allele(ptr[1]);
      }
    }
    bufn = 0;
  } 
    memcpy(c, buf + (bufn * n), n);
    bufn++;
}

void sbfgetcoln(int fd, size_t n, uint8_t *restrict c, size_t nc) {}

void mbfgetcoln(int fd, size_t n, uint8_t *restrict c, size_t nc) {}

#define FGETCOLIW_IMPL(W)                                                      \
  void fgetcoliw##W(void *fd, size_t i, size_t n, uint64_t *restrict c,        \
                    size_t nc) {}                                              \
  void w##W##mrgsi(size_t n, uint64_t const *wc, uint64_t const *wp,           \
                   uint64_t *restrict c, size_t i) {                           \
    uint64_t c1;                                                               \
    for (size_t r = 0; r < n; r++) {                                           \
    }                                                                          \
  }                                                                            \
  void fgetcoliw##W##r(void *fd, size_t i, size_t n, uint64_t *restrict c,     \
                       size_t nc) {}                                           \
  void wr##W##mrgsi(size_t n, uint64_t const *wc, uint64_t const *wp,          \
                    uint64_t *restrict c, size_t i) {}                         \
  void bfgetcolw##W##rn(void *fd, size_t n, uint64_t *restrict c, size_t nc) { \
  }                                                                            \
  void sbfgetcolw##W##rn(int fd, size_t n, uint64_t *restrict c, size_t nc) {} \
  void sbfgetcolw##W##rn_mmap(int fd, size_t n, uint64_t *restrict c,          \
                              size_t nc) {}

FGETCOLIW_IMPL(8)
FGETCOLIW_IMPL(16)
FGETCOLIW_IMPL(32)
FGETCOLIW_IMPL(64)

void fgetcoliwg(void *fd, size_t i, size_t n, uint64_t *restrict c, size_t nc,
                uint8_t w) {}

void fgetcoliwgr(void *fd, size_t i, size_t n, uint64_t *restrict c, size_t nc,
                 uint8_t w) {}

void bfgetcolwgrn(void *fd, size_t n, uint64_t *restrict c, size_t nc,
                  uint8_t w) {}
void sbfgetcolwgrn(int fd, size_t n, uint64_t *restrict c, size_t nc,
                   uint8_t w) {}

void fgetcolwgri(void *fd, size_t i, size_t n, uint64_t *restrict c, size_t nc,
                 uint8_t w) {}

void sfgetcolwgri(int fd, size_t i, size_t n, uint64_t *restrict c, size_t nc,
                  uint8_t w) {}
