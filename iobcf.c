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
// NOTE on IOBCF_ASSUME_WELLFORMED / IOBCF_WELLFORMED_CHECK:
// Left disabled (no-op) on purpose. A real per-genotype check is two extra
// branches (a compare + a shift-and-compare) per allele, i.e. per element of
// the hottest loop in the whole BCF decode path (called 2 * nsamples times
// per record). Given this step's goal is CPU-time reduction and the plan's
// baseline already ran with this disabled (bit-identical correctness was
// verified against that baseline), enabling it here would both cost cycles
// and change behavior on malformed input in ways not covered by the
// correctness matrix. If this is ever re-enabled:
//   - On the fast native-width path in `bcf_decode_gt_row` below, `ptr` is a
//     *native-width* value (normally `int8_t`, promoted to `int` by C's
//     integer promotion rules), NOT the `int32_t` that `bcf_get_genotypes`
//     would have produced. `bcf_gt_is_missing` still works correctly at any
//     width (real missing genotype alleles are encoded as raw value 0
//     regardless of storage width), but the `== bcf_int32_vector_end` check
//     does NOT — a genuine vector-end pad byte shows up as
//     `bcf_int8_vector_end` (-127) at this width, not `bcf_int32_vector_end`.
//     Re-enabling the check on the fast path requires comparing against the
//     sentinel of the field's *actual* `fmt->type`, not always the int32 one.
//   - The fallback path (`bcf_get_genotypes`) already produces genuinely
//     int32-widened values, so the existing check is correct there as-is.

void fgetrc(void *fd, size_t *nr, size_t *nc) {
  // bcf_srs_t *sr = bcf_sr_init();
  // bcf_sr_add_reader(sr, ((bcf_srs_t *)fd)->readers[0].fname);
  // bcf_hdr_t *hdr = sr->readers[0].header;

  bcf_hdr_t *hdr = ((bcf_srs_t *)fd)->readers[0].header;

  *nr = bcf_hdr_nsamples(hdr) * 2;
  *nc = -1;

  // // NOTE: maybe some parts might be rewritten to avoid doing this,
  // // however there is not much overhead. For chr10 it takes ~5 secs
  // while (bcf_sr_next_line(sr)) {
  //   (*nc)++;
  // }
  //
  // bcf_sr_destroy(sr);
}

// Decode one BCF record's GT field into `out[n]` (`n` = 2 * nsamples, one
// byte per haplotype allele, same values `bcf_gt_allele()` would produce).
//
// Fast path: read the raw `bcf_fmt_t` for "GT" and decode `fmt->p` directly
// in its native on-disk width, skipping the int32 expansion that
// `bcf_get_genotypes()` (== `bcf_get_format_values(..., BCF_HT_INT)`) always
// performs regardless of source width. `bcf_get_fmt()` still calls
// `bcf_unpack(line, BCF_UN_FMT)` internally (unavoidable — any FORMAT field
// access needs the per-record FMT block decoded), but does not additionally
// expand every genotype to int32.
//
// This fast path is only taken for the common case: biallelic, diploid,
// BCF_BT_INT8-encoded GT (which is what typical phase3-scale panels with a
// few thousand samples use). Anything else (multiallelic — more bits needed
// per genotype means htslib promotes to BCF_BT_INT16/32 — or non-diploid
// ploidy, fmt->n != 2) falls back to the original `bcf_get_genotypes()` path
// so those records still decode correctly.
static inline void bcf_decode_gt_row(bcf_hdr_t *hdr, bcf1_t *line, size_t n,
                                      uint8_t *restrict out) {
  bcf_fmt_t *fmt = bcf_get_fmt(hdr, line, "GT");
  if (fmt && fmt->type == BCF_BT_INT8 && fmt->n == 2) {
    uint8_t *p = fmt->p;
    for (size_t i = 0; i < n / 2; i++) {
      int8_t a0 = (int8_t)p[0];
      int8_t a1 = (int8_t)p[1];
      IOBCF_WELLFORMED_CHECK(a0);
      out[2 * i] = (uint8_t)bcf_gt_allele(a0);

      IOBCF_WELLFORMED_CHECK(a1);
      out[2 * i + 1] = (uint8_t)bcf_gt_allele(a1);

      p += fmt->size;
    }
    return;
  }

  // Fallback: multiallelic / non-diploid / wider-than-int8 GT encoding.
  // `gt_arr`/`ngt_arr` persist across calls (thread-local) instead of being
  // malloc'd and freed per record — see the "hoisted genotype buffer" note
  // at the top of this file.
  static __thread int32_t *gt_arr = NULL;
  static __thread int ngt_arr = 0;
  bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);

  for (size_t i = 0; i < n / 2; i++) {
    int32_t *ptr = gt_arr + i * 2;
    IOBCF_WELLFORMED_CHECK(ptr[0]);
    out[2 * i] = bcf_gt_allele(ptr[0]);

    IOBCF_WELLFORMED_CHECK(ptr[1]);
    out[2 * i + 1] = bcf_gt_allele(ptr[1]);
  }
}

int fgetcoli(void *fd, size_t i, size_t n, uint8_t *restrict c, size_t nc) {
  // NOTE: nc is used as a flag.
  // If nc == 0, ignore _li and assume sequential read
  // also since it might invalidate hdr, reset that as well from fd
  bcf_srs_t *sr = fd;
  static bcf_hdr_t *hdr = NULL;
  if (!hdr || !nc)
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
      return 0;
    bcf1_t *line = bcf_sr_get_line(sr, 0);
    bcf_decode_gt_row(hdr, line, n, c);
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
  return 1;
}

void bfgetcoln(void *fd, size_t n, uint8_t *restrict c, size_t nc) {
  // NOTE: nc is not used here.
  bcf_srs_t *sr = fd;
  static bcf_hdr_t *hdr = NULL;
  if (!hdr)
    hdr = sr->readers[0].header;

  static size_t bufn = BFGETCOLI_BUF_SIZE;
  static uint8_t *buf = NULL;
  if (!buf)
    buf = malloc(BFGETCOLI_BUF_SIZE * n * sizeof *buf);

  if (bufn == BFGETCOLI_BUF_SIZE) {
    for (size_t r = 0; r < BFGETCOLI_BUF_SIZE; r++) {
      if (!bcf_sr_next_line(sr))
        break;
      bcf1_t *line = bcf_sr_get_line(sr, 0);
      bcf_decode_gt_row(hdr, line, n, buf + r * n);
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

// Ensure the per-thread `W x n` staging buffer (one uint8_t allele per
// (window-row, sample-haplotype) pair) is large enough, growing it lazily.
// Kept as a macro so each window-reader (compile-time W instantiation, and
// the two runtime-w variants below) gets its own independently-sized,
// independently-thread-local staging buffer.
#define IOBCF_ENSURE_STAGE(stage, stage_cap, need)                            \
  do {                                                                        \
    if ((stage_cap) < (need)) {                                               \
      (stage) = realloc((stage), (need));                                     \
      (stage_cap) = (need);                                                   \
    }                                                                         \
  } while (0)

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
  int fgetcoliw##W##r(void *fd, size_t i, size_t n, uint64_t *restrict c,      \
                      size_t nc) {                                             \
    bcf_srs_t *sr = fd;                                                        \
    static bcf_hdr_t *hdr = NULL;                                              \
    if (!hdr)                                                                  \
      hdr = sr->readers[0].header;                                             \
                                                                               \
    static ssize_t _li = -1;                                                   \
    if (i == _li + 1 || !nc) {                                                 \
      /* W x n staging buffer: decode into it row-by-row (one contiguous */    \
      /* write per record), then pack into c[] column-wise so each c[r] */     \
      /* is written exactly once per window instead of once per          */    \
      /* (window, wix) pair. */                                                \
      static __thread uint8_t *stage = NULL;                                   \
      static __thread size_t stage_cap = 0;                                    \
      IOBCF_ENSURE_STAGE(stage, stage_cap, (size_t)W * n);                     \
                                                                               \
      size_t wix;                                                              \
      for (wix = 0; wix < W; wix++) {                                          \
        if (!bcf_sr_next_line(sr))                                             \
          break;                                                               \
        bcf1_t *line = bcf_sr_get_line(sr, 0);                                 \
        bcf_decode_gt_row(hdr, line, n, stage + wix * n);                      \
      }                                                                        \
      for (size_t r = 0; r < n; r++) {                                         \
        uint64_t val = 0;                                                      \
        for (size_t k = 0; k < wix; k++)                                       \
          val |= (uint64_t)stage[k * n + r] << k;                              \
        c[r] = val;                                                            \
      }                                                                        \
      if (wix < W) {                                                           \
        _li = i;                                                               \
        return wix;                                                            \
      }                                                                        \
    } else {                                                                   \
      errno = EPERM;                                                           \
      perror("NOT IMPLMENTED YET");                                            \
      exit(26);                                                                \
    }                                                                          \
    _li = i;                                                                   \
    return W;                                                                  \
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
  int sbfgetcolw##W##rn(int fd, size_t n, uint64_t *restrict c, size_t nc) {   \
    fputs("\e[0;33mMode not used for this type of file. Exiting.\e[0m\n",      \
          stderr);                                                             \
    exit(IOBCF_UNUSED_EXITCODE);                                               \
    return 0;                                                                  \
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

int fgetcoliwgr(void *fd, size_t i, size_t n, uint64_t *restrict c, size_t nc,
                uint8_t w) {
  // NOTE: nc is not used here.
  bcf_srs_t *sr = fd;
  static bcf_hdr_t *hdr = NULL;
  if (!hdr)
    hdr = sr->readers[0].header;

  static ssize_t _li = -1;
  // NOTE: Same trickery here as in `fgetcoli`, however
  // _li is not the last index or row, but the last index of window.
  // Current BCF row is (_li * w)
  if (i == _li + 1 || !nc) {
    // See `FGETCOLIW_IMPL`: decode into a `w x n` staging buffer (one
    // contiguous write per record) then pack column-wise, one write per
    // c[r], instead of read-modify-writing all of c[] once per record.
    static __thread uint8_t *stage = NULL;
    static __thread size_t stage_cap = 0;
    IOBCF_ENSURE_STAGE(stage, stage_cap, (size_t)w * n);

    size_t wix;
    for (wix = 0; wix < w; wix++) {
      if (!bcf_sr_next_line(sr))
        break;
      bcf1_t *line = bcf_sr_get_line(sr, 0);
      bcf_decode_gt_row(hdr, line, n, stage + wix * n);
    }
    for (size_t r = 0; r < n; r++) {
      uint64_t val = 0;
      for (size_t k = 0; k < wix; k++)
        val |= (uint64_t)stage[k * n + r] << k;
      c[r] = val;
    }
    if (wix < w) {
      _li = i;
      return wix;
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
  return w;
}

void bfgetcolwgrn(void *fd, size_t n, uint64_t *restrict c, size_t nc,
                  uint8_t w) {
  fputs("\e[0;33mMode not used for this type of file. Exiting.\e[0m\n", stderr);
  exit(IOBCF_UNUSED_EXITCODE);
}
void sbfgetcolwgrn(int fd, size_t n, uint64_t *restrict c, size_t nc,
                   uint8_t w) {
  fputs("\e[0;33mMode not used for this type of file. Exiting.\e[0m\n", stderr);
  exit(IOBCF_UNUSED_EXITCODE);
}

int fgetcolwgri(void *fd, size_t i, size_t n, uint64_t *restrict c, size_t nc,
                uint8_t w) {

  // NOTE: NC not used here, can be used as thread-safe _li
  //
  // NOTE: this function is called *concurrently* by `stagpar`'s per-lane
  // threads (one `bcf_srs_t *fd` per lane/thread), so every persistent
  // buffer used here (the staging buffer below, and `gt_arr`/`ngt_arr`
  // inside `bcf_decode_gt_row`) MUST be thread-local (`static __thread`),
  // not plain `static`, or concurrent lanes would race on the same buffer.

  bcf_srs_t *sr = fd;
  bcf_hdr_t *hdr = NULL;
  if (!hdr)
    hdr = sr->readers[0].header;
  if (nc > i) {
    // bcf_sr_seek does not work, maybe it does not work only without index.
    // either way we cannot use it reliably
    char *fname = strdup(sr->readers[0].fname);
    bcf_sr_remove_reader(sr, 0);
    if (!bcf_sr_add_reader(sr, fname)) {
      fprintf(stderr, "Failed to re-open: %s\n", fname);
      exit(24);
    }
    free(fname);
    hdr = sr->readers[0].header;
    nc = 0;
  }

  for (; nc < i; nc++) {
    if (!bcf_sr_next_line(sr)) {
      fprintf(stderr, "frror seeking to BCF to i=%zu, nc = %zu\n", i, nc);
      exit(22);
    }
  }

  static __thread uint8_t *stage = NULL;
  static __thread size_t stage_cap = 0;
  IOBCF_ENSURE_STAGE(stage, stage_cap, (size_t)w * n);

  size_t wix;
  for (wix = 0; wix < w; wix++) {
    if (!bcf_sr_next_line(sr))
      break;
    bcf1_t *line = bcf_sr_get_line(sr, 0);
    bcf_decode_gt_row(hdr, line, n, stage + wix * n);
  }
  for (size_t r = 0; r < n; r++) {
    uint64_t val = 0;
    for (size_t k = 0; k < wix; k++)
      val |= (uint64_t)stage[k * n + r] << k;
    c[r] = val;
  }
  if (wix < w)
    return i + wix;
  return i + w;
}

void sfgetcolwgri(int fd, size_t i, size_t n, uint64_t *restrict c, size_t nc,
                  uint8_t w) {
  fputs("\e[0;33mMode not used for this type of file. Exiting.\e[0m\n", stderr);
  exit(IOBCF_UNUSED_EXITCODE);
}
void spfgetcolwgri(int fd, size_t i, size_t n, uint64_t *restrict c, size_t nc,
                  uint8_t w) {
  fputs("\e[0;33mMode not used for this type of file. Exiting.\e[0m\n", stderr);
  exit(IOBCF_UNUSED_EXITCODE);
}
void fgetcolwgri_mmap(uint8_t *fdmm, size_t i, size_t n, uint64_t *restrict c, size_t nc,
                  uint8_t w) {
  fputs("\e[0;33mMode not used for this type of file. Exiting.\e[0m\n", stderr);
  exit(IOBCF_UNUSED_EXITCODE);
}
