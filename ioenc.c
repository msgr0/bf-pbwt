// vim:ft=c
#include "io.h"
#include <assert.h>
#include <string.h>

#define HEADER_BYTES 2

#ifndef IOENC_UNUSED_EXITCODE
#define IOENC_UNUSED_EXITCODE 3
#endif
void fgetrc(void *fd, size_t *nr, size_t *nc) {
  fseek(fd, 0, SEEK_SET);
  *nc = *nr = 0;
  fread(nr, sizeof(uint32_t), 1, fd);
  fread(nc, sizeof(uint32_t), 1, fd);
  size_t nfc = 0; 
  // true column size in the file with padding
  nfc = (*nc) + (((*nc) % 8) ? (8 - (*nc) % 8) : 0);
  // (*nc) = ((8* (*nc))+7) / 8;
  printf("rows: %zu, fileCols: %zu, headerCols: %zu\n", *nr, *nc, nfc );
  (*nc) = nfc;

  // count byte frequency -- possibile alphabet encoding as in syllable?

  // uint64_t freq[256];
  // size_t MM = 0;
  // for(size_t i=0; i<*nc; i=i+8) {
  //   memset(freq, 0, sizeof(freq));
  //   char* col = (char*) malloc(*nr * sizeof(*col));
  //   fread(col, sizeof(uint8_t), *nr, fd);
  //   int M=0;
  //   for(size_t j=0; j<*nr; j++) {
  //     uint8_t i = col[j];
  //     ++freq[i];
  //     M += (freq[i] == 1);
  //     MM = MM < M ? M : MM;
  //   }
  // }

  // printf("%zu\n", MM);
}

int fgetcoli(void *fd, size_t i, size_t n, uint8_t *restrict c, size_t nc) {
  // probably unused!! gets by windows and then shift to the correct bit
  int x;
  size_t ie = i / 8;
  size_t off = i % 8;
  fseek(fd, HEADER_BYTES + ie * nc, SEEK_SET);
  fread(c, sizeof(uint8_t), n, fd);
  for (size_t t = 0; t < n; t++) {
    c[t] = (c[t] >> (7 - off)) & 1;
  }
  return 1;
}

void bfgetcoln(void *fd, size_t n, uint8_t *restrict c, size_t nc) {
  fputs("\e[0;33mMode not used for this type of file. Exiting.\e[0m\n", stderr);
  exit(IOENC_UNUSED_EXITCODE);
}

void sbfgetcoln(int fd, size_t n, uint8_t *restrict c, size_t nc) {
  fputs("\e[0;33mMode not used for this type of file. Exiting.\e[0m\n", stderr);
  exit(IOENC_UNUSED_EXITCODE);
}

void mbfgetcoln(int fd, size_t n, uint8_t *restrict c, size_t nc) {
  fputs("\e[0;33mMode not used for this type of file. Exiting.\e[0m\n", stderr);
  exit(IOENC_UNUSED_EXITCODE);
}

#define FGETCOLIW_IMPL(W)                                                      \
  void fgetcoliw##W(void *fd, size_t i, size_t n, uint64_t *restrict c,        \
                    size_t nc) {                                               \
    int x;                                                                     \
    assert(W % 8 == 0);                                                        \
    size_t mult = W / 8;                                                       \
    do {                                                                       \
      fread(c, sizeof(uint8_t), n, fd);                                        \
    } while (--mult);                                                          \
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
    int x;                                                                     \
    assert(W % 8 == 0);                                                        \
    size_t mult = W / 8;                                                       \
    do {                                                                       \
      fread(c, sizeof(uint8_t), n, fd);                                        \
    } while (--mult);                                                          \
    return 1;                                                                  \
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
    static size_t i = 0;                                                       \
    static size_t bufn = BFGETCOLWR_BUF_SIZE;                                  \
    static uint64_t *buf = NULL;                                               \
    if (!buf)                                                                  \
      buf = malloc(BFGETCOLWR_BUF_SIZE * n * sizeof *buf);                     \
    char rbuf[64 * BFGETCOLWR_BUF_SIZE];                                       \
                                                                               \
    if (bufn == BFGETCOLWR_BUF_SIZE) {                                         \
      uint64_t x;                                                              \
      fseek(fd, i *W, SEEK_SET);                                               \
      for (size_t r = 0; r < n; r++) {                                         \
        fread(&rbuf, 1, 64 * BFGETCOLWR_BUF_SIZE, fd);                         \
        for (size_t s = 0; s < BFGETCOLWR_BUF_SIZE; s++) {                     \
          buf[(n * s) + r] = 0;                                                \
          for (size_t j = 0; j < W; j++) {                                     \
            x = rbuf[W * s + j] - 48;                                          \
            buf[(n * s) + r] |= (x << j);                                      \
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
  }                                                                            \
  int sbfgetcolw##W##rn(int fd, size_t n, uint64_t *restrict c, size_t nc) {   \
    static size_t i = 0;                                                       \
    static size_t bufn = BFGETCOLWR_BUF_SIZE;                                  \
    static uint64_t *buf = NULL;                                               \
    if (!buf)                                                                  \
      buf = malloc(BFGETCOLWR_BUF_SIZE * n * sizeof *buf);                     \
    char rbuf[64 * BFGETCOLWR_BUF_SIZE];                                       \
                                                                               \
    if (bufn == BFGETCOLWR_BUF_SIZE) {                                         \
      uint64_t x;                                                              \
      size_t offset;                                                           \
      offset = i * W;                                                          \
      /*_Pragma("omp parallel for") */                                         \
      for (size_t r = 0; r < n; r++) {                                         \
        pread(fd, &rbuf, 64 * BFGETCOLWR_BUF_SIZE, offset);                    \
        for (size_t s = 0; s < BFGETCOLWR_BUF_SIZE; s++) {                     \
          buf[(n * s) + r] = 0;                                                \
          for (size_t j = 0; j < W; j++) {                                     \
            x = rbuf[W * s + j] - 48;                                          \
            buf[(n * s) + r] |= (x << j);                                      \
          }                                                                    \
        }                                                                      \
        c[r] = buf[r];                                                         \
        offset += nc + 1;                                                      \
      }                                                                        \
      bufn = 1;                                                                \
    } else {                                                                   \
      for (size_t r = 0; r < n; r++) {                                         \
        c[r] = buf[(n * bufn) + r];                                            \
      }                                                                        \
      bufn++;                                                                  \
    }                                                                          \
    i++;                                                                       \
    return 1;                                                                  \
  }                                                                            \
  void sbfgetcolw##W##rn_mmap(int fd, size_t n, uint64_t *restrict c,          \
                              size_t nc) {                                     \
    static size_t i = 0;                                                       \
    static size_t bufn = BFGETCOLWR_BUF_SIZE;                                  \
    static uint64_t *buf = NULL;                                               \
    if (!buf)                                                                  \
      buf = malloc(BFGETCOLWR_BUF_SIZE * n * sizeof *buf);                     \
    char rbuf[64 * BFGETCOLWR_BUF_SIZE];                                       \
    static uint8_t *fdmm = NULL;                                               \
    if (!fdmm) {                                                               \
      struct stat st;                                                          \
      if (fstat(fd, &st) < 0) {                                                \
        perror("fstat");                                                       \
        exit(EXIT_FAILURE);                                                    \
      }                                                                        \
      fdmm = mmap(NULL, st.st_size, PROT_READ, __MMAP_FLAGS, fd, 0);           \
      if (fdmm == MAP_FAILED) {                                                \
        perror("mmap");                                                        \
        exit(EXIT_FAILURE);                                                    \
      }                                                                        \
    }                                                                          \
    if (bufn == BFGETCOLWR_BUF_SIZE) {                                         \
      uint64_t x;                                                              \
      size_t offset;                                                           \
      offset = i * W;                                                          \
      for (size_t r = 0; r < n; r++) {                                         \
        memcpy(&rbuf, &fdmm[offset], 64 * BFGETCOLWR_BUF_SIZE);                \
        for (size_t s = 0; s < BFGETCOLWR_BUF_SIZE; s++) {                     \
          buf[(n * s) + r] = 0;                                                \
          for (size_t j = 0; j < W; j++) {                                     \
            x = rbuf[W * s + j] - 48;                                          \
            buf[(n * s) + r] |= (x << j);                                      \
          }                                                                    \
        }                                                                      \
        c[r] = buf[r];                                                         \
        offset += nc + 1;                                                      \
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

FGETCOLIW_IMPL(8)
FGETCOLIW_IMPL(16)
FGETCOLIW_IMPL(32)
FGETCOLIW_IMPL(64)

void fgetcoliwg(void *fd, size_t i, size_t n, uint64_t *restrict c, size_t nc,
                uint8_t w) {
  
  /// enc file is presented 1byte (8cols) of nrows at time, 
  /// w is the generic window size, must be a multiple of 8
  /// Current i serves for seeking, but currently we can assume no jumps
  /// on the file. 

  size_t wmult = w/8; // how many time we have to grab nrows
  assert( w%8 ==0);
  
  uint8_t *crow = malloc(n * sizeof(uint8_t));
  
  for (size_t t = 0; t < wmult; t++) { 
    for (size_t j = 0; j < n-1; j++) {
    fread (crow, sizeof(uint8_t), n, fd);
      c[j] = (c[j] << 8)  | (crow[j]); 
    }
    fread (crow, sizeof(uint8_t), n, fd);
    c[n-2] |= crow[n-2];
  }

}

int fgetcoliwgr(void *fd, size_t i, size_t n, uint64_t *restrict c, size_t nc,
                uint8_t w) {
  // uint64_t x;
  // fseek(fd, i * w, SEEK_SET);
  // for (size_t r = 0; r < n; r++) {
  //   c[r] = 0;
  //   /*printf("r(");*/
  //   for (size_t s = 0; s < w; s++) {
  //     x = fgetc(fd) - 48;
  //     /*printf("%llu", x);*/
  //     c[r] = (x << s) | c[r];
  //   }
  //   /*printf(")=%llu\n", c[r]);*/
  //   fseek(fd, nc - w + 1, SEEK_CUR);
  // }
  return 1;
}

void bfgetcolwgrn(void *fd, size_t n, uint64_t *restrict c, size_t nc,
                  uint8_t w) {
  //
  // static size_t i = 0;
  // static size_t bufn = BFGETCOLWR_BUF_SIZE;
  // static uint64_t *buf = NULL;
  // if (!buf)
  //   buf = malloc(BFGETCOLWR_BUF_SIZE * n * sizeof *buf);
  // char rbuf[64 * BFGETCOLWR_BUF_SIZE];
  //
  // if (bufn == BFGETCOLWR_BUF_SIZE) {
  //   uint64_t x;
  //   fseek(fd, i * w, SEEK_SET);
  //   for (size_t r = 0; r < n; r++) {
  //     fread(&rbuf, 1, 64 * BFGETCOLWR_BUF_SIZE, fd);
  //     for (size_t s = 0; s < BFGETCOLWR_BUF_SIZE; s++) {
  //       buf[(n * s) + r] = 0;
  //       for (size_t j = 0; j < w; j++) {
  //         /*x = fgetc(fd) - 48;*/
  //         x = rbuf[w * s + j] - 48;
  //         buf[(n * s) + r] = (x << j) | buf[(n * s) + r];
  //       }
  //     }
  //     c[r] = buf[r];
  //     fseek(fd, nc - (BFGETCOLWR_BUF_SIZE * w) + 1, SEEK_CUR);
  //   }
  //   bufn = 1;
  // } else {
  //   for (size_t r = 0; r < n; r++) {
  //     c[r] = buf[(n * bufn) + r];
  //   }
  //   bufn++;
  // }
  // i++;
}
void sbfgetcolwgrn(int fd, size_t n, uint64_t *restrict c, size_t nc,
                   uint8_t w) {
  //
  // static size_t i = 0;
  // static size_t bufn = BFGETCOLWR_BUF_SIZE;
  // static uint64_t *buf = NULL;
  // if (!buf)
  //   buf = malloc(BFGETCOLWR_BUF_SIZE * n * sizeof *buf);
  // char rbuf[64 * BFGETCOLWR_BUF_SIZE];
  //
  // if (bufn == BFGETCOLWR_BUF_SIZE) {
  //   uint64_t x;
  //   size_t offset;
  //   /*lseek(fd, i * w, SEEK_SET);*/
  //   offset = i * w;
  //   for (size_t r = 0; r < n; r++) {
  //     pread(fd, &rbuf, 64 * BFGETCOLWR_BUF_SIZE, offset);
  //     for (size_t s = 0; s < BFGETCOLWR_BUF_SIZE; s++) {
  //       buf[(n * s) + r] = 0;
  //       /*pread(fd, &rbuf, w, offset + w * s);*/
  //       for (size_t j = 0; j < w; j++) {
  //         x = rbuf[w * s + j] - 48;
  //         buf[(n * s) + r] = (x << j) | buf[(n * s) + r];
  //       }
  //     }
  //     c[r] = buf[r];
  //     offset += nc + 1;
  //   }
  //   bufn = 1;
  // } else {
  //   for (size_t r = 0; r < n; r++) {
  //     c[r] = buf[(n * bufn) + r];
  //   }
  //   bufn++;
  // }
  // i++;
}

void fgetcolwgri(void *fd, size_t i, size_t n, uint64_t *restrict c, size_t nc,
                 uint8_t w) {
  // uint64_t x;
  // fseek(fd, i, SEEK_SET);
  // for (size_t r = 0; r < n; r++) {
  //   c[r] = 0;
  //   /*printf("r(");*/
  //   for (size_t s = 0; s < w; s++) {
  //     x = fgetc(fd) - 48;
  //     /*printf("%llu", x);*/
  //     c[r] = (x << s) | c[r];
  //   }
  //   /*printf(")=%llu\n", c[r]);*/
  //   fseek(fd, nc - w + 1, SEEK_CUR);
  // }
}

void sfgetcolwgri(int fd, size_t i, size_t n, uint64_t *restrict c, size_t nc,
                  uint8_t w) {
  // uint64_t x;
  // uint8_t _x;
  // lseek(fd, i, SEEK_SET);
  // for (size_t r = 0; r < n; r++) {
  //   c[r] = 0;
  //   /*printf("r(");*/
  //   for (size_t s = 0; s < w; s++) {
  //     read(fd, &_x, 1);
  //     x = _x - 48;
  //     /*printf("%llu", x);*/
  //     c[r] = (x << s) | c[r];
  //   }
  //   /*printf(")=%llu\n", c[r]);*/
  //   lseek(fd, nc - w + 1, SEEK_CUR);
  // }
}
