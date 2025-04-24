// vim:ft=c
#include "io.h"
#include <string.h>

#define RCBUFSIZE 8192

void fgetrc(void *fd, size_t *nr, size_t *nc) {
  unsigned char buf[RCBUFSIZE];
  size_t bytes;
  fseek(fd, 0, SEEK_SET);

  *nc = *nr = 0;
  int c;
  while ((c = fgetc(fd)) != 0xA)
    (*nc)++;

  (*nr)++;
#if 0
  while ((bytes = fread(buf, 1, RCBUFSIZE, fd)) > 0) {
    for (size_t i = 0; i < bytes; ++i)
      (*nr) += buf[i] == 0xA;
  }
#else
  while (1) {
    char buf[RCBUFSIZE + 1];
    size_t bytes_read = fread(buf, 1, RCBUFSIZE, fd);
    if (bytes_read <= 0)
      return;

    bytes += bytes_read;
    char *end = buf + bytes_read;
    size_t buflines = 0;

    *end = '\n';
    for (char *p = buf; (p = memchr(p, '\n', RCBUFSIZE+1)) < end; p++) {
      buflines++;
    }

    *nr += buflines;
  }
#endif
}

int fgetcoli(void *fd, size_t i, size_t n, uint8_t *restrict c, size_t nc) {
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
  return 1;
}

// get column i from file
// c[n] is a pointer to store the column,
// nc is total number of columns, needed for fseek
void bfgetcoln(void *fd, size_t n, uint8_t *restrict c, size_t nc) {
  // NOTE: this assumes ASCII text file, offset are computed assuming
  // 1-byte size for each character

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

void sbfgetcoln(int fd, size_t n, uint8_t *restrict c, size_t nc) {
  // NOTE: this assumes ASCII text file, offset are computed assuming
  // 1-byte size for each character

  static char rbuf[BFGETCOLI_BUF_SIZE];
  static size_t i = 0;
  static size_t bufn = BFGETCOLI_BUF_SIZE;
  static uint8_t *buf = NULL;
  if (!buf)
    buf = malloc(BFGETCOLI_BUF_SIZE * n * sizeof *buf);

  size_t offset = i;
  if (bufn == BFGETCOLI_BUF_SIZE) {
    /*#pragma omp parallel for*/
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

void mbfgetcoln(int fd, size_t n, uint8_t *restrict c, size_t nc) {
  // NOTE: this assumes ASCII text file, offset are computed assuming
  // 1-byte size for each character

  static char rbuf[BFGETCOLI_BUF_SIZE];
  static size_t i = 0;
  static size_t bufn = BFGETCOLI_BUF_SIZE;
  static uint8_t *buf = NULL;
  if (!buf)
    buf = malloc(BFGETCOLI_BUF_SIZE * n * sizeof *buf);

  static uint8_t *fdmm = NULL;
  if (!fdmm) {
    struct stat st;
    if (fstat(fd, &st) < 0) {
      perror("fstat");
      exit(EXIT_FAILURE);
    }
    fdmm = mmap(NULL, st.st_size, PROT_READ, __MMAP_FLAGS, fd, 0);
    if (fdmm == MAP_FAILED) {
      perror("mmap");
      exit(EXIT_FAILURE);
    }
  }

  size_t offset = i;
  if (bufn == BFGETCOLI_BUF_SIZE) {
    /*#pragma omp parallel for*/
    for (size_t r = 0; r < n; r++) {
      memcpy(&buf[r * BFGETCOLI_BUF_SIZE], &fdmm[offset], BFGETCOLI_BUF_SIZE);
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

#define FGETCOLIW_IMPL(W)                                                      \
  void fgetcoliw##W(void *fd, size_t i, size_t n, uint64_t *restrict c,        \
                    size_t nc) {                                               \
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
  void sbfgetcolw##W##rn(int fd, size_t n, uint64_t *restrict c, size_t nc) {  \
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

void fgetcoliwgr(void *fd, size_t i, size_t n, uint64_t *restrict c, size_t nc,
                 uint8_t w) {
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

void bfgetcolwgrn(void *fd, size_t n, uint64_t *restrict c, size_t nc,
                  uint8_t w) {
  // NOTE: this assumes ASCII text file, offset are computed assuming
  // 1-byte size for each character

  static size_t i = 0;
  static size_t bufn = BFGETCOLWR_BUF_SIZE;
  static uint64_t *buf = NULL;
  if (!buf)
    buf = malloc(BFGETCOLWR_BUF_SIZE * n * sizeof *buf);
  char rbuf[64 * BFGETCOLWR_BUF_SIZE];

  if (bufn == BFGETCOLWR_BUF_SIZE) {
    uint64_t x;
    fseek(fd, i * w, SEEK_SET);
    for (size_t r = 0; r < n; r++) {
      fread(&rbuf, 1, 64 * BFGETCOLWR_BUF_SIZE, fd);
      for (size_t s = 0; s < BFGETCOLWR_BUF_SIZE; s++) {
        buf[(n * s) + r] = 0;
        for (size_t j = 0; j < w; j++) {
          /*x = fgetc(fd) - 48;*/
          x = rbuf[w * s + j] - 48;
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
void sbfgetcolwgrn(int fd, size_t n, uint64_t *restrict c, size_t nc,
                   uint8_t w) {
  // NOTE: this assumes ASCII text file, offset are computed assuming
  // 1-byte size for each character

  static size_t i = 0;
  static size_t bufn = BFGETCOLWR_BUF_SIZE;
  static uint64_t *buf = NULL;
  if (!buf)
    buf = malloc(BFGETCOLWR_BUF_SIZE * n * sizeof *buf);
  char rbuf[64 * BFGETCOLWR_BUF_SIZE];

  if (bufn == BFGETCOLWR_BUF_SIZE) {
    uint64_t x;
    size_t offset;
    /*lseek(fd, i * w, SEEK_SET);*/
    offset = i * w;
    for (size_t r = 0; r < n; r++) {
      pread(fd, &rbuf, 64 * BFGETCOLWR_BUF_SIZE, offset);
      for (size_t s = 0; s < BFGETCOLWR_BUF_SIZE; s++) {
        buf[(n * s) + r] = 0;
        /*pread(fd, &rbuf, w, offset + w * s);*/
        for (size_t j = 0; j < w; j++) {
          x = rbuf[w * s + j] - 48;
          buf[(n * s) + r] = (x << j) | buf[(n * s) + r];
        }
      }
      c[r] = buf[r];
      offset += nc + 1;
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

void fgetcolwgri(void *fd, size_t i, size_t n, uint64_t *restrict c, size_t nc,
                 uint8_t w) {
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

void sfgetcolwgri(int fd, size_t i, size_t n, uint64_t *restrict c, size_t nc,
                  uint8_t w) {
  // NOTE: this assumes ASCII text file, offset are computed assuming
  // 1-byte size for each character
  uint64_t x;
  uint8_t _x;
  lseek(fd, i, SEEK_SET);
  for (size_t r = 0; r < n; r++) {
    c[r] = 0;
    /*printf("r(");*/
    for (size_t s = 0; s < w; s++) {
      read(fd, &_x, 1);
      x = _x - 48;
      /*printf("%llu", x);*/
      c[r] = (x << s) | c[r];
    }
    /*printf(")=%llu\n", c[r]);*/
    lseek(fd, nc - w + 1, SEEK_CUR);
  }
}
