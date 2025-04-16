// vim:ft=c
#ifndef BFPBWT_IO_H
#define BFPBWT_IO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdint.h>

#define BFGETCOLI_BUF_SIZE 128
#define BFGETCOLWR_BUF_SIZE 2

// NOTE: this requires huge pages to be reserved *before* running the tool.
// This program cannot itself reserve huge pages and, on some systems,
// this action requires special privileges;
// hence the presence of the following macro as this feature must be
// intentionally activated and setup beforehand
//
/*#define __USE_MAP_HUGETLB*/
#if defined(MAP_HUGETLB) && defined(__USE_MAP_HUGETLB)
#define __MMAP_FLAGS (MAP_PRIVATE | MAP_HUGETLB)
#else
#define __MMAP_FLAGS (MAP_PRIVATE)
#endif

// get row and column count from file
static inline void fgetrc(FILE *fd, size_t *nr, size_t *nc) {
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

// File GET COLumn I
// get column i from file
// c[n] is a pointer to store the column,
// nc is total number of columns
void fgetcoli(FILE *fd, size_t i, size_t n, uint8_t *c, size_t nc);

// Buffered File GET COLumn Next
// get next column, using buffered freads
// c[n] is a pointer to store the column,
// nc is total number of columns
void bfgetcoln(FILE *fd, size_t n, uint8_t *c, size_t nc);

// Syscall Buffered File GET COLumn Next
// get next column, using buffered syscall preads
// c[n] is a pointer to store the column,
// nc is total number of columns
void sbfgetcoln(int fd, size_t n, uint8_t *c, size_t nc);

// Mmap Buffered File GET COLumn Next
// get next column, using buffered mmap access
// c[n] is a pointer to store the column,
// nc is total number of columns
void mbfgetcoln(int fd, size_t n, uint8_t *c, size_t nc);

#define FGETCOLIW_DECLARE(W)                                                   \
  void fgetcoliw##W(FILE *fd, size_t i, size_t n, uint64_t *c, size_t nc);     \
  void w##W##mrgsi(size_t n, uint64_t const *wc, uint64_t const *wp,           \
                   uint64_t *c, size_t i);                                     \
  void fgetcoliw##W##r(FILE *fd, size_t i, size_t n, uint64_t *c, size_t nc);  \
  void wr##W##mrgsi(size_t n, uint64_t const *wc, uint64_t const *wp,          \
                    uint64_t *c, size_t i);                                    \
  void bfgetcolw##W##rn(FILE *fd, size_t n, uint64_t *c, size_t nc);           \
  void sbfgetcolw##W##rn(int fd, size_t n, uint64_t *c, size_t nc);            \
  void sbfgetcolw##W##rn_mmap(int fd, size_t n, uint64_t *c, size_t nc);

FGETCOLIW_DECLARE(8)
FGETCOLIW_DECLARE(16)
FGETCOLIW_DECLARE(32)
FGETCOLIW_DECLARE(64)

// File GET COLumn I Window General-length
// get column i*w from file,
// returing the bit-packed version of lenght w
// c[n] is a pointer to store the column,
// nc is total number of columns
void fgetcoliwg(FILE *fd, size_t i, size_t n, uint64_t *c, size_t nc,
                uint8_t w);

// File GET COLumn I Window General-length Reversed
// get column i*w from file,
// returing the bit-packed reversed version of lenght w
// c[n] is a pointer to store the column,
// nc is total number of columns
void fgetcoliwgr(FILE *fd, size_t i, size_t n, uint64_t *c, size_t nc,
                 uint8_t w);

// Buffered File GET COLumn Window General-length Reversed Next
// get next i*w column, using buffered freads
// returing the bit-packed reversed version of lenght w
// c[n] is a pointer to store the column,
// nc is total number of columns
void bfgetcolwgrn(FILE *fd, size_t n, uint64_t *c, size_t nc, uint8_t w);

// Syscall Buffered File GET COLumn Window General-length Reversed Next
// get next i*w column, using buffered syscall preads
// returing the bit-packed reversed version of lenght w
// c[n] is a pointer to store the column,
// nc is total number of columns
void sbfgetcolwgrn(int fd, size_t n, uint64_t *c, size_t nc, uint8_t w);

// File GET COLumn Window General-length starting at I
// get w-length bit-packed column starting at position i, using freads
// c[n] is a pointer to store the column,
// nc is total number of columns
//
// This version does not use window size to move in the file.
// It starts reading from column `i` as-is without window-offset computation
void fgetcolwgri(FILE *fd, size_t i, size_t n, uint64_t *c, size_t nc,
                  uint8_t w);

// Syscall File GET COLumn Window General-length starting at I
// get w-length bit-packed column starting at position i, using syscall preads
// c[n] is a pointer to store the column,
// nc is total number of columns
//
// This version does not use window size to move in the file.
// It starts reading from column `i` as-is without window-offset computation
void sfgetcolwgri(int fd, size_t i, size_t n, uint64_t *c, size_t nc,
                   uint8_t w);

// Window MeRGe Shift I:
// merge windows `wc[n]` and `wp[n]`, storing in `c`
// with `n` number of rows, shifting `-i` positions
static inline void wmrgsi(size_t n, uint64_t const *wc, uint64_t const *wp,
                   uint64_t *c, size_t i, uint8_t w) {
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
// with `n` number of rows, shifting `-i` positions
static inline void wrmrgsi(size_t n, uint64_t const *wc, uint64_t const *wp,
                    uint64_t *c, size_t i, uint8_t w) {
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

#endif // !BFPBWT_IO_H
