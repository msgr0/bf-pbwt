#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BITS 8

int main(int argc, char *argv[]) {
  if (argc < 4) {
    fprintf(stderr, "usage: %s enc|txt ROWS COLS [err]\n", argv[0]);
    return EXIT_FAILURE;
  }
  size_t rows = atoi(argv[2]);
  size_t cols = atoi(argv[3]);

  float err = 0.0;
  srand(42);

  if (argc > 4 && strcmp(argv[1], "enc") == 0) {

    // printf("doing ENC\n");
    err = atof(argv[4]);

    size_t padding = cols % 8;
    size_t nw = (cols / 8) + (padding > 0);
    size_t fcols = cols + padding;

    char *frow = (char *)malloc(fcols + sizeof(char));

    for (size_t i = 0; i < cols; i++) {
      frow[i] = (rand() & 1U);
      //printf("%c", frow[i] + 48);
    }
    for (size_t i = 0; i < padding; i++) {
      frow[i + cols] = 0;
    }
    //printf("b-nw:%zu-b\n", nw);

    char *ftcol_pack = (char *)malloc(rows + sizeof(char));
    size_t b = 0;
    for (size_t j = 0; j < nw; j++) { // per ogni w-sized colonna
      for (size_t i = 0; i < rows; i++) { // per ogni riga
        for (size_t b = 8; b >0; b--) { // per ogni bit della w-sized colonna
          ftcol_pack[i] = (ftcol_pack[i] << 1) | ( rand() > (int)(err * INT32_MAX) ? (frow[(j)*8 + b-1] & 1U) : (rand() & 1U) );
        }
        putc(ftcol_pack[i], stdout);
      }
      puts("");
    }


//    char *frow_pack = (char *)malloc(cols + sizeof(char));
//    size_t b = 0;
//    for (size_t i = 0; i < cols; i++) {
//      // printf("curr i %zu\n", i);
//
//      for (size_t j = 0; j < rows; j++) {
//        if (rand() > (int)(err * INT32_MAX)) {
//          frow_pack[j] = ((frow[j] & 1U) << j) | (frow_pack[j]);
//        } else {
//          frow_pack[j] = ((rand() & 1U) << j) | (frow_pack[j]);
//        }
//      }
//
//      if (b % 8 == 7) {
//        // printf("printing encoded %zu \n", (i*8)+b);
//        for (size_t j = 0; j < cols; j++) {
//          // printf("%c", frow_pack[j]);
//          putc(frow_pack[j], stdout);
//        }
//        puts("");
//      }
//      b++;
//    }
//
//    if (b % 8 != 0) {
//      for (size_t j = 0; j < cols; j++) {
//        putc(frow_pack[j], stdout);
//      }
//      puts("");
//    }
//    // for (size_t j = 0; j < cols; j++) {
//    //   printf("%c", frow[j] + 48);
//    // }
//    // puts("");
//    free(frow_pack);
//    free(frow);

    return EXIT_SUCCESS;
  }

  if (argc > 4 && strcmp(argv[1], "txt") == 0) {
    printf("doing TXT\n");
    err = atof(argv[4]);
    char *frow = (char *)malloc(cols * sizeof(char));
    for (size_t i = 0; i < cols; i++) {
      frow[i] = (rand() & 1U) + 48;
      putc(frow[i], stdout);
    }
    puts("");

    for (size_t i = 1; i < rows; i++) {
      for (size_t j = 0; j < cols; j++) {
        if (rand() > (int)(err * UINT32_MAX)) {
          putc(frow[j], stdout);
        } else {
          putc(((rand() & 1U) + 48), stdout);
        }
      }
      puts("");
    }
  } else {

    for (size_t i = 0; i < rows; i++) {
      for (size_t j = 0; j < cols; j++) {
        putc((rand() & 1U) + 48, stdout);
      }
      puts("");
    }
  }

  return EXIT_SUCCESS;
}
