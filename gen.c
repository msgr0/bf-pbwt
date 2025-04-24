#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BITS 8

int main(int argc, char *argv[]) {
  if (argc < 5) {
    fprintf(stderr, "usage: %s enc|txt ROWS COLS ERR [file] \n", argv[0]);
    return EXIT_FAILURE;
  }
  size_t rows = atoi(argv[2]);
  size_t cols = atoi(argv[3]);

  float err = 0.0;
  srand(42);

  if (argc >= 5 && strcmp(argv[1], "enc") == 0) {
    FILE of;
    if (argc == 6) {
      of = *fopen(argv[5], "wb");
    }
    else {
      of = *stdout;
    }

    // printf("doing ENC\n");
    err = atof(argv[4]);

    size_t padding = cols % 8;
    size_t nw = (cols / 8) + (padding > 0);
    size_t fcols = cols + 8 - padding;

    // we'll write rows, cols as an header
    // if cols % 8 != 0, then we know we have to add 
    // a padding to reach the next_int % 8 == 0

    fwrite(&rows, sizeof(uint32_t), 1, &of);
    fwrite(&cols, sizeof(uint32_t), 1, &of);

    char *frow = (char *)malloc(fcols * sizeof(char));

    for (size_t i = 0; i < cols; i++) {
      frow[i] = (rand() & 1U);
      //printf("%c", frow[i] + 48);
    }
    for (size_t i = 0; i < 8 - padding; i++) {
      frow[i + cols] = 0;
    }
    //printf("b-nw:%zu-b\n", nw);

    char *ftcol_pack = (char *)malloc(rows * sizeof(char));
    size_t b = 0;
    for (size_t j = 0; j < nw; j++) { // per ogni w-sized colonna
      for (size_t i = 0; i < rows; i++) { // per ogni riga
        for (size_t b = 8; b >0; b--) { // per ogni bit della w-sized colonna
          ftcol_pack[i] = (ftcol_pack[i] << 1) | ( rand() > (int)(err * INT32_MAX) ? (frow[(j)*8 + b-1] & 1U) : (rand() & 1U) );
        }
      }
      fwrite(ftcol_pack, sizeof(char), rows, &of);
    }

    free(frow);
    free(ftcol_pack);
    fclose(&of);

    return EXIT_SUCCESS;
  }

  if (argc > 4 && strcmp(argv[1], "txt") == 0) {
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
