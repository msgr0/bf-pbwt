#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int main(int argc, char *argv[]) {
  if (argc < 3) {
    fprintf(stderr, "usage: %s ROWS COLS [err]\n", argv[0]);
    return EXIT_FAILURE;
  }
  size_t rows = atoi(argv[1]);
  size_t cols = atoi(argv[2]);

  float err = 0.0;
  srand(21);

  if (argc > 3) {
    err = atof(argv[3]);
    char* frow = (char*) malloc(cols * sizeof(char));
    for (size_t i = 0; i < cols; i++) {
      frow[i] = (rand() & 1U) + 48; 
      putc(frow[i], stdout);   	  
    }
    puts("");


    for (size_t i = 1; i < rows; i++) {
      for (size_t j = 0; j < cols; j++) {
        if (rand() > (int)( err * UINT32_MAX)) {
          putc(frow[j], stdout);
        }
        else {
          putc(((rand() & 1U) + 48), stdout);
        }
      }
      puts("");
    }         
  }
  else {

    for (size_t i = 0; i < rows; i++) {
      for (size_t j = 0; j < cols; j++) {
        putc((rand() & 1U) + 48, stdout);
      }
      puts("");
    }

  }

  return EXIT_SUCCESS;
}
