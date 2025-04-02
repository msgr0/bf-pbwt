#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int main(int argc, char *argv[]) {
  if (argc < 4) {
    fprintf(stderr, "usage: %s ROWS COLS ERR\n", argv[0]);
    return EXIT_FAILURE;
  }

  size_t rows = atoi(argv[1]);
  size_t cols = atoi(argv[2]);
  float err = atof(argv[3]);

  srand(21);

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

  return EXIT_SUCCESS;
}
