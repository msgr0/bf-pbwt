#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  if (argc < 3) {
    fprintf(stderr, "usage: %s ROWS COLS\n", argv[0]);
    return EXIT_FAILURE;
  }

  size_t rows = atoi(argv[1]);
  size_t cols = atoi(argv[2]);

  srand(21);

  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      putc((rand() & 1U) + 48, stdout);
    }
    puts("");
  }

  return EXIT_SUCCESS;
}
