// vim:ft=c
#include "htslib/synced_bcf_reader.h"
#include "htslib/vcf.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#define BFGETCOLI_BUF_SIZE 128
#define BFGETCOLWR_BUF_SIZE 2

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

#define W 8

typedef uint8_t uw_t;
int writebcfwg(void *fd, void *hdr, size_t nr, size_t nc, FILE *of) {
  bcf_srs_t *sr = fd;
  static size_t c = 0;
  
  uw_t *ftcol_pack = (uw_t *)malloc(nr * sizeof(uw_t));
  // iterating over bcf rows, packing them by W
  int i = 0;
  int cont = 0;
  for (i = 0;!cont; i++) { // iterating over full windows
    printf("\rwriting window: %zu", c);
    for (size_t w = 0; w < W; w++) { // reading W bcf lines, each line appending to ftcol_pack;
      if (!bcf_sr_next_line(sr)) {
        cont++;
      }
      else {

      bcf1_t *line = bcf_sr_get_line(sr, 0);

      int32_t *gt_arr = NULL;
      int32_t ngt_arr = 0;

      int ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);
      for (size_t r = 0; r < nr / 2; r++) {
        int32_t *ptr = gt_arr + r * 2;
        IOBCF_WELLFORMED_CHECK(ptr[0]);
        ftcol_pack[2 * r] |= ((uint64_t)bcf_gt_allele(ptr[0]) << w);
        IOBCF_WELLFORMED_CHECK(ptr[1]);
        ftcol_pack[2 * r + 1] |= ((uint8_t)bcf_gt_allele(ptr[1]) << w);
      }
      free(gt_arr);
      }
    }
    fwrite(ftcol_pack, sizeof(uint8_t), nr, of);
    c++;
    memset(ftcol_pack, 0, sizeof(*ftcol_pack) * nr);
  }
  
  printf("\n written a total of %zu windows", c);
  free(ftcol_pack);
  return W * (i-1) + (W-cont);
}

int main(int argc, char *argv[]) {

  bcf_srs_t *sr = bcf_sr_init();
  bcf_sr_add_reader(sr, argv[1]);

  FILE of;
  if (argc > 2) {
    of = *fopen(argv[2], "wb");
  } else {
    of = *stdout;
  }

  size_t nrow, ncol = 0;
  bcf_hdr_t *hdr = ((bcf_srs_t *)sr)->readers[0].header;
  nrow = bcf_hdr_nsamples(hdr) * 2;

  fwrite(&nrow, sizeof(uint32_t), 1, &of);
  fwrite(&ncol, sizeof(uint32_t), 1, &of);

  ncol = writebcfwg(sr, hdr, nrow, ncol, &of);
  fseek(&of,0,0);
  fwrite(&nrow, sizeof(uint32_t), 1, &of);
  fwrite(&ncol, sizeof(uint32_t), 1, &of);
  // problem with bcf_sr when seeking back to 0,
  // solved by counting the column during writebcfwg
  // room for improvement.
  uint32_t x = (ncol / W) + (ncol % W > 0); 
  printf("\n exited; written %zu rows times %u windows (bytes) representing %zu columns\n", nrow, x, ncol);
  fseek(&of, 0, 0);
  fclose(&of);

  bcf_sr_destroy(sr);
  return EXIT_SUCCESS;
}
