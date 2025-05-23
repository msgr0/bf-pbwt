// vim:ft=c
#include "htslib/synced_bcf_reader.h"
#include "htslib/vcf.h"
#include <string.h>
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
int writebcfwg(void *fd, size_t nr, size_t nc, FILE *of) {
  bcf_srs_t *sr = fd;
  static bcf_hdr_t *hdr = NULL;
  if (!hdr)
    hdr = sr->readers[0].header;
  
  size_t padding = nc % W;
  size_t nw = (nc / W) + (padding > 0); // 1byte windows written to the file
  size_t fcol = nc +  W - padding; // columns actually written to the file, not used.

  
  uw_t *ftcol_pack = (uw_t *)malloc(nr * sizeof(uw_t));

  // iterating over bcf rows, packing them by W
  for (size_t i = 0; i < nc / W; i++) {
    for (size_t w = 0; w < W; w++) {
       if (!bcf_sr_next_line(sr)) {
        return w;
      } 
      bcf1_t *line = bcf_sr_get_line(sr, 0);
      int32_t *gt_arr = NULL;
      int32_t ngt_arr = 0;

      int ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);
      for (size_t r = 0; r < nr /2; r++) {
        int32_t *ptr = gt_arr + r * 2;
        IOBCF_WELLFORMED_CHECK(ptr[0]);
        ftcol_pack[2 * r] |= ((uint64_t)bcf_gt_allele(ptr[0]) << w);
        IOBCF_WELLFORMED_CHECK(ptr[1]);
        ftcol_pack[2 * r + 1] |= ((uint64_t)bcf_gt_allele(ptr[1]) << w);
      }
      free(gt_arr);
    }
    printf("\rsaving window %zu", i);
    fwrite(ftcol_pack, sizeof(uint8_t), nr, of);
    memset(ftcol_pack, 0, sizeof(uint8_t)*nr);
  }
  for (size_t w = 0; w < W; w++) {
    if ((int64_t)(nc % W - w) > 0 ) {
      if (!bcf_sr_next_line(sr)) {
        return (nc/W + 1)* W + w ;
      }
           
      bcf1_t *line = bcf_sr_get_line(sr, 0);
      int32_t *gt_arr = NULL;
      int32_t ngt_arr = 0;

      int ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);
      for (size_t r = 0; r < nr /2; r++) {
        int32_t *ptr = gt_arr + r  *   2;
        IOBCF_WELLFORMED_CHECK(ptr[0]);
        ftcol_pack[2 * r] |= ((uint64_t)bcf_gt_allele(ptr[0]) << w);
        IOBCF_WELLFORMED_CHECK(ptr[1]);
        ftcol_pack[2 * r + 1] |= ((uint64_t)bcf_gt_allele(ptr[1]) << w);
      }
      free(gt_arr);

        
    } 
    else {
      for (size_t r = 0; r < nr; r++) {
        ftcol_pack[r] |= ((uint64_t)0 << w); 
      }
    }
    printf("\rsaving last window %zu", (nc/W + 1));
    fwrite(ftcol_pack, sizeof(uint8_t), nr, of);
  }
  free(ftcol_pack);
  return 0;

}

int main(int argc, char *argv[]) {
  

  bcf_srs_t *sr = bcf_sr_init();
  bcf_sr_add_reader(sr, argv[1]);

  FILE of;
  if (argc >2) {
    of = *fopen(argv[2], "wb");
  } else {
    of = *stdout;
  }

  size_t nrow, ncol = 0;
  bcf_hdr_t *hdr = ((bcf_srs_t *)sr)->readers[0].header;
  nrow = bcf_hdr_nsamples(hdr) * 2;
  while (bcf_sr_next_line(sr)) {
    (ncol)++;
  }

  // bcf_sr_remove_reader(sr, 0);
  
  // we'll write rows, cols as an header
  // if cols % 8 != 0, then we know we have to add
  // a padding to reach the next_int % 8 == 0

  fwrite(&nrow, sizeof(uint32_t), 1, &of);
  fwrite(&ncol, sizeof(uint32_t), 1, &of);
  printf("encoding %zu rows and %zu columns in %zu windows, the last one of %zu bits (then padded at 1byte)\n", nrow, ncol, ncol/W +1, ncol%W);
  // fclose(&of);

  bcf_sr_seek(sr, NULL, 0);
  uint64_t x = writebcfwg(sr, nrow, ncol, &of);
  printf("\n exited at column %llu of %zu", x, ncol);
  fclose(&of);


  bcf_sr_destroy(sr);
  return EXIT_SUCCESS;
}
