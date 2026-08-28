#include "htslib/synced_bcf_reader.h"
#include "htslib/vcf.h"
#include <stdio.h>
#include <stdlib.h>

// Simple program to extract first N records from a BCF file
int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "usage: %s <input.bcf> <output.bcf> <max_records>\n", argv[0]);
        return EXIT_FAILURE;
    }

    const char *infile = argv[1];
    const char *outfile = argv[2];
    int max_records = atoi(argv[3]);

    htsFile *in = hts_open(infile, "r");
    if (!in) {
        fprintf(stderr, "Failed to open input file\n");
        return EXIT_FAILURE;
    }

    bcf_hdr_t *hdr = bcf_hdr_read(in);
    if (!hdr) {
        fprintf(stderr, "Failed to read BCF header\n");
        hts_close(in);
        return EXIT_FAILURE;
    }

    htsFile *out = hts_open(outfile, "wb");
    if (!out) {
        fprintf(stderr, "Failed to open output file\n");
        bcf_hdr_destroy(hdr);
        hts_close(in);
        return EXIT_FAILURE;
    }

    if (bcf_hdr_write(out, hdr) < 0) {
        fprintf(stderr, "Failed to write header\n");
        hts_close(out);
        bcf_hdr_destroy(hdr);
        hts_close(in);
        return EXIT_FAILURE;
    }

    bcf1_t *rec = bcf_init();
    int nrecs = 0;

    while (bcf_read(in, hdr, rec) >= 0 && nrecs < max_records) {
        if (bcf_write(out, hdr, rec) < 0) {
            fprintf(stderr, "Failed to write record\n");
            break;
        }
        nrecs++;
    }

    fprintf(stderr, "Extracted %d records\n", nrecs);

    bcf_destroy(rec);
    hts_close(out);
    bcf_hdr_destroy(hdr);
    hts_close(in);

    return EXIT_SUCCESS;
}
