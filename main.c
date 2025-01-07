// To compile:
//   gcc -g -O2 example.c libminimap2.a -lz

#include "kseq.h"
#include "minimap.h"
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <zlib.h>
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[]) {
  mm_idxopt_t iopt;
  mm_mapopt_t mopt;
  int n_threads = 3;

  mm_verbose = 2; // disable message output to stderr
  mm_set_opt(0, &iopt, &mopt);
  mopt.flag |= MM_F_CIGAR; // perform alignment

  if (argc < 3) {
    fprintf(stderr, "Usage: mm2-backrack <target.fa> <query.fa>\n");
    return 1;
  }

  // open query file for reading; you may use your favorite FASTA/Q parser
  gzFile f = gzopen(argv[2], "r");
  assert(f);
  kseq_t *ks = kseq_init(f);

  // open index reader
  mm_idx_reader_t *r = mm_idx_reader_open(argv[1], &iopt, 0);
  mm_idx_t *mi;
  while ((mi = mm_idx_reader_read(r, n_threads)) !=
         0) { // traverse each part of the index
    mm_mapopt_update(&mopt,
                     mi); // this sets the maximum minimizer occurrence; TODO:
                          // set a better default in mm_mapopt_init()!
    mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading,
                                      // allocate one tbuf for each thread
    gzrewind(f);
    kseq_rewind(ks);
    while (kseq_read(ks) >=
           0) { // each kseq_read() call reads one query sequence
      uint64_t *u;
      int j, i, n_reg;
      // get output from mg_backtrack for comparison
      mm_map_bt(mi, &u, ks->seq.l, ks->seq.s, &n_reg, tbuf, &mopt, 0);
      printf("%lu", u[0]);
      free(u);
    }
    mm_tbuf_destroy(tbuf);
    mm_idx_destroy(mi);
  }
  mm_idx_reader_close(r); // close the index reader
  kseq_destroy(ks);       // close the query file
  gzclose(f);
  return 0;
}
