#include "bseq.h"
#include "kalloc.h"
#include "khash.h"
#include "kthread.h"
#include "kvec.h"
#include "mmpriv.h"
#include "sdust.h"
#include <assert.h>
#include <errno.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

struct mm_tbuf_s {
  void *km;
  int rep_len, frag_gap; // updated per read.
  double timers[MM_N_THR_TIMERS];
}; // per thread

mm_tbuf_t *mm_tbuf_init(void) {
  mm_tbuf_t *b;
  b = (mm_tbuf_t *)calloc(1, sizeof(mm_tbuf_t));
  if (!(mm_dbg_flag & 1))
    b->km = km_init();
  return b;
}

void mm_tbuf_destroy(mm_tbuf_t *b) {
  if (b == 0)
    return;
  km_destroy(b->km);
  free(b);
}

void *mm_tbuf_get_km(mm_tbuf_t *b) { return b->km; }

static int mm_dust_minier(void *km, int n, mm128_t *a, int l_seq,
                          const char *seq, int sdust_thres) {
  int n_dreg, j, k, u = 0;
  const uint64_t *dreg;
  sdust_buf_t *sdb;
  if (sdust_thres <= 0)
    return n;
  sdb = sdust_buf_init(km);
  dreg = sdust_core((const uint8_t *)seq, l_seq, sdust_thres, 64, &n_dreg, sdb);
  for (j = k = 0; j < n;
       ++j) { // squeeze out minimizers that significantly overlap with LCRs
    int32_t qpos = (uint32_t)a[j].y >> 1, span = a[j].x & 0xff;
    int32_t s = qpos - (span - 1), e = s + span;
    while (u < n_dreg && (int32_t)dreg[u] <= s)
      ++u;
    if (u < n_dreg && (int32_t)(dreg[u] >> 32) < e) {
      int v, l = 0;
      for (v = u; v < n_dreg && (int32_t)(dreg[v] >> 32) < e;
           ++v) { // iterate over LCRs overlapping this minimizer
        int ss = s > (int32_t)(dreg[v] >> 32) ? s : dreg[v] >> 32;
        int ee = e < (int32_t)dreg[v] ? e : (uint32_t)dreg[v];
        l += ee - ss;
      }
      if (l <= span >> 1)
        a[k++] = a[j]; // keep the minimizer if less than half of it falls in
                       // masked region
    } else
      a[k++] = a[j];
  }
  sdust_buf_destroy(sdb);
  return k; // the new size
}

static void collect_minimizers(void *km, const mm_mapopt_t *opt,
                               const mm_idx_t *mi, int n_segs, const int *qlens,
                               const char **seqs, mm128_v *mv) {
  int i, n, sum = 0;
  mv->n = 0;
  for (i = n = 0; i < n_segs; ++i) {
    size_t j;
    mm_sketch(km, seqs[i], qlens[i], mi->w, mi->k, i, mi->flag & MM_I_HPC, mv);
    for (j = n; j < mv->n; ++j)
      mv->a[j].y += sum << 1;
    if (opt->sdust_thres > 0) // mask low-complexity minimizers
      mv->n = n + mm_dust_minier(km, mv->n - n, mv->a + n, qlens[i], seqs[i],
                                 opt->sdust_thres);
    sum += qlens[i], n = mv->n;
  }
}

#include "ksort.h"
#define heap_lt(a, b) ((a).x > (b).x)
KSORT_INIT(heap, mm128_t, heap_lt)

static inline int skip_seed(int flag, uint64_t r, const mm_seed_t *q,
                            const char *qname, int qlen, const mm_idx_t *mi,
                            int *is_self) {
  *is_self = 0;
  if (qname && (flag & (MM_F_NO_DIAG | MM_F_NO_DUAL))) {
    const mm_idx_seq_t *s = &mi->seq[r >> 32];
    int cmp;
    cmp = strcmp(qname, s->name);
    if ((flag & MM_F_NO_DIAG) && cmp == 0 && (int)s->len == qlen) {
      if ((uint32_t)r >> 1 == (q->q_pos >> 1))
        return 1; // avoid the diagnonal anchors
      if ((r & 1) == (q->q_pos & 1))
        *is_self =
            1; // this flag is used to avoid spurious extension on self chain
    }
    if ((flag & MM_F_NO_DUAL) && cmp > 0) // all-vs-all mode: map once
      return 1;
  }
  if (flag & (MM_F_FOR_ONLY | MM_F_REV_ONLY)) {
    if ((r & 1) == (q->q_pos & 1)) { // forward strand
      if (flag & MM_F_REV_ONLY)
        return 1;
    } else {
      if (flag & MM_F_FOR_ONLY)
        return 1;
    }
  }
  return 0;
}

static mm128_t *collect_seed_hits_heap(void *km, const mm_mapopt_t *opt,
                                       int max_occ, const mm_idx_t *mi,
                                       const char *qname, const mm128_v *mv,
                                       int qlen, int64_t *n_a, int *rep_len,
                                       int *n_mini_pos, uint64_t **mini_pos) {
  int i, n_m, heap_size = 0;
  int64_t j, n_for = 0, n_rev = 0;
  mm_seed_t *m;
  mm128_t *a, *heap;

  m = mm_collect_matches(km, &n_m, qlen, max_occ, opt->max_max_occ,
                         opt->occ_dist, mi, mv, n_a, rep_len, n_mini_pos,
                         mini_pos);

  heap = (mm128_t *)kmalloc(km, n_m * sizeof(mm128_t));
  a = (mm128_t *)kmalloc(km, *n_a * sizeof(mm128_t));

  for (i = 0, heap_size = 0; i < n_m; ++i) {
    if (m[i].n > 0) {
      heap[heap_size].x = m[i].cr[0];
      heap[heap_size].y = (uint64_t)i << 32;
      ++heap_size;
    }
  }
  ks_heapmake_heap(heap_size, heap);
  while (heap_size > 0) {
    mm_seed_t *q = &m[heap->y >> 32];
    mm128_t *p;
    uint64_t r = heap->x;
    int32_t is_self, rpos = (uint32_t)r >> 1;
    if (!skip_seed(opt->flag, r, q, qname, qlen, mi, &is_self)) {
      if ((r & 1) == (q->q_pos & 1)) { // forward strand
        p = &a[n_for++];
        p->x = (r & 0xffffffff00000000ULL) | rpos;
        p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;
      } else { // reverse strand
        p = &a[(*n_a) - (++n_rev)];
        p->x = 1ULL << 63 | (r & 0xffffffff00000000ULL) | rpos;
        p->y = (uint64_t)q->q_span << 32 |
               (qlen - ((q->q_pos >> 1) + 1 - q->q_span) - 1);
      }
      p->y |= (uint64_t)q->seg_id << MM_SEED_SEG_SHIFT;
      if (q->is_tandem)
        p->y |= MM_SEED_TANDEM;
      if (is_self)
        p->y |= MM_SEED_SELF;
    }
    // update the heap
    if ((uint32_t)heap->y < q->n - 1) {
      ++heap[0].y;
      heap[0].x = m[heap[0].y >> 32].cr[(uint32_t)heap[0].y];
    } else {
      heap[0] = heap[heap_size - 1];
      --heap_size;
    }
    ks_heapdown_heap(0, heap_size, heap);
  }
  kfree(km, m);
  kfree(km, heap);

  // reverse anchors on the reverse strand, as they are in the descending order
  for (j = 0; j < n_rev >> 1; ++j) {
    mm128_t t = a[(*n_a) - 1 - j];
    a[(*n_a) - 1 - j] = a[(*n_a) - (n_rev - j)];
    a[(*n_a) - (n_rev - j)] = t;
  }
  if (*n_a > n_for + n_rev) {
    memmove(a + n_for, a + (*n_a) - n_rev, n_rev * sizeof(mm128_t));
    *n_a = n_for + n_rev;
  }
  return a;
}

static mm128_t *collect_seed_hits(void *km, const mm_mapopt_t *opt, int max_occ,
                                  const mm_idx_t *mi, const char *qname,
                                  const mm128_v *mv, int qlen, int64_t *n_a,
                                  int *rep_len, int *n_mini_pos,
                                  uint64_t **mini_pos) {
  int i, n_m;
  mm_seed_t *m;
  mm128_t *a;
  m = mm_collect_matches(km, &n_m, qlen, max_occ, opt->max_max_occ,
                         opt->occ_dist, mi, mv, n_a, rep_len, n_mini_pos,
                         mini_pos);
  a = (mm128_t *)kmalloc(km, *n_a * sizeof(mm128_t));
  for (i = 0, *n_a = 0; i < n_m; ++i) {
    mm_seed_t *q = &m[i];
    const uint64_t *r = q->cr;
    uint32_t k;
    for (k = 0; k < q->n; ++k) {
      int32_t is_self, rpos = (uint32_t)r[k] >> 1;
      mm128_t *p;
      if (skip_seed(opt->flag, r[k], q, qname, qlen, mi, &is_self))
        continue;
      p = &a[(*n_a)++];
      if ((r[k] & 1) == (q->q_pos & 1)) { // forward strand
        p->x = (r[k] & 0xffffffff00000000ULL) | rpos;
        p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;
      } else if (!(opt->flag & MM_F_QSTRAND)) { // reverse strand and not in the
                                                // query-strand mode
        p->x = 1ULL << 63 | (r[k] & 0xffffffff00000000ULL) | rpos;
        p->y = (uint64_t)q->q_span << 32 |
               (qlen - ((q->q_pos >> 1) + 1 - q->q_span) - 1);
      } else { // reverse strand; query-strand
        int32_t len = mi->seq[r[k] >> 32].len;
        p->x = 1ULL << 63 | (r[k] & 0xffffffff00000000ULL) |
               (len - (rpos + 1 - q->q_span) -
                1); // coordinate only accurate for non-HPC seeds
        p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;
      }
      p->y |= (uint64_t)q->seg_id << MM_SEED_SEG_SHIFT;
      if (q->is_tandem)
        p->y |= MM_SEED_TANDEM;
      if (is_self)
        p->y |= MM_SEED_SELF;
    }
  }
  kfree(km, m);
  radix_sort_128x(a, a + (*n_a));
  return a;
}

void mm_map_frag_bt(const mm_idx_t *mi, int n_segs, const int *qlens,
                    const char **seqs, int *n_regs, uint64_t **_u, mm_tbuf_t *b,
                    const mm_mapopt_t *opt, const char *qname) {
  int i, j, rep_len, qlen_sum, n_regs0, n_mini_pos;
  int max_chain_gap_qry, max_chain_gap_ref,
      is_splice = !!(opt->flag & MM_F_SPLICE), is_sr = !!(opt->flag & MM_F_SR);
  uint32_t hash;
  int64_t n_a;
  uint64_t *u;
  uint64_t *mini_pos;
  mm128_t *a;
  mm128_v mv = {0, 0, 0};
  mm_reg1_t *regs0;
  km_stat_t kmst;
  float chn_pen_gap, chn_pen_skip;
  double *timers = b->timers;
  double t1 = realtime();

  for (i = 0, qlen_sum = 0; i < n_segs; ++i)
    qlen_sum += qlens[i], n_regs[i] = 0;

  if (qlen_sum == 0 || n_segs <= 0 || n_segs > MM_MAX_SEG)
    return;
  if (opt->max_qlen > 0 && qlen_sum > opt->max_qlen)
    return;

  // ====================== SEEDING =====================
  hash = qname && !(opt->flag & MM_F_NO_HASH_NAME) ? __ac_X31_hash_string(qname)
                                                   : 0;
  hash ^= __ac_Wang_hash(qlen_sum) + __ac_Wang_hash(opt->seed);
  hash = __ac_Wang_hash(hash);

  collect_minimizers(b->km, opt, mi, n_segs, qlens, seqs, &mv);
  if (opt->q_occ_frac > 0.0f)
    mm_seed_mz_flt(b->km, &mv, opt->mid_occ, opt->q_occ_frac);
  if (opt->flag & MM_F_HEAP_SORT)
    a = collect_seed_hits_heap(b->km, opt, opt->mid_occ, mi, qname, &mv,
                               qlen_sum, &n_a, &rep_len, &n_mini_pos,
                               &mini_pos);
  else
    a = collect_seed_hits(b->km, opt, opt->mid_occ, mi, qname, &mv, qlen_sum,
                          &n_a, &rep_len, &n_mini_pos, &mini_pos);

  if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
    fprintf(stderr, "RS\t%d\n", rep_len);
    for (i = 0; i < n_a; ++i)
      fprintf(stderr, "SD\t%s\t%d\t%c\t%d\t%d\t%d\n",
              mi->seq[a[i].x << 1 >> 33].name, (int32_t)a[i].x,
              "+-"[a[i].x >> 63], (int32_t)a[i].y,
              (int32_t)(a[i].y >> 32 & 0xff),
              i == 0 ? 0
                     : ((int32_t)a[i].y - (int32_t)a[i - 1].y) -
                           ((int32_t)a[i].x - (int32_t)a[i - 1].x));
  }
  timers[MM_TIME_SEED] += realtime() - t1;

  // ==================== CHAINING ===========================
  t1 = realtime();
  // set max chaining gap on the query and the reference sequence
  if (is_sr)
    max_chain_gap_qry = qlen_sum > opt->max_gap ? qlen_sum : opt->max_gap;
  else
    max_chain_gap_qry = opt->max_gap;
  if (opt->max_gap_ref > 0) {
    max_chain_gap_ref =
        opt->max_gap_ref; // always honor mm_mapopt_t::max_gap_ref if set
  } else if (opt->max_frag_len > 0) {
    max_chain_gap_ref = opt->max_frag_len - qlen_sum;
    if (max_chain_gap_ref < opt->max_gap)
      max_chain_gap_ref = opt->max_gap;
  } else
    max_chain_gap_ref = opt->max_gap;

  chn_pen_gap = opt->chain_gap_scale * 0.01 * mi->k;
  chn_pen_skip = opt->chain_skip_scale * 0.01 * mi->k;
  if (opt->flag & MM_F_RMQ) {
    a = mg_lchain_rmq(opt->max_gap, opt->rmq_inner_dist, opt->bw,
                      opt->max_chain_skip, opt->rmq_size_cap, opt->min_cnt,
                      opt->min_chain_score, chn_pen_gap, chn_pen_skip, n_a, a,
                      &n_regs0, &u, b->km);
  } else {
    a = mg_lchain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw,
                     opt->max_chain_skip, opt->max_chain_iter, opt->min_cnt,
                     opt->min_chain_score, chn_pen_gap, chn_pen_skip, is_splice,
                     n_segs, n_a, a, &n_regs0, &u, b->km);
  }

  if (opt->bw_long > opt->bw &&
      (opt->flag & (MM_F_SPLICE | MM_F_SR | MM_F_NO_LJOIN)) == 0 &&
      n_segs == 1 && n_regs0 > 1) { // re-chain/long-join for long sequences
    int32_t st = (int32_t)a[0].y, en = (int32_t)a[(int32_t)u[0] - 1].y;
    if (qlen_sum - (en - st) > opt->rmq_rescue_size ||
        en - st > qlen_sum * opt->rmq_rescue_ratio) {
      int32_t i;
      for (i = 0, n_a = 0; i < n_regs0; ++i)
        n_a += (int32_t)u[i];
      kfree(b->km, u);
      radix_sort_128x(a, a + n_a);
      a = mg_lchain_rmq(opt->max_gap, opt->rmq_inner_dist, opt->bw_long,
                        opt->max_chain_skip, opt->rmq_size_cap, opt->min_cnt,
                        opt->min_chain_score, chn_pen_gap, chn_pen_skip, n_a, a,
                        &n_regs0, &u, b->km);
    }
  } else if (opt->max_occ > opt->mid_occ && rep_len > 0 &&
             !(opt->flag & MM_F_RMQ)) { // re-chain, mostly for short reads
    int rechain = 0;
    if (n_regs0 > 0) { // test if the best chain has all the segments
      int n_chained_segs = 1, max = 0, max_i = -1, max_off = -1, off = 0;
      for (i = 0; i < n_regs0; ++i) { // find the best chain
        if (max < (int)(u[i] >> 32))
          max = u[i] >> 32, max_i = i, max_off = off;
        off += (uint32_t)u[i];
      }
      for (i = 1; i < (int32_t)u[max_i];
           ++i) // count the number of segments in the best chain
        if ((a[max_off + i].y & MM_SEED_SEG_MASK) !=
            (a[max_off + i - 1].y & MM_SEED_SEG_MASK))
          ++n_chained_segs;
      if (n_chained_segs < n_segs)
        rechain = 1;
    } else
      rechain = 1;
    if (rechain) { // redo chaining with a higher max_occ threshold
      kfree(b->km, a);
      // kfree(b->km, u);
      // free(u);
      kfree(b->km, mini_pos);
      if (opt->flag & MM_F_HEAP_SORT)
        a = collect_seed_hits_heap(b->km, opt, opt->max_occ, mi, qname, &mv,
                                   qlen_sum, &n_a, &rep_len, &n_mini_pos,
                                   &mini_pos);
      else
        a = collect_seed_hits(b->km, opt, opt->max_occ, mi, qname, &mv,
                              qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
      a = mg_lchain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw,
                       opt->max_chain_skip, opt->max_chain_iter, opt->min_cnt,
                       opt->min_chain_score, chn_pen_gap, chn_pen_skip,
                       is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);
    }
  }
  b->frag_gap = max_chain_gap_ref;
  b->rep_len = rep_len;
  timers[MM_TIME_CHAIN] += realtime() - t1;

  kfree(b->km, mv.a);
  kfree(b->km, a);
  kfree(b->km, mini_pos);

  *_u = malloc(n_regs0 * 8);
  memcpy(*_u, u, n_regs0 * 8);
  kfree(b->km, u);

  if (b->km) {
    km_stat(b->km, &kmst);
    if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
      fprintf(stderr, "QM\t%s\t%d\tcap=%ld,nCore=%ld,largest=%ld\n", qname,
              qlen_sum, kmst.capacity, kmst.n_cores, kmst.largest);
    assert(kmst.n_blocks == kmst.n_cores); // otherwise, there is a memory leak
    if (kmst.largest > 1U << 28 ||
        (opt->cap_kalloc > 0 && kmst.capacity > opt->cap_kalloc)) {
      if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
        fprintf(stderr, "[W::%s] reset thread-local memory after read %s\n",
                __func__, qname);
      km_destroy(b->km);
      b->km = km_init();
    }
  }
}

void mm_map_bt(const mm_idx_t *mi, uint64_t **u, int qlen, const char *seq,
               int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt,
               const char *qname) {
  mm_map_frag_bt(mi, 1, &qlen, &seq, n_regs, u, b, opt, qname);
}

/**************************
 * Multi-threaded mapping *
 **************************/

typedef struct {
  int n_processed, n_threads, n_fp;
  int64_t mini_batch_size;
  const mm_mapopt_t *opt;
  mm_bseq_file_t **fp;
  const mm_idx_t *mi;
  kstring_t str;

  int n_parts;
  uint32_t *rid_shift;
  FILE *fp_split, **fp_parts;
} pipeline_t;

typedef struct {
  const pipeline_t *p;
  int n_seq, n_frag;
  mm_bseq1_t *seq;
  int *n_reg, *seg_off, *n_seg, *rep_len, *frag_gap;
  mm_reg1_t **reg;
  mm_tbuf_t **buf;
} step_t;

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

// consolidate timers from worker threads
void mm_consolidate_timers(step_t *s, pipeline_t *p) {
  // TODO: disabled for sysbio submission
  return;
  mm_time_seed_min = 0;
  mm_time_chain_min = 0;
  mm_time_align_min = 0;
  mm_time_seed_max = 0;
  mm_time_chain_max = 0;
  mm_time_align_max = 0;
  mm_time_seed_avg = 0;
  mm_time_chain_avg = 0;
  mm_time_align_avg = 0;
  for (int i = 0; i < p->n_threads; ++i) {
    mm_time_seed_min = MIN(mm_time_seed_min, s->buf[i]->timers[MM_TIME_SEED]);
    mm_time_chain_min =
        MIN(mm_time_chain_min, s->buf[i]->timers[MM_TIME_CHAIN]);
    mm_time_align_min =
        MIN(mm_time_align_min, s->buf[i]->timers[MM_TIME_ALIGN]);
    mm_time_seed_max = MAX(mm_time_seed_max, s->buf[i]->timers[MM_TIME_SEED]);
    mm_time_chain_max =
        MAX(mm_time_chain_max, s->buf[i]->timers[MM_TIME_CHAIN]);
    mm_time_align_max =
        MAX(mm_time_align_max, s->buf[i]->timers[MM_TIME_ALIGN]);
    mm_time_seed_avg += s->buf[i]->timers[MM_TIME_SEED];
    mm_time_chain_avg += s->buf[i]->timers[MM_TIME_CHAIN];
    mm_time_align_avg += s->buf[i]->timers[MM_TIME_ALIGN];
    mm_time_seed_sum += s->buf[i]->timers[MM_TIME_SEED];
    mm_time_chain_sum += s->buf[i]->timers[MM_TIME_CHAIN];
    mm_time_align_sum += s->buf[i]->timers[MM_TIME_ALIGN];
  }
  mm_time_seed_avg /= p->n_threads;
  mm_time_chain_avg /= p->n_threads;
  mm_time_align_avg /= p->n_threads;

  fprintf(stderr, "----------------------------------------------------\n");
  fprintf(stderr, "              Min (sec)  Max (sec)  Avg (sec)  \n");
  fprintf(stderr, "----------------------------------------------------\n");
  fprintf(stderr, "Seed    = %11.3f %11.3f %11.3f\n", mm_time_seed_min,
          mm_time_seed_max, mm_time_seed_avg);
  fprintf(stderr, "Chain   = %11.3f %11.3f %11.3f\n", mm_time_chain_min,
          mm_time_chain_max, mm_time_chain_avg);
  fprintf(stderr, "Align   = %11.3f %11.3f %11.3f\n", mm_time_align_min,
          mm_time_align_max, mm_time_align_avg);
  fprintf(stderr, "----------------------------------------------------\n");
  fprintf(stderr, "Avg (seed + chain + align) per thread = %.3f secs\n",
          (mm_time_seed_avg + mm_time_chain_avg + mm_time_align_avg));
  fprintf(stderr,
          "Total (seed + chain + align) (all batches) for %d thread(s) = %.3f "
          "secs\n",
          p->n_threads,
          (mm_time_seed_sum + mm_time_chain_sum + mm_time_align_sum));
}

static mm_bseq_file_t **open_bseqs(int n, const char **fn) {
  mm_bseq_file_t **fp;
  int i, j;
  fp = (mm_bseq_file_t **)calloc(n, sizeof(mm_bseq_file_t *));
  for (i = 0; i < n; ++i) {
    if ((fp[i] = mm_bseq_open(fn[i])) == 0) {
      if (mm_verbose >= 1)
        fprintf(stderr, "ERROR: failed to open file '%s': %s\n", fn[i],
                strerror(errno));
      for (j = 0; j < i; ++j)
        mm_bseq_close(fp[j]);
      free(fp);
      return 0;
    }
  }
  return fp;
}
