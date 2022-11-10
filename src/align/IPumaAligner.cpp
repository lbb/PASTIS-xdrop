#include "../../inc/align/IPumaAligner.hpp"

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <set>
#include <utility>
#include <string_view>

#include "../../inc/util.hpp"
#include "../../ipuma-lib/src/swatlib/matrix.h"

// #include "ipumultidriver.hpp"
// #include "batch.hpp"

using std::max;
using std::min;
using std::string;
using std::to_string;
using std::vector;

extern shared_ptr<pastis::ParallelOps> parops;

namespace pastis {

void IPumaAligner::construct_seqs(std::shared_ptr<DistFastaData> dfd) {
  uint64_t rseq_cnt = dfd->get_nrseqs();
  uint64_t cseq_cnt = dfd->get_ncseqs();

  rseqs_.resize(rseq_cnt);
  cseqs_.resize(cseq_cnt);

  int cur_nbr = 0;
  uint64_t cur_rseq = 0;
  uint64_t cur_cseq = 0;
  for (auto &nbr : dfd->my_nbrs) {
    uint64_t nbr_seqs_count = nbr.nbr_seq_end_idx - nbr.nbr_seq_start_idx + 1;
    FastaData *curfd = dfd->fd;
    uint64_t sidx = nbr.nbr_seq_start_idx;
    if (nbr.nbr_rank != parops->g_rank)  // not myself
    {
      curfd = dfd->recv_fds_[cur_nbr];
      sidx = 0;
      ++cur_nbr;
    }

    vector<string> &cur_seqs = (nbr.rc_flag == 1 ? rseqs_ : cseqs_);
    uint64_t seq_beg = (nbr.rc_flag == 1 ? cur_rseq : cur_cseq);

#pragma omp parallel
    {
      ushort len;
      uint64_t start_offset, end_offset_inclusive;

#pragma omp for
      for (uint64_t i = 0; i < nbr_seqs_count; ++i) {
        char *buff = curfd->get_sequence(sidx + i, len, start_offset, end_offset_inclusive);
        cur_seqs[seq_beg + i] = std::move(string(buff + start_offset, len));
      }
    }

    if (nbr.rc_flag == 1)
      cur_rseq += nbr_seqs_count;
    else
      cur_cseq += nbr_seqs_count;
  }

  // diag cell row and col seqs are the same
  if (parops->grid->GetRankInProcRow() == parops->grid->GetRankInProcCol()) {
    assert(cseqs_.empty());
    cseqs_.assign(rseqs_.begin(), rseqs_.end());
  }
}

void IPumaAligner::aln_batch(
    std::tuple<uint64_t, uint64_t, CommonKmerLoc *> *mattuples,
    uint64_t beg,
    uint64_t end,
    uint64_t bl_roffset,
    uint64_t bl_coffset,
    const params_t &params) {
  parops->tp->start_timer("sim:align_pre");

  uint64_t npairs = end - beg;
  vector<string> seqs_V(npairs);    // queries - shorter seqs
  vector<string> seqs_H(npairs);    // refs - longer seqs

  // std::vector<ipu::MultiComparison> comparisons(npairs);

  int numThreads = 1;
#ifdef THREADED
#pragma omp parallel
  {
    numThreads = omp_get_num_threads();
    omp_set_num_threads(numThreads);
  }
#endif

  PLOGI.printf("seed_cnt = %d", seed_cnt);
  PLOGI.printf("seed_len = %d", seed_len);
  PLOGI.printf("npairs = %d", npairs);
  if (seed_cnt != NSEEDS) {
    std::cerr << "unsupported seed coutn set " << std::endl;
  }

  PLOGW << "Seed Length" << seed_len;

  std::vector<uint64_t> colstartrange;
  int lastcol = -1;
  size_t totalDim = 0;
  for (auto i = 0; i < npairs; i++) {
    // uint64_t lr = std::get<0>(mattuples[beg + i]);
    uint64_t lr = std::get<0>(mattuples[beg + i]);
    uint64_t lc = std::get<1>(mattuples[beg + i]);
    totalDim = max(totalDim, max(lc, lr));
    if (lc != lastcol) {
      colstartrange.push_back(i);
      lastcol = lc; 
    }
    // colstartrange.push_back(i);
  }
  totalDim += 1;



  std::vector<std::string_view> sequences_fake(npairs * 2);

  std::vector<std::string_view> sequences(totalDim * 2);
  for (int i = 0; i < totalDim; i++) {
    sequences[i] = rseqs_[i + bl_roffset];
    sequences[totalDim+i] = cseqs_[i + bl_coffset];
  }


  // Prepare sequences

 auto getRseqI = [&](int i) {
   return i;
 };

 auto getCseqI = [&](int i) {
  return totalDim + i;
 };


// const auto maxbatch_size = ALGOCONFIG.maxComparisonsPerVertex/2 - 1;
const auto maxbatch_size = 1;

PLOGF << "maxbatch size: " << maxbatch_size;
std::vector<int> hist(maxbatch_size, 0);
std::vector<ipu::MultiComparison> comparisons;
auto cnt_comparisons = 0;

swatlib::TickTock multibatch_timer;
multibatch_timer.tick();
// #pragma omp parallel for 
  // for (uint64_t i = 0; i < npairs; ++i) {
  for (uint64_t colix = 0; colix < colstartrange.size(); colix++) {
    
    auto col = std::get<1>(mattuples[beg + colstartrange[colix]]);
    std::vector<ipu::Comparison> tmpcmps;
    auto seqstore = 0;

    for (uint64_t i = colstartrange[colix]; i < npairs; ++i) {
      // Start here for ocls 
      CommonKmerLoc *ckl = std::get<2>(mattuples[beg + i]);
      uint64_t lr = std::get<0>(mattuples[beg + i]);
      uint64_t lc = std::get<1>(mattuples[beg + i]);
      if (lc != col) {break;}
      // assert(lr + bl_roffset < rseqs_.size() && lc + bl_coffset < cseqs_.size());


      int seeds[2][2];
      for (int cnt = 0; cnt < seed_cnt; ++cnt) {
        ushort l_row_seed_start_offset = (cnt == 0) ? ckl->first.first : ckl->second.first;
        ushort l_col_seed_start_offset = (cnt == 0) ? ckl->first.second : ckl->second.second;
        seeds[cnt][0] = l_col_seed_start_offset;
        seeds[cnt][1] = l_row_seed_start_offset;
      }

      if (sequences[getRseqI(lr)].size() + seqstore + sequences[getCseqI(lc)].size() >= ALGOCONFIG.vertexBufferSize || tmpcmps.size() >= maxbatch_size) {
          comparisons.emplace_back(tmpcmps, seed_len);
          hist[tmpcmps.size() - 1] += 1;
          tmpcmps.resize(0);
          seqstore = 0;
      }

      seqstore += sequences[getCseqI(lc)].size();
      cnt_comparisons += NSEEDS;

      sequences_fake[i * 2] = sequences[getCseqI(lc)];
      sequences_fake[i * 2 + 1] = sequences[getRseqI(lr)];
      tmpcmps.push_back({
        i,
        i * 2,
        sequences[getCseqI(lc)].size(), 
        i * 2 + 1,
        sequences[getRseqI(lr)].size(), 
        {{{seeds[0][0], seeds[0][1]}, {seeds[1][0], seeds[1][1]}}},
        0
      });

      // tmpcmps.push_back({
        // i,
        // getCseqI(lc),
        // sequences[getCseqI(lc)].size(), 
        // getRseqI(lr),
        // sequences[getRseqI(lr)].size(), 
        // {{{seeds[0][0], seeds[0][1]}, {seeds[1][0], seeds[1][1]}}},
        // 0
      // });
    }
    if (tmpcmps.size() != 0) {
      comparisons.emplace_back(tmpcmps, seed_len);
      hist[tmpcmps.size() - 1] += 1;
      tmpcmps.resize(0);
    }
  }
  multibatch_timer.tock();
  PLOGI << "MultiBatching took [ms]: " << multibatch_timer.duration();
  PLOGI << "Sigular comparison count: " << cnt_comparisons;
  PLOGD << "We got MultiComparison: " << comparisons.size();
  for (size_t i = 0; i < hist.size(); i++) {
    PLOGD.printf("Bucket %5d has %7d entries.", i + 1, hist[i]);
  }
  parops->info("prepare batches");
  parops->info("prepare batches");

  // std::vector<ipu::partition::BatchMapping> mappings;
  auto &cmps = comparisons;
  auto &seqs = sequences_fake;
  // {
  //   ofstream myfile;
  //   myfile.open ("cmps.txt");
  //   myfile << json(cmps).dump() << std::endl;
  //   myfile.close();
  // }
  // {
  //   ofstream myfile;
  //   myfile.open ("seqs.txt");
  //   myfile << json(seqs).dump() << std::endl;
  //   myfile.close();
  // }

  /////////////////////////////////////////////////////////////////////////////
  // RUN BATCHES
  /////////////////////////////////////////////////////////////////////////////


  // std::vector<ipu::Batch> batches;
  swatlib::TickTock mapping_timer;
  mapping_timer.tick();
  PLOGI << "Use multicomparisons for creating batches";
  auto mappings = ipu::partition::mapBatches(ALGOCONFIG, seqs, cmps);
  // batches = ipu::create_batches(seqs, mcmps, config.ipuconfig, config.swconfig);
  mapping_timer.tock();
  PLOGI << "Mapping took [ms]: " << mapping_timer.duration();

  std::vector<ipu::batchaffine::Job *> jobs(mappings.size());
  std::vector<ipu::Batch> batches(mappings.size());
  bool lazyGeneration = true;

  if (mappings.size() > 750 || lazyGeneration) {
    lazyGeneration = true;
    PLOGE << "Generate batches lazy, as batches " << mappings.size();
  }

  if (!lazyGeneration) {
    swatlib::TickTock batch_timer;
    batch_timer.tick();
#pragma omp parallel for
    for (int i = 0; i < mappings.size(); ++i) {
      batches[i] = ipu::create_batch(mappings[i], seqs, ALGOCONFIG, SW_CONFIGURATION);
    }
    batch_timer.tock();
    PLOGI << "Batch copies took [ms]: " << batch_timer.duration();
  }

  vector<int> scorepairs(npairs, 0);
  swatlib::TickTock align_timer;
  int64_t gcells = 0;
  int progress = 0;

  align_timer.tick();
#pragma omp parallel for
  for (int i = 0; i < mappings.size(); ++i) {
    ipu::Batch batch;
    if (lazyGeneration) {
      PLOGW << "lazy generated batch";
      batch = ipu::create_batch(mappings[i], seqs, ALGOCONFIG, SW_CONFIGURATION);
    } else {
      batch = std::move(batches[i]);
    }
    PLOGW << batch.inputs.size();
    jobs[i] = driver_algo->async_submit(&batch);
    PLOGW << "submitted job";
    driver_algo->blocking_join(*jobs[i]);

#pragma omp critical
    {
      progress++;
      gcells += batch.cellCount;
      PLOGI << "Received batch " << progress << " / " << jobs.size();
    }

    auto result = batch.get_result();
    delete jobs[i];
    // for (int ii = 0; ii < batch.origin_comparison_index.size(); ++ii) {
    //   auto [orig_i, orig_seed] = ipu::unpackOriginIndex(batch.origin_comparison_index[ii]);
    //   if (orig_i >= 0) {
    //     int lpartScore = result.a_range_result[ii];
    //     int rpartScore = result.b_range_result[ii];
    //     // PLOGF << orig_i << ":" << batch.origin_comparison_index[i] << " " << lpartScore << " " << rpartScore;
    //     scorepairs[orig_i] = std::max(lpartScore + rpartScore + SW_CONFIGURATION.seedLength, scorepairs[orig_i]);
    //   }
    // }
  }

  align_timer.tock();
  PLOGI << "GCells: " << gcells;
  PLOGI << "Align time took [ms]: " << align_timer.duration();

  /////////////////////////////////////////////////////////////////////////////
  // END RUN BATCHES
  /////////////////////////////////////////////////////////////////////////////

  // parops->tp->start_timer
  parops->tp->stop_timer("sim:align");
  parops->tp->start_timer("sim:align_post");

#pragma omp for
  for (uint64_t i = 0; i < npairs; i++) {
    // int len_r_left = seqs_H[i * 2].size();
    // int len_r_right = seqs_H[i * 2 + 1].size();
    // int len_q_left = seqs_V[i * 2].size();
    // int len_q_right = seqs_V[i * 2 + 1].size();

    // int id_left = (double)scorepairs[i * 2] / min(len_r_left, len_q_left);
    // int id_right = (double)scorepairs[i * 2 + 1] / min(len_r_right, len_q_right);
    CommonKmerLoc *ckl = std::get<2>(mattuples[beg + i]);

    ckl->score_aln = 0;  // (float)(id_left + id_right) / 2.0;
  }
  parops->tp->stop_timer("sim:align_post");
}

void IPumaAligner::construct_seqs_bl(std::shared_ptr<DistFastaData> dfd) {
  parops->logger->log("construct_seqs_bl not implemented");
  // exit(1);
}

}  // namespace pastis
