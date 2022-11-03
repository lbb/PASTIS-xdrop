// Created by Saliya Ekanayake on 2019-09-03.

#include "../../inc/align/IPumaAligner.hpp"


#include <iostream>
#include <fstream>


#include <algorithm>
#include <utility>

#include <set>

#include "../../inc/util.hpp"
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

std::tuple<int16_t, int16_t> convertPackedRange(int32_t packed) {
  int16_t begin = packed & 0xffff;
  int16_t end = packed >> 16;
  return {begin, end};
}

struct Batch {
        int startIndex;
        int endIndex;
        int numCmps;
        int totalSize;
        std::vector<int32_t> inputBuffer;
        std::vector<int32_t> resultBuffer;
        std::vector<int32_t> mappingBuffer;
        ipu::batchaffine::Job* job = nullptr;
        ~Batch() {
          if (job != nullptr) {
                  delete job;
          }
        }
};

using Batches = std::vector<Batch>;

Batches createBatches(std::vector<std::string>& A, std::vector<std::string>& B, const int batchCmpLimit, const int batchDataLimit) {
        if (A.size() != B.size()) {
                throw std::runtime_error("Sizes of A and B are not equal");
        }
        Batches batches;
        int numCmps = 0;
        int totalSize = 0;
        int batchBegin = 0;
        for (int i = 0; i < A.size(); ++i) {
                const auto alen = A[i].size();
                const auto blen = B[i].size();

                if ((numCmps + 1 < batchCmpLimit) && (totalSize + alen + blen < batchDataLimit)) {
                        numCmps++;
                        totalSize += alen + blen;
                } else {
                        batches.push_back({.startIndex = batchBegin, .endIndex = batchBegin + numCmps, .numCmps = numCmps, .totalSize = totalSize});

                        batchBegin += numCmps;
                        numCmps = 1;
                        totalSize = alen + blen;
                }
        }
        if (numCmps > 0) {
                batches.push_back({.startIndex = batchBegin, .endIndex = batchBegin + numCmps, .numCmps = numCmps, .totalSize = totalSize});
        }

        // PLOGD.printf("Created %d batches for a total of %d comparisons.", batches.size(), A.size());
        return batches;
}

void IPumaAligner::aln_batch(std::tuple<uint64_t, uint64_t, CommonKmerLight *> *mattuples, uint64_t beg, uint64_t end, uint64_t bl_roffset, uint64_t bl_coffset, const params_t &params) {
  parops->tp->start_timer("sim:align_pre");

  uint64_t npairs = end - beg;
  vector<string> seqs_q(npairs);  // queries - shorter seqs
  vector<string> seqs_r(npairs);  // refs - longer seqs
  uint64_t max_rlen = 0;
  uint64_t max_qlen = 0;

  int numThreads = 1;
#ifdef THREADED
  #pragma omp parallel
  {
    numThreads = omp_get_num_threads();
    omp_set_num_threads(numThreads);
  }
#endif

#pragma omp parallel for reduction(max : max_rlen, max_qlen)
  for (uint64_t i = beg; i < end; ++i) {
    uint64_t lr = std::get<0>(mattuples[i]);
    uint64_t lc = std::get<1>(mattuples[i]);
    assert(lr + bl_roffset < rseqs_.size() && lc + bl_coffset < cseqs_.size());
    string &rseq = rseqs_[lr + bl_roffset];
    string &cseq = cseqs_[lc + bl_coffset];
    if (rseq.size() < cseq.size()) {
      seqs_q[i - beg] = rseq;
      seqs_r[i - beg] = cseq;
    } else {
      seqs_q[i - beg] = cseq;
      seqs_r[i - beg] = rseq;
    }

    max_rlen = max(max_rlen, seqs_r[i - beg].size());
    max_qlen = max(max_qlen, seqs_q[i - beg].size());
  }

  parops->info("prepare batches");
  // if (driver_algo == nullptr) {
  //   init.join();
  // }
  const auto inputBufferSize = driver_algo->algoconfig.getInputBufferSize32b();
  const auto mappingBufferSize = driver_algo->algoconfig.getTotalNumberOfComparisons();
  const auto resultBufferSize = driver_algo->algoconfig.getTotalNumberOfComparisons() * 3;

  const int batchCmpLimit = driver_algo->algoconfig.getTotalNumberOfComparisons() - driver_algo->algoconfig.maxBatches;
  const int batchDataLimit = driver_algo->algoconfig.getTotalBufsize32b() * 4 - driver_algo->algoconfig.bufsize * 100;

  auto bs = createBatches(seqs_q, seqs_r, batchCmpLimit, batchDataLimit);

  auto prepare_batch = [&](Batch& batch) {
          batch.inputBuffer.resize(inputBufferSize);
          batch.mappingBuffer.resize(mappingBufferSize);
          batch.resultBuffer.resize(resultBufferSize);
          driver_algo->prepare_remote(driver_algo->config, driver_algo->algoconfig, 
            {seqs_q.begin()+batch.startIndex, seqs_q.begin()+batch.endIndex},
            {seqs_r.begin()+batch.startIndex, seqs_r.begin()+batch.endIndex},
            &*batch.inputBuffer.begin(),
            &*batch.inputBuffer.end(),
            batch.mappingBuffer.data()
          );
  };

  parops->info("prepare batches parallel loop");
  #pragma omp parallel for
  for (size_t i = 0; i < bs.size(); i++) {
      prepare_batch(bs[i]);
  }
  

  parops->tp->stop_timer("sim:align_pre");

// {
//   ofstream myfile;
//   myfile.open ("As.txt", ios::app);
//   for (auto &&i : seqs_r) {
//     myfile << i << '\n';
//   }
//   myfile.close();
// }
// {
//   ofstream myfile;
//   myfile.open ("Bs.txt", ios::app);
//   for (auto &&i : seqs_q) {
//     myfile << i << '\n';
//   }
//   myfile.close();
// }

  parops->info("start algo");
  parops->tp->start_timer("sim:align");
  for (size_t i = 0; i < bs.size(); i++) {
    bs[i].job = driver_algo->async_submit_prepared_remote_compare(
      &*bs[i].inputBuffer.begin(),
      &*bs[i].inputBuffer.end(),
      &*bs[i].resultBuffer.begin(),
      &*bs[i].resultBuffer.end()
      );
  }
  parops->info("join algo");
  for (size_t i = 0; i < bs.size(); i++) {
    driver_algo->blocking_join_prepared_remote_compare(*bs[i].job);
  }
  // // driver_algo->compare_local(seqs_q, seqs_r);
  // auto [scores, r_range_result, q_range_result] = driver_algo->get_result();

  // parops->tp->start_timer
  parops->tp->stop_timer("sim:align");

  // This seems to be the same block ignore for now

  parops->tp->start_timer("sim:align_post");


  parops->info("iter results");
  for (size_t j = 0; j < bs.size(); j++) {
    vector<int32_t> q_range_result(driver_algo->algoconfig.getTotalNumberOfComparisons());
    vector<int32_t> r_range_result(driver_algo->algoconfig.getTotalNumberOfComparisons());
    vector<int32_t> scores(driver_algo->algoconfig.getTotalNumberOfComparisons());
      driver_algo->transferResults(
        &*bs[j].resultBuffer.begin(),
        &*bs[j].resultBuffer.end(),
        &*bs[j].mappingBuffer.begin(),
        &*bs[j].mappingBuffer.end(),
        &*scores.begin(),
        &*scores.end(),
        &*q_range_result.begin(),
        &*q_range_result.end(),
        &*r_range_result.begin(),
        &*r_range_result.end()
      );

#pragma omp for
    for (uint64_t i = 0; i < bs[j].numCmps; ++i) {
      int len_r = seqs_r[ bs[j].startIndex+i].size();
      int len_q = seqs_q[ bs[j].startIndex+i].size();
      auto [ref_begin, ref_end] = convertPackedRange(r_range_result[i]);
      auto [query_begin, query_end] = convertPackedRange(q_range_result[i]);

      int alen_r = ref_end - ref_begin;
      int alen_q = query_end - query_begin;

      double cov_r = (double)(alen_r) / len_r;
      double cov_q = (double)(alen_q) / len_q;

      if (max(cov_r, cov_q) >= params.aln_cov_thr)  // coverage constraint
      {
        CommonKmerLight *ckl = std::get<2>(mattuples[beg + bs[j].startIndex + i]);
        ckl->score_aln = (float)(scores[i]) /
                         (float)min(len_r, len_q);
      }
    }
  }

  parops->tp->stop_timer("sim:align_post");
}

void IPumaAligner::construct_seqs_bl(std::shared_ptr<DistFastaData> dfd) {
  parops->logger->log("construct_seqs_bl not implemented");
  exit(1);
}

}  // namespace pastis
