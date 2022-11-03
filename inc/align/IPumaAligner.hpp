#pragma once

#include <string>
#include <tuple>
#include <vector>
#include <thread>

#include "PWAlign.hpp"
#include "../../ipuma-lib/src/driver.hpp"
#include "../../ipuma-lib/src/swatlib/vector.hpp"

namespace pastis {

// uses light representation
class
    IPumaAligner : public PWAlign {
 private:
  std::vector<std::string> rseqs_;
  std::vector<std::string> cseqs_;
  std::tuple<int, int> gaps_;
  // std::thread init;
  ipu::batchaffine::SWAlgorithm *driver_algo = nullptr;
  uint32_t g_batch_sz_;

 public:
  IPumaAligner(int gap_open, int gap_ext, uint32_t bsz = 1e6) : PWAlign(), g_batch_sz_(bsz) {
    this->gaps_ = std::make_tuple<>(gap_open, gap_ext);

    // init_single_ipu(SW_CONFIGURATION, ALGO_CONFIGURATION);

    // init = std::thread([&](){
      swatlib::TickTock t;
      t.tick();
      const ipu::SWConfig SW_CONFIGURATION = {
            .gapInit = gap_open,
            .gapExtend = gap_ext,
            .matchValue = 1,
            .mismatchValue = -1,
            .ambiguityValue = -1,
            .similarity = swatlib::Similarity::blosum62,
            .datatype = swatlib::DataType::aminoAcid,
      };

      const ipu::IPUAlgoConfig ALGOCONFIG = {
          .tilesUsed = 1472,         // number of active vertices
          .maxAB = 1028,        // maximum length of a single comparison
          .maxBatches = 260,  // maximum number of comparisons in a single batch
          .bufsize = 160000,         // total size of buffer for A and B individually
          .vtype = ipu::VertexType::multiasm,
          .fillAlgo = ipu::Algorithm::greedy,
          .forwardOnly = false,  // do not calculate the start position of a match, this should approx 2x performance, as no reverse
                                 // pass is needed
          .useRemoteBuffer = false,
          .transmissionPrograms = 1,  // number of separate transmission programs, use only with remote!
          .ioTiles = 0,
      };
      const auto ipus=1;
      this->driver_algo = new ipu::batchaffine::SWAlgorithm(SW_CONFIGURATION, ALGOCONFIG, 0, 1, ipus);
      t.tock();
      std::cout << "Init took " << t.duration() << " ms" << std::endl;
    // });
  }

  ~IPumaAligner() {
    delete this->driver_algo;
  }

  void
  construct_seqs(std::shared_ptr<DistFastaData> dfd) override;

  void
  construct_seqs_bl(std::shared_ptr<DistFastaData> dfd) override;

  void
  aln_batch(std::tuple<uint64_t, uint64_t, CommonKmerLight *> *mattuples,
            uint64_t beg, uint64_t end,
            uint64_t bl_roffset, uint64_t bl_coffset,
            const params_t &params) override;

  size_t
  rseq_len(uint64_t i) {
    return rseqs_[i].size();
  }

  size_t
  cseq_len(uint64_t i) {
    return cseqs_[i].size();
  }
};

}  // namespace pastis
