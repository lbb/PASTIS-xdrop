#pragma once

#include <plog/Appenders/ColorConsoleAppender.h>
#include <plog/Formatters/TxtFormatter.h>
#include <plog/Initializers/RollingFileInitializer.h>
#include <plog/Log.h>

#include <string>
#include <thread>
#include <tuple>
#include <vector>

#include "../../ipuma-lib/src/driver.hpp"
#include "../../ipuma-lib/src/swatlib/vector.hpp"
#include "PWAlign.hpp"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

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
  int seed_len;
  int seed_cnt;

 public:
  IPumaAligner(int gap_open, int gap_ext, int seed_len, int seed_cnt, uint32_t bsz = 1e6) : PWAlign(), g_batch_sz_(bsz), seed_len(seed_len), seed_cnt(seed_cnt) {
    static plog::ColorConsoleAppender<plog::TxtFormatter> consoleAppender;
    plog::init(plog::verbose, &consoleAppender);

    this->gaps_ = std::make_tuple<>(gap_open, gap_ext);

    // init_single_ipu(SW_CONFIGURATION, ALGO_CONFIGURATION);

    // init = std::thread([&](){
    swatlib::TickTock t;
    t.tick();
     SW_CONFIGURATION = {
        .gapInit = gap_open,
        .gapExtend = gap_ext,
        .matchValue = 1,
        .mismatchValue = -1,
        .ambiguityValue = -1,
        .similarity = swatlib::Similarity::blosum62,
        .datatype = swatlib::DataType::aminoAcid,
        .seedLength = seed_len,
        .xDrop = 49,
    };

     ALGOCONFIG = {
      .numVertices = 1472,
      .maxSequenceLength = 19295,
      .maxComparisonsPerVertex = 600,
      .vertexBufferSize = 160000,
      .vtype = ipu::VertexType::xdroprestrictedseedextend,
      .fillAlgo = ipu::Algorithm::fillFirst,
      .complexityAlgo = ipu::Complexity::xdrop,
      .partitionadd = ipu::PartitionAdd::alternating,
      .partitioningSortComparisons = true,
      .forwardOnly = false,
      .ioTiles = 0,
      .bandPercentageXDrop = 0.42
    };

    PLOGI << " SW_CONFIG: " << json(SW_CONFIGURATION).dump() << std::endl;
    PLOGI << "ALGOCONFIG: " << json(ALGOCONFIG).dump() << std::endl;

    const auto ipus = 1;
    this->driver_algo = new ipu::batchaffine::SWAlgorithm(SW_CONFIGURATION, ALGOCONFIG, 0, 1, ipus);
    t.tock();
    std::cout << "Init took " << t.duration() << " ms" << std::endl;
    // });
  }

  ~IPumaAligner() {
    std::cout << "delete IPU" << std::endl; 
    delete this->driver_algo;
  }

  void
  construct_seqs(std::shared_ptr<DistFastaData> dfd) override;

  void
  construct_seqs_bl(std::shared_ptr<DistFastaData> dfd) override;

  void
  aln_batch(std::tuple<uint64_t, uint64_t, CommonKmerLoc *> *mattuples,
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
  ipu::SWConfig SW_CONFIGURATION;
  ipu::IPUAlgoConfig ALGOCONFIG;
};

}  // namespace pastis
