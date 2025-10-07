/*!
 * \file MvtxFakeHitRate.h
 * \brief Module to calculate the fake hit rate for the Mvtx
 * \author Tanner Mengel <tmengel@bnl.gov>
 * \version $Version: 1.0.1 $
 * \date $Date: 05/23/2025.
 */

#ifndef MVTXCALIB_MVTXFAKEHITRATE_H
#define MVTXCALIB_MVTXFAKEHITRATE_H

#include <fun4all/SubsysReco.h>

#include <cstdint>
#include <string>
#include <vector>

class PHCompositeNode;
class MvtxHitMap;
class MvtxPixelMask;
class MvtxRawEvtHeader;
class MvtxRawHitContainer;
class TTree;
class TH1;

class MvtxFakeHitRate : public SubsysReco
{
 public:
  MvtxFakeHitRate(const std::string& name = "MvtxFakeHitRate")
    : SubsysReco(name)
  {}
  
  ~MvtxFakeHitRate() override;

  // standard Fun4All functions
  int InitRun(PHCompositeNode* /*topNode*/) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* /*topNode*/) override;

  void SetOutputfile(const std::string& name) { m_outputfile = name; }
  void SetMaxMaskedPixels(int n) { m_max_masked_pixels = n; }
  void SelectLayer(int layer = 0) { m_target_layer = layer; }
  void StartFromCDB(bool start = true) { m_load_from_cdb = start; }

 private:
  // optional output
  std::string m_outputfile{"mvtx_fhr.root"};

  int m_max_masked_pixels{1000};
  bool m_masked_pixels_in_file{false};
  uint64_t m_last_strobe{0};

  int m_target_layer{-1};

  // mask hot pixels
  bool m_load_from_cdb{false};
  MvtxPixelMask* m_hot_pixel_mask{nullptr};

  // hit map
  MvtxHitMap* m_hit_map{nullptr};

  // raw event header and hit container
  MvtxRawEvtHeader* m_mvtx_raw_event_header{nullptr};
  MvtxRawHitContainer* m_mvtx_raw_hit_container{nullptr};

  TTree* m_tree{nullptr};
  int m_num_strobes{0};
  int m_num_masked_pixels{0};
  double m_noise_threshold{0.0};
  std::vector<uint64_t> m_masked_pixels{};

  TTree* m_current_mask{nullptr};
  int m_current_nmasked{0};
  double m_current_threshold{0.0};
  std::vector<uint64_t> m_current_masked_pixels{};

  TH1* m_threshold_vs_nmasked{nullptr};

  int get_nodes(PHCompositeNode* topNode);
  int FillCurrentMaskTree();
  int CalcFHR();

  double calc_threshold(int nhits) const;
};

#endif
