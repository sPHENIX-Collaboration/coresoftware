#ifndef G4TPC_PHG4TPCPADPLANEREADOUT_H
#define G4TPC_PHG4TPCPADPLANEREADOUT_H

#include "PHG4TpcPadPlane.h"
#include "TpcClusterBuilder.h"

#include <g4main/PHG4HitContainer.h>

#include <gsl/gsl_rng.h>

#include <array>
#include <climits>
#include <cmath>
#include <string>  // for string
#include <vector>
#include <map>

typedef std::map<TrkrDefs::hitsetkey, std::vector<TrkrDefs::hitkey>> hitMaskTpc;

class PHCompositeNode;
class PHG4TpcGeomContainer;
class PHG4TpcGeom;
class TH2;
class TF1;
class TNtuple;
class TrkrHitSetContainer;
class TrkrHitTruthAssoc;

class PHG4TpcPadPlaneReadout : public PHG4TpcPadPlane
{
 public:
  PHG4TpcPadPlaneReadout(const std::string &name = "PHG4TpcPadPlaneReadout");

  ~PHG4TpcPadPlaneReadout() override;

  int InitRun(PHCompositeNode *topNode) override;

  void UseGain(const int flagToUseGain);
  void SetUseModuleGainWeights(const int flag) { m_use_module_gain_weights = flag; }
  void SetModuleGainWeightsFileName(const std::string &name) { m_tpc_module_gain_weights_file = name; }
  void ReadGain();
  void SetUsePolyaGEMGain(const int flagPolya) { m_usePolya = flagPolya; }
  void SetUseLangauGEMGain(const int flagLangau) { m_useLangau = flagLangau; }
  void SetLangauParsFileName(const std::string &name) { m_tpc_langau_pars_file = name; }

  // otherwise warning of inconsistent overload since only one MapToPadPlane methow is overridden
  using PHG4TpcPadPlane::MapToPadPlane;

  void MapToPadPlane(TpcClusterBuilder &tpc_truth_clusterer, TrkrHitSetContainer *single_hitsetcontainer, TrkrHitSetContainer *hitsetcontainer, TrkrHitTruthAssoc * /*hittruthassoc*/, const double x_gem, const double y_gem, const double t_gem, const unsigned int side, PHG4HitContainer::ConstIterator hiter, TNtuple * /*ntpad*/, TNtuple * /*nthit*/) override;

  void SetDefaultParameters() override;
  void UpdateInternalParameters() override;

  void SetDeadChannelMapName(const std::string& dcmap) 
  {
    m_maskDeadChannels = true;
    m_deadChannelMapName = dcmap;
  }
  void SetHotChannelMapName(const std::string& hmap) 
  {
    m_maskHotChannels = true;
    m_hotChannelMapName = hmap;
  }

 private:
  //  void populate_rectangular_phibins(const unsigned int layernum, const double phi, const double cloud_sig_rp, std::vector<int> &pad_phibin, std::vector<double> &pad_phibin_share);
  void populate_zigzag_phibins(const unsigned int side, const unsigned int layernum, const double phi, const double cloud_sig_rp, std::vector<int> &phibin_pad, std::vector<double> &phibin_pad_share);

  void sampaTimeDistribution(double tzero,  std::vector<int> &adc_tbin, std::vector<double> &adc_tbin_share);
  double sampaShapingResponseFunction(double tzero, double t) const;
  
  double check_phi(const unsigned int side, const double phi, const double radius);

  void makeChannelMask(hitMaskTpc& aMask, const std::string& dbName, const std::string& totalChannelsToMask);

  PHG4TpcGeomContainer *GeomContainer = nullptr;
  PHG4TpcGeom *LayerGeom = nullptr;

  double neffelectrons_threshold {std::numeric_limits<double>::quiet_NaN()};

  std::array<double, 3> MinRadius{};
  std::array<double, 3> MaxRadius{};

  static constexpr int NSides {2};
  static constexpr int NSectors {12};
  static const int NRSectors {3};

  double sigmaT {std::numeric_limits<double>::quiet_NaN()};
  std::array<double, 2> sigmaL{};
  double phi_bin_width{};

  int NTBins {std::numeric_limits<int>::max()};
  int m_NHits {0};
  // Using Gain maps is turned off by default
  int m_flagToUseGain {0};

  // Optionally apply a module-by-module weight to the GEM gain
  // Weights are input from a file for all 72 TPC modules
  bool m_use_module_gain_weights {false};
  std::string m_tpc_module_gain_weights_file;

  // gaussian sampling
  static constexpr double _nsigmas {5};

  double Ts {55.0}; // SAMPA v5 peaking time

  double averageGEMGain {std::numeric_limits<double>::quiet_NaN()};
  double polyaTheta {std::numeric_limits<double>::quiet_NaN()};

  std::array<std::vector<double>, NSides> sector_min_Phi;
  std::array<std::vector<double>, NSides> sector_max_Phi;

  // return random distribution of number of electrons after amplification of GEM for each initial ionizing electron
  double getSingleEGEMAmplification();
  double getSingleEGEMAmplification(double weight);
  static double getSingleEGEMAmplification(TF1 *f);
  bool m_usePolya {false};

  bool m_useLangau {false};
  std::string m_tpc_langau_pars_file;

  gsl_rng *RandomGenerator {nullptr};

  std::array<TH2 *, 2> h_gain{nullptr};

  double m_module_gain_weight[2][3][12] {
      {{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
       {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
       {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
      {{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
       {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
       {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}}};

  TF1 *flangau[2][3][12] {{{nullptr}}};

  hitMaskTpc m_deadChannelMap;
  hitMaskTpc m_hotChannelMap; 

  bool m_maskDeadChannels {false};
  bool m_maskHotChannels {false};
  std::string m_deadChannelMapName; 
  std::string m_hotChannelMapName; 
};

#endif
