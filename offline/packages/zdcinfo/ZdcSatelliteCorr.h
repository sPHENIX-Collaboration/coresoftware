// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef ZDCINFO_ZDCSATELLITECORR_H
#define ZDCINFO_ZDCSATELLITECORR_H

//===========================================================
/// \brief Correction of the ZDC arm energies for the bunch (satellite) and the run period.
///
/// The ZDC pulse is narrower than the 17.8 ns sample spacing, so the reconstructed energy depends on where the samples fall on the
/// pulse: on the run period (sampling phase of the digitizer clock) and, inside a run, on the collision time of the bunch
/// (main bunch and the satellites at -5 / +5 ns). This module divides the Zdcinfo arm energy by
///      G(run) x f(period, t0)
/// G     : main-bunch energy scale of the run (fallback: of the period),
/// f     : energy at time t0 over the main-bunch energy at the same centrality, in bins of MBD t0 (fallback: one value per bunch),
/// t0    : MBD t0 from MbdReco; f was measured with t0 centered on each run's main peak, so the MBD timing calibration must center t0.
/// Run after MbdReco and ZdcReco. Events without an MBD time, outside the binned t0 range, or in a run range flagged unusable
/// (status 0) are left unchanged. The constants are valid for the FermiExp (FUNCFIT type 2) ZDC tower energies.
///
/// CDB object (CDBTTree, one payload per run period), default name "ZDC_SatelliteCorr":
///   key 0 (south), 1 (north): "gain", "f_m5", "f_p5", "status", "phi0" (diagnostic), "fb00" ... "fbNN"
///   single values: "nt0bin", "t0bin_lo", "t0bin_w", "nruns_t0", "run_first", "run_last"
///   keys 100 + i (i < nruns_t0): "run", "gain_run_s", "gain_run_n"
//===========================================================

#include <fun4all/SubsysReco.h>

#include <array>
#include <string>
#include <vector>

class CDBTTree;
class PHCompositeNode;

class ZdcSatelliteCorr : public SubsysReco
{
 public:
  explicit ZdcSatelliteCorr(const std::string &name = "ZdcSatelliteCorr");
  ~ZdcSatelliteCorr() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  /// CDB payload name
  void set_calibName(const std::string &name) { m_calibName = name; }
  /// testing: read the constants from this local CDBTTree file (any path, need not be in the CDB) instead of the CDB lookup.
  /// The file is used for every run; the run aborts if it does not exist, and a warning is printed if the run is outside its run_first..run_last.
  void set_calibFile(const std::string &file) { m_calibFile = file; }
  /// same as set_calibFile (name used by CaloTowerCalib)
  void set_directURL(const std::string &url) { m_calibFile = url; }
  /// true: write the corrected energies into "Zdcinfo" (default); false: into a new node "ZdcinfoCorr"
  void set_overwrite(const bool b) { m_overwrite = b; }
  /// use the run's own main-bunch scale when the payload has it (default true)
  void set_perRunGain(const bool b) { m_perRunGain = b; }
  /// use the t0-binned f when the payload has it (default true); false: one value per bunch
  void set_useT0Bins(const bool b) { m_useT0Bins = b; }

 private:
  int createNodes(PHCompositeNode *topNode);
  float factor(int arm, float t0) const;

  CDBTTree *m_cdb{nullptr};
  std::string m_calibName{"ZDC_SatelliteCorr"};
  std::string m_calibFile;
  bool m_overwrite{true};
  bool m_perRunGain{true};
  bool m_useT0Bins{true};
  bool m_have{false};

  float m_edge1{2.5};  // main bunch |t0| < m_edge1 [ns]
  float m_edge2{7.5};  // satellites up to |t0| < m_edge2 [ns]
  std::array<float, 2> m_gain{1., 1.};
  std::array<float, 2> m_fm5{1., 1.};
  std::array<float, 2> m_fp5{1., 1.};
  std::array<int, 2> m_status{1, 1};
  int m_nt0bin{0};
  float m_t0lo{-7.5};
  float m_t0w{0.5};
  std::array<std::vector<float>, 2> m_fb;

  unsigned long m_nevt{0};
  unsigned long m_ncorr{0};
};

#endif
