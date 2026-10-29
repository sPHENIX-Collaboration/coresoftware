// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef ZDCINFO_ZDCRECO_H
#define ZDCINFO_ZDCRECO_H

//===========================================================
/// \author Ejiro Umaka
///
/// Satellite-bunch and run-period correction of the arm energies (applied at the end of process_event when a payload exists):
/// the ZDC pulse is narrower than the 17.6 ns sample spacing, so the energy depends on where the samples fall on the pulse: on the
/// run period (sampling phase of the digitizer clock) and, inside a run, on the collision time (main bunch, satellites at -5 / +5 ns).
/// The arm energy is divided by G(run) x f(period, t0), with t0 the MBD t0 (MbdReco must run before ZdcReco). Events without an
/// MBD time, or with |t0| beyond the satellites, are corrected as main-bunch events (t0 = 0): G(run) x f(t0 = 0).
/// No payload for the run (e.g. simulation): energies unchanged. Payloads are made for FermiExp (FUNCFIT type 2) ZDC tower energies.
/// CDB object "ZDC_SatelliteCorr" (CDBTTree, one per run period): key 0 / 1 (south / north): "gain", "f_m5", "f_p5", "status",
/// "fb00".."fbNN" (t0 bins); single: "nt0bin", "t0bin_lo", "t0bin_w", "nruns_t0", "run_first", "run_last";
/// keys 100 + i: "run", "gain_run_s", "gain_run_n" (per-run scale).
//===========================================================

#include <cdbobjects/CDBTTree.h>

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_const_cgsm.h>

#include <array>
#include <string>  // for string
#include <vector>  // for vector

class PHCompositeNode;
class Zdcinfo;

class ZdcReco : public SubsysReco
{
 public:
  ZdcReco(const std::string &name = "ZdcReco");
  ~ZdcReco() override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  void ResetMe();
  void set_zdc1_cut(const float f) { _zdc1_e = f; }
  void set_zdc2_cut(const float g) { _zdc2_e = g; }
  void set_smd_threshold(const float d) { _smd_e = d; }
  /// satellite-bunch / run-period correction of the arm energies (default on; it does nothing when the run has no payload)
  void set_doSatelliteCorr(const bool b) { m_doSatCorr = b; }
  /// CDB name of the satellite-correction payload (default "ZDC_SatelliteCorr")
  void set_satCalibName(const std::string &name) { m_satCalibName = name; }
  /// testing: read the satellite-correction constants from this local CDBTTree file instead of the CDB (run aborts if the file is missing)
  void set_satCalibFile(const std::string &file) { m_satCalibFile = file; }
  /// use the run's own scale when the payload has it (default true)
  void set_satPerRunGain(const bool b) { m_satPerRunGain = b; }
  /// use the t0-binned bunch factor when the payload has it (default true); false: one value per bunch
  void set_satUseT0Bins(const bool b) { m_satUseT0Bins = b; }

 private:
  void CompSmdPos();
  int LoadSatelliteCorr(PHCompositeNode *topNode);
  float SatFactor(int arm, float t0) const;
  CDBTTree *cdbttree{nullptr};
  // satellite-bunch / run-period correction
  CDBTTree *m_satcdb{nullptr};
  std::string m_satCalibName{"ZDC_SatelliteCorr"};
  std::string m_satCalibFile;
  bool m_doSatCorr{true};
  bool m_satHave{false};
  bool m_satPerRunGain{true};
  bool m_satUseT0Bins{true};
  bool m_satWarnNoMbd{false};
  float m_satEdge1{2.5};  // main bunch |t0| < 2.5 ns
  float m_satEdge2{7.5};  // satellites up to |t0| < 7.5 ns
  std::array<float, 2> m_satGain{1., 1.};
  std::array<float, 2> m_satFm5{1., 1.};
  std::array<float, 2> m_satFp5{1., 1.};
  std::array<int, 2> m_satStatus{1, 1};
  int m_satNt0bin{0};
  float m_satT0lo{-7.5};
  float m_satT0w{0.5};
  std::array<std::vector<float>, 2> m_satFb;
  long m_satNcorr{0};
  long m_satNnoT0{0};
  Zdcinfo *m_zdcinfo{nullptr};
  std::string m_Detector{"ZDC"};
  std::string m_fieldname;
  std::string m_calibName;
  bool m_overrideCalibName{false};
  bool m_overrideFieldName{false};
  bool smd_north_fired{false};
  bool smd_south_fired{false};
  std::vector<float> vsmdadc;
  std::vector<float> vsmdtime;
  std::vector<float> vzdcadc;
  std::vector<float> vzdctime;
  float smd_adc[32]{0.0};
  float smd_pos[4]{0.0};
  float _radius_south{0.};
  float _radius_north{0.};
  const double _t{17.623};                                // convert to ns
  const double _c{GSL_CONST_CGSM_SPEED_OF_LIGHT / 1.e9};  // speed of light in cm/ns
  float _z_vertex{0.};
  float _zdc1_e{65.0};
  float _zdc2_e{20.0};
  float _smd_e{5.0};
  float _sumS{0.};
  float _sumN{0.};
  float _sumSden{0.};
  float _sumNden{0.};
  float _sumSt{0.};
  float _sumNt{0.};
  int _nhor{0};
  int _nver{0};
  int _shor{0};
  int _sver{0};
};

#endif  // ZDCRECO_H
