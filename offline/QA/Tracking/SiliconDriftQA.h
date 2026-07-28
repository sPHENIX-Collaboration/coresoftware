#ifndef QA_TRACKING_SILICONDRIFTQA_H
#define QA_TRACKING_SILICONDRIFTQA_H

/*
 * QA version of SiliconDriftEvaluator (B. Sayki, LANL).
 *
 * Monitors the TPC drift velocity calibration by comparing the z position of
 * the TPC seed and the silicon seed at the beam line, following the standard
 * sPHENIX QA module conventions
 * Claude code was used in formatting and debugging of this module
 * v_new = v_in / (1 + slope)
 * t0    = (offset_pos - offset_neg) / (2 v_new)
 */

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class TH1;
class TH2;

class SiliconDriftQA : public SubsysReco
{
 public:
  explicit SiliconDriftQA(const std::string& name = "SiliconDriftQA");

  ~SiliconDriftQA() override = default;

  //! run initialization: create and register histograms
  int InitRun(PHCompositeNode* topNode) override;

  //! event processing: fill histograms
  int process_event(PHCompositeNode* topNode) override;

  //! end of processing: fit accumulated distributions, fill summary histogram
  int End(PHCompositeNode* topNode) override;

  //! track map name
  void set_trackmapname(const std::string& value) { m_trackmapname = value; }

  //! initial drift velocity (cm/ns); starting point for the fit and used in the crossing correction
  void set_drift_velocity(double value) { m_drift_velocity = value; }

  //! bunch-crossing interval in ns (default: 106.65237 ns)
  void set_crossing_interval(double value) { m_crossing_interval = value; }

  //! minimum pT cut on tracks (GeV)
  void set_min_pt(double value) { m_min_pt = value; }

  //! minimum number of TPC clusters required
  void set_min_nclusters_tpc(unsigned int value) { m_min_nclusters_tpc = value; }

  //! minimum number of MVTX clusters required
  void set_min_nclusters_mvtx(unsigned int value) { m_min_nclusters_mvtx = value; }

  //! minimum number of INTT clusters required
  void set_min_nclusters_intt(unsigned int value) { m_min_nclusters_intt = value; }

  //! maximum abs(eta) of TPC seed accepted
  void set_max_eta(double value) { m_max_eta = value; }

  //! half-range of the z_si histogram axis (cm)
  void set_max_z(double value) { m_max_z = value; }

  //! half-range of the dz histogram axis (cm)
  void set_max_dz(double value) { m_max_dz = value; }

  //! minimum entries per z slice required by FitSlicesY
  void set_min_slice_entries(int value) { m_min_slice_entries = value; }

 private:
  void createHistos();
  std::string getHistoPrefix() const;

  //!@name histograms (owned by the QA histogram manager)
  //@{

  //! z_si vs dz, one per eta bin (0: eta<0, 1: eta>=0)
  TH2* h_zsi_dz[2]{nullptr, nullptr};

  //! crossing-corrected dz = z_tpc_corr - z_si (cm)
  TH1* h_dz{nullptr};

  //! number of accepted tracks per event
  TH1* h_ntracks{nullptr};

  //! drift velocity fit summary, filled in End()
  TH1* h_driftSummary{nullptr};

  //@}

  //! track map name
  std::string m_trackmapname{"SvtxTrackMap"};

  //! initial drift velocity (cm/ns)
  double m_drift_velocity{0.00749};

  //! bunch-crossing interval (ns)
  double m_crossing_interval{106.65237};

  //!@name track selection cuts
  //@{
  double m_min_pt{0.5};
  unsigned int m_min_nclusters_tpc{20};
  unsigned int m_min_nclusters_mvtx{3};
  unsigned int m_min_nclusters_intt{2};
  double m_max_eta{0.9};
  //@}

  //! histogram ranges
  double m_max_z{20.0};
  double m_max_dz{10.0};

  //! minimum entries per z slice for FitSlicesY
  int m_min_slice_entries{10};
};

#endif  // QA_TRACKING_SILICONDRIFTQA_H