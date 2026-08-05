// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef JETBASE_JETCALIB_H
#define JETBASE_JETCALIB_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class TF1;
class TFile;
class TH2;

// Jet energy scale calibration:
//   1) first pass: calibrated (truth-equivalent) pT read from a 2D map versus
//      (reco pT, EMCal energy fraction), by clamped bilinear interpolation;
//   2) residual correction: multiplicative scale factor versus (signed z-vertex,
//      fractional eta position f within the z-corrected acceptance window);
//      events without a reconstructed z-vertex instead use a 1D correction
//      function versus f, with f computed at z = 0.

class JetCalib : public SubsysReco
{
 public:
  explicit JetCalib(const std::string &name = "JetCalib");
  ~JetCalib() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int CreateNodeTree(PHCompositeNode *topNode);

  // Setters.
  void set_InputNode(const std::string &inputNode) { m_inputNode = inputNode; }
  void set_OutputNode(const std::string &outputNode) { m_outputNode = outputNode; }
  void set_JetRadius(float radius) { jet_radius = radius; }
  void set_ApplyResidualCalib(bool apply) { m_ApplyResidualCalib = apply; }
  void set_CalibFile(const std::string &path) { m_calibFileOverride = path; }  // local file override (default: CDB)

 private:
  // Functions.
  static std::string fetchCalibDir(const char *calibType);
  // z-corrected eta acceptance window (intersection of EMCal/IHCal/OHCal reach,
  // padded inward by the jet radius); must match the derivation exactly.
  static void getEtaAcceptance(float zvrtx, float radius, float &minlimit, float &maxlimit);
  // Bilinear interpolation with both coordinates clamped to the bin-center range.
  static float interpolateClamped(TH2 *h2, double x, double y);
  float getFirstPassPt(float jetPt, float emFrac) const;
  float getResidualScale(float zvrtx, bool hasVertex, float eta) const;

  // Input.
  std::string m_inputNode{"AntiKt_Tower_r04"};
  std::string m_outputNode{"AntiKt_Tower_r04_Calib"};
  std::string m_calibFileOverride;      // if set, read this file instead of the CDB payload
  float jet_radius{0.4};                // Jet radius.
  bool m_ApplyResidualCalib{true};      // Apply the residual (z-vertex, eta) correction.

  // Calibration objects (for the configured radius).
  TFile *m_JetCalibFile{nullptr};
  TH2 *m_JesCalibMap{nullptr};   // h2d_jes_calib_r0R: truth-equivalent pT vs (reco pT, EMCal fraction)
  TH2 *m_ZetaCorrMap{nullptr};   // h2_zeta_corr_r0R: scale factor vs (signed z-vertex, f)
  TF1 *m_NozCorrFunc{nullptr};   // f_zeta_noz_corr_r0R: scale factor vs f (f at z = 0), no-vertex events
  bool m_warnedNoEmFrac{false};
};

#endif  // JETBASE_JETCALIB_H
