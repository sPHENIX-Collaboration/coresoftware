// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef NNSCALEMAP_H
#define NNSCALEMAP_H

#include <fun4all/SubsysReco.h>

#include <array>
#include <string>
#include <vector>

class PHCompositeNode;
class SvtxTrackMap;

// Runs after track reconstruction. Reads the tracks in an existing
// SvtxTrackMap, applies a per-track pT scale factor kappa(q, pT, eta, phi)
// taken from a trilinearly-interpolated lookup table (as produced by the
// MomentumScaleStudies ml_momentum_calibration_reso*.py training scripts,
// kappa_lookup.csv), and writes the corrected tracks into a new
// SvtxTrackMap node. The input track map and node tree are left untouched.
//
// By default (setUseCDB(true)) the lookup CSV is located from the
// Calibration Database, keyed by run number/timestamp via CDBInterface,
// under the configured domain name -- with the fixed local path from
// setKappaLookupFile(), if any, used as the fallback if the CDB has no
// payload for this run. Call setUseCDB(false) to always use the fixed
// local path instead. CDB resolution happens once per run in InitRun(),
// so a job spanning multiple runs picks up the correct map for each one
// automatically.
//
// The CSV needs at least columns named "q", "pT", "eta", "phi", "kappa"
// (located by header name, so extra diagnostic columns -- e.g. "curv",
// "kappa_curv", "eps", "delta" from a curvature-space training run -- are
// ignored and column order doesn't matter). "kappa" must already be the
// pT multiplier: pT_corrected = kappa * pT_reco.
class NNScaleMap : public SubsysReco
{
 public:
  explicit NNScaleMap(const std::string& name = "NNScaleMap");
  ~NNScaleMap() override = default;

  int Init(PHCompositeNode* topNode) override;
  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  void Print(const std::string& what = "ALL") const override;

  //! local path to a kappa_lookup.csv. Used directly if setUseCDB(false); used
  //! as the CDBInterface::getUrl() fallback filename if setUseCDB(true) (the
  //! default) and the CDB has no payload for this run.
  void setKappaLookupFile(const std::string& path) { m_lookupFile = path; }

  //! look up the lookup CSV from the Calibration Database per-run instead of
  //! (or as an override of) the fixed local path (default true)
  void setUseCDB(bool doIt) { m_useCDB = doIt; }

  //! CDB domain/payload name to query when setUseCDB(true) (default "TRACKING_MOMENTUM_SCALE")
  void setCdbDomainName(const std::string& domain) { m_cdbDomainName = domain; }

  //! node name of the input track map to read (default "SvtxTrackMap")
  void setInputTrackMapName(const std::string& name) { m_inputTrackMapName = name; }

  //! node name of the corrected track map this module creates (default "SvtxTrackMapMomentumScaleCorrected")
  void setOutputTrackMapName(const std::string& name) { m_outputTrackMapName = name; }

  //! also rescale the momentum block of the track covariance matrix (default true)
  void setScaleCovariance(bool doIt) { m_scaleCovariance = doIt; }

 private:
  int CreateNodes(PHCompositeNode* topNode);
  int GetNodes(PHCompositeNode* topNode);
  std::string ResolveLookupFile() const;
  bool LoadLookupTable(const std::string& path);
  float GetKappa(int charge, float pt, float eta, float phi) const;

  std::string m_lookupFile;
  bool m_useCDB = true;
  std::string m_cdbDomainName = "TRACKING_MOMENTUM_SCALE";
  std::string m_loadedFile;  // path actually loaded, so InitRun can skip a needless reparse

  std::string m_inputTrackMapName = "SvtxTrackMap";
  std::string m_outputTrackMapName = "SvtxTrackMapMomentumScaleCorrected";
  bool m_scaleCovariance = true;

  // lookup grid axes, read directly from the CSV (not assumed a priori;
  // a curvature-space training run yields a grid that is non-uniform in
  // pT -- dense below ~3 GeV, sparse above)
  std::vector<double> m_ptEdges;
  std::vector<double> m_etaEdges;
  std::vector<double> m_phiEdges;

  // kappa[chargeIdx][(ipt*nEta + ieta)*nPhi + iphi], chargeIdx: 0 -> q<0, 1 -> q>0
  std::array<std::vector<float>, 2> m_kappaGrid;

  SvtxTrackMap* m_inputTrackMap = nullptr;
  SvtxTrackMap* m_outputTrackMap = nullptr;

  unsigned int m_nEvents = 0;
  unsigned long m_nTracksSeen = 0;
  unsigned long m_nTracksCorrected = 0;
  unsigned long m_nTracksSkippedZeroPt = 0;
};

#endif  // NNSCALEMAP_H
