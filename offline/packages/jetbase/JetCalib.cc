#include "JetCalib.h"

#include "Jet.h"
#include "JetContainerv1.h"

<<<<<<< Updated upstream
=======
#include <cdbobjects/CDBTF.h>  // for CDBTF1 (legacy method)

>>>>>>> Stashed changes
#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TAxis.h>
#include <TF1.h>
#include <TFile.h>
#include <TH2.h>
#include <TString.h>  // for Form

#include <algorithm>  // for std::max, std::min
#include <cmath>
#include <cstdlib>   // for exit
#include <iostream>  // for operator<<, basic_ostream

JetCalib::JetCalib(const std::string &name)
  : SubsysReco(name)
{
  if (Verbosity() > 0)
  {
    std::cout << "JetCalib::JetCalib(const std::string &name) Calling ctor." << std::endl;
  }
}

JetCalib::~JetCalib()
{
  if (m_JetCalibFile)
  {
    m_JetCalibFile->Close();
    delete m_JetCalibFile;
  }
  if (Verbosity() > 0)
  {
    std::cout << "JetCalib::~JetCalib() : Calling dtor." << std::endl;
  }
}

int JetCalib::InitRun(PHCompositeNode *topNode)
{
  // Create calib jet node.
  CreateNodeTree(topNode);

<<<<<<< Updated upstream
  // Locate the calibration file: local override first, otherwise the CDB payload.
  std::string calibFile = m_calibFileOverride.empty() ? fetchCalibDir("Default") : m_calibFileOverride;
  if (calibFile.empty())
  {
    std::cout << "JetCalib::InitRun(PHCompositeNode *topNode) : No default calibration available! Will apply calib factor 1." << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_JetCalibFile = TFile::Open(calibFile.c_str(), "READ");
  if (!m_JetCalibFile || m_JetCalibFile->IsZombie())
  {
    std::cout << "JetCalib::InitRun(PHCompositeNode *topNode) : Could not open calibration file " << calibFile << "!" << std::endl;
=======
  int ret = m_useEMfracCalib ? initEMfracCalib() : initLegacyCalib();
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  if (Verbosity() > 0)
  {
    std::cout << "JetCalib::InitRun(PHCompositeNode *topNode) - InitRun with inputNode: " << m_inputNode
              << " outputNode: " << m_outputNode
              << " method: " << (m_useEMfracCalib ? "EMfrac" : "legacy") << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

// ---------------------------------------------------------------------------
// EMfrac method initialization: CDB payload "JES_Calib_EMfrac" (or the local
// file set with set_CalibFile), holding per-radius h2d_jes_calib_r0R,
// h2_zeta_corr_r0R and f_zeta_noz_corr_r0R.
// ---------------------------------------------------------------------------
int JetCalib::initEMfracCalib()
{
  std::string calibFile = m_calibFileOverride.empty() ? fetchCalibDir("EMfrac") : m_calibFileOverride;
  if (calibFile.empty())
  {
    std::cout << "JetCalib::initEMfracCalib() : No EMfrac calibration available! Will apply calib factor 1." << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_JetCalibFile = TFile::Open(calibFile.c_str(), "READ");
  if (!m_JetCalibFile || m_JetCalibFile->IsZombie())
  {
    std::cout << "JetCalib::initEMfracCalib() : Could not open calibration file " << calibFile << "!" << std::endl;
>>>>>>> Stashed changes
    exit(1);
  }

  int r_index = (int) (jet_radius * 10 + 0.1);

  // First-pass map (required).
  m_JesCalibMap = dynamic_cast<TH2 *>(m_JetCalibFile->Get(Form("h2d_jes_calib_r0%d", r_index)));
  if (!m_JesCalibMap)
  {
<<<<<<< Updated upstream
    std::cout << "JetCalib::InitRun(PHCompositeNode *topNode) : Could not find first-pass calibration map h2d_jes_calib_r0" << r_index << "!" << std::endl;
=======
    std::cout << "JetCalib::initEMfracCalib() : Could not find first-pass calibration map h2d_jes_calib_r0" << r_index << "!" << std::endl;
>>>>>>> Stashed changes
    exit(1);
  }

  // Residual (z-vertex, eta) correction (map + no-vertex function).
  if (m_ApplyResidualCalib)
  {
    m_ZetaCorrMap = dynamic_cast<TH2 *>(m_JetCalibFile->Get(Form("h2_zeta_corr_r0%d", r_index)));
    m_NozCorrFunc = dynamic_cast<TF1 *>(m_JetCalibFile->Get(Form("f_zeta_noz_corr_r0%d", r_index)));
    if (!m_ZetaCorrMap || !m_NozCorrFunc)
    {
<<<<<<< Updated upstream
      std::cout << "JetCalib::InitRun(PHCompositeNode *topNode) : Could not find residual correction ("
=======
      std::cout << "JetCalib::initEMfracCalib() : Could not find residual correction ("
>>>>>>> Stashed changes
                << (m_ZetaCorrMap ? "" : Form("h2_zeta_corr_r0%d ", r_index))
                << (m_NozCorrFunc ? "" : Form("f_zeta_noz_corr_r0%d", r_index))
                << "); residual correction disabled." << std::endl;
      m_ZetaCorrMap = nullptr;
      m_NozCorrFunc = nullptr;
      m_ApplyResidualCalib = false;
    }
  }

  if (Verbosity() > 0)
  {
<<<<<<< Updated upstream
    std::cout << "JetCalib::InitRun(PHCompositeNode *topNode) - InitRun with inputNode: " << m_inputNode
              << " outputNode: " << m_outputNode << " calibFile: " << calibFile
=======
    std::cout << "JetCalib::initEMfracCalib() : calibFile: " << calibFile
>>>>>>> Stashed changes
              << " residual: " << (m_ApplyResidualCalib ? "on" : "off") << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

// ---------------------------------------------------------------------------
// Legacy method initialization: CDB payload "JES_Calib_Default", holding
// per-(z bin, eta bin) TF1 corrections versus pT.
// ---------------------------------------------------------------------------
int JetCalib::initLegacyCalib()
{
  if (fetchCalibDir("Default").empty())
  {
    std::cout << "JetCalib::initLegacyCalib() : No default calibration available! Will apply calib factor 1." << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_LegacyCalibFile = new CDBTF(fetchCalibDir("Default"));
  if (!m_LegacyCalibFile)
  {
    std::cout << "JetCalib::initLegacyCalib() : Could not open calibration file!" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_LegacyCalibFile->LoadCalibrations();
  std::string JetCalibFunc_name = "JES_Calib_Func_R0" + std::to_string((int) (jet_radius * 10));
  if (!ApplyZvrtxDependentCalib && !ApplyEtaDependentCalib)
  {
    int nZvrtxBins = 1;
    int nEtaBins = 1;
    m_JetCalibFunc.resize(nZvrtxBins);
    for (auto &row : m_JetCalibFunc)
    {
      row.resize(nEtaBins, nullptr);
    }
    TF1 *func_temp = m_LegacyCalibFile->getTF(JetCalibFunc_name + "_Default");
    if (!func_temp)
    {
      std::cout << "JetCalib::initLegacyCalib() : Could not find calibration function: " << JetCalibFunc_name << "_Default" << std::endl;
      exit(1);
    }
    m_JetCalibFunc[0][0] = func_temp;
  }
  else if (ApplyZvrtxDependentCalib && !ApplyEtaDependentCalib)
  {
    int nZvrtxBins = 5;
    int nEtaBins = 1;
    m_JetCalibFunc.resize(nZvrtxBins);
    for (auto &row : m_JetCalibFunc)
    {
      row.resize(nEtaBins, nullptr);
    }
    for (int iz = 0; iz < nZvrtxBins; ++iz)
    {
      TF1 *func_temp = m_LegacyCalibFile->getTF(JetCalibFunc_name + "_Z" + std::to_string(iz));
      if (!func_temp)
      {
        std::cout << "JetCalib::initLegacyCalib() : Could not find calibration function: " << JetCalibFunc_name << "_Z" << iz << std::endl;
        exit(1);
      }
      m_JetCalibFunc[iz][0] = func_temp;
    }
  }
  else if (ApplyEtaDependentCalib)
  {
    if (!ApplyZvrtxDependentCalib)
    {
      std::cout << "JetCalib::initLegacyCalib() : Must apply Zvrtx dependent calibration to apply eta dependent calibration. Applying Zvrtx + eta dependent calibration." << std::endl;
    }
    int nZvrtxBins = 5;
    int nEtaBins = 4;
    m_JetCalibFunc.resize(nZvrtxBins);
    for (auto &row : m_JetCalibFunc)
    {
      row.resize(nEtaBins, nullptr);
    }
    for (int iz = 0; iz < nZvrtxBins; ++iz)
    {
      if (iz < 3 && jet_radius <= 0.4)
      {
        for (int ieta = 0; ieta < nEtaBins; ++ieta)
        {
          TF1 *func_temp = m_LegacyCalibFile->getTF(JetCalibFunc_name + "_Z" + std::to_string(iz) + "_Eta" + std::to_string(ieta));
          if (!func_temp)
          {
            std::cout << "JetCalib::initLegacyCalib() : Could not find calibration function: " << JetCalibFunc_name << "_Z" << iz << "_Eta" << ieta << std::endl;
            exit(1);
          }
          m_JetCalibFunc[iz][ieta] = func_temp;
        }
      }
      else if (iz < 3 && jet_radius > 0.4)
      {
        for (int ieta = 0; ieta < nEtaBins; ++ieta)
        {
          TF1 *func_temp = m_LegacyCalibFile->getTF(JetCalibFunc_name + "_Z" + std::to_string(iz));
          if (!func_temp)
          {
            std::cout << "JetCalib::initLegacyCalib() : Could not find calibration function: " << JetCalibFunc_name << "_Z" << iz << std::endl;
            exit(1);
          }
          m_JetCalibFunc[iz][ieta] = func_temp;
        }
      }
      else
      {
        TF1 *func_temp = m_LegacyCalibFile->getTF(JetCalibFunc_name + "_Z" + std::to_string(iz));
        if (!func_temp)
        {
          std::cout << "JetCalib::initLegacyCalib() : Could not find calibration function: " << JetCalibFunc_name << "_Z" << iz << std::endl;
          exit(1);
        }
        m_JetCalibFunc[iz][0] = func_temp;
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int JetCalib::process_event(PHCompositeNode *topNode)
{
  // Get raw and calib jet nodes.
  JetContainer *_raw_jets = findNode::getClass<JetContainer>(topNode, m_inputNode);
  if (!_raw_jets)
  {
    std::cout << "JetCalib::process_event(PHCompositeNode *topNode) : Raw jet node not found! Cannot apply calibration." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  JetContainer *_calib_jets = findNode::getClass<JetContainer>(topNode, m_outputNode);
  if (!_calib_jets)
  {
    std::cout << "JetCalib::process_event(PHCompositeNode *topNode) : Calib jet node not found! Cannot apply calibration." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

<<<<<<< Updated upstream
  // z-vertex for the residual correction: taken from the raw jet container, i.e.
  // the vertex that was ACTUALLY USED in the jet reconstruction (recorded from the
  // jet inputs) -- not re-fetched from the GlobalVertexMap, which may differ from
  // what the jets were built with. Events whose jets carry no z-vertex
  // (has_zvertex() false) take the no-vertex correction path; get_vertex_z()
  // alone cannot distinguish them, since the jet reco falls back to z = 0.
  float zvertex = 0;
  bool hasVertex = false;
  if (m_ApplyResidualCalib)
=======
  float zvertex = m_useEMfracCalib ? 0 : -999.0;
  bool hasVertex = false;
  if ((m_useEMfracCalib && m_ApplyResidualCalib) || (!m_useEMfracCalib && ApplyZvrtxDependentCalib))
>>>>>>> Stashed changes
  {
    hasVertex = _raw_jets->has_zvertex();
    if (hasVertex)
    {
      zvertex = _raw_jets->get_vertex_z();
      if (!std::isfinite(zvertex))
      {
        hasVertex = false;
<<<<<<< Updated upstream
        zvertex = 0;
=======
        zvertex = m_useEMfracCalib ? 0 : -999.0;
>>>>>>> Stashed changes
      }
    }
    if (!hasVertex && Verbosity() > 0)
    {
<<<<<<< Updated upstream
      std::cout << "JetCalib::process_event(PHCompositeNode *topNode) : jets reconstructed without a z-vertex; applying the no-vertex correction." << std::endl;
=======
      std::cout << "JetCalib::process_event(PHCompositeNode *topNode) : jets reconstructed without a z-vertex; applying the no-vertex treatment." << std::endl;
>>>>>>> Stashed changes
    }
  }

  // Apply calibration to jets and add to calib jet node.
  int ijet = 0;
  for (auto *jet : *_raw_jets)
  {
    float pt = jet->get_pt();
    float eta = jet->get_eta();
    float phi = jet->get_phi();
    float emfrac = jet->get_emcal_frac();

    float calib_pt = pt;
    if (m_JesCalibMap && std::isfinite(emfrac))
    {
      calib_pt = getFirstPassPt(pt, emfrac);
      if (m_ApplyResidualCalib)
      {
        calib_pt *= getResidualScale(zvertex, hasVertex, eta);
      }
    }
    else if (m_JesCalibMap && !m_warnedNoEmFrac)
    {
      std::cout << "JetCalib::process_event(PHCompositeNode *topNode) : Jet has no EMCal fraction stored"
                << " (Jet::get_emcal_frac() is NaN; enable the calo fractions in the jet reco)."
                << " Applying calib factor 1 to such jets." << std::endl;
      m_warnedNoEmFrac = true;
    }

<<<<<<< Updated upstream
=======
    float calib_pt = pt;
    if (m_useEMfracCalib)
    {
      float emfrac = jet->get_emcal_frac();
      if (m_JesCalibMap && std::isfinite(emfrac))
      {
        calib_pt = getFirstPassPt(pt, emfrac);
        if (m_ApplyResidualCalib)
        {
          calib_pt *= getResidualScale(zvertex, hasVertex, eta);
        }
      }
      else if (m_JesCalibMap && !m_warnedNoEmFrac)
      {
        std::cout << "JetCalib::process_event(PHCompositeNode *topNode) : Jet has no EMCal fraction stored"
                  << " (Jet::get_emcal_frac() is NaN; enable the calo fractions in the jet reco)."
                  << " Applying calib factor 1 to such jets." << std::endl;
        m_warnedNoEmFrac = true;
      }
    }
    else if (!m_JetCalibFunc.empty())
    {
      calib_pt = doCalibration(m_JetCalibFunc, pt, zvertex, eta);
    }

>>>>>>> Stashed changes
    auto *calib_jet = _calib_jets->add_jet();
    calib_jet->set_px(calib_pt * std::cos(phi));
    calib_jet->set_py(calib_pt * std::sin(phi));
    calib_jet->set_pz(calib_pt * std::sinh(eta));
    calib_jet->set_id(ijet);
    calib_jet->insert_comp(jet->get_comp_vec(), true);
    calib_jet->set_isCalib(1);
    ijet++;
  }

  if (Verbosity() > 0)
  {
    std::cout << "JetCalib::process_event(PHCompositeNode *topNode) : nRawJets: " << _raw_jets->size() << " nCalibJets: " << _calib_jets->size() << std::endl;
    if (_calib_jets->size() != _raw_jets->size())
    {
      std::cout << "JetCalib::process_event(PHCompositeNode *topNode) : different number of raw jets vs. calib jets! Something is amiss! " << std::endl;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int JetCalib::CreateNodeTree(PHCompositeNode *topNode)
{
  // Check nodes.
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "JetCalib::CreateNodeTree : DST Node missing, aborting!" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHCompositeNode *antiktNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "ANTIKT"));
  if (!antiktNode)
  {
    std::cout << PHWHERE << "JetCalib::CreateNodeTree : ANTIKT node missing, aborting!" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHCompositeNode *towerNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "TOWER"));
  if (!towerNode)
  {
    std::cout << PHWHERE << "TOWER node not found, aborting!" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Create calib jet node if it doesn't exist.
  JetContainer *test_jets = findNode::getClass<JetContainer>(topNode, m_outputNode);
  if (!test_jets)
  {
    JetContainer *calib_jets = new JetContainerv1();
    PHIODataNode<PHObject> *calibjetNode;
    if (Verbosity() > 0)
    {
      std::cout << "JetCalib::CreateNode : creating " << m_outputNode << std::endl;
    }
    calibjetNode = new PHIODataNode<PHObject>(calib_jets, m_outputNode, "PHObject");
    towerNode->addNode(calibjetNode);
  }
  else
  {
    std::cout << "JetCalib::CreateNode : " << m_outputNode << " already exists! " << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

std::string JetCalib::fetchCalibDir(const char *calibType)
{
  std::string calibName = std::string("JES_Calib_") + calibType;
  return CDBInterface::instance()->getUrl(calibName);
}

<<<<<<< Updated upstream
// z-corrected eta acceptance window: intersection of the EMCal/IHCal/OHCal reach as
// seen from the vertex, padded inward by the jet radius. Must match the derivation.
void JetCalib::getEtaAcceptance(float zvrtx, float radius, float &minlimit, float &maxlimit)
=======
// ---------------------------------------------------------------------------
// EMfrac method helpers.
// ---------------------------------------------------------------------------

// z-corrected eta acceptance window: intersection of the EMCal/IHCal/OHCal reach as
// seen from the vertex, padded inward by the jet radius. Must match the derivation.
void JetCalib::getEtaAcceptance(float zvrtx, float radius, float &minlimit, float &maxlimit)
{
  // Calorimeter half-lengths [cm] and radii [cm].
  const float minz_EM = -130.23, maxz_EM = 130.23, radius_EM = 93.5;
  const float minz_IH = -170.299, maxz_IH = 170.299, radius_IH = 127.503;
  const float minz_OH = -301.683, maxz_OH = 301.683, radius_OH = 225.87;

  float emcal_min = std::asinh((minz_EM - zvrtx) / radius_EM);
  float emcal_max = std::asinh((maxz_EM - zvrtx) / radius_EM);
  float ihcal_min = std::asinh((minz_IH - zvrtx) / radius_IH);
  float ihcal_max = std::asinh((maxz_IH - zvrtx) / radius_IH);
  float ohcal_min = std::asinh((minz_OH - zvrtx) / radius_OH);
  float ohcal_max = std::asinh((maxz_OH - zvrtx) / radius_OH);

  minlimit = std::max({emcal_min, ihcal_min, ohcal_min}) + radius;
  maxlimit = std::min({emcal_max, ihcal_max, ohcal_max}) - radius;
}

// Bilinear interpolation with both coordinates clamped to the bin-center range
// (TH2::Interpolate is undefined outside it).
float JetCalib::interpolateClamped(TH2 *h2, double x, double y)
{
  TAxis *xa = h2->GetXaxis();
  TAxis *ya = h2->GetYaxis();
  double x_lo = xa->GetBinCenter(1);
  double x_hi = xa->GetBinCenter(xa->GetNbins());
  double y_lo = ya->GetBinCenter(1);
  double y_hi = ya->GetBinCenter(ya->GetNbins());
  if (x < x_lo)
  {
    x = x_lo;
  }
  if (x > x_hi)
  {
    x = x_hi;
  }
  if (y < y_lo)
  {
    y = y_lo;
  }
  if (y > y_hi)
  {
    y = y_hi;
  }
  return h2->Interpolate(x, y);
}

// First pass: calibrated (truth-equivalent) pT from the inverted (reco pT, EMCal
// fraction) map.
float JetCalib::getFirstPassPt(float jetPt, float emFrac) const
{
  if (!m_JesCalibMap)
  {
    return jetPt;
  }
  float calibPt = interpolateClamped(m_JesCalibMap, jetPt, emFrac);
  if (!std::isfinite(calibPt) || calibPt <= 0)
  {
    return jetPt;  // unpopulated map region (no inversion solution): leave uncalibrated
  }
  return calibPt;
}

// Residual (z-vertex, eta) scale factor. With a vertex: the (signed z, f) map;
// without: the no-vertex function evaluated at f computed with z = 0.
float JetCalib::getResidualScale(float zvrtx, bool hasVertex, float eta) const
{
  if (hasVertex)
  {
    if (!m_ZetaCorrMap)
    {
      return 1.0;
    }
    float minlimit = 0;
    float maxlimit = 0;
    getEtaAcceptance(zvrtx, jet_radius, minlimit, maxlimit);
    if (maxlimit <= minlimit)
    {
      return 1.0;  // acceptance closed at this z
    }
    double f = (eta - minlimit) / (double) (maxlimit - minlimit);
    float scale = interpolateClamped(m_ZetaCorrMap, zvrtx, f);  // SIGNED z, no mirror
    if (!std::isfinite(scale) || scale <= 0)
    {
      return 1.0;
    }
    return scale;
  }

  // No reconstructed z-vertex: f at z = 0, 1D correction function.
  if (!m_NozCorrFunc)
  {
    return 1.0;
  }
  float minlimit = 0;
  float maxlimit = 0;
  getEtaAcceptance(0.0, jet_radius, minlimit, maxlimit);
  if (maxlimit <= minlimit)
  {
    return 1.0;
  }
  double f = (eta - minlimit) / (double) (maxlimit - minlimit);
  if (f < 0 || f > 1)
  {
    return 1.0;  // outside the acceptance window
  }
  double scale = m_NozCorrFunc->Eval(f);
  if (!std::isfinite(scale) || scale < 0.5 || scale > 2.0)
  {
    return 1.0;  // sanity guard
  }
  return (float) scale;
}

// ---------------------------------------------------------------------------
// Legacy method helpers (per-(z bin, eta bin) TF1 corrections versus pT).
// ---------------------------------------------------------------------------

static int getZvrtxBin(float zvrtx)
>>>>>>> Stashed changes
{
  // Calorimeter half-lengths [cm] and radii [cm].
  const float minz_EM = -130.23, maxz_EM = 130.23, radius_EM = 93.5;
  const float minz_IH = -170.299, maxz_IH = 170.299, radius_IH = 127.503;
  const float minz_OH = -301.683, maxz_OH = 301.683, radius_OH = 225.87;

  float emcal_min = std::asinh((minz_EM - zvrtx) / radius_EM);
  float emcal_max = std::asinh((maxz_EM - zvrtx) / radius_EM);
  float ihcal_min = std::asinh((minz_IH - zvrtx) / radius_IH);
  float ihcal_max = std::asinh((maxz_IH - zvrtx) / radius_IH);
  float ohcal_min = std::asinh((minz_OH - zvrtx) / radius_OH);
  float ohcal_max = std::asinh((maxz_OH - zvrtx) / radius_OH);

  minlimit = std::max({emcal_min, ihcal_min, ohcal_min}) + radius;
  maxlimit = std::min({emcal_max, ihcal_max, ohcal_max}) - radius;
}

<<<<<<< Updated upstream
// Bilinear interpolation with both coordinates clamped to the bin-center range
// (TH2::Interpolate is undefined outside it).
float JetCalib::interpolateClamped(TH2 *h2, double x, double y)
=======
static int getEtaBin(int zvrtxbin, float eta, float jet_radius)
>>>>>>> Stashed changes
{
  TAxis *xa = h2->GetXaxis();
  TAxis *ya = h2->GetYaxis();
  double x_lo = xa->GetBinCenter(1);
  double x_hi = xa->GetBinCenter(xa->GetNbins());
  double y_lo = ya->GetBinCenter(1);
  double y_hi = ya->GetBinCenter(ya->GetNbins());
  if (x < x_lo)
  {
    x = x_lo;
  }
  if (x > x_hi)
  {
    x = x_hi;
  }
  if (y < y_lo)
  {
    y = y_lo;
  }
  if (y > y_hi)
  {
    y = y_hi;
  }
  return h2->Interpolate(x, y);
}

// First pass: calibrated (truth-equivalent) pT from the inverted (reco pT, EMCal
// fraction) map.
float JetCalib::getFirstPassPt(float jetPt, float emFrac) const
{
  if (!m_JesCalibMap)
  {
    return jetPt;
  }
  float calibPt = interpolateClamped(m_JesCalibMap, jetPt, emFrac);
  if (!std::isfinite(calibPt) || calibPt <= 0)
  {
    return jetPt;  // unpopulated map region (no inversion solution): leave uncalibrated
  }
  return calibPt;
}

// Residual (z-vertex, eta) scale factor. With a vertex: the (signed z, f) map;
// without: the no-vertex function evaluated at f computed with z = 0.
float JetCalib::getResidualScale(float zvrtx, bool hasVertex, float eta) const
{
  if (hasVertex)
  {
    if (!m_ZetaCorrMap)
    {
      return 1.0;
    }
    float minlimit = 0;
    float maxlimit = 0;
    getEtaAcceptance(zvrtx, jet_radius, minlimit, maxlimit);
    if (maxlimit <= minlimit)
    {
      return 1.0;  // acceptance closed at this z
    }
    double f = (eta - minlimit) / (double) (maxlimit - minlimit);
    float scale = interpolateClamped(m_ZetaCorrMap, zvrtx, f);  // SIGNED z, no mirror
    if (!std::isfinite(scale) || scale <= 0)
    {
      return 1.0;
    }
    return scale;
  }

  // No reconstructed z-vertex: f at z = 0, 1D correction function.
  if (!m_NozCorrFunc)
  {
    return 1.0;
  }
  float minlimit = 0;
  float maxlimit = 0;
  getEtaAcceptance(0.0, jet_radius, minlimit, maxlimit);
  if (maxlimit <= minlimit)
  {
    return 1.0;
  }
  double f = (eta - minlimit) / (double) (maxlimit - minlimit);
  if (f < 0 || f > 1)
  {
    return 1.0;  // outside the acceptance window
  }
  double scale = m_NozCorrFunc->Eval(f);
  if (!std::isfinite(scale) || scale < 0.8 || scale > 1.2)
  {
    return 1.0;  // sanity guard
  }
  return (float) scale;
}
