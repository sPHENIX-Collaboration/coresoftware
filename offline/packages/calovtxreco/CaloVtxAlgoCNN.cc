#include "CaloVtxAlgoCNN.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <onnxruntime_cxx_api.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <utility>

namespace
{
  const std::string TowerNode[CaloVtxAlgoCNN::kNLayer] = {
      "TOWERINFO_CALIB_CEMC", "TOWERINFO_CALIB_HCALIN", "TOWERINFO_CALIB_HCALOUT"};
  const std::string GeomNodeEmc = "TOWERGEOM_CEMC";
  const std::string GeomNodeIhc = "TOWERGEOM_HCALIN";
}  // namespace

// onnxruntime session (pImpl, keeps Ort types out of the header)
struct CaloVtxAlgoCNN::OnnxSession
{
  //NOLINTBEGIN(misc-non-private-member-variables-in-classes)
  Ort::Env env{ORT_LOGGING_LEVEL_WARNING, "CaloVtxAlgoCNN"};
  Ort::SessionOptions opts;
  std::unique_ptr<Ort::Session> session;
  Ort::MemoryInfo memInfo{
      Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault)};
// NOLINTEND(misc-non-private-member-variables-in-classes)

  explicit OnnxSession(const std::string &path)
  {
    opts.SetIntraOpNumThreads(1);
    opts.SetGraphOptimizationLevel(ORT_ENABLE_ALL);
    session = std::make_unique<Ort::Session>(env, path.c_str(), opts);
  }
};

CaloVtxAlgoCNN::~CaloVtxAlgoCNN() = default;

int CaloVtxAlgoCNN::Init(PHCompositeNode * /*topNode*/)
{
  try
  {
    m_onnx = std::make_unique<OnnxSession>(m_modelFile);
  }
  catch (const std::exception &e)
  {
    std::cout << PHWHERE << "failed to load model " << m_modelFile << ": " << e.what() << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloVtxAlgoCNN::CalculateVertex(PHCompositeNode *topNode, float &zvtx)
{
  zvtx = std::numeric_limits<float>::quiet_NaN();

  if (fillTowerImage(topNode) != 0)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  double etot = 0.;
  for (auto & layer : m_image)
  {
    for (auto & ieta : layer)
    {
      for (float iphi : ieta)
      {
        etot += iphi;
      }
    }
  }
  if (etot <= m_minTotalEnergy)
  {
    // empty image: leave zvtx NaN, do not run the network
    return Fun4AllReturnCodes::EVENT_OK;
  }

  float z = std::numeric_limits<float>::quiet_NaN();
  if (!predict(z))
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  zvtx = z;
  return Fun4AllReturnCodes::EVENT_OK;
}

bool CaloVtxAlgoCNN::predict(float &z)
{
  const int64_t shape[4] = {1, kNLayer, kNEtaImg, kNPhiImg};
  try
  {
    Ort::Value input = Ort::Value::CreateTensor<float>(m_onnx->memInfo, &m_image[0][0][0], static_cast<size_t>(kNLayer) * kNEtaImg * kNPhiImg, shape, 4);
    const char *inNames[] = {"raw_image"};
    const char *outNames[] = {"z_cal_cm"};
    auto outs = m_onnx->session->Run(Ort::RunOptions{nullptr}, inNames, &input, 1, outNames, 1);
    z = outs[0].GetTensorData<float>()[0];
    return true;
  }
  catch (const std::exception &e)
  {
    if (!m_warnedPredict)
    {
      std::cout << PHWHERE << "model evaluation failed: " << e.what() << std::endl;
      m_warnedPredict = true;
    }
    return false;
  }
}

// pre-processing: fill the 3-layer tower image from the nodes, including retowering of the EMCAL
int CaloVtxAlgoCNN::fillTowerImage(PHCompositeNode *topNode)
{
  for (auto & layer : m_image)
  {
    for (auto & ieta : layer)
    {
      for (float & iphi : ieta)
      {
        iphi = 0.;
      }
    }
  }
  // HCal layers are already on the common 24x64 grid.
  if (fillHcalLayer(topNode, kIHC) != 0)
  {
    return -1;
  }
  if (fillHcalLayer(topNode, kOHC) != 0)
  {
    return -1;
  }
  return fillEmcRetower(topNode);
}

int CaloVtxAlgoCNN::fillHcalLayer(PHCompositeNode *topNode, int layer)
{
  TowerInfoContainer *towers = findNode::getClass<TowerInfoContainer>(topNode, TowerNode[layer]);
  if (!towers)
  {
    if (!m_warnedMissingTowers)
    {
      std::cout << PHWHERE << "tower node missing: " << TowerNode[layer] << std::endl;
      m_warnedMissingTowers = true;
    }
    return -1;
  }

  const unsigned int ntow = towers->size();
  for (unsigned int ch = 0; ch < ntow; ++ch)
  {
    TowerInfo *tower = towers->get_tower_at_channel(ch);
    if (!tower || !tower->get_isGood())
    {
      continue;
    }
    const float e = tower->get_energy();
    if (e < m_towerEMin[layer])
    {
      continue;
    }
    const unsigned int key = towers->encode_key(ch);
    const int ieta = towers->getTowerEtaBin(key);
    const int iphi = towers->getTowerPhiBin(key);
    if (ieta < 0 || ieta >= kNEtaImg || iphi < 0 || iphi >= kNPhiImg)
    {
      continue;
    }
    m_image[layer][ieta][iphi] += e;
  }
  return 0;
}

int CaloVtxAlgoCNN::fillEmcRetower(PHCompositeNode *topNode)
{
  if (buildEmcRetowerMap(topNode) != 0)
  {
    return -1;
  }

  TowerInfoContainer *towers = findNode::getClass<TowerInfoContainer>(topNode, TowerNode[kEMC]);
  if (!towers)
  {
    if (!m_warnedMissingTowers)
    {
      std::cout << PHWHERE << "tower node missing: " << TowerNode[kEMC] << std::endl;
      m_warnedMissingTowers = true;
    }
    return -1;
  }

  for (auto &row : m_rawEmcFine)
  {
    for (double &e : row)
    {
      e = 0.;
    }
  }

  const unsigned int ntow = towers->size();
  for (unsigned int ch = 0; ch < ntow; ++ch)
  {
    TowerInfo *tower = towers->get_tower_at_channel(ch);
    if (!tower || !tower->get_isGood())
    {
      continue;
    }
    const float e = tower->get_energy();
    if (e < m_towerEMin[kEMC])
    {
      continue;
    }
    const unsigned int key = towers->encode_key(ch);
    const int ieta = towers->getTowerEtaBin(key);
    const int iphi = towers->getTowerPhiBin(key);
    if (ieta < 0 || ieta >= kNEtaEmcFine || iphi < 0 || iphi >= kNPhiEmcFine)
    {
      continue;
    }
    m_rawEmcFine[ieta][iphi] = e;
  }

  // eta-fraction / phi-grouping sums, as in RetowerCEMC
  for (int ietaHcal = 0; ietaHcal < kNEtaImg; ++ietaHcal)
  {
    for (int iphiHcal = 0; iphiHcal < kNPhiImg; ++iphiHcal)
    {
      double retowerE = 0.;
      for (int ietaEmc = m_retowerLowerEta[ietaHcal]; ietaEmc <= m_retowerUpperEta[ietaHcal]; ++ietaEmc)
      {
        double fraction = 1.;
        if (ietaEmc == m_retowerLowerEta[ietaHcal])
        {
          fraction = m_retowerLowerFrac[ietaHcal];
        }
        else if (ietaEmc == m_retowerUpperEta[ietaHcal])
        {
          fraction = m_retowerUpperFrac[ietaHcal];
        }
        for (int iphiEmc = m_retowerPhiOffset + iphiHcal * 4; iphiEmc < m_retowerPhiOffset + iphiHcal * 4 + 4; ++iphiEmc)
        {
          int iphiEmcWrap = iphiEmc;
          if (iphiEmcWrap > kNPhiEmcFine - 1)
          {
            iphiEmcWrap -= kNPhiEmcFine;
          }
          retowerE += m_rawEmcFine[ietaEmc][iphiEmcWrap] * fraction;
        }
      }
      m_image[kEMC][ietaHcal][iphiHcal] = retowerE;
    }
  }
  return 0;
}

int CaloVtxAlgoCNN::buildEmcRetowerMap(PHCompositeNode *topNode)
{
  if (m_retowerMapReady)
  {
    return 0;
  }

  RawTowerGeomContainer *geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, GeomNodeEmc);
  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, GeomNodeIhc);
  if (!geomEM || !geomIH)
  {
    if (!m_warnedRetowerMap)
    {
      std::cout << PHWHERE << "cannot build EMCal retower map, missing " << GeomNodeEmc << " or " << GeomNodeIhc << std::endl;
      m_warnedRetowerMap = true;
    }
    return -1;
  }

  // first fine EMCal phi bin belonging to HCal phi 0, cf. RetowerCEMC::get_first_phi_index()
  bool foundFirstLowerBound = false;
  int iphiEmc = 0;
  while (iphiEmc < kNPhiEmcFine)
  {
    const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CEMC, 0, iphiEmc);
    RawTowerGeom *towerGeom = geomEM->get_tower_geometry(key);
    if (towerGeom && geomIH->get_phibin(towerGeom->get_phi()) == 0)
    {
      foundFirstLowerBound = true;
      break;
    }
    ++iphiEmc;
  }

  if (foundFirstLowerBound && iphiEmc == 0)
  {
    bool outOfRange = false;
    int iphiEmcTemp = kNPhiEmcFine - 1;
    while (iphiEmcTemp > iphiEmc)
    {
      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CEMC, 0, iphiEmcTemp);
      RawTowerGeom *towerGeom = geomEM->get_tower_geometry(key);
      if (towerGeom && geomIH->get_phibin(towerGeom->get_phi()) == kNPhiImg - 1)
      {
        outOfRange = true;
        break;
      }
      --iphiEmcTemp;
    }
    if (!outOfRange)
    {
      if (!m_warnedRetowerMap)
      {
        std::cout << PHWHERE << "cannot build EMCal retower map, no wrap-around " << "phi match" << std::endl;
        m_warnedRetowerMap = true;
      }
      return -1;
    }
    m_retowerPhiOffset = (iphiEmcTemp + 1 == kNPhiEmcFine) ? 0 : iphiEmcTemp + 1;
  }
  else if (!foundFirstLowerBound)
  {
    if (!m_warnedRetowerMap)
    {
      std::cout << PHWHERE << "cannot build EMCal retower map, no EMCal phi bin " << "maps to HCal phi 0" << std::endl;
      m_warnedRetowerMap = true;
    }
    return -1;
  }
  else
  {
    m_retowerPhiOffset = iphiEmc;
  }

  // eta-bound overlaps (edge bins fractional), cf. RetowerCEMC::get_weighted_fraction()
  int ietaEmc = 0;
  for (int ietaHcal = 0; ietaHcal < kNEtaImg; ++ietaHcal)
  {
    const std::pair<double, double> rangeHcal = geomIH->get_etabounds(ietaHcal);
    const double hcalLower = rangeHcal.first;
    const double hcalUpper = rangeHcal.second;
    bool foundLower = false;
    bool foundUpper = false;

    while ((!foundLower || !foundUpper) && ietaEmc < kNEtaEmcFine)
    {
      const std::pair<double, double> rangeEmc = geomEM->get_etabounds(ietaEmc);
      const double emcLower = rangeEmc.first;
      const double emcUpper = rangeEmc.second;

      if (!foundLower)
      {
        if (emcUpper > hcalLower && emcLower <= hcalLower)
        {
          m_retowerLowerEta[ietaHcal] = ietaEmc;
          m_retowerLowerFrac[ietaHcal] = (emcUpper - hcalLower) / (emcUpper - emcLower);
          foundLower = true;
        }
        if (emcUpper > hcalLower && emcLower > hcalLower)
        {
          m_retowerLowerEta[ietaHcal] = ietaEmc;
          m_retowerLowerFrac[ietaHcal] = 1.;
          foundLower = true;
        }
      }
      else
      {
        if (emcUpper >= hcalUpper && emcLower < hcalUpper)
        {
          m_retowerUpperEta[ietaHcal] = ietaEmc;
          m_retowerUpperFrac[ietaHcal] = (hcalUpper - emcLower) / (emcUpper - emcLower);
          foundUpper = true;
        }
        if (emcUpper > hcalUpper && emcLower > hcalUpper)
        {
          --ietaEmc;
          m_retowerUpperEta[ietaHcal] = ietaEmc;
          m_retowerUpperFrac[ietaHcal] = 1.;
          foundUpper = true;
        }
      }

      if (!(foundLower && foundUpper))
      {
        ++ietaEmc;
      }
    }

    if (!foundLower || !foundUpper)
    {
      if (!m_warnedRetowerMap)
      {
        std::cout << PHWHERE << "cannot build EMCal retower map, missing " << (foundLower ? "upper" : "lower") << " eta overlap for HCal ieta " << ietaHcal << std::endl;
        m_warnedRetowerMap = true;
      }
      return -1;
    }
  }

  m_retowerMapReady = true;
  return 0;
}
