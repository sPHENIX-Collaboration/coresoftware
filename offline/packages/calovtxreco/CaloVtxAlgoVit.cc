#include "CaloVtxAlgoVit.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <onnxruntime_cxx_api.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

// onnxruntime session (pImpl, keeps Ort types out of the header)
struct CaloVtxAlgoVit::OnnxSession
{
  Ort::Env env{ORT_LOGGING_LEVEL_WARNING, "CaloVtxAlgoVit"};
  Ort::SessionOptions opts;
  std::unique_ptr<Ort::Session> session;
  Ort::MemoryInfo memInfo{
      Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault)};

  explicit OnnxSession(const std::string &path)
  {
    // pin to one compute thread: on shared (Condor) nodes onnxruntime must
    // not spawn one thread per node core
    opts.SetIntraOpNumThreads(1);
    opts.SetInterOpNumThreads(1);
    opts.SetGraphOptimizationLevel(ORT_ENABLE_ALL);
    session = std::make_unique<Ort::Session>(env, path.c_str(), opts);
  }
};

CaloVtxAlgoVit::CaloVtxAlgoVit() = default;
CaloVtxAlgoVit::~CaloVtxAlgoVit() = default;

int CaloVtxAlgoVit::Init(PHCompositeNode * /*topNode*/)
{
  for (int calo = 0; calo < kNCalo; ++calo)
  {
    m_input[calo].assign(static_cast<size_t>(kNChan) * kNEta[calo] * kNPhi[calo], 0.);
  }

  try
  {
    m_onnx = std::make_unique<OnnxSession>(m_modelFile);
  }
  catch (const std::exception &e)
  {
    std::cout << PHWHERE << "failed to load model " << m_modelFile << ": " << e.what() << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  const size_t nIn = m_onnx->session->GetInputCount();
  if (nIn != kNCalo)
  {
    std::cout << PHWHERE << "model " << m_modelFile << " has " << nIn
              << " inputs, expected " << kNCalo << " (emcal, ihcal, ohcal)" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloVtxAlgoVit::CalculateVertex(PHCompositeNode *topNode, float &zvtx)
{
  zvtx = std::numeric_limits<float>::quiet_NaN();

  if (fillInputs(topNode) != 0)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (m_etot <= m_minTotalEnergy)
  {
    // empty event: leave zvtx NaN, do not run the network
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

int CaloVtxAlgoVit::fillInputs(PHCompositeNode *topNode)
{
  m_etot = 0.;
  for (int calo = 0; calo < kNCalo; ++calo)
  {
    std::fill(m_input[calo].begin(), m_input[calo].end(), 0.);
    if (fillCalo(topNode, calo) != 0)
    {
      return -1;
    }
  }
  return 0;
}

int CaloVtxAlgoVit::fillCalo(PHCompositeNode *topNode, int calo)
{
  TowerInfoContainer *towers = findNode::getClass<TowerInfoContainer>(topNode, m_towerNode[calo]);
  if (!towers)
  {
    if (!m_warnedMissingTowers)
    {
      std::cout << PHWHERE << "tower node missing: " << m_towerNode[calo] << std::endl;
      m_warnedMissingTowers = true;
    }
    return -1;
  }

  const int nEta = kNEta[calo];
  const int nPhi = kNPhi[calo];
  float *eChan = m_input[calo].data();                                      // channel 0: raw energy
  float *tChan = m_input[calo].data() + static_cast<size_t>(nEta) * nPhi;  // channel 1: time

  const unsigned int ntow = towers->size();
  for (unsigned int ch = 0; ch < ntow; ++ch)
  {
    TowerInfo *tower = towers->get_tower_at_channel(ch);
    if (!tower)
    {
      continue;
    }
    if (m_useGoodTowersOnly && !tower->get_isGood())
    {
      continue;
    }
    const unsigned int key = towers->encode_key(ch);
    const int ieta = towers->getTowerEtaBin(key);
    const int iphi = towers->getTowerPhiBin(key);
    if (ieta < 0 || ieta >= nEta || iphi < 0 || iphi >= nPhi)
    {
      continue;
    }
    // raw energy floored at 0, NO log1p (the graph applies it); keep the time
    const float e = std::max(tower->get_energy(), 0.f);
    const size_t idx = static_cast<size_t>(ieta) * nPhi + iphi;
    eChan[idx] = e;
    tChan[idx] = tower->get_time();
    m_etot += e;
  }
  return 0;
}

bool CaloVtxAlgoVit::predict(float &z)
{
  try
  {
    std::vector<Ort::Value> inputs;
    inputs.reserve(kNCalo);
    const char *inNames[kNCalo];
    for (int calo = 0; calo < kNCalo; ++calo)
    {
      const int64_t shape[4] = {1, kNChan, kNEta[calo], kNPhi[calo]};
      inputs.push_back(Ort::Value::CreateTensor<float>(
          m_onnx->memInfo, m_input[calo].data(), m_input[calo].size(), shape, 4));
      inNames[calo] = m_inputName[calo].c_str();
    }
    const char *outNames[] = {m_outputName.c_str()};
    auto outs = m_onnx->session->Run(Ort::RunOptions{nullptr}, inNames, inputs.data(), kNCalo, outNames, 1);
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
