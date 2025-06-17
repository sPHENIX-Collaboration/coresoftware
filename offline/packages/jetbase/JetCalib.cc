#include "JetCalib.h"

#include "JetContainerv1.h"

#include <cdbobjects/CDBTF.h>  // for CDBTF1

#include <ffamodules/CDBInterface.h>
#include <ffaobjects/EventHeader.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <TF1.h>
#include <TSystem.h>

#include <boost/format.hpp>

#include <cmath>
#include <cstdlib>    // for exit
#include <exception>  // for exception
#include <iostream>   // for operator<<, basic_ostream
#include <stdexcept>  // for runtime_error

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
  if (Verbosity() > 0)
  {
    std::cout << "JetCalib::~JetCalib() : Calling dtor." << std::endl;
  }
}

int JetCalib::InitRun(PHCompositeNode *topNode)
{
  // Create calib jet node.
  CreateNodeTree(topNode);

  // Load calibration file and function.
  if (fetchCalibDir("Default").empty())
  {
    std::cout << "JetCalib::InitRun(PHCompositeNode *topNode) : No default calibration available! Will apply calib factor 1." << std::endl;
  }
  else
  {
    m_JetCalibFile = new CDBTF(fetchCalibDir("Default"));
    if (m_JetCalibFile)
    {
      m_JetCalibFile->LoadCalibrations();
      std::string JetCalibFunc_name = "JES_Calib_Func_R0" + std::to_string(jet_radius * 10);
      if (!ApplyZvrtxDependentCalib && !ApplyEtaDependentCalib)
      {
        int nZvrtxBins = 1;
        int nEtaBins = 1;
        m_JetCalibFunc.resize(nZvrtxBins);
        for (auto &row : m_JetCalibFunc)
        {
          row.resize(nEtaBins, nullptr);
        }
        TF1 *func_temp = m_JetCalibFile->getTF(JetCalibFunc_name + "_Default");
        if (!func_temp)
        {
          std::cout << "JetCalib::InitRun(PHCompositeNode *topNode) : Could not find calibration function: " << JetCalibFunc_name << "_Default" << std::endl;
          exit(1);
        }
        m_JetCalibFunc[0][0] = func_temp;
      }
      else if (ApplyZvrtxDependentCalib && !ApplyEtaDependentCalib)
      {
        int nZvrtxBins = 3;
        int nEtaBins = 1;
        m_JetCalibFunc.resize(nZvrtxBins);
        for (auto &row : m_JetCalibFunc)
        {
          row.resize(nEtaBins, nullptr);
        }
        for (int iz = 0; iz < nZvrtxBins; ++iz)
        {
          TF1 *func_temp = m_JetCalibFile->getTF(JetCalibFunc_name + "_Zvrtx" + std::to_string(iz));
          if (!func_temp)
          {
            std::cout << "JetCalib::InitRun(PHCompositeNode *topNode) : Could not find calibration function: " << JetCalibFunc_name << "_Zvrtx" << iz << std::endl;
            exit(1);
          }
          m_JetCalibFunc[iz][0] = func_temp;
        }
      }
      else if (ApplyEtaDependentCalib)
      {
        if (!ApplyZvrtxDependentCalib)
        {
          std::cout << "JetCalib::InitRun(PHCompositeNode *topNode) : Must apply Zvrtx dependent calibration to apply eta dependent calibration. Applying Zvrtx + eta dependent calibration." << std::endl;
        }
        int nZvrtxBins = 3;
        int nEtaBins = 4;
        m_JetCalibFunc.resize(nZvrtxBins);
        for (auto &row : m_JetCalibFunc)
        {
          row.resize(nEtaBins, nullptr);
        }
        for (int iz = 0; iz < nZvrtxBins; ++iz)
        {
          for (int ieta = 0; ieta < nEtaBins; ++ieta)
          {
            TF1 *func_temp = m_JetCalibFile->getTF(JetCalibFunc_name + "_Zvrtx" + std::to_string(iz) + "_Eta" + std::to_string(ieta));
            if (!func_temp)
            {
              std::cout << "JetCalib::InitRun(PHCompositeNode *topNode) : Could not find calibration function: " << JetCalibFunc_name << "_Zvrtx" << iz << "_Eta" << ieta << std::endl;
              exit(1);
            }
            m_JetCalibFunc[iz][ieta] = func_temp;
          }
        }
      }
    }
    else
    {
      std::cout << "JetCalib::InitRun(PHCompositeNode *topNode) : Could not open calibration file!" << std::endl;
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << "JetCalib::InitRun(PHCompositeNode *topNode) - InitRun with inputNode: " << m_inputNode << " outputNode: " << m_outputNode << std::endl;
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
  GlobalVertexMap *vertexmap = nullptr;
  if (ApplyZvrtxDependentCalib)
  {
    vertexmap = findNode::getClass<GlobalVertexMap>(topNode, m_zvrtxNode);
    if (!vertexmap)
    {
      std::cout << "JetCalib::process_event(PHCompositeNode *topNode) : GlobalVertexMap node not found! Cannot apply Z-vertex dependent calibration." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }
  // Apply calibration to jets and add to calib jet node.
  int ijet = 0;
  for (auto *jet : *_raw_jets)
  {
    float pt = jet->get_pt();
    float eta = jet->get_eta();
    float phi = jet->get_phi();
    float zvertex = 0.0;
    if (ApplyZvrtxDependentCalib)
    {
      if (!vertexmap->empty())
      {
        GlobalVertex *vtx = vertexmap->begin()->second;
        if (vtx)
        {
          zvertex = vtx->get_z();
        }
      }
      else
      {
        std::cout << "JetCalib::process_event(PHCompositeNode *topNode) : GlobalVertexMap is empty. Assign zvertex = 0." << std::endl;
      }
    }
    float calib_pt = doCalibration(m_JetCalibFunc, pt, zvertex, eta);
    auto *calib_jet = _calib_jets->add_jet();
    calib_jet->set_px(calib_pt * std::cos(phi));
    calib_jet->set_py(calib_pt * std::sin(phi));
    calib_jet->set_pz(calib_pt * std::sinh(eta));
    calib_jet->set_id(ijet);
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

int getZvrtxBin(float zvrtx)
{
  if (zvrtx < -30.0)
  {
    return 1;  // -inf to -30
  }
  if (zvrtx < 30.0)
  {
    return 0;  // -30 to 30
  }

  return 2;  // 30 to inf
}

int getEtaBin(int zvrtxbin, float eta, float jet_radius)
{
  float eta_low = 0;
  float eta_high = 0;
  if (zvrtxbin == 0)  // -30 < zvrtx < 30
  {
    eta_low = -1.2 + jet_radius;
    eta_high = 1.2 - jet_radius;
  }
  else if (zvrtxbin == 1)  // -60 < zvrtx < -30
  {
    eta_low = -0.95 + jet_radius;
    eta_high = 1.25 - jet_radius;
  }
  else if (zvrtxbin == 2)  // 30 < zvrtx < inf
  {
    eta_low = -1.25 + jet_radius;
    eta_high = 0.95 - jet_radius;
  }

  float threshold1 = eta_low + ((eta_high - eta_low) / 4.0);
  float threshold2 = eta_low + ((eta_high - eta_low) / 2.0);
  float threshold3 = eta_low + (3.0 * (eta_high - eta_low) / 4.0);

  if (eta < threshold1)
  {
    return 0;  // -inf to threshold1
  }
  if (eta < threshold2)
  {
    return 1;  // threshold1 to threshold2
  }
  if (eta < threshold3)
  {
    return 2;  // threshold2 to threshold3
  }

  return 3;  // threshold3 to inf
}

float JetCalib::doCalibration(const std::vector<std::vector<TF1 *>> &JetCalibFunc, float jetPt, float zvrtx, float eta) const
{
  float calib = 1;
  int zvrtxbin = 0;
  int etabin = 0;
  if (ApplyZvrtxDependentCalib)
  {
    zvrtxbin = getZvrtxBin(zvrtx);
    if (ApplyEtaDependentCalib)
    {
      etabin = getEtaBin(zvrtxbin, eta, jet_radius);
    }
  }
  calib = JetCalibFunc[zvrtxbin][etabin]->Eval(jetPt);
  return calib;
}
