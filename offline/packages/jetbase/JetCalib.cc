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
  delete m_JetCalibFile;
  delete m_JetCalibFunc;
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
      m_JetCalibFunc = m_JetCalibFile->getTF("JES_Calib_Default_Func");
      if (!m_JetCalibFunc)
      {
        std::cout << "JetCalib::InitRun(PHCompositeNode *topNode) : Calibration function not found!" << std::endl;
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
  JetContainer *_calib_jets = findNode::getClass<JetContainer>(topNode, m_outputNode);

  // Apply calibration to jets and add to calib jet node.
  int ijet = 0;
  for (auto jet : *_raw_jets)
  {
    float pt = jet->get_pt();
    float calib_pt = getTotalCorrFactor(m_JetCalibFunc, pt);
    auto calib_jet = _calib_jets->add_jet();
    calib_jet->set_px(calib_pt * cos(jet->get_phi()));
    calib_jet->set_py(calib_pt * sin(jet->get_phi()));
    calib_jet->set_pz(calib_pt * sinh(jet->get_eta()));
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

float JetCalib::getTotalCorrFactor(TF1 *JetCalibFunc, float jetPt)
{
  float calib = 1;
  calib = JetCalibFunc->Eval(jetPt);
  return calib;
}
