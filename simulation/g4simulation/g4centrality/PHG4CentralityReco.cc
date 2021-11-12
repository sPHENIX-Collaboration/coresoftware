#include "PHG4CentralityReco.h"

#include <centrality/CentralityInfo.h>    // for CentralityInfo, CentralityI...
#include <centrality/CentralityInfov1.h>  // for CentralityInfov1

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <iostream>   // for operator<<, basic_ostream
#include <map>        // for _Rb_tree_const_iterator
#include <stdexcept>  // for runtime_error
#include <utility>    // for pair

PHG4CentralityReco::PHG4CentralityReco(const std::string &name)
  : SubsysReco(name)
{
}

int PHG4CentralityReco::InitRun(PHCompositeNode *topNode)
{
  if (Verbosity() >= 1)
    std::cout << " PHG4CentralityReco::InitRun : enter " << std::endl;

  try
  {
    CreateNode(topNode);
  }
  catch (std::exception &e)
  {
    std::cout << PHWHERE << ": " << e.what() << std::endl;
    throw;
  }

  auto bhits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_BBC");
  if (!bhits)
    std::cout << "PHG4CentralityReco::InitRun : cannot find G4HIT_BBC, will not use MBD centrality";

  auto ehits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_EPD");
  if (!ehits)
    std::cout << "PHG4CentralityReco::InitRun : cannot find G4HIT_EPD, will not use sEPD centrality";

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4CentralityReco::process_event(PHCompositeNode *topNode)
{
  std::cout << "PHG4CentralityReco::process_event -- heartbeat" << std::endl;

  _mbd_N = 0;
  _mbd_S = 0;
  _mbd_NS = 0;

  auto bhits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_BBC");

  if (bhits)
  {
    auto brange = bhits->getHits();
    for (auto it = brange.first; it != brange.second; ++it)
    {
      if ((it->second->get_t(0) > -50) && (it->second->get_t(1) < 50))
      {
        _mbd_NS += it->second->get_edep();
        int id = it->second->get_layer();
        if ((id & 0x40) == 0)
          _mbd_N += it->second->get_edep();
        else
          _mbd_S += it->second->get_edep();
      }
    }

    if (Verbosity() >= 5)
      std::cout << "PHG4CentralityReco::process_event : MBD Sum Charge N / S / N+S = " << _mbd_N << " / " << _mbd_S << " / " << _mbd_NS << std::endl;
  }
  else
  {
    if (Verbosity() >= 5)
      std::cout << "PHG4CentralityReco::process_event : No MBD info, setting all Sum Charges = 0" << std::endl;
  }

  _epd_N = 0;
  _epd_S = 0;
  _epd_NS = 0;

  auto ehits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_EPD");

  if (ehits)
  {
    auto erange = ehits->getHits();
    for (auto it = erange.first; it != erange.second; ++it)
      if ((it->second->get_t(0) > -50) && (it->second->get_t(1) < 50))
      {
        _epd_NS += it->second->get_edep();
        int id = it->second->get_scint_id();
        if ((id & 0x200) == 0)
          _epd_N += it->second->get_edep();
        else
          _epd_S += it->second->get_edep();
      }

    if (Verbosity() >= 5)
      std::cout << "PHG4CentralityReco::process_event : sEPD Sum Energy N / S / N+S = " << _epd_N << " / " << _epd_S << " / " << _epd_NS << std::endl;
  }
  else
  {
    if (Verbosity() >= 5)
      std::cout << "PHG4CentralityReco::process_event : No sEPD info, setting all Sum Energies = 0" << std::endl;
  }

  if (Verbosity() >= 1)
    std::cout << "PHG4CentralityReco::process_event : summary MBD (N, S, N+S) = (" << _mbd_N << ", " << _mbd_S << ", " << _mbd_NS << "), sEPD (N, S, N+S) = (" << _epd_N << ", " << _epd_S << ", " << _epd_NS << ")" << std::endl;

  FillNode(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4CentralityReco::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4CentralityReco::FillNode(PHCompositeNode *topNode)
{
  CentralityInfo *cent = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
  if (!cent)
  {
    std::cout << " ERROR -- can't find CentralityInfo node after it should have been created" << std::endl;
    return;
  }
  else
  {
    cent->set_quantity(CentralityInfo::PROP::mbd_N, _mbd_N);
    cent->set_quantity(CentralityInfo::PROP::mbd_S, _mbd_S);
    cent->set_quantity(CentralityInfo::PROP::mbd_NS, _mbd_NS);
    cent->set_quantity(CentralityInfo::PROP::epd_N, _epd_N);
    cent->set_quantity(CentralityInfo::PROP::epd_S, _epd_S);
    cent->set_quantity(CentralityInfo::PROP::epd_NS, _epd_NS);
  }
}

void PHG4CentralityReco::CreateNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find DST node in PHG4CentralityReco::CreateNode");
  }

  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "GLOBAL"));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode("GLOBAL");
    dstNode->addNode(DetNode);
  }

  CentralityInfo *cent = new CentralityInfov1();

  PHIODataNode<PHObject> *centNode = new PHIODataNode<PHObject>(cent, "CentralityInfo", "PHObject");
  DetNode->addNode(centNode);
}
