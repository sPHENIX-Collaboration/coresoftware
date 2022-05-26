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

#include <ffaobjects/EventHeaderv1.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <iostream>  // for operator<<, basic_ostream
#include <map>       // for _Rb_tree_const_iterator
#include <sstream>
#include <stdexcept>  // for runtime_error
#include <utility>    // for pair

PHG4CentralityReco::PHG4CentralityReco(const std::string &name)
  : SubsysReco(name)
  , _centrality_calibration_params(name)
{
}

int PHG4CentralityReco::InitRun(PHCompositeNode *topNode)
{
  if (Verbosity() >= 1)
    std::cout << "PHG4CentralityReco::InitRun : enter " << std::endl;

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
    std::cout << "PHG4CentralityReco::InitRun : cannot find G4HIT_BBC, will not use MBD centrality" << std::endl;

  auto ehits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_EPD");
  if (!ehits)
    std::cout << "PHG4CentralityReco::InitRun : cannot find G4HIT_EPD, will not use sEPD centrality" << std::endl;

  if (_centrality_calibration_params.exist_string_param("description"))
  {
    if (Verbosity() >= 1)
    {
      std::cout << "PHG4CentralityReco::InitRun : Centrality calibration description : " << std::endl
                << "    ";
      std::cout << _centrality_calibration_params.get_string_param("description") << std::endl;
    }
    // search for possible centile definitions
    for (int n = 0; n < 101; n++)
    {
      std::ostringstream s1;
      s1 << "epd_centile_" << n;
      if (_centrality_calibration_params.exist_double_param(s1.str().c_str()))
      {
        _cent_cal_epd[_centrality_calibration_params.get_double_param(s1.str().c_str())] = n;
        if (Verbosity() >= 2)
          std::cout << "PHG4CentralityReco::InitRun : sEPD centrality calibration, centile " << n << "% is " << _centrality_calibration_params.get_double_param(s1.str().c_str()) << std::endl;
      }
    }
    for (int n = 0; n < 101; n++)
    {
      std::ostringstream s2;
      s2 << "mbd_centile_" << n;
      if (_centrality_calibration_params.exist_double_param(s2.str().c_str()))
      {
        _cent_cal_mbd[_centrality_calibration_params.get_double_param(s2.str().c_str())] = n;
        if (Verbosity() >= 2)
          std::cout << "PHG4CentralityReco::InitRun : MBD centrality calibration, centile " << n << "% is " << _centrality_calibration_params.get_double_param(s2.str().c_str()) << std::endl;
      }
    }
    for (int n = 0; n < 101; n++)
    {
      std::ostringstream s3;
      s3 << "bimp_centile_" << n;
      if (_centrality_calibration_params.exist_double_param(s3.str().c_str()))
      {
        _cent_cal_bimp[_centrality_calibration_params.get_double_param(s3.str().c_str())] = n;
        if (Verbosity() >= 2)
          std::cout << "PHG4CentralityReco::InitRun : b (impact parameter) centrality calibration, centile " << n << "% is " << _centrality_calibration_params.get_double_param(s3.str().c_str()) << std::endl;
      }
    }
  }
  else
  {
    std::cout << "PHG4CentralityReco::InitRun : no centrality calibration found!" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4CentralityReco::process_event(PHCompositeNode *topNode)
{
  _bimp = 101;
  auto event_header = findNode::getClass<EventHeaderv1>(topNode, "EventHeader");
  if (event_header)
  {
    _bimp = event_header->get_floatval("bimp");

    if (Verbosity() >= 5)
      std::cout << "PHG4CentralityReco::process_event : Hijing impact parameter b = " << _bimp << std::endl;
  }
  else
  {
    if (Verbosity() >= 5)
      std::cout << "PHG4CentralityReco::process_event : No Hijing impact parameter info, setting b = 101" << std::endl;
  }

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

  if (_do_centrality_calibration)
  {
    // sEPD centrality
    float low_epd_val = -10000;
    float high_epd_val = 10000;
    int low_epd_centile = -1;
    int high_epd_centile = -1;

    for (std::map<float, int>::iterator it = _cent_cal_epd.begin(); it != _cent_cal_epd.end(); ++it)
    {
      float signal = it->first;
      int cent = it->second;

      if (signal < _epd_NS && signal > low_epd_val)
      {
        low_epd_val = signal;
        low_epd_centile = cent;
      }
      if (signal > _epd_NS && signal < high_epd_val)
      {
        high_epd_val = signal;
        high_epd_centile = cent;
      }

    }  // close iterate through sEPD cuts

    if (low_epd_centile >= 0 && high_epd_centile >= 0)
    {
      _epd_cent = (low_epd_centile + high_epd_centile) / 2.0;
      if (Verbosity() >= 10)
        std::cout << "PHG4CentralityReco::process_event : lower EPD value is " << low_epd_val << " (" << low_epd_centile << "%), higher is " << high_epd_val << " (" << high_epd_centile << "%), assigning " << _epd_cent << "%" << std::endl;
    }
    else
    {
      _epd_cent = 101;
      if (Verbosity() >= 5)
        std::cout << "PHG4CentralityReco::process_event : not able to map EPD value to a centrality. debug info = " << low_epd_val << "/" << low_epd_centile << "/" << high_epd_val << "/" << high_epd_centile << std::endl;
    }

    // MBD centrality
    float low_mbd_val = -10000;
    float high_mbd_val = 10000;
    int low_mbd_centile = -1;
    int high_mbd_centile = -1;

    for (std::map<float, int>::iterator it = _cent_cal_mbd.begin(); it != _cent_cal_mbd.end(); ++it)
    {
      float signal = it->first;
      int cent = it->second;

      if (signal < _mbd_NS && signal > low_mbd_val)
      {
        low_mbd_val = signal;
        low_mbd_centile = cent;
      }
      if (signal > _mbd_NS && signal < high_mbd_val)
      {
        high_mbd_val = signal;
        high_mbd_centile = cent;
      }

    }  // close iterate through MBD cuts

    if (low_mbd_centile >= 0 && high_mbd_centile >= 0)
    {
      _mbd_cent = (low_mbd_centile + high_mbd_centile) / 2.0;
      if (Verbosity() >= 10)
        std::cout << "PHG4CentralityReco::process_event : lower MBD value is " << low_mbd_val << " (" << low_mbd_centile << "%), higher is " << high_mbd_val << " (" << high_mbd_centile << "%), assigning " << _mbd_cent << "%" << std::endl;
    }
    else
    {
      _mbd_cent = 101;
      if (Verbosity() >= 5)
        std::cout << "PHG4CentralityReco::process_event : not able to map MBD value to a centrality. debug info = " << low_mbd_val << "/" << low_mbd_centile << "/" << high_mbd_val << "/" << high_mbd_centile << std::endl;
    }

    // b (impact parameter) centrality
    float low_bimp_val = -1;
    float high_bimp_val = 10000;
    int low_bimp_centile = -1;
    int high_bimp_centile = -1;

    for (std::map<float, int>::iterator it = _cent_cal_bimp.begin(); it != _cent_cal_bimp.end(); ++it)
    {
      float signal = it->first;
      int cent = it->second;

      if (signal < _bimp && signal > low_bimp_val)
      {
        low_bimp_val = signal;
        low_bimp_centile = cent;
      }
      if (signal > _bimp && signal < high_bimp_val)
      {
        high_bimp_val = signal;
        high_bimp_centile = cent;
      }

    }  // close iterate through bimp cuts

    if (low_bimp_centile >= 0 && high_bimp_centile >= 0)
    {
      _bimp_cent = (low_bimp_centile + high_bimp_centile) / 2.0;
      if (Verbosity() >= 10)
        std::cout << "PHG4CentralityReco::process_event : lower b value is " << low_bimp_val << " (" << low_bimp_centile << "%), higher is " << high_bimp_val << " (" << high_bimp_centile << "%), assigning " << _bimp_cent << "%" << std::endl;
    }
    else
    {
      _bimp_cent = 101;
      if (Verbosity() >= 5)
        std::cout << "PHG4CentralityReco::process_event : not able to map b value to a centrality. debug info = " << low_bimp_val << "/" << low_bimp_centile << "/" << high_bimp_val << "/" << high_bimp_centile << std::endl;
    }

  }  // close centrality calibration

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
    cent->set_quantity(CentralityInfo::PROP::bimp, _bimp);

    cent->set_centile(CentralityInfo::PROP::epd_NS, _epd_cent);
    cent->set_centile(CentralityInfo::PROP::mbd_NS, _mbd_cent);
    cent->set_centile(CentralityInfo::PROP::bimp, _bimp_cent);
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
