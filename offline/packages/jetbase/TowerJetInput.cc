#include "TowerJetInput.h"

#include "Jet.h"
#include "Jetv2.h"

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/RawTowerDefs.h>           // for encode_towerid
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <phool/getClass.h>

#include <cassert>
#include <cmath>  // for asinh, atan2, cos, cosh
#include <iostream>
#include <map>      // for _Rb_tree_const_iterator
#include <utility>  // for pair
#include <vector>

TowerJetInput::TowerJetInput(Jet::SRC input)
  : m_input(input)
{
}

void TowerJetInput::identify(std::ostream &os)
{
  os << "   TowerJetInput: ";
  if (m_input == Jet::CEMC_TOWER)
  {
    os << "TOWER_CEMC to Jet::CEMC_TOWER";
  }
  else if (m_input == Jet::EEMC_TOWER)
  {
    os << "TOWER_EEMC to Jet::EEMC_TOWER";
  }
  else if (m_input == Jet::HCALIN_TOWER)
  {
    os << "TOWER_HCALIN to Jet::HCALIN_TOWER";
  }
  else if (m_input == Jet::HCALOUT_TOWER)
  {
    os << "TOWER_HCALOUT to Jet::HCALOUT_TOWER";
  }
  else if (m_input == Jet::FEMC_TOWER)
  {
    os << "TOWER_FEMC to Jet::FEMC_TOWER";
  }
  else if (m_input == Jet::FHCAL_TOWER)
  {
    os << "TOWER_FHCAL to Jet::FHCAL_TOWER";
  }
  os << std::endl;
}

std::vector<Jet *> TowerJetInput::get_input(PHCompositeNode *topNode)
{
  if (Verbosity() > 0) std::cout << "TowerJetInput::process_event -- entered" << std::endl;

  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexmap)
  {
    std::cout << "TowerJetInput::get_input - Fatal Error - GlobalVertexMap node is missing. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << std::endl;
    assert(vertexmap);  // force quit

    return std::vector<Jet *>();
  }

  if (vertexmap->empty())
  {
    std::cout << "TowerJetInput::get_input - Fatal Error - GlobalVertexMap node is empty. Please turn on the do_bbc or tracking reco flags in the main macro in order to reconstruct the global vertex." << std::endl;
    return std::vector<Jet *>();
  }
  m_use_towerinfo = false;

  /* std::string name =(m_input == Jet::CEMC_TOWER ? "CEMC_TOWER" */
  /*                        : m_input == Jet::CEMC_TOWERINFO ? "CEMC_TOWERINFO" */
  /*                        : m_input == Jet::EEMC_TOWER ? "Jet::EEMC_TOWER" */ 
  /*                        : m_input == Jet::HCALIN_TOWER ? "Jet::HCALIN_TOWER" */ 
  /*                        : m_input == Jet::HCALIN_TOWERINFO ? "Jet::HCALIN_TOWERINFO" */ 
  /*                        : m_input == Jet::HCALOUT_TOWER ? "Jet::HCALOUT_TOWER" */ 
  /*                        : m_input == Jet::HCALOUT_TOWERINFO ? "Jet::HCALOUT_TOWERINFO" */ 
  /*                        : m_input == Jet::FEMC_TOWER ? "Jet::FEMC_TOWER" */ 
  /*                        : m_input == Jet::FHCAL_TOWER ? "Jet::FHCAL_TOWER" */ 
  /*                        : m_input == Jet::CEMC_TOWER_RETOWER ? "Jet::CEMC_TOWER_RETOWER" */ 
  /*                        : m_input == Jet::CEMC_TOWERINFO_RETOWER ? "Jet::CEMC_TOWERINFO_RETOWER" */ 
  /*                        : m_input == Jet::CEMC_TOWER_SUB1 ? "Jet::CEMC_TOWER_SUB1" */ 
  /*                        : m_input == Jet::CEMC_TOWERINFO_SUB1 ? "Jet::CEMC_TOWERINFO_SUB1" */ 
  /*                        : m_input == Jet::HCALIN_TOWER_SUB1 ? "Jet::HCALIN_TOWER_SUB1" */ 
  /*                        : m_input == Jet::HCALIN_TOWERINFO_SUB1 ? "Jet::HCALIN_TOWERINFO_SUB1" */ 
  /*                        : m_input == Jet::HCALOUT_TOWER_SUB1 ? "Jet::HCALOUT_TOWER_SUB1" */ 
  /*                        : m_input == Jet::HCALOUT_TOWERINFO_SUB1 ? "Jet::HCALOUT_TOWERINFO_SUB1" */ 
  /*                        : m_input == Jet::CEMC_TOWER_SUB1CS ? "Jet::CEMC_TOWER_SUB1CS" */ 
  /*                        : m_input == Jet::HCALIN_TOWER_SUB1CS ? "Jet::HCALIN_TOWER_SUB1CS" */ 
  /*                        : m_input == Jet::HCALOUT_TOWER_SUB1CS ? "Jet::HCALOUT_TOWER_SUB1CS" */
  /*                        : "NO NAME"); */
  /* std::cout << " TowerJetInput (" << name << ")" << std::endl; */

  RawTowerContainer *towers = nullptr;
  TowerInfoContainer *towerinfos = nullptr; 
  RawTowerGeomContainer *geom = nullptr;
  if (m_input == Jet::CEMC_TOWER)
  {
    towers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC");
    geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
    if ((!towers) || !geom)
    {
      return std::vector<Jet *>();
    }
  }
  else if (m_input == Jet::CEMC_TOWERINFO)
  {
    m_use_towerinfo = true;
    towerinfos = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
    geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
    geocaloid = RawTowerDefs::CalorimeterId::CEMC;
    if ((!towerinfos) || !geom)
      {
	return std::vector<Jet *>();
      }
  }
  else if (m_input == Jet::CEMC_TOWERINFO_EMBED)
  {
    m_use_towerinfo = true;
    towerinfos = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_EMBED_CEMC");
    geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
    geocaloid = RawTowerDefs::CalorimeterId::CEMC;
    if ((!towerinfos) || !geom)
      {
	return std::vector<Jet *>();
      }
  }
  else if (m_input == Jet::CEMC_TOWERINFO_SIM)
  {
    m_use_towerinfo = true;
    towerinfos = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_SIM_CEMC");
    geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
    geocaloid = RawTowerDefs::CalorimeterId::CEMC;
    if ((!towerinfos) || !geom)
      {
	return std::vector<Jet *>();
      }
  }
  else if (m_input == Jet::EEMC_TOWER)
    {
      towers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_EEMC");
      geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_EEMC");
      if ((!towers && !towerinfos) || !geom)
	{
	  return std::vector<Jet *>();
	}
    }
  else if (m_input == Jet::HCALIN_TOWER)
    {
      towers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALIN");
      geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
      if ((!towers) || !geom)
	{
	  return std::vector<Jet *>();
	}
    }
  else if (m_input == Jet::HCALIN_TOWERINFO)
    {
      m_use_towerinfo = true;
      towerinfos = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
      geocaloid = RawTowerDefs::CalorimeterId::HCALIN;
      geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
      if ((!towerinfos) || !geom)
	{
	  return std::vector<Jet *>();
	}
    }
  else if (m_input == Jet::HCALIN_TOWERINFO_EMBED)
  {
    m_use_towerinfo = true;
    towerinfos = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_EMBED_HCALIN");
    geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
    geocaloid = RawTowerDefs::CalorimeterId::HCALIN;
    if ((!towerinfos) || !geom)
      {
	return std::vector<Jet *>();
      }
  }
  else if (m_input == Jet::HCALIN_TOWERINFO_SIM)
  {
    m_use_towerinfo = true;
    towerinfos = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_SIM_HCALIN");
    geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
    geocaloid = RawTowerDefs::CalorimeterId::HCALIN;
    if ((!towerinfos) || !geom)
      {
	return std::vector<Jet *>();
      }
  }
  else if (m_input == Jet::HCALOUT_TOWER)
    {
      towers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALOUT");
      geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
      if ((!towers) || !geom)
	{
	  return std::vector<Jet *>();
	}
    }
  else if (m_input == Jet::HCALOUT_TOWERINFO)
    {
      m_use_towerinfo = true;
      towerinfos = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
      geocaloid = RawTowerDefs::CalorimeterId::HCALOUT;
      geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
      if ((!towerinfos) || !geom)
	{
	  return std::vector<Jet *>();
	}
    }
  else if (m_input == Jet::HCALOUT_TOWERINFO_EMBED)
  {
    m_use_towerinfo = true;
    towerinfos = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_EMBED_HCALOUT");
    geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
    geocaloid = RawTowerDefs::CalorimeterId::HCALOUT;
    if ((!towerinfos) || !geom)
      {
	return std::vector<Jet *>();
      }
  }
  else if (m_input == Jet::HCALOUT_TOWERINFO_SIM)
  {
    m_use_towerinfo = true;
    towerinfos = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_SIM_HCALOUT");
    geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
    geocaloid = RawTowerDefs::CalorimeterId::HCALOUT;
    if ((!towerinfos) || !geom)
      {
	return std::vector<Jet *>();
      }
  }


  else if (m_input == Jet::FEMC_TOWER)
    {
      towers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_FEMC");
      geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_FEMC");
      if ((!towers) || !geom)
	{
	  return std::vector<Jet *>();
	}
    }
  else if (m_input == Jet::FHCAL_TOWER)
  {
    towers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_FHCAL");
    geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_FHCAL");
    if ((!towers) || !geom)
    {
      return std::vector<Jet *>();
    }
  }
  else if (m_input == Jet::CEMC_TOWER_RETOWER)
  {
    towers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC_RETOWER");
    geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
    if ((!towers) || !geom)
    {
      return std::vector<Jet *>();
    }
  }
  else if (m_input == Jet::CEMC_TOWERINFO_RETOWER)
    {
      m_use_towerinfo = true;
      towerinfos = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
      geocaloid = RawTowerDefs::CalorimeterId::HCALIN;
      geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
      if ((!towerinfos) || !geom)
	{
	  return std::vector<Jet *>();
	}
    }
  else if (m_input == Jet::CEMC_TOWER_SUB1)
  {
    towers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC_RETOWER_SUB1");
    geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
    if ((!towers) || !geom)
    {
      return std::vector<Jet *>();
    }
  }
  else if (m_input == Jet::CEMC_TOWERINFO_SUB1)
    {
      m_use_towerinfo = true;
      towerinfos = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER_SUB1");
      geocaloid = RawTowerDefs::CalorimeterId::HCALIN;
      geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
      if ((!towerinfos) || !geom)
	{
	  return std::vector<Jet *>();
	}
    }
  else if (m_input == Jet::HCALIN_TOWER_SUB1)
  {  
    towers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALIN_SUB1");
    geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
    if ((!towers) || !geom)
    {
      return std::vector<Jet *>();
    }
  }
  else if (m_input == Jet::HCALIN_TOWERINFO_SUB1)
    {
      m_use_towerinfo = true;
      towerinfos = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN_SUB1");
      geocaloid = RawTowerDefs::CalorimeterId::HCALIN;
      geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
      if ((!towerinfos) || !geom)
	{
	  return std::vector<Jet *>();
	}
    }
  else if (m_input == Jet::HCALOUT_TOWER_SUB1)
  {
    towers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALOUT_SUB1");
    geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
    if ((!towers) || !geom)
    {
      return std::vector<Jet *>();
    }
  }
  else if (m_input == Jet::HCALOUT_TOWERINFO_SUB1)
    {
      m_use_towerinfo = true;
      towerinfos = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT_SUB1");
      geocaloid = RawTowerDefs::CalorimeterId::HCALOUT;
      geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
      if ((!towerinfos) || !geom)
	{
	  return std::vector<Jet *>();
	}
    }
  else if (m_input == Jet::CEMC_TOWER_SUB1CS)
    {
      towers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC_RETOWER_SUB1CS");
      geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
      if ((!towers) || !geom)
	{
	  return std::vector<Jet *>();
	}
    }
  else if (m_input == Jet::HCALIN_TOWER_SUB1CS)
    {
      towers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALIN_SUB1CS");
      geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
      if ((!towers) || !geom)
	{
	  return std::vector<Jet *>();
	}
    }
  else if (m_input == Jet::HCALOUT_TOWER_SUB1CS)
    {
      towers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALOUT_SUB1CS");
      geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
      if ((!towers) || !geom)
	{
	  return std::vector<Jet *>();
	}
    }
  else
    {
      return std::vector<Jet *>();
    }

  // first grab the event vertex or bail
  GlobalVertex *vtx = vertexmap->begin()->second;
  float vtxz = NAN;
  if (vtx)
  {
    vtxz = vtx->get_z();
  }
  else
  {
    return std::vector<Jet *>();
  }

  if (std::isnan(vtxz))
  {
    static bool once = true;
    if (once)
    {
      once = false;

      std::cout << "TowerJetInput::get_input - WARNING - vertex is NAN. Drop all tower inputs (further NAN-vertex warning will be suppressed)." << std::endl;
    }

    return std::vector<Jet *>();
  }

  std::vector<Jet *> pseudojets;
  if (m_use_towerinfo)
    {
      if (!towerinfos)
	{
	  return std::vector<Jet *>();
	}

      unsigned int nchannels = towerinfos->size();
      for (unsigned int channel = 0; channel < nchannels;channel++)
	{
	  TowerInfo *tower = towerinfos->get_tower_at_channel(channel);
	  assert(tower);

	  unsigned int calokey = towerinfos->encode_key(channel);
	  int ieta = towerinfos->getTowerEtaBin(calokey);
	  int iphi = towerinfos->getTowerPhiBin(calokey);
	  const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(geocaloid, ieta, iphi);
	  //skip masked towers
	  if (tower->get_time()==-10 || tower->get_time()==-11) continue;
	  RawTowerGeom *tower_geom = geom->get_tower_geometry(key);
	  assert(tower_geom);

	  double r = tower_geom->get_center_radius();
	  double phi = atan2(tower_geom->get_center_y(), tower_geom->get_center_x());
	  double z0 = tower_geom->get_center_z();
	  double z = z0 - vtxz;
	  double eta = asinh(z / r);  // eta after shift from vertex
	  double pt = tower->get_energy() / cosh(eta);
	  if (tower->get_energy() == NAN) {pt = 0/cosh(eta);}
	  double px = pt * cos(phi);
	  double py = pt * sin(phi);
	  double pz = pt * sinh(eta);

	  Jet *jet = new Jetv2();
	  jet->set_px(px);
	  jet->set_py(py);
	  jet->set_pz(pz);
	  jet->set_e(tower->get_energy());
	  jet->insert_comp(m_input,channel);
	  pseudojets.push_back(jet);
	}
    }
  else
    {
      RawTowerContainer::ConstRange begin_end = towers->getTowers();
      RawTowerContainer::ConstIterator rtiter;
      for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
	{
	  RawTower *tower = rtiter->second;
	  
	  RawTowerGeom *tower_geom = geom->get_tower_geometry(tower->get_key());
	  assert(tower_geom);
	  
	  double r = tower_geom->get_center_radius();
	  double phi = atan2(tower_geom->get_center_y(), tower_geom->get_center_x());
	  double z0 = tower_geom->get_center_z();
	  double z = z0 - vtxz;
	  double eta = asinh(z / r);  // eta after shift from vertex
	  double pt = tower->get_energy() / cosh(eta);
	  double px = pt * cos(phi);
	  double py = pt * sin(phi);
	  double pz = pt * sinh(eta);
	  
	  Jet *jet = new Jetv2();
	  jet->set_px(px);
	  jet->set_py(py);
	  jet->set_pz(pz);
	  jet->set_e(tower->get_energy());
	  jet->insert_comp(m_input, tower->get_id());
	  pseudojets.push_back(jet);
	}
    }
  if (Verbosity() > 0) std::cout << "TowerJetInput::process_event -- exited" << std::endl;

  return pseudojets;
}
