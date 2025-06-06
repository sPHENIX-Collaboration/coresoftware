#include "TowerJetInput.h"

#include "Jet.h"
#include "Jetv2.h"

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>  // for encode_towerid
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <phool/getClass.h>

#include <cassert>
#include <cmath>  // for asinh, atan2, cos, cosh
#include <iostream>
#include <map>      // for _Rb_tree_const_iterator
#include <utility>  // for pair
#include <vector>

TowerJetInput::TowerJetInput(Jet::SRC input, const std::string &prefix)
  : m_input(input)
  , m_towerNodePrefix(prefix)
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
  if (Verbosity() > 0)
  {
    std::cout << "TowerJetInput::process_event -- entered" << std::endl;
  }
  float vtxz = 0;  // default to 0
  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexmap)
  {
    std::cout << "TowerJetInput::get_input - Fatal Error - GlobalVertexMap node is missing. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << std::endl;
    assert(vertexmap);  // force quit

    return std::vector<Jet *>();
  }
  if (vertexmap->empty())
  {
    if (Verbosity() > 0)
    {
      std::cout << "TowerJetInput::get_input - empty vertex map, continuing as if zvtx = 0" << std::endl;
    }
  }
  else
  {
    GlobalVertex *vtx = vertexmap->begin()->second;
    if (vtx)
    {
      if (m_use_vertextype)
      {
        auto typeStartIter = vtx->find_vertexes(m_vertex_type);
        auto typeEndIter = vtx->end_vertexes();
        for (auto iter = typeStartIter; iter != typeEndIter; ++iter)
        {
          const auto &[type, vertexVec] = *iter;
          if (type != m_vertex_type)
          {
            continue;
          }
          for (const auto *vertex : vertexVec)
          {
            if (!vertex)
            {
              continue;
            }
            vtxz = vertex->get_z();
          }
        }
      }
      else
      {
        vtxz = vtx->get_z();
      }
    }
  }
  if (std::isnan(vtxz))
  {
    static bool once = true;
    if (once)
    {
      once = false;
      std::cout << "TowerJetInput::get_input - WARNING - vertex is NAN. Continue with zvtx = 0 (further vertex warning will be suppressed)." << std::endl;
    }
    vtxz = 0;
  }

  if (std::abs(vtxz) > 1e3)  // code crashes with very large z vertex, so skip these events
  {
    static bool once = true;
    if (once)
    {
      once = false;

      std::cout << "TowerJetInput::get_input - WARNING - vertex is " << vtxz << ". Set vtxz = 0 (further vertex warning will be suppressed)." << std::endl;
    }
    vtxz = 0;
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
  RawTowerGeomContainer *EMCal_geom = nullptr;

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
    towerName = m_towerNodePrefix + "_CEMC";
    towerinfos = findNode::getClass<TowerInfoContainer>(topNode, towerName);
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
    towerName = m_towerNodePrefix + "_EMBED_CEMC";
    towerinfos = findNode::getClass<TowerInfoContainer>(topNode, towerName);
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
    towerName = m_towerNodePrefix + "_SIM_CEMC";
    towerinfos = findNode::getClass<TowerInfoContainer>(topNode, towerName);
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
    towerName = m_towerNodePrefix + "_HCALIN";
    towerinfos = findNode::getClass<TowerInfoContainer>(topNode, towerName);
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
    towerName = m_towerNodePrefix + "_EMBED_HCALIN";
    towerinfos = findNode::getClass<TowerInfoContainer>(topNode, towerName);
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
    towerName = m_towerNodePrefix + "_SIM_HCALIN";
    towerinfos = findNode::getClass<TowerInfoContainer>(topNode, towerName);
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
    towerName = m_towerNodePrefix + "_HCALOUT";
    towerinfos = findNode::getClass<TowerInfoContainer>(topNode, towerName);
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
    towerName = m_towerNodePrefix + "_EMBED_HCALOUT";
    towerinfos = findNode::getClass<TowerInfoContainer>(topNode, towerName);
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
    towerName = m_towerNodePrefix + "_SIM_HCALOUT";
    towerinfos = findNode::getClass<TowerInfoContainer>(topNode, towerName);
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
    towerName = m_towerNodePrefix + "_CEMC_RETOWER";
    towerinfos = findNode::getClass<TowerInfoContainer>(topNode, towerName);
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
    towerName = m_towerNodePrefix + "_CEMC_RETOWER_SUB1";
    towerinfos = findNode::getClass<TowerInfoContainer>(topNode, towerName);
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
    towerName = m_towerNodePrefix + "_HCALIN_SUB1";
    towerinfos = findNode::getClass<TowerInfoContainer>(topNode, towerName);
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
    towerName = m_towerNodePrefix + "_HCALOUT_SUB1";
    towerinfos = findNode::getClass<TowerInfoContainer>(topNode, towerName);
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

  // for those cases we need to use the EMCal R and IHCal eta phi to calculate the vertex correction
  if (m_input == Jet::CEMC_TOWER_RETOWER || m_input == Jet::CEMC_TOWERINFO_RETOWER || m_input == Jet::CEMC_TOWER_SUB1 || m_input == Jet::CEMC_TOWERINFO_SUB1 || m_input == Jet::CEMC_TOWER_SUB1CS)
  {
    EMCal_geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
    if (!EMCal_geom)
    {
      return std::vector<Jet *>();
    }
  }

  // first grab the event vertex or bail

  std::vector<Jet *> pseudojets;
  if (m_use_towerinfo)
  {
    if (!towerinfos)
    {
      return std::vector<Jet *>();
    }

    unsigned int nchannels = towerinfos->size();
    for (unsigned int channel = 0; channel < nchannels; channel++)
    {
      TowerInfo *tower = towerinfos->get_tower_at_channel(channel);
      assert(tower);

      unsigned int calokey = towerinfos->encode_key(channel);
      int ieta = towerinfos->getTowerEtaBin(calokey);
      int iphi = towerinfos->getTowerPhiBin(calokey);
      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(geocaloid, ieta, iphi);
      // skip masked towers
      if (tower->get_isHot() || tower->get_isNoCalib() || tower->get_isNotInstr() || tower->get_isBadChi2())
      {
        continue;
      }
      if (std::isnan(tower->get_energy()))
      {
        continue;
      }
      RawTowerGeom *tower_geom = geom->get_tower_geometry(key);
      assert(tower_geom);

      double r = tower_geom->get_center_radius();
      if (m_input == Jet::CEMC_TOWER_RETOWER || m_input == Jet::CEMC_TOWERINFO_RETOWER || m_input == Jet::CEMC_TOWER_SUB1 || m_input == Jet::CEMC_TOWERINFO_SUB1 || m_input == Jet::CEMC_TOWER_SUB1CS)
      {
        const RawTowerDefs::keytype EMCal_key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, 0, 0);
        RawTowerGeom *EMCal_tower_geom = EMCal_geom->get_tower_geometry(EMCal_key);
        assert(EMCal_tower_geom);
        r = EMCal_tower_geom->get_center_radius();
      }
      double phi = atan2(tower_geom->get_center_y(), tower_geom->get_center_x());
      double towereta = tower_geom->get_eta();
      double z0 = sinh(towereta) * r;
      double z = z0 - vtxz;
      double eta = asinh(z / r);  // eta after shift from vertex
      double pt = tower->get_energy() / cosh(eta);
      double e = tower->get_energy();
      double px = pt * cos(phi);
      double py = pt * sin(phi);
      double pz = pt * sinh(eta);

      Jet *jet = new Jetv2();
      jet->set_px(px);
      jet->set_py(py);
      jet->set_pz(pz);
      jet->set_e(e);
      jet->insert_comp(m_input, channel);
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
      if (m_input == Jet::CEMC_TOWER_RETOWER || m_input == Jet::CEMC_TOWERINFO_RETOWER || m_input == Jet::CEMC_TOWER_SUB1 || m_input == Jet::CEMC_TOWERINFO_SUB1 || m_input == Jet::CEMC_TOWER_SUB1CS)
      {
        const RawTowerDefs::keytype EMCal_key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, 0, 0);
        RawTowerGeom *EMCal_tower_geom = EMCal_geom->get_tower_geometry(EMCal_key);
        assert(EMCal_tower_geom);
        r = EMCal_tower_geom->get_center_radius();
      }
      double phi = atan2(tower_geom->get_center_y(), tower_geom->get_center_x());
      double towereta = tower_geom->get_eta();
      double z0 = sinh(towereta) * r;
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
  if (Verbosity() > 0)
  {
    std::cout << "TowerJetInput::process_event -- exited" << std::endl;
  }
  return pseudojets;
}
