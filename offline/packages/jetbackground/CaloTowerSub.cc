#include "CaloTowerSub.h"

#include "TowerBackground.h"
#include "TowerBackgroundv1.h"

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <eventplaneinfo/Eventplaneinfo.h>
#include <eventplaneinfo/EventplaneinfoMap.h>

#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TLorentzVector.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <utility>
#include <vector>
#include <set>

CaloTowerSub::CaloTowerSub(const std::string &name)
  : SubsysReco(name)
{
  _UE.resize(3, std::vector<float>(1, 0));
}

CaloTowerSub::~CaloTowerSub()
{

  for (auto & input : m_inputs)
  {
    delete input;
  }
  m_inputs.clear();
  
  if (m_seed_algo)
  {
    delete m_seed_algo;
  }

}

int CaloTowerSub::InitRun( PHCompositeNode *topNode )
{
  return CreateNode(topNode);
}

int CaloTowerSub::getZvertex( PHCompositeNode *topNode )
{

  if ( Verbosity() )
  {
    std::cout << PHWHERE \
      << "Getting Z vertex from GlobalVertexMap node" << std::endl;
  }

  m_zvrtx = 0.0;
  auto vrtxmap = findNode::getClass< GlobalVertexMap >( topNode, "GlobalVertexMap");
  if ( !vrtxmap )
  {
    std::cerr << PHWHERE \
      << "GlobalVertexMap node not found, fatal error, exiting." << std::endl;
    exit(1);
  }

  if ( vrtxmap->empty() && Verbosity() )
  {
    std::cout << PHWHERE \
      << "empty vertex map, continuing as if zvtx = 0" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;

  }

  auto vrtx = vrtxmap->begin()->second;
  if ( !vrtx && Verbosity() )
  {
    std::cout << PHWHERE \
      << "GlobalVertex is null, continuing as if zvtx = 0" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // if you're getting picky
  if ( m_use_vertextype )
  {
    // loop over the container of vertices within the GlobalVertex container within the GlobalVertexMap within ...
    auto start_here = vrtx -> find_vertexes( m_vertex_type );
    for ( auto iter = start_here; 
          iter != vrtx -> end_vertexes(); 
          ++iter )
    {
      const auto& [type, it_never_ends] = *iter;
      if ( type != m_vertex_type )
      {
        continue;
      }
      
      for ( const auto * i_told_you : it_never_ends )
      {
        if ( !i_told_you )
        {
          continue;
        }
        // at last
        m_zvtx = i_told_you -> get_z();
      }
    } // end loop over vertices of specified type
    else
    {
      // this was always an option
      m_zvtx = vrtx -> get_z();
    }
   
  } // you finally found the z vertex!

  // you though we were done?
  if ( std::isnan(m_zvtx)  || std::abs(m_zvtx) > 1e3 )
  {
    static bool once = true;
    if ( once )
    { // we don't evn care enough to log multiple times
      once = false;
      std::cout << PHWHERE \
        << "vertex is either NAN or too large so I'm setting it to 0 and I'll never tell you if it happens again." \
        << "Continue with zvtx = 0 (further vertex warning will be suppressed)." \
        << std::endl;
    }
    m_zvtx = 0;
  }


  // see ya
  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloTowerSub::genSeedMask( PHCompositeNode *topNode )
{
  if ( Verbosity() > 0)
  {
    std::cout << "CaloTowerSub::genSeedMask: generating seed mask with _seed_type = " << _seed_type << std::endl;
  
    if (_seed_type == 0)
    {
      std::cout << "CaloTowerSub::genSeedMask: 1st iteration params: \n";
      std::cout << "  D = " << _seed_jet_D << std::endl;
    }
    else if (_seed_type == 1)
    {
      std::cout << "CaloTowerSub::genSeedMask: using pT = " << _seed_jet_pt << std::endl;
    }
    else
    {
      std::cout << "CaloTowerSub::genSeedMask: UNDEFINED seed behavior! Exiting. " << std::endl;
      exit(1);
    }
  
  
  }

  std::string jetnode = "AntiKt_TowerInfo_HIRecoSeedsRaw_r02";
  if (_seed_type == 1)
  {
    jetnode = "AntiKt_TowerInfo_HIRecoSeedsSub_r02";
  }

  
  auto jets = findNode::getClass<JetContainer>(topNode, jetnode);
  if (! jets )
  {
    std::cerr << PHWHERE \
      << "CaloTowerSub::genSeedMask: cannot find jet node " << jetnode << ", exiting." \
      << std::endl;
    exit(1);
  }

  _index_SeedD = jets->property_index(Jet::PROPERTY::prop_SeedD);
  _index_SeedItr = jets->property_index(Jet::PROPERTY::prop_SeedItr);

  // reset seed vectors
  _seed_eta.resize(0);
  _seed_phi.resize(0);

  for ( auto * this_jet : * jets )
  {

    float this_pt = this_jet->get_pt();
    float this_phi = this_jet->get_phi();
    float this_eta = this_jet->get_eta();

    if ( this_pt < _seed_jet_pt[_seed_type] )
    {

      // if this is the first iteration, setting this to 0 will
      // mark that this jet was not selected as a seed (and did not have D determined)
      // if this is the second iteration, setting this to 0 will
      // mark that this jet was considered but not used as a seed

      this_jet->set_property(_index_SeedItr, 0.0);
      if ( _seed_type == 0 )
      {
        // only zero this out if we haven't calculated D yet
        this_jet->set_property(_index_SeedD, 0);
      }

      continue;
    }

    if ( _seed_type == 1 )
    {
      // seed type 1 is the set of those jets above which, when their
      // kinematics are updated for the first background subtraction, have
      // pT > _seed_jet_pt[1]
      _seed_eta.push_back(this_eta);
      _seed_phi.push_back(this_phi);
      this_jet->set_property(_index_SeedItr, 2.0); // mark as seed used in 2nd iteration

      if ( Verbosity() > 1 )
      {
        std::cout << PHWHERE \
          << "Adding seed jet with pT = " << this_pt \
          << ", eta = " << this_eta \
          << ", phi = " << this_phi \
        << std::endl;
      }
      
      // done with second iteration seeds
      continue;
    }

    // we will only get here in seed type 0
    if ( _seed_type != 0 )
    {
      std::cerr << PHWHERE \
        << "Seeding iteration 1 somehow got into seed type 0 loop, exiting." \
        << std::endl;
      exit(1);
    }

    std::map<int, double> constituent_eT_sum {}; // map of constituent tower ET sums, key = (ieta*1000 + iphi)
    for ( const auto & comp : this_jet->get_comp_vec() )
    {

      TowerInfoContainer * towerinfo = nullptr;
      RawTowerGeomContainer * geom = nullptr;
      RawTowerDefs::CalorimeterId geocaloid { RawTowerDefs::CalorimeterId::NONE };
      bool is_retowered = false;
  
      if ( comp.first == Jet::SRC::HCALIN_TOWER || comp.first == Jet::SRC::HCALIN_TOWERINFO )
      {
        towerinfo = towerinfosIH3;
        geom = geomIH;
        geocaloid = RawTowerDefs::CalorimeterId::HCALIN;
      }
      else if ( comp.first == Jet::SRC::HCALOUT_TOWER || comp.first == Jet::SRC::HCALOUT_TOWERINFO )
      {
        towerinfo = towerinfosOH3;
        geom = geomOH;
        geocaloid = RawTowerDefs::CalorimeterId::HCALOUT;
      }
      else if ( comp.first == Jet::SRC::CEMC_TOWER_RETOWER || comp.first == Jet::SRC::CEMC_TOWERINFO_RETOWER )
      {
        towerinfo = towerinfosEM3;
        geom = geomIH;
        geocaloid = RawTowerDefs::CalorimeterId::HCALIN;
        is_retowered = true;
      }
      
      if ( !towerinfo  )
      {
        if ( Verbosity() > 0 )
        {
          std::cout << PHWHERE \ 
            << "Unrecognized jet constituent source " << comp.first << ", skipping constituent." \ 
            << std::endl;
        }

        continue;

      }
      if ( !geom )
      {
        if ( Verbosity() > 0 )
        {
          std::cout << PHWHERE \ 
            << "No geometry found for constituent source " << comp.first \ 
            << ", skipping constituent." \ 
            << std::endl;
        }

        continue;

      }

      auto tower = towerinfo -> get_tower_at_channel( comp.second );
      if ( !tower ) 
      {
        if ( Verbosity() > 0 )
        {
          std::cout << PHWHERE \ 
            << "No tower found for constituent source " << comp.first \ 
            << " at channel " << comp.second << ", skipping constituent." \ 
            << std::endl;
        }

        continue;

      }

      // check for bad towers
      if ( tower -> get_isHot() 
          || tower -> get_isNoCalib() 
          || tower -> get_isNotInstr() 
          || tower -> get_isBadChi2() 
      )
      { // all for not
        if ( Verbosity() > 4 )
        {
          std::cout << PHWHERE \ 
            << "Constituent source " << comp.first \ 
            << " at channel " << comp.second \ 
            << " is marked bad (hot/no calib/not instr/bad chi2), skipping constituent." \ 
          << std::endl;
        }
        continue;
      }
      
      // get ieta and iphi
      unsigned int towerkey = towerinfo -> encode_key( comp.second );
      int ieta = towerinfo -> getTowerEtaBin( towerkey );
      int iphi = towerinfo -> getTowerPhiBin( towerkey );
      int comp_supkey = (1000 * ieta) + iphi;

      if ( Verbosity() > 4 )
      {
        std::cout << PHWHERE \ 
          << "Constituent source " << comp.first \ 
          << " at channel " << comp.second \ 
          << " has ieta / iphi = " << ieta << " / " << iphi \ 
          << " and superkey = " << comp_supkey << "." \ 
        << std::endl;
      }

      auto tower_geom = geom -> get_tower_geometry( 
        RawTowerDefs::encode_towerid(
          geocaloid, 
          ieta, 
          iphi
        )
      );
      if ( !tower_geom ) 
      {
        if ( Verbosity() > 0 )
        {
          std::cout << PHWHERE \ 
            << "No tower geometry found for constituent source " << comp.first \ 
            << " at ieta / iphi = " << ieta << " / " << iphi << ", skipping constituent." \ 
            << std::endl;
        }

        continue;
      }

      double this_comp_eta = tower_geom->get_eta();
      double calo_radius = tower_geom->get_center_radius();
      if ( is_retowered ) // tricky tricky
      {
        auto cemc_geo = geomEM ->get_tower_geometry( 
          RawTowerDefs::encode_towerid(
            RawTowerDefs::CalorimeterId::CEMC, 
            0, 
            0
          )
        );
        assert(cemc_geo);
        calo_radius = cemc_geo->get_center_radius();
      }
      
      // correct eta for vertex
      correct_calo_eta(this_comp_eta, calo_radius);

      double comp_eT = tower->get_energy() / cosh(this_comp_eta);
      constituent_eT_sum[comp_supkey] += comp_eT;
      
      // add to constituent ET sum map
      if ( Verbosity() > 4 )
      {
        std::cout << PHWHERE \ 
          << "Added eT = " << comp_eT \
          << " to constituent ET sum for key " << comp_supkey \
          << ", new sum = " << constituent_eT_sum[comp_supkey] \
        << std::endl;
      }
      
      

    } // end loop over jet constituents
      
    // now iterate over constituent_ET sums to find maximum and mean
    float max_super_eT = 0, sum_super_eT = 0;
    int n_super_towers = 0;
    for ( const auto & map_iter : constituent_eT_sum )
    {
      n_super_towers++;
      sum_super_eT += map_iter.second;
      if ( map_iter.second > max_super_eT )
      {
        max_super_eT = map_iter.second;
      }
    } 
    if ( n_super_towers == 0 )
    {
      if ( Verbosity() > 0 )
      {
        std::cout << PHWHERE \ 
          << "Jet at eta / phi = " << this_eta << " / " << this_phi \ 
          << " has no valid constituents, skipping D calculation." \ 
          << std::endl;
      }
      
      this_jet->set_property(_index_SeedItr, 0);
      this_jet->set_property(_index_SeedD, 0);
      
      continue;
    }

    float avg_super_eT = sum_super_eT / n_super_towers;
    float seed_D = max_super_eT / avg_super_eT;

    this_jet->set_property(_index_SeedD, seed_D);
    if (Verbosity() > 3)
    {
      std::cout << PHWHERE \ 
        << "Jet at eta / phi = " << this_eta << " / " << this_phi \ 
        << " has max eT = " << max_super_eT \ 
        << ", avg eT = " << avg_super_eT \ 
        << ", and D = " << seed_D \
      << std::endl;
    }

    if ( seed_D > _seed_jet_D 
      && max_super_eT > _seed_max_const // this will be max_super_eT > 0 if not set
    ) 
    { // accept seed

      _seed_eta.push_back(this_eta);
      _seed_phi.push_back(this_phi);

      this_jet->set_property(_index_SeedItr, 1.0);

      if (Verbosity() > 1)
      {
        std::cout << PHWHERE << \
          "CaloTowerSub::process_event: --> adding seed at eta / phi = " << this_eta << " / " << this_phi \ 
          << " ( R=0.2 jet with pt = " << this_pt << ", D = " << seed_D << " ) " \
        << std::endl;
      }

    }
    else
    {
      // mark that this jet was considered but not used as a seed
      this_jet->set_property(_index_SeedItr, 0.0);

      if (Verbosity() > 3)
      {
        std::cout << PHWHERE << \
          "CaloTowerSub::process_event: --> discarding potential seed at eta / phi = " << this_eta << " / " << this_phi \ 
          << " ( R=0.2 jet with pt = " << this_pt << ", D = " << seed_D << " ) " \
        << std::endl;
      }

    } 

  } // end loop over jets

  // eta phi of accepted seeds are now in _seed_eta and _seed_phi corresponding to whichever _seed_type this is

  // fill the masks to exclude R = 0.4 around each seed in subsequent background calculations
  if ( Verbosity() > 0 )
  {
    std::cout << PHWHERE << \
      " --> found " << _seed_eta.size() \
      << " seeds for seed type " \
      << _seed_type \
    << std::endl;
  }


  // clear calo tower seed masks
  _CALO_TOWER_MASK.assign( _N_LAYERS, std::vector<std::vector<int>>( _HCAL_NETA, std::vector<int>( _HCAL_NPHI, 0 ) ) );
  // loop over seeds to fill masks
  for ( size_t iseed = 0; iseed < _seed_eta.size(); ++iseed )
  {
    float this_eta = _seed_eta.at(iseed);
    float this_phi = _seed_phi.at(iseed);

    // loop over all towers to find those within R = 0.4 of this seed
    for ( int ilayer = 0; ilayer < _N_LAYERS; ++ilayer )
    {
      for ( int ieta = 0; ieta < _HCAL_NETA; ++ieta )
      {
        for ( int iphi = 0; iphi < _HCAL_NPHI; ++iphi )
        {
           
          // mask identified bad towers as well
          bool is_bad = getCaloStatus( ilayer, ieta, iphi );
          if ( is_bad )
          {
            _CALO_TOWER_MASK[ilayer][ieta][iphi] = 1;
            continue; // easy peasy
          }

          float tower_eta = getCaloEta( ilayer, ieta, iphi );
          float tower_phi = getCaloPhi( ilayer, ieta, iphi );

          float dR = deltaR( this_eta, this_phi, tower_eta, tower_phi );
          if ( dR < _seed_exclusion_dR )
          {
            // mark tower as excluded
            _CALO_TOWER_MASK[ilayer][ieta][iphi] = 1;

            if ( Verbosity() > 4 )
            {
              std::cout << PHWHERE << \
                " --> excluding tower at layer / ieta / iphi = " << ilayer << " / " << ieta << " / " << iphi \
                << " with eta / phi = " << tower_eta << " / " << tower_phi \
                << " ( deltaR = " << dR << " ) " \
              << std::endl;
            }

          }
        } // end loop over iphi
      } // end loop over ieta
    } // end loop over layers
  
  } // end loop over seeds

  return Fun4AllReturnCodes::EVENT_OK;

}

int CaloTowerSub::fillCaloArrays( PHCompositeNode *topNode )
{

  if ( Verbosity() > 0 )
  {
    std::cout << PHWHERE << \
      "fillCaloArrays: filling calo tower arrays from tower infos."
      << std::endl;
  }

  // get the binning from the geometry (different for 1D vs 2D...)
  if ( _HCAL_NETA < 0 ) // fisrt event
  {
    _HCAL_NETA = geomIH->get_etabins();
    _HCAL_NPHI = geomIH->get_phibins();
    
    // resize UE density and energy vectors
    _UE.resize(3 , std::vector<float>(_HCAL_NETA, 0));

    _CALO_E.resize( _N_LAYERS, std::vector<std::vector<float>>( _HCAL_NETA, std::vector<float>( _HCAL_NPHI, 0 ) ) );
    _CALO_ISBAD.resize( _N_LAYERS, std::vector<std::vector<int>>( _HCAL_NETA, std::vector<int>( _HCAL_NPHI, 0 ) ) );
    _CALO_POS.resize( _N_LAYERS, std::vector<std::vector<std::pair<float,float>>>( _HCAL_NETA, std::vector<std::pair<float,float>>( _HCAL_NPHI, std::make_pair(0.0,0.0) ) ) );

    // defualt set weights to 1.0 for all phi bins
    _CALO_WEIGHTS.resize( _N_LAYERS, std::vector<std::vector<float>>( _HCAL_NETA, std::vector<float>( _HCAL_NPHI, 1.0 ) ) );
    _CALO_FLOW_WEIGHTS.resize( _N_LAYERS, std::vector<std::vector<float>>( _HCAL_NETA, std::vector<float>( _HCAL_NPHI, 1.0 ) ) );

    // for flow determination, build up a 1-D phi distribution of
    // energies from all layers summed together, populated only from eta
    // strips which do not have any excluded phi towers
    _FULLCALOFLOW_PHI_E.resize(_HCAL_NPHI, 0);
    _FULLCALOFLOW_PHI_VAL.resize(_HCAL_NPHI, 0);

    if (Verbosity() > 0)
    {
      std::cout << PHWHERE << \
        "initialized calo tower arrays with _HCAL_NETA = " << _HCAL_NETA \
        << " and _HCAL_NPHI = " << _HCAL_NPHI << std::endl;    
    }
  }

  // reset all maps map
  _UE.assign(3, std::vector<float>(_HCAL_NETA, 0));
  _CALO_E.assign( _N_LAYERS, std::vector<std::vector<float>>( _HCAL_NETA, std::vector<float>( _HCAL_NPHI, 0 ) ) );
  _CALO_ISBAD.assign( _N_LAYERS, std::vector<std::vector<int>>( _HCAL_NETA, std::vector<int>( _HCAL_NPHI, 0 ) ) );
  _CALO_POS.assign( _N_LAYERS, std::vector<std::vector<std::pair<float,float>>>( _HCAL_NETA, std::vector<std::pair<float,float>>( _HCAL_NPHI, std::make_pair(0.0,0.0) ) ) );

  // fill calo energy and status arrays
  for ( int ilayer = 0; ilayer < _N_LAYERS; ++ilayer )
  {
    TowerInfoContainer * towerinfos = nullptr;
    RawTowerGeomContainer * geom = nullptr;
    RawTowerDefs::CalorimeterId geocaloid { RawTowerDefs::CalorimeterId::NONE };
    double calo_radius = 0.0;
    if ( ilayer == 0 )
    {
      towerinfos = towerinfosEM3;
      geom = geomIH;
      geocaloid = RawTowerDefs::CalorimeterId::HCALIN;
      // get EMCal radius from EMCal geometry
      auto tower_geom = geomEM -> get_tower_geometry( 
        RawTowerDefs::encode_towerid(
          RawTowerDefs::CalorimeterId::CEMC, 
          0, 
          0
        )
      );
      calo_radius = tower_geom->get_center_radius();
    }
    else if ( ilayer == 1 )
    {
      towerinfos = towerinfosIH3;
      geom = geomIH;
      geocaloid = RawTowerDefs::CalorimeterId::HCALIN;
      auto tower_geom = geom -> get_tower_geometry( 
        RawTowerDefs::encode_towerid(
          geocaloid, 
          0, 
          0
        )
      );
      calo_radius = tower_geom->get_center_radius();
    }
    else if ( ilayer == 2 )
    {
      towerinfos = towerinfosOH3;
      geom = geomOH;
      geocaloid = RawTowerDefs::CalorimeterId::HCALOUT;
      auto tower_geom = geom -> get_tower_geometry( 
        RawTowerDefs::encode_towerid(
          geocaloid, 
          0, 
          0
        )
      );
      calo_radius = tower_geom->get_center_radius();
    }
    else
    {
      std::cerr << PHWHERE \
        << "fillCaloArrays: invalid layer " << ilayer << ", exiting." \
        << std::endl;
      exit(1);
    }

    if ( !towerinfos )
    {
      std::cerr << PHWHERE \
        << "fillCaloArrays: null towerinfos for layer " << ilayer << ", exiting." \
        << std::endl;
      exit(1);
    }
    if ( !geom )
    {
      std::cerr << PHWHERE \
        << "fillCaloArrays: null geometry for layer " << ilayer << ", exiting." \
        << std::endl;
      exit(1);
    }


    for (unsigned int channel = 0; 
        channel < towerinfos->size(); 
        channel++
    )
    {
      
      unsigned int key = towerinfos->encode_key(channel);
      int this_etabin = towerinfos->getTowerEtaBin(key);
      int this_phibin = towerinfos->getTowerPhiBin(key);
      auto tower = towerinfos->get_tower_at_channel(channel);
      assert(tower);

      // get tower geometry
      auto tower_geom = geom -> get_tower_geometry( 
        RawTowerDefs::encode_towerid(
          geocaloid, 
          this_etabin, 
          this_phibin
        )
      );
      assert(tower_geom);
      float tower_eta = tower_geom->get_eta();
      float tower_phi = tower_geom->get_phi();
      
      // correct eta for vertex
      correct_calo_eta(tower_eta, calo_radius);
      _CALO_POS[ilayer][this_etabin][this_phibin] = std::make_pair(tower_eta, tower_phi);

      int this_isBad = tower->get_isHot() || tower->get_isNoCalib() || tower->get_isNotInstr() || tower->get_isBadChi2();
      _CALO_ISBAD[ilayer][this_etabin][this_phibin] = this_isBad;

      if ( !this_isBad )
      {
        float this_E = tower->get_energy();
        _CALO_E[ilayer][this_etabin][this_phibin] = this_E;
      }
      
    } // end loop over channels

  } // end loop over layers

}

int CaloTowerSub::getCaloNodes( PHCompositeNode *topNode )
{
  
  if ( Verbosity() > 0 )
  {
    std::cout << PHWHERE << \
      "getCaloNodes: getting calo tower info and geometry nodes."
      << std::endl;
  }

  // pull out the tower containers and geometry objects at the start
  TowerInfoContainer *towerinfosEM3 = findNode::getClass<TowerInfoContainer>(topNode, EMTowerName);
  TowerInfoContainer *towerinfosIH3 = findNode::getClass<TowerInfoContainer>(topNode, IHTowerName);
  TowerInfoContainer *towerinfosOH3 = findNode::getClass<TowerInfoContainer>(topNode, OHTowerName);
  if ( !towerinfosEM3 )
  {
    std::cerr << PHWHERE \
      << "CaloTowerSub::process_event: Cannot find node " << EMTowerName << ", exiting." \
      << std::endl;
    exit(1);
  }
  if ( !towerinfosIH3 )
  {
    std::cerr << PHWHERE \
      << "CaloTowerSub::process_event: Cannot find node " << IHTowerName << ", exiting." \
      << std::endl;
    exit(1);
  }
  if ( !towerinfosOH3 )
  {
    std::cerr << PHWHERE \
      << "CaloTowerSub::process_event: Cannot find node " << OHTowerName << ", exiting." \
      << std::endl;
    exit(1);
  }

  RawTowerGeomContainer *geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  RawTowerGeomContainer *geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  if ( !geomEM )
  {
    std::cerr << PHWHERE \
      << "CaloTowerSub::process_event: Cannot find node TOWERGEOM_CEMC, exiting." \
      << std::endl;
    exit(1);
  }
  if ( !geomIH )
  {
    std::cerr << PHWHERE \
      << "CaloTowerSub::process_event: Cannot find node TOWERGEOM_HCALIN, exiting." \
      << std::endl;
    exit(1);
  }
  if ( !geomOH )
  {
    std::cerr << PHWHERE \
      << "CaloTowerSub::process_event: Cannot find node TOWERGEOM_HCALOUT, exiting." \
      << std::endl;
    exit(1);
  }

  return Fun4AllReturnCodes::EVENT_OK;

}


int CaloTowerSub::process_event( PHCompositeNode *topNode )
{

  getZvertex( topNode );
  
  getCaloNodes( topNode ); 

  fillCaloArrays( topNode );

  loadCBDWeights( topNode ); // load avg E weights


  for ( int iter = 0; iter < 2; iter++ )
  {

    genSeedMask( topNode );

    calcFlow( topNode );

    calcEnergyDensity( topNode );

    if (iter == 0 )
    {
      subJets( topNode );
    }
    else if ( iter == 1 )
    {
      subTowers( topNode );
    }
  
  }

    return Fun4AllReturnCodes::EVENT_OK;
}

int CaloTowerSub::recoSeeds( PHCompositeNode *topNode )
{

  if ( Verbosity() > 0 )
  {
    std::cout << PHWHERE << \
      "recoSeeds: reconstructing seeds for iteration " << iter \
    << std::endl;
  }

  // grab tower jet input objects
  std::vector< Jet* > inputs {};
  for ( auto & input : m_inputs )
  {

    auto towerjets = input -> get_input( topNode );

    m_zvrtx = input -> get_zvertex();
    auto this_node_name = input -> get_input_node_name();
    auto this_caloid = input -> get_input_caloId();
    Jet::SRC this_source = input -> get_src();

    auto calo_radius = getCaloRadius( this_caloid , this_source );
    auto ilayer = getCaloLayer( this_caloid , this_source );

    auto this_geom = findNode::getClass<RawTowerGeomContainer>( 
      topNode, 
      getGeomName( this_caloid )
    ); 
    if ( !this_geom )
    {
      std::cerr << PHWHERE \
        << "can't find geometry node " << geom_name << ", exiting." \
      << std::endl;
      exit(1);
    }

    auto this_towerinfo = findNode::getClass<TowerInfoContainer>( 
      topNode, 
      this_node_name
    );
    if ( !this_towerinfo )
    {
      std::cerr << PHWHERE \
        << "can't find towerinfo node " << this_node_name << ", exiting." \
      << std::endl;
      exit(1);
    }


    for ( auto & p : towerjets )
    {
      inputs.push_back( p );
      inputs.back() -> set_id( inputs.size() - 1 ); // unique ids ensured
    }
  }

  // get jetcontainer to fill
  auto jetcon = findNode::getClass<JetContainer>( topNode, m_raw_seed_output );
  if ( !jetcon )
  {
    std::cerr << PHWHERE \
      << "can't find JetContainer " << m_raw_seed_output << \
      << " shoudl have been created by upstream module, exiting." \
    << std::endl;
    exit(1);
  }

  // reset jet container
  jetcon -> Reset();

  m_jetalgo -> cluster_and_fill( inputs, jetcon ); // fills the jet
  for ( auto & input : m_inputs )
  {
    jetcon -> insert_src( input -> get_src() );
  }
  
  if ( Verbosity() > 7 )
  {
    std::cout << PHWHERE << \
      "jets in container " << m_raw_seed_output \
    << std::endl;
    jetcon -> print_jets();
  }

  // clean up inputs
  for ( auto & input : inputs )
  {
    delete input;
  }
  inputs.clear();

  if ( Verbosity() > 0 )
  {
    std::cout << PHWHERE << \
      "recoSeeds: found " << jetcon -> size() \
      << " jets in container " << m_raw_seed_output \
    << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CopyAndSubtractJets::subJets( PHCompositeNode *topNode )
{
  if (Verbosity() > 0)
  {
    std::cout << "CopyAndSubtractJets::process_event: entering, with _use_flow_modulation = " << _use_flow_modulation << std::endl;
  }

  
  TowerInfoContainer *towerinfosEM3 = nullptr;
  TowerInfoContainer *towerinfosIH3 = nullptr;
  TowerInfoContainer *towerinfosOH3 = nullptr;
  EMTowerName = m_towerNodePrefix + "_CEMC_RETOWER";
  IHTowerName = m_towerNodePrefix + "_HCALIN";
  OHTowerName = m_towerNodePrefix + "_HCALOUT";
  towerinfosEM3 = findNode::getClass<TowerInfoContainer>(topNode, EMTowerName);
  towerinfosIH3 = findNode::getClass<TowerInfoContainer>(topNode, IHTowerName);
  towerinfosOH3 = findNode::getClass<TowerInfoContainer>(topNode, OHTowerName);
  if (!towerinfosEM3)
  {
    std::cout << "CopyAndSubtractJets::process_event: Cannot find node " << EMTowerName << std::endl;
    exit(1);
  }
  if (!towerinfosIH3)
  {
    std::cout << "CopyAndSubtractJets::process_event: Cannot find node " << IHTowerName << std::endl;
    exit(1);
  }
  if (!towerinfosOH3)
  {
    std::cout << "CopyAndSubtractJets::process_event: Cannot find node " << OHTowerName << std::endl;
    exit(1);
  }
  

  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  RawTowerGeomContainer *geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");

  // pull out jets and background
  JetContainer *unsub_jets;
  JetContainer *sub_jets;
  TowerBackground *background;
  if (m_use_towerinfo)
  {
    unsub_jets = findNode::getClass<JetContainer>(topNode, "AntiKt_TowerInfo_HIRecoSeedsRaw_r02");
    sub_jets = findNode::getClass<JetContainer>(topNode, "AntiKt_TowerInfo_HIRecoSeedsSub_r02");
    background = findNode::getClass<TowerBackground>(topNode, "TowerInfoBackground_Sub1");
  }

  std::vector<float> background_UE_0 = background->get_UE(0);
  std::vector<float> background_UE_1 = background->get_UE(1);
  std::vector<float> background_UE_2 = background->get_UE(2);

  float background_v2 = background->get_v2();
  float background_Psi2 = background->get_Psi2();

  if (Verbosity() > 0)
  {
    std::cout << "CopyAndSubtractJets::process_event: entering with # unsubtracted jets = " << unsub_jets->size() << std::endl;
    std::cout << "CopyAndSubtractJets::process_event: entering with # subtracted jets = " << sub_jets->size() << std::endl;
  }

  // iterate over old jets
  int ijet = 0;
  for (auto *this_jet : *unsub_jets)
  {
    float this_pt = this_jet->get_pt();
    float this_phi = this_jet->get_phi();
    float this_eta = this_jet->get_eta();

    float new_total_px = 0;
    float new_total_py = 0;
    float new_total_pz = 0;
    float new_total_e = 0;

    // if (this_jet->get_pt() < 5) continue;

    if (Verbosity() > 1 && this_jet->get_pt() > 5)
    {
      std::cout << "CopyAndSubtractJets::process_event: unsubtracted jet with pt / eta / phi = " << this_pt << " / " << this_eta << " / " << this_phi << std::endl;
    }

    for (const auto &comp : this_jet->get_comp_vec())
    {
      RawTower *tower = nullptr;
      RawTowerGeom *tower_geom = nullptr;
      TowerInfo *towerinfo = nullptr;

      double comp_e = 0;
      double comp_eta = 0;
      double comp_phi = 0;

      int comp_ieta = 0;

      double comp_background = 0;

      if (m_use_towerinfo)
      {
        if (comp.first == 5 || comp.first == 26)
        {
          towerinfo = towerinfosIH3->get_tower_at_channel(comp.second);
          unsigned int towerkey = towerinfosIH3->encode_key(comp.second);
          comp_ieta = towerinfosIH3->getTowerEtaBin(towerkey);
          int comp_iphi = towerinfosIH3->getTowerPhiBin(towerkey);
          const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, comp_ieta, comp_iphi);

          tower_geom = geomIH->get_tower_geometry(key);
          comp_background = background_UE_1.at(comp_ieta);
        }
        else if (comp.first == 7 || comp.first == 27)
        {
          towerinfo = towerinfosOH3->get_tower_at_channel(comp.second);
          unsigned int towerkey = towerinfosOH3->encode_key(comp.second);
          comp_ieta = towerinfosOH3->getTowerEtaBin(towerkey);
          int comp_iphi = towerinfosOH3->getTowerPhiBin(towerkey);
          const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, comp_ieta, comp_iphi);
          tower_geom = geomOH->get_tower_geometry(key);
          comp_background = background_UE_2.at(comp_ieta);
        }
        else if (comp.first == 13 || comp.first == 28)
        {
          towerinfo = towerinfosEM3->get_tower_at_channel(comp.second);
          unsigned int towerkey = towerinfosEM3->encode_key(comp.second);
          comp_ieta = towerinfosEM3->getTowerEtaBin(towerkey);
          int comp_iphi = towerinfosEM3->getTowerPhiBin(towerkey);
          const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, comp_ieta, comp_iphi);

          tower_geom = geomIH->get_tower_geometry(key);
          comp_background = background_UE_0.at(comp_ieta);
        }
        if (towerinfo)
        {
          comp_e = towerinfo->get_energy();
        }
      }
      else
      {
        if (comp.first == 5)
        {
          tower = towersIH3->getTower(comp.second);
          tower_geom = geomIH->get_tower_geometry(tower->get_key());

          comp_ieta = geomIH->get_etabin(tower_geom->get_eta());
          comp_background = background_UE_1.at(comp_ieta);
        }
        else if (comp.first == 7)
        {
          tower = towersOH3->getTower(comp.second);
          tower_geom = geomOH->get_tower_geometry(tower->get_key());

          comp_ieta = geomOH->get_etabin(tower_geom->get_eta());
          comp_background = background_UE_2.at(comp_ieta);
        }
        else if (comp.first == 13)
        {
          tower = towersEM3->getTower(comp.second);
          tower_geom = geomIH->get_tower_geometry(tower->get_key());

          comp_ieta = geomIH->get_etabin(tower_geom->get_eta());
          comp_background = background_UE_0.at(comp_ieta);
        }
        if (tower)
        {
          comp_e = tower->get_energy();
        }
      }

      if (tower_geom)
      {
        comp_eta = tower_geom->get_eta();
        comp_phi = tower_geom->get_phi();
      }

      if (Verbosity() > 4 && this_jet->get_pt() > 5)
      {
        std::cout << "CopyAndSubtractJets::process_event: --> constituent in layer " << comp.first << ", has unsub E = " << comp_e << ", is at ieta #" << comp_ieta << ", and has UE = " << comp_background << std::endl;
      }

      // flow modulate background if turned on
      if (_use_flow_modulation)
      {
        comp_background = comp_background * (1 + 2 * background_v2 * cos(2 * (comp_phi - background_Psi2)));
        if (Verbosity() > 4 && this_jet->get_pt() > 5)
        {
          std::cout << "CopyAndSubtractJets::process_event: --> --> flow mod, at phi = " << comp_phi << ", v2 and Psi2 are = " << background_v2 << " , " << background_Psi2 << ", UE after modulation = " << comp_background << std::endl;
        }
      }

      // update constituent energy based on the background
      double comp_sub_e = comp_e - comp_background;

      // now define new kinematics

      double comp_px = comp_sub_e / cosh(comp_eta) * cos(comp_phi);
      double comp_py = comp_sub_e / cosh(comp_eta) * sin(comp_phi);
      double comp_pz = comp_sub_e * tanh(comp_eta);

      new_total_px += comp_px;
      new_total_py += comp_py;
      new_total_pz += comp_pz;
      new_total_e += comp_sub_e;
    }

    auto *new_jet = sub_jets->add_jet();  // returns a new Jet_v2

    new_jet->set_px(new_total_px);
    new_jet->set_py(new_total_py);
    new_jet->set_pz(new_total_pz);
    new_jet->set_e(new_total_e);
    new_jet->set_id(ijet);

    if (Verbosity() > 1 && this_pt > 5)
    {
      std::cout << "CopyAndSubtractJets::process_event: old jet #" << ijet << ", old px / py / pz / e = " << this_jet->get_px() << " / " << this_jet->get_py() << " / " << this_jet->get_pz() << " / " << this_jet->get_e() << std::endl;
      std::cout << "CopyAndSubtractJets::process_event: new jet #" << ijet << ", new px / py / pz / e = " << new_jet->get_px() << " / " << new_jet->get_py() << " / " << new_jet->get_pz() << " / " << new_jet->get_e() << std::endl;
    }

    ijet++;
  }

  if (Verbosity() > 0)
  {
    std::cout << "CopyAndSubtractJets::process_event: exiting with # subtracted jets = " << sub_jets->size() << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloTowerSub::CreateNode( PHCompositeNode *topNode )
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // store the jet background stuff under a sub-node directory
  PHCompositeNode *bkgNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "JETBACKGROUND"));
  if (!bkgNode)
  {
    bkgNode = new PHCompositeNode("JETBACKGROUND");
    dstNode->addNode(bkgNode);
  }

  // create the TowerBackground node...
  TowerBackground *towerbackground = findNode::getClass<TowerBackground>(topNode, _backgroundName);
  if (!towerbackground)
  {
    towerbackground = new TowerBackgroundv1();
    PHIODataNode<PHObject> *bkgDataNode = new PHIODataNode<PHObject>(towerbackground, _backgroundName, "PHObject");
    bkgNode->addNode(bkgDataNode);
  }
  else
  {
    std::cout << PHWHERE << "::ERROR - " << _backgroundName << " pre-exists, but should not" << std::endl;
    exit(-1);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloTowerSub::FillNode( PHCompositeNode *topNode )
{
  TowerBackground *towerbackground = findNode::getClass<TowerBackground>(topNode, _backgroundName);
  if (!towerbackground)
  {
    std::cout << " ERROR -- can't find TowerBackground node after it should have been created" << std::endl;
    return;
  }

  towerbackground->set_UE(0, _UE[0]);
  towerbackground->set_UE(1, _UE[1]);
  towerbackground->set_UE(2, _UE[2]);

  towerbackground->set_v2(_v2);

  towerbackground->set_Psi2(_Psi2);

  towerbackground->set_nStripsUsedForFlow(_nStrips);

  towerbackground->set_nTowersUsedForBkg(_nTowers);

  towerbackground->set_flow_failure_flag(_is_flow_failure);

  return;
}






