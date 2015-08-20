
#include "CaloEvaluator.h"

#include "CaloTruthEval.h"
#include "CaloRawTowerEval.h"
#include "CaloRawClusterEval.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/getClass.h>
#include <phool/PHCompositeNode.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4Particle.h>
#include <g4hough/SvtxVertexMap.h>
#include <g4cemc/RawTowerContainer.h>
#include <g4cemc/RawTower.h>
#include <g4cemc/RawClusterContainer.h>
#include <g4cemc/RawCluster.h>

#include <TNtuple.h>
#include <TFile.h>

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

CaloEvaluator::CaloEvaluator(const string &name, const string &caloname, const string &filename) 
  : SubSystemReco(name),
    _caloname(caloname),
    _filename(filename) {
  verbosity = 0;
}

int CaloEvaluator::Init(PHCompositeNode *topNode)
{
  _ievent = 0;

  _gshower_list.clear();
  _g4hit_gshower_map.clear();
  _g4hitid_gshower_map.clear();
  _tfile = new TFile(_filename.c_str(), "RECREATE");

  _ntp_event = new TNtuple("ntp_event","event-wise ntuple",
                           "event:vx:vy:vz:"
                           "gvx:gvy:gvz:ngshowers:"
                           "ng4hits:ntowers:nclusters");

  _ntp_gshower = new TNtuple("ntp_gshower","shower-wise ntuple",
  			     "event:particleID:flavor:"
			     "px:py:pz:e:vx:vy:vz:nhits:mrad:edep:eabs:embed");
  
  _ntp_tower = new TNtuple("ntp_tower","tower-wise ntuple",
			   "event:ieta:iphi:eta:phi:e:"
			   "gparticleID:gflavor:geta:gphi:ge:gpt:gmrad:gedep:gembed:"
			   "epurity");

  _ntp_cluster = new TNtuple("ntp_cluster","cluster-wise ntuple",
			     "event:clusterID:ntowers:eta:phi:e:"
			     "gparticleID:gflavor:geta:gphi:ge:gpt:gmrad:gedep:gembed:"
			     "epurity");

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloEvaluator::process_event(PHCompositeNode *topNode)
{
  // pull the g4 truth
  _truth_info_container = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if(!_truth_info_container) 
    {
      cout << PHWHERE << " WARNING: Can't find PHG4TruthInfoContainer." << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  // pull the g4hits
  string hitnodename = "G4HIT_" + detector;
  _g4hitList = findNode::getClass<PHG4HitContainer>(topNode,hitnodename.c_str());    
  if(!_g4hitList) 
    {
      cerr << PHWHERE << " ERROR: Can't find " << hitnodename << " node, aborting!" << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  // pull the g4hits in the absorber if present
  string hitabsorbernodename = "G4HIT_ABSORBER_" + detector;
  _g4hitAbsorberList = findNode::getClass<PHG4HitContainer>(topNode,hitabsorbernodename.c_str());    

  // pull the cells
  string cellnodename = "G4CELL_" + detector;
  _cyl_cel_container = findNode::getClass<PHG4CylinderCellContainer>(topNode,cellnodename.c_str());
  if(!_cyl_cel_container) 
    {
      cerr << PHWHERE << " ERROR: Can't find " << cellnodename << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  // pull the towers
  string towernodename = "TOWER_" + detector;
  _towerList = findNode::getClass<RawTowerContainer>(topNode,towernodename.c_str());
  if(!_towerList) 
    {
      cerr << PHWHERE << " ERROR: Can't find node " << towernodename << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  // pull the clusters
  string clusternodename = "CLUSTER_" + detector;
  _clusterList = findNode::getClass<RawClusterContainer>(topNode,clusternodename.c_str());
  if(!_clusterList) 
    {
      cerr << PHWHERE << " ERROR: Can't find node " << clusternodename << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
 
  // pull the tracks
  _trackList = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
  if(!_trackList) 
    {
      static bool firsttrackwarn = true;
      if (firsttrackwarn) {
	cout << PHWHERE << " WARNING: Can't find SvtxTrackMap. The reconstructed tracks are inaccessible." << endl;
	firsttrackwarn = false;
      }
    }
  
  _vertexList = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMap");
  if (!_vertexList) 
    {
      static bool firstvertexwarn = true;
      if (firstvertexwarn) {
	cout << PHWHERE << " WARNING: Can't find SvtxVertexMap. The reconstructed vertex is inaccessible." << endl;
	firstvertexwarn = false;
      }
    }
  
  string towergeomnode = "TOWERGEOM_" + detector;
  towergeom = findNode::getClass<RawTowerGeom>(topNode,towergeomnode.c_str());
  if (!towergeom)
    {
      cout << "cannot find node " << towergeomnode << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  //-----------------------------------
  // print what is coming into the code
  //-----------------------------------
  
  printInputInfo();
  
  if(_trace_truth)
    {
      //-------------------------------
      // fill the Gshower storage vector
      //-------------------------------
  
      fillGshowerObjects();

      //--------------------------------
      // fill the truth association maps
      //--------------------------------

      fillTowerToGshowerMap();
      fillClusterToGshowerMap();
    }

  //---------------------------
  // fill the Evaluator NTuples
  //---------------------------

  fillOutputNtuples();

  //--------------------------------------------------
  // Print out the ancestry information for this event
  //--------------------------------------------------
  
  printOutputInfo();
  
  _ievent++;

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloEvaluator::End(PHCompositeNode *topNode)
{
  _tfile->cd();

  _ntp_event->Write();
  _ntp_gshower->Write();
  _ntp_tower->Write();
  _ntp_cluster->Write();
  
  _tfile->Close();

  delete _tfile;

  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloEvaluator::fillGshowerObjects()
{
  if(verbosity > 1) cout << "CaloEvaluator::fillGshowerObjects() entered" << endl;

  // fill a vector of truth particles that know which g4hits they made

  // reset the vector
  _gshower_list.clear();
  _g4hit_gshower_map.clear();
  _g4hitid_gshower_map.clear();

  PHG4HitContainer::ConstIterator hiter;
  PHG4HitContainer::ConstRange hit_begin_end = _g4hitList->getHits();
  for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; hiter++)
    {
      PHG4Hit *g4hit = hiter->second;

      int particle_id = g4hit->get_trkid();

      PHG4Particle *particle = _truth_info_container->GetHit( particle_id );
      if(!particle)
	{
	  cout << "CaloEvaluator::fillGshowerObjects() - ERROR, Corrupt truth information, skip hit!" << endl;
	  cout << "   particle_id = " << particle_id << " particle pointer = " << particle << endl;
	  continue;
	}

      // do we already have this gshower already?
      bool isNewGshower = true;
      
      // loop to see if we have this primary gshower already
      for(unsigned int igshower = 0; igshower < _gshower_list.size(); igshower++)
	{
	  // if we do, add this ghit to the list
	  if (particle->get_primary_id() == _gshower_list[igshower].get_particle_id())
	    {
	      isNewGshower = false;
	      
	      CalGshower *gshower = &_gshower_list[igshower];
		      
	      gshower->add_g4hit(g4hit);
	    }
	}

      if(isNewGshower)
	{
	  // get the primary truth particle
	  if(particle->get_primary_id() != (int)(0xFFFFFFFF))
	    {
	      particle = _truth_info_container->GetHit( particle->get_primary_id() );
	    }

	  particle_id = particle->get_track_id();
	  
	  // get the vertex
	  PHG4VtxPoint* vertex = _truth_info_container->GetVtx( particle->get_vtx_id() );

	  float vx = NAN;
	  float vy = NAN;
	  float vz = NAN;

	  if(vertex)
	    {
	      vx = vertex->get_x();
	      vy = vertex->get_y();
	      vz = vertex->get_z();
	    }	  

	  CalGshower gshower;
	  
	  gshower.set_particle_id(particle_id);
	  gshower.set_flavor(particle->get_pid());
	  
	  gshower.set_px(particle->get_px());
	  gshower.set_py(particle->get_py());
	  gshower.set_pz(particle->get_pz());
		  
	  gshower.set_e(particle->get_e());
	  
	  gshower.set_vx(vx);
	  gshower.set_vy(vy);
	  gshower.set_vz(vz);

	  gshower.set_embed(_truth_info_container->isEmbeded(particle_id));	  
	  
	  gshower.add_g4hit(g4hit);

	  if (_g4hitAbsorberList) {
	    float eabs = 0.0;

	    PHG4HitContainer::ConstIterator hiter;
	    PHG4HitContainer::ConstRange hit_begin_end = _g4hitAbsorberList->getHits();
	    for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; hiter++)
	      {
		PHG4Hit *g4hit = hiter->second;
		int particle_id = g4hit->get_trkid();

		PHG4Particle *particle = _truth_info_container->GetHit( particle_id );
		if (!particle)
		  {
		    cout << "CaloEvaluator::fillGshowerObjects() - ERROR, Corrupt truth information, skip hit!" << endl;
		    cout << "   particle_id = " << particle_id << " particle pointer = " << particle << endl;
		    continue;
		  }

		if (particle->get_primary_id() != particle_id) continue;
		
		eabs += g4hit->get_edep();
	      }
	    
	    gshower.set_eabs(eabs);
	  }
	  
	  _gshower_list.push_back(gshower);
	}//isnewtrack
      // } // if truth particle iterator
    } // loop over all gh4its
      
  // fill the moliere radius (no account for magnetic bend)
  for(unsigned int igshower = 0; igshower < _gshower_list.size(); igshower++)
    {
      CalGshower *gshower = &_gshower_list[igshower];

      std::multimap<float,float> radii_energy_mmap;
      for(unsigned int ig4hit = 0; ig4hit < gshower->get_ng4hits(); ig4hit++)
	{
	  PHG4Hit *g4hit = gshower->get_g4hit(ig4hit);
	      
	  // momentum vector
	  float p_x = gshower->get_px();
	  float p_y = gshower->get_py();
	  float p_z = gshower->get_pz();
	  float p   = sqrt(pow(p_x,2)+pow(p_y,2)+pow(p_z,2));

	  // relative position vector (vertex-to-ghit)
	  float d_x = gshower->get_vx() - g4hit->get_avg_x();
	  float d_y = gshower->get_vy() - g4hit->get_avg_y();
	  float d_z = gshower->get_vz() - g4hit->get_avg_z();
	  float d   = sqrt(pow(d_x,2)+pow(d_y,2)+pow(d_z,2));

	  // angle between them
	  float phi = acos( (p_x*d_x+p_y*d_y+p_z*d_z)/p/d );

	  // distance between them at ghit
	  float r = d*sin(phi); 
	  float edep = g4hit->get_edep();

	  radii_energy_mmap.insert(make_pair(r,edep));
	}

      float sum_e = 0.0;
      float frac_e = 0.0;

      float r_in = 0.0;
      float r_out = 0.0;

      for(std::multimap<float,float>::iterator iter = radii_energy_mmap.begin();
	  iter != radii_energy_mmap.end();
	  iter++)
	{
	  r_out = iter->first;
	  sum_e = sum_e + iter->second;
	  frac_e = sum_e / gshower->get_edep();

	  if(frac_e > 0.90) break;

	  r_in = r_out;
	}
	      
      gshower->set_moliere_radius(0.5*(r_in+r_out));
    }

  // now loop over all gshowers and fill a reverse lookup map from the
  // g4hit back to the gshower
  for(unsigned int igshower = 0; igshower < _gshower_list.size(); igshower++)
    {
      CalGshower *gshower = &_gshower_list[igshower];
	  
      for(unsigned int ig4hit = 0; ig4hit < gshower->get_ng4hits(); ig4hit++)
	{
	  PHG4Hit *g4hit = gshower->get_g4hit(ig4hit);
	      
	  // fill map
	  _g4hit_gshower_map.insert(make_pair(g4hit,gshower));
	  _g4hitid_gshower_map.insert(make_pair(g4hit->get_hit_id(),gshower));
	}
    }
}

void CaloEvaluator::printInputInfo()
{
  if(verbosity > 1) cout << "CaloEvaluator::printInputInfo() entered" << endl;

  // print out the truth container

  // if(verbosity > 0)
  //   {
  //     cout << "PHG4TruthInfoContainer contents: " << endl; 

  //     PHG4TruthInfoContainer::Range truthrange = _truth_info_container->GetHitRange();
  //     for(PHG4TruthInfoContainer::Iterator truthiter = truthrange.first;
  // 	  truthiter != truthrange.second;
  // 	  truthiter++)
  // 	{
  // 	  PHG4Particle *particle = truthiter->second;

  // 	  cout << truthiter->first << " => pid: " << particle->get_pid() << " pt: " << sqrt(pow(particle->get_px(),2)+pow(particle->get_py(),2)) << endl;
  // 	}
  //   }

  return;
}

void CaloEvaluator::printOutputInfo()
{
  if(verbosity > 1) cout << "CaloEvaluator::printOutputInfo() entered" << endl;

  //==========================================
  // print out some useful stuff for debugging
  //==========================================

  if(verbosity)
    {
      // event information
      cout << endl;
      cout << PHWHERE << "   NEW OUTPUT FOR EVENT " << _ievent << endl;
      cout << endl;

      PHG4VtxPoint *gvertex = _truth_info_container->GetPrimaryVtx( _truth_info_container->GetPrimaryVertexIndex() );
      float gvx = gvertex->get_x();
      float gvy = gvertex->get_y();
      float gvz = gvertex->get_z();

      float vx = NAN;
      float vy = NAN;
      float vz = NAN;
      if (_vertexList) {
	if (!_vertexList->empty()) {
	  SvtxVertex* vertex = &(_vertexList->begin()->second);
	
	  vx = vertex->get_x();
	  vy = vertex->get_y();
	  vz = vertex->get_z();
	}
      }

      cout << "vtrue = (" << gvx << "," << gvy << "," << gvz << ") => vreco = (" << vx << "," << vy << "," << vz << ")" << endl;
    
      float ngshowers = _gshower_list.size();
      float ng4hits  = _g4hitList->size();
      float ntowers  = _towerList->size();
      float nclusters = _clusterList->size();

      cout << "nGshowers = " << ngshowers << endl;
      cout << " => nGhits = " << ng4hits << endl;
      cout << " => nTowers = " << ntowers << endl;
      cout << " => nClusters = " << nclusters << endl;

      if(verbosity > 1)
	{
	  for(unsigned igshower = 0; igshower < _gshower_list.size(); igshower++)
	    {
	      CalGshower *gshower = &_gshower_list[igshower];
	      
	      // track-wise information
	      cout << endl;
      
	      cout << "===CalGshower===================================================" << endl;
	      cout << " CalGshower id = " << gshower->get_particle_id() << endl;
	      cout << " flavor = " << gshower->get_flavor() << endl;
	      cout << " ptrue = (";
	      cout.width(5); cout << gshower->get_px();
	      cout << ",";
	      cout.width(5); cout << gshower->get_py();
	      cout << ",";
	      cout.width(5); cout << gshower->get_pz();
	      cout << ")" << endl;
	      cout << " vtrue = (";
	      cout.width(5); cout << gshower->get_vx();
	      cout << ",";
	      cout.width(5); cout << gshower->get_vy();
	      cout << ",";
	      cout.width(5); cout << gshower->get_vz();
	      cout << ")" << endl;
	      cout << " ---Associated-PHG4Hits-----------------------------------------" << endl;
	      
	      for(unsigned int ig4hit = 0; ig4hit < gshower->get_ng4hits(); ig4hit++)
		{
		  if((ig4hit > 5)&&(ig4hit < gshower->get_ng4hits() - 5)) continue;

		  PHG4Hit *g4hit = gshower->get_g4hit(ig4hit);
		  
		  float x = 0.5*(g4hit->get_x(1)+g4hit->get_x(0));
		  float y = 0.5*(g4hit->get_y(1)+g4hit->get_y(0));
		  float z = 0.5*(g4hit->get_z(1)+g4hit->get_z(0));
		  
		  cout << " #" << ig4hit << " xtrue = (";
		  cout.width(5); cout << x;
		  cout << ",";
		  cout.width(5); cout << y;
		  cout << ",";
		  cout.width(5); cout << z;
		  cout << ")";
		  cout << " e = " << g4hit->get_edep();
		  
		  /*
		    typedef multimap<PHG4Hit*,SvtxCluster*>::iterator mapiter2;
		    typedef pair<mapiter2,mapiter2> maprange2;
		    maprange2 therange2 = _g4hit_cluster_mmap.equal_range( g4hit );
		    for(mapiter2 theiter2=therange2.first; theiter2!=therange2.second; theiter2++) 
		    {
		    SvtxCluster *cluster = theiter2->second;
		    
		    float x = cluster->getHitPosition(0);
		    float y = cluster->getHitPosition(1);
		    float z = cluster->getHitPosition(2);
	    
		    cout << " => #" << cluster->getClusterID() << " xreco = (";
		    cout.width(5); cout << x;
		    cout << ",";
		    cout.width(5); cout << y;
		    cout << ",";
		    cout.width(5); cout << z;
		    cout << ")";
		    }
		  */

		  cout << endl;
		}
	    }      
	}
    }

  return;
}

void CaloEvaluator::fillTowerToGshowerMap()
{
  if(verbosity > 1) cout << "CaloEvaluator::fillTowerToGshowerMap() entered" << endl;

  // clear maps
  _tower_gshower_map.clear();
  _tower_ghit_mmap.clear();
  _tower_epurity_map.clear();

  // energy contribution maps
  std::map <CalGshower*, float> tower_energy_contributors;

  // for every tower
  RawTowerContainer::ConstRange begin_end = _towerList->getTowers();
  RawTowerContainer::ConstIterator rtiter;
  for (rtiter = begin_end.first; rtiter !=  begin_end.second; ++rtiter)
    {
      RawTower *tower = rtiter->second;

	  tower_energy_contributors.clear();

	  // loop over all cells within those towers
	  std::pair< std::map<unsigned int,float>::const_iterator, std::map<unsigned int,float>::const_iterator > cell_range = tower->get_g4cells();
	  for(std::map<unsigned int, float>::const_iterator cell_iter = cell_range.first; cell_iter != cell_range.second; cell_iter++)
	    {
	      unsigned int cell_id = cell_iter->first;

	      PHG4CylinderCell *cell = _cyl_cel_container->findCylinderCell( cell_id );
	      if(!cell)
		{
		  cout << PHWHERE << "Error, Invalid cell id, cell " << cell_id << "is not in cell container" << endl;
		  continue;
		}

	      // loop over all ghits within those cells
	      PHG4CylinderCell::EdepConstRange g4hits = cell->get_g4hits();
	      PHG4CylinderCell::EdepConstIterator g4iter = g4hits.first;
	      for (g4iter = g4hits.first; g4iter != g4hits.second; g4iter++)
		{
		  // ask each ghit how much energy it deposited and where that energy came from
		  PHG4HitDefs::keytype g4hit_id = g4iter->first;
		  float g4hit_energy = g4iter->second;
		  CalGshower *gshower = _g4hitid_gshower_map[g4hit_id];
		  PHG4Hit *ghit = _g4hitList->findHit( g4hit_id );
		  if(!ghit)
		    {
		      cout << PHWHERE << "Error, Invalid ghit id, ghit " << g4hit_id << "is not in ghit container" << endl;
		      continue;
		    }

		  // fill quick loop up map
		  _tower_ghit_mmap.insert( make_pair(tower,ghit) );

		  // see if this gshower already exists in the tower map
		  std::map<CalGshower*,float>::iterator entry_iter = tower_energy_contributors.find( gshower );
		  if(entry_iter != tower_energy_contributors.end())
		    {
		      // if so, add the energy to the prexisting gshower entry
		      entry_iter->second = entry_iter->second + g4hit_energy;
		    }
		  else
		    {
		      // if not, create a new entry with this energy
		      tower_energy_contributors.insert( make_pair(gshower,g4hit_energy) );
		    }
		  
		} // g4hit-loop
	      //	    } // cell-search loop
	} // cellid-loop
	  
	  float max_energy = -9999.0;
	  CalGshower *tower_truth = NULL;
	  
	  for(std::map<CalGshower*,float>::iterator tower_iter = tower_energy_contributors.begin();
	      tower_iter != tower_energy_contributors.end();
	      tower_iter++)
	    {
	      // if this truth gives a higher energy
	      if(max_energy < tower_iter->second)
		{
		  tower_truth = tower_iter->first;
		  max_energy = tower_iter->second;
		}
	    }
	  
	  _tower_gshower_map.insert( make_pair( tower, tower_truth ) );	  
      
	  float epurity = max_energy/tower->get_energy();
	  _tower_epurity_map.insert( make_pair( tower, epurity ) );
    }
  return;
}

void CaloEvaluator::fillClusterToGshowerMap()
{
  if(verbosity > 1) cout << "CaloEvaluator::fillClusterToGshowerMap() entered" << endl;

  // clear maps
  _cluster_gshower_map.clear();
  _cluster_epurity_map.clear();
  _gshower_cluster_map.clear();

  // energy contribution maps
  std::map <CalGshower*, float> cluster_energy_contributors;

  // loop over all clusters
  for(unsigned int icluster = 0; icluster < _clusterList->size(); icluster++) 
    {
      cluster_energy_contributors.clear();

      RawCluster *cluster = _clusterList->getCluster(icluster);

      // loop over all towers with the cluster
      for(unsigned int itower = 0; itower < cluster->getNTowers(); itower++)
	{
	  std::pair<int,int> tower_pos = cluster->getTowerBin(itower);

	  RawTower *tower = _towerList->getTower(tower_pos.first,tower_pos.second);

	  // loop over all the ghits within those towers
	  std::pair< std::multimap<RawTower*,PHG4Hit*>::iterator, std::multimap<RawTower*,PHG4Hit*>::iterator > range = _tower_ghit_mmap.equal_range(tower);
	  for(std::multimap<RawTower*,PHG4Hit*>::iterator ghit_iter = range.first;
	      ghit_iter != range.second;
	      ghit_iter++)
	    {
	      PHG4Hit *ghit = ghit_iter->second;

	      PHG4HitDefs::keytype g4hit_id = ghit->get_hit_id();
	      float g4hit_energy = ghit->get_edep();

	      CalGshower *gshower = _g4hitid_gshower_map[g4hit_id];

	      // see if this gshower already exists in the tower map
	      std::map<CalGshower*,float>::iterator entry2_iter = cluster_energy_contributors.find( gshower );
	      if(entry2_iter != cluster_energy_contributors.end())
		{
		  // if so, add the energy to the prexisting gshower entry
		  entry2_iter->second = entry2_iter->second + g4hit_energy;
		}
	      else
		{
		  // if not, create a new entry with this energy
		  cluster_energy_contributors.insert( make_pair(gshower,g4hit_energy) );
		}
	      
	    } // g4hit-loop
	} // tower-loop

      float max_energy = -9999.0;
      CalGshower *cluster_truth = NULL;

      for(std::map<CalGshower*,float>::iterator cluster_iter = cluster_energy_contributors.begin();
	  cluster_iter != cluster_energy_contributors.end();
	  cluster_iter++)
	{
	  // if this truth gives a higher energy
	  if(max_energy < cluster_iter->second)
	    {
	      cluster_truth = cluster_iter->first;
	      max_energy = cluster_iter->second;
	    }
	}

      _cluster_gshower_map.insert( make_pair( cluster, cluster_truth ) );

      float epurity = max_energy/cluster->get_energy();
      _cluster_epurity_map.insert( make_pair( cluster, epurity ));
    } // cluster-loop


  // loop over all gshower objects
  // for each gshower object determine the largest energy cluster


  return;
}

void CaloEvaluator::fillOutputNtuples()
{
  if(verbosity > 1) cout << "CaloEvaluator::fillOutputNtuples() entered" << endl;

  //----------------------
  // fill the Event NTuple
  //----------------------

  PHG4VtxPoint *gvertex = _truth_info_container->GetPrimaryVtx( _truth_info_container->GetPrimaryVertexIndex() );
  float gvx = gvertex->get_x();
  float gvy = gvertex->get_y();
  float gvz = gvertex->get_z();

  float vx = NAN;
  float vy = NAN;
  float vz = NAN;
  if (_vertexList) {
    if (!_vertexList->empty()) {
      SvtxVertex* vertex = &(_vertexList->begin()->second);
      
      vx = vertex->get_x();
      vy = vertex->get_y();
      vz = vertex->get_z();
    }
  }

  float ngshowers = _gshower_list.size();
  float ng4hits  = _g4hitList->size();
  float ntowers  = _towerList->size();
  float nclusters = _clusterList->size();

  float event_data[12] = {_ievent,
			  vx,
			  vy,
			  vz,
			  gvx,
			  gvy,
			  gvz,
			  ngshowers,
			  ng4hits,
			  ntowers,
			  nclusters};

  _ntp_event->Fill(event_data);
  
  //------------------------
  // fill the Gshower NTuple
  //------------------------
  
  if(verbosity > 1) cout << "CaloEvaluator::filling gshower ntuple..." << endl;

  for(unsigned int ishower = 0; ishower < _gshower_list.size(); ishower++)
    {
      float particleID = _gshower_list[ishower].get_particle_id();
      float flavor     = _gshower_list[ishower].get_flavor();

      float px         = _gshower_list[ishower].get_px();
      float py         = _gshower_list[ishower].get_py();
      float pz         = _gshower_list[ishower].get_pz();
      float e          = _gshower_list[ishower].get_e();

      float vx         = _gshower_list[ishower].get_vx();
      float vy         = _gshower_list[ishower].get_vy();
      float vz         = _gshower_list[ishower].get_vz();

      float nhits      = (float)_gshower_list[ishower].get_ng4hits(); 

      float mrad       = _gshower_list[ishower].get_moliere_radius();

      float edep       = _gshower_list[ishower].get_edep();
      float eabs       = _gshower_list[ishower].get_eabs();

      float embed       = _gshower_list[ishower].get_embed();

      float shower_data[15] = {_ievent,
			       particleID,
			       flavor,
			       px,
			       py,
			       pz,
			       e,
			       vx,
			       vy,
			       vz,
			       nhits,
			       mrad,
			       edep,
			       eabs,
			       embed
      };

      _ntp_gshower->Fill(shower_data);
    }


  //----------------------
  // fill the Tower NTuple
  //----------------------
  
  if(verbosity > 1) cout << "CaloEvaluator::filling tower ntuple..." << endl;

  // for every tower
  RawTowerContainer::ConstRange begin_end = _towerList->getTowers();
  RawTowerContainer::ConstIterator rtiter;
  for (rtiter = begin_end.first; rtiter !=  begin_end.second; ++rtiter)
    {
	  RawTower *tower = rtiter->second;
	  if(!tower) continue;

	  CalGshower *gshower = _tower_gshower_map[ tower ];
	  
	  float ieta    = tower->get_bineta();
	  float iphi    = tower->get_binphi();
	  float eta     = towergeom->get_etacenter(tower->get_bineta());
	  float phi     = towergeom->get_phicenter(tower->get_binphi());
	  float e       = tower->get_energy();

	  float gparticleID = -9999.0;
	  float gflavor     = -9999.0;
	  float gphi        = -9999.0;
	  float ge          = -9999.0;
	  float gpt         = -9999.0;
	  float geta        = -9999.0;
	  float gmrad       = -9999.0;
	  float gedep       = -9999.0;
	  float gembed      = -9999.0;
	  float epurity     = -9999.0;

	  if(gshower)
	    {
	      gparticleID = gshower->get_particle_id();
	      gflavor     = gshower->get_flavor();
	      gphi        = atan2(gshower->get_py(),gshower->get_px());
	      ge          = gshower->get_e();
	      gpt         = sqrt(pow(gshower->get_px(),2)+pow(gshower->get_py(),2));
	      geta        = asinh(gshower->get_pz()/gpt);
	      gmrad       = gshower->get_moliere_radius();
	      gedep       = gshower->get_edep();
	      gembed      = gshower->get_embed();

	      epurity     = _tower_epurity_map[ tower ];
	    }

	  float tower_data[16] = {_ievent,
				  ieta,
				  iphi,
				  eta,
				  phi,
				  e,
				  gparticleID,
				  gflavor,
				  geta,
				  gphi,
				  ge,
				  gpt,
				  gmrad,
				  gedep,
				  gembed,
				  epurity};

	  _ntp_tower->Fill(tower_data);
	
    }

  //------------------------
  // fill the Cluster NTuple
  //------------------------

  if(verbosity > 1) cout << "CaloEvaluator::filling gcluster ntuple..." << endl;
  
  // for every cluster
  for(unsigned int icluster = 0; icluster < _clusterList->size(); icluster++) 
    {
      RawCluster *cluster = _clusterList->getCluster(icluster);

      CalGshower *gshower = _cluster_gshower_map[ cluster ];

      float clusterID = icluster;
      float ntowers   = cluster->getNTowers();
      float eta       = cluster->get_eta();
      float phi       = cluster->get_phi();
      float e         = cluster->get_energy();

      float gparticleID = -9999.0;
      float gflavor     = -9999.0;
      float gphi        = -9999.0;
      float ge          = -9999.0;
      float gpt         = -9999.0;
      float geta        = -9999.0;
      float gmrad       = -9999.0;
      float gedep       = -9999.0;
      float gembed      = -9999.0;
      
      float epurity     = -9999.0;

      if(gshower)
	{
	  gparticleID = gshower->get_particle_id();
	  gflavor     = gshower->get_flavor();
	  gphi        = atan2(gshower->get_py(),gshower->get_px());
	  ge          = gshower->get_e();
	  gpt         = sqrt(pow(gshower->get_px(),2)+pow(gshower->get_py(),2));
	  geta        = asinh(gshower->get_pz()/gpt);
	  gmrad       = gshower->get_moliere_radius();
	  gedep       = gshower->get_edep();
	  gembed      = gshower->get_embed();

	  epurity     = _cluster_epurity_map[ cluster ];
	}

      float cluster_data[16] = {_ievent,
				clusterID,
				ntowers,
				eta,
				phi,
				e,
				gparticleID,
				gflavor,
				geta,
				gphi,
				ge,
				gpt,
				gmrad,
				gedep,
				gembed,
				epurity
      };

      _ntp_cluster->Fill(cluster_data);
    }

  return;
}


double
CalGshower::get_edep() const
{
    double energy = 0.0;
    for(unsigned int ihit = 0; ihit < get_ng4hits(); ihit++)
      {
	if (_g4hits[ihit]->get_edep() > 0)
	  {
	    energy += _g4hits[ihit]->get_edep();
	  }
      }
    return energy;
  }
