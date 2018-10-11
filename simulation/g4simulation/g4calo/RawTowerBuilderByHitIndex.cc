#include "RawTowerBuilderByHitIndex.h"

#include <calobase/RawTowerv1.h>
#include <calobase/RawTowerContainer.h>

#include <calobase/RawTowerGeomv3.h>
#include <calobase/RawTowerGeomContainerv1.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <TRotation.h>
#include <TVector3.h>

#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <map>

using namespace std;

RawTowerBuilderByHitIndex::RawTowerBuilderByHitIndex(const std::string& name):
  SubsysReco(name),
  towers_(NULL),
  geoms_(NULL),
  detector_("NONE"),
  node_name_hits_("DEFAULT"),
  node_name_towers_("DEFAULT"),
  node_name_tower_geometries_("DEFAULT"),
  mapping_tower_file_("default.txt"),
  calo_id_( RawTowerDefs::NONE ),
  global_place_in_x_(0),
  global_place_in_y_(0),
  global_place_in_z_(0),
  rot_in_x_(0),
  rot_in_y_(0),
  rot_in_z_(0),
  emin_(1e-6),
  timer_( PHTimeServer::get()->insert_new(name) )
{}

int
RawTowerBuilderByHitIndex::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
      exit(1);
    }

  try
    {
      CreateNodes(topNode);
    }
  catch (std::exception& e)
    {
      std::cout << e.what() << std::endl;
      //exit(1);
    }

  try
    {
      ReadGeometryFromTable();
    }
  catch (std::exception& e)
    {
      std::cout << e.what() << std::endl;
      //exit(1);
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int
RawTowerBuilderByHitIndex::process_event(PHCompositeNode *topNode)
{
  // get hits
  node_name_hits_ = "G4HIT_" + detector_;
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, node_name_hits_.c_str());
  if (!g4hit)
    {
      cout << "Could not locate g4 hit node " << node_name_hits_ << endl;
      exit(1);
    }

  // loop over all hits in the event
  PHG4HitContainer::ConstIterator hiter;
  PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits();

  for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; hiter++)
    {
      PHG4Hit* g4hit_i =  hiter->second ;
      
      // Don't include hits with zero energy
      if(g4hit_i->get_edep()<=0) continue; 

      /* encode CaloTowerID from j, k index of tower / hit and calorimeter ID */
      RawTowerDefs::keytype calotowerid = RawTowerDefs::encode_towerid( calo_id_ ,
							       g4hit_i->get_index_j() ,
							       g4hit_i->get_index_k() );

      /* add the energy to the corresponding tower */
      RawTowerv1 *tower = dynamic_cast<RawTowerv1 *> (towers_->getTower( calotowerid ));
      if (! tower)
        {
          tower = new RawTowerv1( calotowerid );
	  tower->set_energy( 0 );
          towers_->AddTower( tower->get_id() , tower );
        }

      tower->add_ecell( (g4hit_i->get_index_j()<<16) + g4hit_i->get_index_k(), g4hit_i->get_light_yield());
      tower->set_energy( tower->get_energy() + g4hit_i->get_light_yield() );
      tower->add_eshower(g4hit_i->get_shower_id(), g4hit_i->get_edep());           
    }

  float towerE = 0.;

  if (Verbosity())
    {
      towerE = towers_->getTotalEdep();
    }

  towers_->compress(emin_);
  if (Verbosity())
    {
      cout << "Energy lost by dropping towers with less than " << emin_
	   << " energy, lost energy: "  << towerE - towers_->getTotalEdep() << endl;
      towers_->identify();
      RawTowerContainer::ConstRange begin_end = towers_->getTowers();
      RawTowerContainer::ConstIterator iter;
      for (iter =  begin_end.first; iter != begin_end.second; ++iter)
        {
          iter->second->identify();
        }
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int
RawTowerBuilderByHitIndex::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}


void
RawTowerBuilderByHitIndex::Detector( const std::string &d )
{
  detector_ = d;
  calo_id_ = RawTowerDefs::convert_name_to_caloid( detector_ );
}


void
RawTowerBuilderByHitIndex::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
    {
      std::cerr << PHWHERE << "Run Node missing, doing nothing." << std::endl;
      throw std::runtime_error("Failed to find Run node in RawTowerBuilderByHitIndex::CreateNodes");
    }

  PHCompositeNode *dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cerr << PHWHERE << "DST Node missing, doing nothing." << std::endl;
      throw std::runtime_error("Failed to find DST node in RawTowerBuilderByHitIndex::CreateNodes");
    }

  // Create the tower geometry node on the tree
  geoms_ = new RawTowerGeomContainerv1( RawTowerDefs::convert_name_to_caloid( detector_ ) );
  node_name_tower_geometries_ = "TOWERGEOM_" + detector_;

  PHIODataNode<PHObject> *geomNode = new PHIODataNode<PHObject>(geoms_, node_name_tower_geometries_.c_str(), "PHObject");
  runNode->addNode(geomNode);

  // Find detector node (or create new one if not found)
  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode*>(dstiter.findFirst(
      "PHCompositeNode", detector_));
  if (!DetNode)
    {
      DetNode = new PHCompositeNode(detector_);
      dstNode->addNode(DetNode);
    }

  // Create the tower nodes on the tree
  towers_ = new RawTowerContainer( RawTowerDefs::convert_name_to_caloid( detector_ ) );

  if ( sim_tower_node_prefix_.length() == 0 )
    {
      // no prefix, consistent with older convension
      node_name_towers_ = "TOWER_" + detector_;
    }
  else
    {
      node_name_towers_ = "TOWER_" + sim_tower_node_prefix_ + "_" + detector_;
    }

  PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(towers_, node_name_towers_.c_str(), "PHObject");
  DetNode->addNode(towerNode);

  return;
}


bool RawTowerBuilderByHitIndex::ReadGeometryFromTable() {

  /* Stream to read table from file */
  ifstream istream_mapping;

  /* Open the datafile, if it won't open return an error */
  if (!istream_mapping.is_open())
    {
      istream_mapping.open( mapping_tower_file_.c_str() );
      if(!istream_mapping)
	{
	  cerr << "CaloTowerGeomManager::ReadGeometryFromTable - ERROR Failed to open mapping file " << mapping_tower_file_ << endl;
	  exit(1);
	}
    }

  string line_mapping;

  while ( getline( istream_mapping, line_mapping ) )
    {
      /* Skip lines starting with / including a '#' */
      if ( line_mapping.find("#") != string::npos )
	{
	  if ( Verbosity() > 0 )
	    {
	      cout << "RawTowerBuilderByHitIndex: SKIPPING line in mapping file: " << line_mapping << endl;
	    }
	  continue;
	}

      istringstream iss(line_mapping);

      /* If line starts with keyword Tower, add to tower positions */
      if ( line_mapping.find("Tower ") != string::npos )
	{
	  unsigned idx_j, idx_k, idx_l;
	  double pos_x, pos_y, pos_z;
	  double size_x, size_y, size_z;
	  double rot_x, rot_y, rot_z;
	  double type; 
	  string dummys;

	  /* read string- break if error */
	  if ( !( iss >> dummys >> type >> idx_j >> idx_k >> idx_l >> pos_x >> pos_y >> pos_z >> size_x >> size_y >> size_z >> rot_x >> rot_y >> rot_z ) )
	    {
	      cerr << "ERROR in RawTowerBuilderByHitIndex: Failed to read line in mapping file " << mapping_tower_file_ << endl;
	      exit(1);
	    }

	  /* Construct unique Tower ID */
	  unsigned int temp_id = RawTowerDefs::encode_towerid( calo_id_ , idx_j , idx_k );

	  /* Create tower geometry object */
	  RawTowerGeom* temp_geo = new RawTowerGeomv3( temp_id );
	  temp_geo->set_center_x( pos_x );
	  temp_geo->set_center_y( pos_y );
	  temp_geo->set_center_z( pos_z );
	  temp_geo->set_size_x( size_x );
	  temp_geo->set_size_y( size_y );
	  temp_geo->set_size_z( size_z );
	  temp_geo->set_tower_type( (int) type );

	  /* Insert this tower into position map */
	  geoms_->add_tower_geometry( temp_geo );
	}

      /* If line does not start with keyword Tower, read as parameter */
      else
	{
	  /* If this line is not a comment and not a tower, save parameter as string / value. */
	  string parname;
	  double parval;

	  /* read string- break if error */
	  if ( !( iss >> parname >> parval ) )
	    {
	      cerr << "ERROR in RawTowerBuilderByHitIndex: Failed to read line in mapping file " << mapping_tower_file_ << endl;
	      exit(1);
	    }

	  map_global_parameter_.insert( make_pair( parname , parval ) );

	}
    }

  /* Update member variables for global parameters based on parsed parameter file */
  std::map<string,double>::iterator parit;

  parit = map_global_parameter_.find("Gx0");
  if (parit != map_global_parameter_.end())
    global_place_in_x_ = parit->second;

  parit = map_global_parameter_.find("Gy0");
  if (parit != map_global_parameter_.end())
    global_place_in_y_ = parit->second;

  parit = map_global_parameter_.find("Gz0");
  if (parit != map_global_parameter_.end())
    global_place_in_z_ = parit->second;

  parit = map_global_parameter_.find("Grot_x");
  if (parit != map_global_parameter_.end())
    rot_in_x_ = parit->second;

  parit = map_global_parameter_.find("Grot_y");
  if (parit != map_global_parameter_.end())
    rot_in_y_ = parit->second;

  parit = map_global_parameter_.find("Grot_z");
  if (parit != map_global_parameter_.end())
    rot_in_z_ = parit->second;


  /* Correct tower geometries for global calorimter translation / rotation 
   * after reading parameters from file */
  RawTowerGeomContainer::ConstRange all_towers = geoms_->get_tower_geometries();

  for ( RawTowerGeomContainer::ConstIterator it = all_towers.first;
	it != all_towers.second; ++it)
    {
      double x_temp = it->second->get_center_x();
      double y_temp = it->second->get_center_y();
      double z_temp = it->second->get_center_z();

      /* Rotation */
      TRotation rot;
      rot.RotateX( rot_in_x_ );
      rot.RotateY( rot_in_y_ );
      rot.RotateZ( rot_in_z_ );

      TVector3 v_temp_r( x_temp , y_temp , z_temp );
      v_temp_r.Transform( rot );

      /* Translation */
      double x_temp_rt = v_temp_r.X() + global_place_in_x_;
      double y_temp_rt = v_temp_r.Y() + global_place_in_y_;
      double z_temp_rt = v_temp_r.Z() + global_place_in_z_;

      /* Update tower geometry object */
      it->second->set_center_x( x_temp_rt );
      it->second->set_center_y( y_temp_rt );
      it->second->set_center_z( z_temp_rt );

      if ( Verbosity() > 2 )
	{
	  cout << "* Local tower x y z : " << x_temp << " " << y_temp << " " << z_temp << endl;
	  cout << "* Globl tower x y z : " << x_temp_rt << " " << y_temp_rt << " " << z_temp_rt << endl;
	}
    }

  return true;

}
