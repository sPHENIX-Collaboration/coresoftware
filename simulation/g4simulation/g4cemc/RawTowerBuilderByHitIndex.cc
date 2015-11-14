#include "RawTowerBuilderByHitIndex.h"

#include "RawTowerv1.h"
#include "RawTowerContainer.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <map>

using namespace std;

RawTowerBuilderByHitIndex::RawTowerBuilderByHitIndex(const std::string& name):
  SubsysReco(name),
  towers_(NULL),
  detector_("NONE"),
  calo_id_( RawTowerDefs::NONE ),
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

      /* encode CaloTowerID from j, k index of tower / hit and calorimeter ID */
      unsigned int calotowerid = RawTowerDefs::encode_towerid( calo_id_ ,
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
      tower->add_ecell(g4hit_i->get_trkid(), g4hit_i->get_edep());
      tower->set_energy( tower->get_energy() + g4hit_i->get_edep() );
    }

  float towerE = 0.;

  if (verbosity)
    {
      towerE = towers_->getTotalEdep();
    }

  towers_->compress(emin_);
  if (verbosity)
    {
      cout << "Energy lost by dropping towers with less than "
           << emin_ << " energy, lost energy: "  << towerE - towers_->getTotalEdep() << endl;
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

  // Create the tower nodes on the tree
  towers_ = new RawTowerContainer();
  node_name_towers_ = "TOWER_" + detector_;

  PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(towers_, node_name_towers_.c_str(), "PHObject");
  dstNode->addNode(towerNode);

  return;
}
