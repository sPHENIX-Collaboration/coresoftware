#include "RawTowerBuilderCone.h"
#include "RawTowerContainer.h"
#include "RawTowerGeomv1.h"
#include "RawTowerv2.h"
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4detectors/PHG4CylinderCell.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/getClass.h>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <map>

using namespace std;

RawTowerBuilderCone::RawTowerBuilderCone(const std::string& name):
  SubsysReco(name),
  _towers(NULL),
  rawtowergeom(NULL),
  detector("NONE"),
  emin(1e-6),
  _nlayers(-1),
  _nphibins(36),
  _netabins(30),
  _zmin(350.),
  _zmax(450.),
  _etamin(1.15),
  _etamax(5.),
  _phimin(0),
  _phimax(2*M_PI),
  _thetamin(0),
  _thetamax(0),
  _thetastep(NAN),
  _phistep(NAN),
  _timer( PHTimeServer::get()->insert_new(name) )
{}

int
RawTowerBuilderCone::InitRun(PHCompositeNode *topNode)
{

  _phistep = (_phimax - _phimin ) / static_cast<double>( _nphibins );

  if ( _etamin > 0 )
  {
  _thetamin = eta2theta( _etamax );
  _thetamax = eta2theta( _etamin );
  }
  else
  {
  _thetamin = eta2theta( _etamin );
  _thetamax = eta2theta( _etamax );  
  }
  _thetastep = (_thetamax - _thetamin ) / static_cast<double>( _netabins );

  if (verbosity)
    {
      std::cout << "Eta:   " << _netabins << " from " << _etamin << " to " << _etamax << std::endl;
      std::cout << "Theta: " << _netabins << " from " << _thetamin << " to " << _thetamax << " , width " << _thetastep << std::endl;
      std::cout << "Phi:   " << _nphibins << " from " << _phimin << " to " << _phimax << " , width " << _phistep << std::endl;
    }

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
RawTowerBuilderCone::process_event(PHCompositeNode *topNode)
{
  if (verbosity)
    {
      std::cout << PHWHERE << "Process event entered" << std::endl;
    }

  // get hits
  hitnodename = "G4HIT_" + detector;
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  if (!g4hit)
    {
      cout << "Could not locate g4 hit node " << hitnodename << endl;
      exit(1);
    }

  // loop over all hits in the event
  PHG4HitContainer::ConstIterator hiter;
  PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits();

  for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; hiter++)
    {
      pair<double, double> etaphi;
      double phibin;
      double etabin;

      etaphi = get_etaphi( ( hiter->second->get_x(0) + hiter->second->get_x(1) ) / 2.0,
                           ( hiter->second->get_y(0) + hiter->second->get_y(1) ) / 2.0,
                           ( hiter->second->get_z(0) + hiter->second->get_z(1) ) / 2.0 );

      //cout << "Hit coordinates theta - phi: " << eta2theta( etaphi.first ) << " - " << etaphi.second << endl;

      etabin = get_thetabin( eta2theta( etaphi.first ) );
      phibin = get_phibin( etaphi.second );

      // check bin range
      if (etabin < 0 || etabin >= _netabins)
        {
          continue;
        }
      if (phibin < 0 || phibin >= _nphibins)
        {
          continue;
        }

      // add the energy to the corresponding tower
      RawTowerv2 *tower = dynamic_cast<RawTowerv2 *> (_towers->getTower( etabin, phibin ));
      if (! tower)
        {
          tower = new RawTowerv2( etabin, phibin );
	  tower->set_edep( 0.0 );
          tower->set_thetaMin(_thetamin + etabin * _thetastep);
          tower->set_thetaSize(_thetastep);
          tower->set_phiMin(_phimin + phibin * _phistep);
          tower->set_phiSize(_phistep);
          tower->set_zMin(_zmin);
          tower->set_zSize(_zmax - _zmin);

          _towers->AddTower( etabin, phibin, tower);
        }
      int cellid = 0;
      tower->add_ecell(cellid, hiter->second->get_edep());
      tower->set_edep( tower->get_edep() + hiter->second->get_edep() );

      //cout << "Tower coordinates theta - phi: (" << tower->get_thetaMin() << "-" << tower->get_thetaMin()+tower->get_thetaSize()
      // << ") - (" << tower->get_phiMin() << "-" <<  tower->get_phiMin()+tower->get_phiSize()<< ")" << endl;

      //cout << "Tower eta bin = " << etabin << " , phi bin = " << phibin << ": add energy of " << hiter->second->get_edep() << endl;
    }

  float towerE = 0.;

  if (verbosity)
    {
      towerE = _towers->getTotalEdep();
    }

  _towers->compress(emin);
  if (verbosity)
    {
      cout << "Energy lost by dropping towers with less than "
           << emin << " energy, lost energy: "  << towerE - _towers->getTotalEdep() << endl;
      _towers->identify();
      RawTowerContainer::ConstRange begin_end = _towers->getTowers();
      RawTowerContainer::ConstIterator iter;
      for (iter =  begin_end.first; iter != begin_end.second; ++iter)
        {
          iter->second->identify();
        }
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int
RawTowerBuilderCone::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void
RawTowerBuilderCone::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
    {
      std::cerr << PHWHERE << "Run Node missing, doing nothing." << std::endl;
      throw std::runtime_error("Failed to find Run node in RawTowerBuilderCone::CreateNodes");
    }

  PHCompositeNode *dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cerr << PHWHERE << "DST Node missing, doing nothing." << std::endl;
      throw std::runtime_error("Failed to find DST node in RawTowerBuilderCone::CreateNodes");
    }

  // Create the tower nodes on the tree
  _towers = new RawTowerContainer();
  TowerNodeName = "TOWER_" + detector;
  if (GroupID.length())
    {
      TowerNodeName = TowerNodeName + "_"+GroupID;
    }

  PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(_towers, TowerNodeName.c_str(), "PHObject");
  dstNode->addNode(towerNode);

  return;
}

pair<double, double>
RawTowerBuilderCone::get_etaphi(const double x, const double y, const double z)
{
  double eta;
  double phi;
  double radius;
  double theta;
  radius = sqrt(x * x + y * y);
  phi = atan2(y, x);
  theta = atan2(radius, z);
  eta = -log(tan(theta / 2.));
  return make_pair(eta, phi);
}


double
RawTowerBuilderCone::eta2theta(const double eta)
{
  double theta = 2.0 * atan( exp( -1.0 * eta ) );
  return theta;
}

double
RawTowerBuilderCone::theta2eta(const double theta)
{
  double eta = -log(tan(theta / 2.));
  return eta;
}


double
RawTowerBuilderCone::get_thetabin(const double theta)
{
  if ( theta < _thetamin || theta > _thetamax )
    {
//      cout << "Asking for bin for theta outside of theta range: " << theta << endl;
      return -1;
    }
  return floor( (theta-_thetamin)/_thetastep );
}


double
RawTowerBuilderCone::get_phibin(const double phi)
{
  double norm_phi = phi;
  if(phi < _phimin || phi > _phimax)
    {
      int nwraparound = -floor((phi-_phimin) * 0.5 / M_PI);
      norm_phi += 2*M_PI*nwraparound;
    }
  return floor( (norm_phi-_phimin)/_phistep );
}


