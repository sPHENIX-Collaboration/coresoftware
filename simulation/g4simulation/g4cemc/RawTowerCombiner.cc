#include "RawTowerCombiner.h"
#include "RawTowerContainer.h"
#include "RawTowerGeomv2.h"
#include "RawTowerv1.h"
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <boost/foreach.hpp>

#include <iostream>
#include <stdexcept>
#include <map>

using namespace std;

RawTowerCombiner::RawTowerCombiner(const std::string& name):
  SubsysReco(name),
  chkenergyconservation(0),
  iphibins(0),
  ietabins(0),
  _timer( PHTimeServer::get()->insert_new(name) )
{}

int
RawTowerCombiner::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  PHCompositeNode *runNode;
  dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
      exit(1);
    }
  runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
    {
      std::cout << PHWHERE << "RUN Node missing, doing nothing." << std::endl;
      exit(1);
    }
  TowerGeomNodeName = "TOWERGEOM_" + detector;
  TowerNodeName = "TOWER_"  + detector;
  // loop over geometry nodes to make sure they are compatible
  int ifirst = 1;
  double minradius = 1.e7;
  double thickness = 0.;
  RawTowerGeom *reftowergeom = NULL;
  BOOST_FOREACH(auto &p, dettuple)
    {
      string det = p.first;
      cout << det << endl;
      std::string geonodename = "TOWERGEOM_" + det;
      RawTowerGeom *towergeo = findNode::getClass<RawTowerGeom>(topNode, geonodename.c_str());
      if (!towergeo)
	{
	  std::cerr << PHWHERE << " " << geonodename << " Node missing, doing nothing." << std::endl;
	  throw std::runtime_error("Failed to find " + geonodename + " node in RawTowerCombiner::InitRun");
	}
      if (ifirst)
	{
	  reftowergeom = towergeo;
	  ifirst = 0;
	}
      else
	{
	  if (!CompareGeometries(reftowergeom,towergeo))
	    {
	      cout << "ERROR in comparing geometry objects for " << det << endl;
	      exit(1);
	    }
	}
      if (minradius > towergeo->get_radius())
	{
	  minradius =  towergeo->get_radius();
	}
      thickness += towergeo->get_thickness();
    }
  if (reftowergeom)
    {
      RawTowerGeom *rawtowergeom = new RawTowerGeomv2(reftowergeom);
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(rawtowergeom, TowerGeomNodeName.c_str(), "PHObject");
      runNode->addNode(newNode);
      // overwrite changed radius and thickness
      rawtowergeom->set_radius(minradius);
      rawtowergeom->set_thickness(thickness);
      ietabins = rawtowergeom->get_etabins();
      iphibins = rawtowergeom->get_phibins();
    }
  else
    {
      cout << "no input towers found" << endl;
      exit(1);
    }
  RawTowerContainer *towers = new RawTowerContainer( RawTowerDefs::convert_name_to_caloid( detector ) );
  PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(towers, TowerNodeName.c_str(), "PHObject");
  dstNode->addNode(towerNode);
  // stash the pointers to our input tower nodes
  BOOST_FOREACH(auto &p, dettuple)
    {
      string nodename = "TOWER_" + p.first;
      towers = findNode::getClass<RawTowerContainer>(topNode,nodename.c_str());
      boost::get<2>(p.second) = towers; 
      //	     towercontainers.push_back(towers);
    }
  BOOST_FOREACH(auto &p, dettuple)
    {
      cout << boost::get<0>(p.second) << endl;
      cout << boost::get<1>(p.second) << endl;
      cout << boost::get<2>(p.second) << endl;
      cout << boost::get<3>(p.second) << endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}



int 
RawTowerCombiner::process_event(PHCompositeNode *topNode)
{
  RawTowerContainer *outtowercontainer = findNode::getClass<RawTowerContainer>(topNode,TowerNodeName.c_str());
  // this is brute force - just loop over all phi/eta bins and get the towers and add them up
  for (int i=0; i< ietabins; i++)
    {
      for (int j=0; j<iphibins; j++)
	{
          RawTower *outtower = NULL;
	  BOOST_FOREACH(auto &p, dettuple)
	    {
	      RawTowerContainer *twrcont =  boost::get<2>(p.second);
	      RawTower *twr = twrcont->getTower(i,j);
	      if (twr)
		{
		  if (!outtower)
		    {
		      outtower = new RawTowerv1(i,j);
		    }
		  outtower->add_ecell(boost::get<1>(p.second),twr->get_energy()/boost::get<3>(p.second));
		}
	    }
	  if (outtower)
	    {
	      	      outtowercontainer->AddTower(i,j,outtower);
	    }

	}
    }
  if (chkenergyconservation)
    {
      float esum = 0;
      BOOST_FOREACH(auto &p, dettuple)
	{
	  RawTowerContainer *twrcont =  boost::get<2>(p.second);
	  esum += twrcont->getTotalEdep()/boost::get<3>(p.second);
	  cout << boost::get<0>(p.second) << " energy before: " 
	       << twrcont->getTotalEdep() << " after sampling correction: " 
	       << twrcont->getTotalEdep()/boost::get<3>(p.second) << endl;
	}
      cout << "added: " << esum << " saved: " << outtowercontainer->getTotalEdep() << endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

int 
RawTowerCombiner::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

bool
RawTowerCombiner::CompareGeometries(RawTowerGeom *geo1, RawTowerGeom *geo2)
{
  if (geo1->get_phibins() != geo2->get_phibins())
    {
      cout << "diff in phibins " << geo1->get_phibins() << ", " << geo2->get_phibins() << endl;
      return false;
    }
  if (geo1->get_phistep() != geo2->get_phistep())
    {
      cout << "diff in phistep " << geo1->get_phistep() << ", " << geo2->get_phistep() << endl;
      return false;
    }
  if (geo1->get_phimin() != geo2->get_phimin())
    {
      cout << "diff in phimin " << geo1->get_phimin() << ", " << geo2->get_phimin() << endl;
      return false;
    }
  if (geo1->get_etabins() != geo2->get_etabins())
    {
      cout << "diff in etabins " << geo1->get_etabins() << ", " << geo2->get_etabins() << endl;
      return false;
    }
  for (int i=0; i<geo2->get_etabins(); i++)
    {
      if (geo1->get_etabounds(i) != geo2->get_etabounds(i))
	{
	  cout << "diff in eta bounds for bin " << i << endl; 
      return false;
	}
    }
  return true;
}

void
RawTowerCombiner::AddInputDetector(const std::string &d, const int index, const double sf)
{
  dettuple[d] = boost::make_tuple(d,index,(RawTowerContainer *)NULL,sf);
  return;
}
