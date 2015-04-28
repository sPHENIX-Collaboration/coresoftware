#include "PHG4CylinderEventAction.h"
#include <g4main/PHG4HitContainer.h>

#include <fun4all/getClass.h>

#include <Geant4/G4Event.hh>

#include <boost/foreach.hpp>

using namespace std;

//___________________________________________________
PHG4CylinderEventAction::PHG4CylinderEventAction( PHCompositeNode *node, const string &name ):
  topNode(node)
{
  nodename_set.insert(name);
}


void
PHG4CylinderEventAction::AddNode(const std::string &name)
{
  nodename_set.insert(name);
}

//___________________________________________________
void PHG4CylinderEventAction::EndOfEventAction(const G4Event* evt)
{
  BOOST_FOREACH(string node, nodename_set)
    {
      PHG4HitContainer* cylinder_hits =  findNode::getClass<PHG4HitContainer>( topNode , node.c_str());
      if (cylinder_hits)
	{
	  cylinder_hits->RemoveZeroEDep();
	}
    }
}
