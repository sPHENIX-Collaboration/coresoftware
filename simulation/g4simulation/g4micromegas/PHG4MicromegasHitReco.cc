// this is the new trackbase version

/*!
 * \file PHG4MicromegasHitReco.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "PHG4MicromegasHitReco.h"

#include <micromegas/CylinderGeomMicromegas.h>
#include <micromegas/MicromegasDefs.h>
#include <micromegas/MicromegasTile.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitTruthAssoc.h>

#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Utils.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/phool.h>

//___________________________________________________________________________
PHG4MicromegasHitReco::PHG4MicromegasHitReco(const std::string &name, const std::string& detector)
  : SubsysReco(name)
  , PHParameterInterface(name)
  , m_detector(detector)
{ SetDefaultParameters(); }

//___________________________________________________________________________
int PHG4MicromegasHitReco::InitRun(PHCompositeNode *topNode)
{

  PHNodeIterator iter(topNode);

  // setup tiles
  setup_tiles( topNode );

  // setup timing window
  m_tmin = get_double_param("tmin");
  m_tmax = get_double_param("tmax");

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4MicromegasHitReco::process_event(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________
void PHG4MicromegasHitReco::SetDefaultParameters()
{
  std::cout << "PHG4MicromegasHitReco: Setting timing window defaults to tmin = -5000 and  tmax = 5000 ns" << std::endl;
  set_default_double_param("tmin", -5000.0);
  set_default_double_param("tmax", 5000.0);
  return;
}

//___________________________________________________________________________
void PHG4MicromegasHitReco::setup_tiles(PHCompositeNode* topNode)
{

  // get geometry
  std::string geonodename = "CYLINDERGEOM_" + m_detector;
  auto geonode = findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename);
  if (!geonode)
  {
    std::cout << PHWHERE << "Could not locate geometry node " << geonodename << std::endl;
    exit(1);
  }

  // get cylinders
  PHG4CylinderGeomContainer::ConstRange range = geonode->get_begin_end();
  for( auto iter = range.first; iter != range.second; ++iter )
  {
    std::cout << "PHG4MicromegasHitReco::setup_tiles - processing layer " << iter->first << std::endl;

    // for now, put a single tile at 0,0, with 50cm along z and 25 cm along phi
    // TODO: allow tiles to be setup from the macro and propagated to the geometry here, rather than hardcoded
    static_cast<CylinderGeomMicromegas*>(iter->second)->set_tiles( { MicromegasTile( 0, 0, 25, 50 )} );
  }

}
