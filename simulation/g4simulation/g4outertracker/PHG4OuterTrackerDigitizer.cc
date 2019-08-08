// This is the new trackbase container version

#include "PHG4OuterTrackerDigitizer.h"
#include <outertracker/CylinderGeomOuterTracker.h>
#include <outertracker/OuterTrackerDefs.h>

// Move to new storage containers
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>                      // for TrkrHit
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                     // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHNode.h>                           // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>                            // for PHWHERE

#include <gsl/gsl_rng.h>                            // for gsl_rng_alloc

#include <cstdlib>                                 // for exit
#include <iostream>

using namespace std;

PHG4OuterTrackerDigitizer::PHG4OuterTrackerDigitizer(const string &name)
  : SubsysReco(name)
{
  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  cout << Name() << " random seed: " << seed << endl;
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(RandomGenerator, seed);

  if (Verbosity() > 0)
    cout << "Creating PHG4OuterTrackerDigitizer with name = " << name << endl;
}

PHG4OuterTrackerDigitizer::~PHG4OuterTrackerDigitizer()
{
  gsl_rng_free(RandomGenerator);
}

int PHG4OuterTrackerDigitizer::InitRun(PHCompositeNode *topNode)
{
  //-------------
  // Add Hit Node
  //-------------
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  CalculateOuterTrackerADCScale(topNode);
  
  //----------------
  // Report Settings
  //----------------

  if (Verbosity() > 0)
  {
    cout << "====================== PHG4OuterTrackerDigitizer::InitRun() =====================" << endl;
    cout << " Max ADC " << _max_adc << " energy scale " << _energy_scale << endl;
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4OuterTrackerDigitizer::process_event(PHCompositeNode *topNode)
{
  // This code now only does the OuterTracker
  DigitizeOuterTracker(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}


void PHG4OuterTrackerDigitizer::CalculateOuterTrackerADCScale(PHCompositeNode *topNode)
{
  // defaults to 8-bit ADC, short-axis MIP placed at 1/4 dynamic range

  PHG4CylinderGeomContainer *geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_OuterTracker");
  if (!geom_container) 
    {
      cout << "Di not find CYLINDERGEOM_OuterTracker node" << endl;
      return;
    }

  if (Verbosity())
    cout << "Found CYLINDERGEOM_OuterTracker node" << endl;

PHG4CylinderGeomContainer::ConstRange layerrange =  geom_container->get_begin_end();
  for (PHG4CylinderGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter)
  {
    int layer = layeriter->first;

    CylinderGeomOuterTracker *layergeom = dynamic_cast<CylinderGeomOuterTracker *>(geom_container->GetLayerGeom(layer));

    double thickness = layergeom->get_outer_radius() - layergeom->get_inner_radius();

    double mip_e = 0.003876 * thickness;

    if (Verbosity())
      cout << "mip_e = " << mip_e << endl;
    _max_adc = 255;
    _energy_scale = mip_e / 64;
  }

  return;
}


void PHG4OuterTrackerDigitizer::DigitizeOuterTracker(PHCompositeNode *topNode)
{
  //----------
  // Get Nodes
  //----------

  // new containers
  //=============
  // Get the TrkrHitSetContainer node
  TrkrHitSetContainer *trkrhitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!trkrhitsetcontainer)
  {
    cout << "Could not locate TRKR_HITSET node, quit! " << endl;
    exit(1);
  }

  // Digitization

  // We want all hitsets for the OuterTracker
  TrkrHitSetContainer::ConstRange hitset_range = trkrhitsetcontainer->getHitSets(TrkrDefs::TrkrId::outertrackerId);
  for (TrkrHitSetContainer::ConstIterator hitset_iter = hitset_range.first;
       hitset_iter != hitset_range.second;
       ++hitset_iter)
  {
    // we have an itrator to one TrkrHitSet for the otrack from the trkrHitSetContainer
    // get the hitset key so we can find the layer
    TrkrDefs::hitsetkey hitsetkey = hitset_iter->first;
    int layer = TrkrDefs::getLayer(hitsetkey);
    if (Verbosity() > 1) cout << "PHG4OuterTrackerDigitizer: found hitset with key: " << hitsetkey << " in layer " << layer << endl;

    // get all of the hits from this hitset
    TrkrHitSet *hitset = hitset_iter->second;
    TrkrHitSet::ConstRange hit_range = hitset->getHits();
    for (TrkrHitSet::ConstIterator hit_iter = hit_range.first;
         hit_iter != hit_range.second;
         ++hit_iter)
    {
      TrkrHit *hit = hit_iter->second;

      // Convert the signal value to an ADC value and write that to the hit
      unsigned int adc = hit->getEnergy() / _energy_scale;
      if (adc > _max_adc) adc = _max_adc;
      hit->setAdc(adc);

      if (Verbosity() > 0) cout << "    PHG4OuterTrackerDigitizer: found hit with key: " << hit_iter->first << " and signal " << hit->getEnergy() << " and adc " << adc << endl;
    }
  }

  // end new containers
  //===============

  return;
}
