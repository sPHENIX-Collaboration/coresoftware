// This is the new trackbase container version

/*!
 * \file PHG4MicromegasDigitizer.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "PHG4MicromegasDigitizer.h"

// Move to new storage containers
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <cassert>
#include <set>

//____________________________________________________________________________
PHG4MicromegasDigitizer::PHG4MicromegasDigitizer(const std::string &name)
  : SubsysReco(name)
{
  // fixed seed is handled in this funtcion
  const unsigned int seed = PHRandomSeed();
  std::cout << Name() << " random seed: " << seed << std::endl;

  // initialize rng
  m_rng.reset( gsl_rng_alloc(gsl_rng_mt19937) );
  gsl_rng_set( m_rng.get(), seed );
}

//____________________________________________________________________________
int PHG4MicromegasDigitizer::InitRun(PHCompositeNode *topNode)
{
  // TODO: set default values for m_max_adc, energy_scale and threshold
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________
int PHG4MicromegasDigitizer::process_event(PHCompositeNode *topNode)
{
  // Get Nodes
  // Get the TrkrHitSetContainer node
  auto trkrhitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  assert( trkrhitsetcontainer );

  // Digitization
  // We want all hitsets for the Mvtx
  const auto hitset_range = trkrhitsetcontainer->getHitSets(TrkrDefs::TrkrId::micromegasId);

  for( auto hitset_it = hitset_range.first; hitset_it != hitset_range.second; ++hitset_it )
  {

    // get key and layer
    const TrkrDefs::hitsetkey hitsetkey = hitset_it->first;
    const int layer = TrkrDefs::getLayer(hitsetkey);
    if (Verbosity() > 1) std::cout << "PHG4MicromegasDigitizer::process_event - hitsetkey: " << hitsetkey << " layer:" << layer << std::endl;

    // get all of the hits from this hitset
    TrkrHitSet* hitset = hitset_it->second;
    TrkrHitSet::ConstRange hit_range = hitset->getHits();

    // keep track of hits to be removed
    std::set<TrkrDefs::hitkey> removed_keys;

    // loop over hits
    for( auto hit_it = hit_range.first; hit_it != hit_range.second; ++hit_it )
    {

      TrkrHit *hit = hit_it->second;

      // Convert the signal value to an ADC value and write that to the hit
      unsigned int adc = hit->getEnergy()/m_energy_scale;

      // TODO: add noise
      // compare to max adc
      if( m_max_adc > 0 && adc > m_max_adc) adc = m_max_adc;

      // assign to hit
      hit->setAdc(adc);

      // Remove the hits with energy under threshold
      if( m_energy_threshold > 0 && hit->getEnergy() < m_energy_threshold )
      { removed_keys.insert( hit_it->first ); }

      // TODO: should also add ADC threshold ?

    }

    // remove hits
    for( const auto& key:removed_keys )
    { hitset->removeHit(key); }

  }
  return Fun4AllReturnCodes::EVENT_OK;
}
