/*!
 * \file PHG4MicromegasDigitizer.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "PHG4MicromegasDigitizer.h"

// Move to new storage containers
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitTruthAssoc.h>

#include <phparameter/PHParameterInterface.h>  // for PHParameterInterface

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>                       // for gsl_rng_alloc, gsl_rng...

#include <cassert>
#include <iostream>                            // for operator<<, basic_ostream
#include <set>
#include <utility>                             // for pair

namespace
{

  // local version of std::clamp, which is only available for c++17
  template<class T>
    constexpr const T& clamp( const T& v, const T& lo, const T& hi )
  { return (v < lo) ? lo : (hi < v) ? hi : v; }

}
//____________________________________________________________________________
PHG4MicromegasDigitizer::PHG4MicromegasDigitizer(const std::string &name)
  : SubsysReco(name)
  , PHParameterInterface(name)
{
  // initialize rng
  const unsigned int seed = PHRandomSeed();
  m_rng.reset( gsl_rng_alloc(gsl_rng_mt19937) );
  gsl_rng_set( m_rng.get(), seed );

  InitializeParameters();
}

//____________________________________________________________________________
int PHG4MicromegasDigitizer::InitRun(PHCompositeNode */*topNode*/)
{

  UpdateParametersWithMacro();

  // load parameters
  m_adc_threshold = get_double_param( "micromegas_adc_threshold" );
  m_enc = get_double_param( "micromegas_enc" );
  m_pedestal = get_double_param( "micromegas_pedestal" );
  m_volts_per_charge = get_double_param( "micromegas_volts_per_charge" );

  // printout
  std::cout
    << "PHG4MicromegasDigitizer::InitRun\n"
    << " m_adc_threshold: " << m_adc_threshold << " electrons\n"
    << " m_enc: " << m_enc << " electrons\n"
    << " m_pedestal: " << m_pedestal << " electrons\n"
    << " m_volts_per_charge: " << m_volts_per_charge << " mV/fC\n"
    << std::endl;

  // threshold is effectively applied on top of pedestal
  m_adc_threshold += m_pedestal;

  /*
   * Factor that convertes charge in a voltage in each z bin
   * the scale up factor of 2.4 is meant to account for shaping time (80ns)
   * it only applies to the signal
   * see: simulations/g4simulations/g4tpc/PHG4TpcDigitizer::InitRun
   */
  m_volt_per_electron_signal = m_volts_per_charge * 1.602e-4 * 2.4;
  m_volt_per_electron_noise = m_volts_per_charge * 1.602e-4;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________
int PHG4MicromegasDigitizer::process_event(PHCompositeNode *topNode)
{
  // Get Nodes
  // Get the TrkrHitSetContainer node
  auto trkrhitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  assert( trkrhitsetcontainer );

  // Get the TrkrHitTruthAssoc node
  auto hittruthassoc = findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");

  // get all micromegas hitsets
  const auto hitset_range = trkrhitsetcontainer->getHitSets(TrkrDefs::TrkrId::micromegasId);
  for( auto hitset_it = hitset_range.first; hitset_it != hitset_range.second; ++hitset_it )
  {

    // get key
    const TrkrDefs::hitsetkey hitsetkey = hitset_it->first;

    // get all of the hits from this hitset
    TrkrHitSet* hitset = hitset_it->second;
    TrkrHitSet::ConstRange hit_range = hitset->getHits();

    // keep track of hits to be removed
    std::set<TrkrDefs::hitkey> removed_keys;

    // loop over hits
    for( auto hit_it = hit_range.first; hit_it != hit_range.second; ++hit_it )
    {
      // store key and hit
      const TrkrDefs::hitkey& key = hit_it->first;
      TrkrHit *hit = hit_it->second;

      // get energy (electrons)
      const double signal = hit->getEnergy();
      const double noise = add_noise();

      // convert to mV
      const double voltage = (m_pedestal + noise)*m_volt_per_electron_noise + signal*m_volt_per_electron_signal;

      // compare to threshold
      if( voltage > m_adc_threshold*m_volt_per_electron_noise )
      {
        // keep hit, update adc
        hit->setAdc( clamp<uint>( voltage*m_adc_per_volt, 0, 1023 ) );
      } else {
        // mark hit as removable
        removed_keys.insert( key );
      }

    }

    // remove hits
    for( const auto& key:removed_keys )
    {
      hitset->removeHit(key);
      if( hittruthassoc ) hittruthassoc->removeAssoc(hitsetkey, key);
    }

  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________
void PHG4MicromegasDigitizer::SetDefaultParameters()
{
  // all values taken from TPC sampa chips (simulations/g4simulations/g4tpc/PHG4TpcDigitizer)
  set_default_double_param( "micromegas_enc", 670 );
  set_default_double_param( "micromegas_adc_threshold", 2680 );
  set_default_double_param( "micromegas_pedestal", 50000 );
  set_default_double_param( "micromegas_volts_per_charge", 20 );
}

//___________________________________________________________________________
double PHG4MicromegasDigitizer::add_noise() const
{ return gsl_ran_gaussian( m_rng.get(), m_enc); }
