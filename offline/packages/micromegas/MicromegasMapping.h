#ifndef MICROMEGAS_MICROMEGASMAPPING_H
#define MICROMEGAS_MICROMEGASMAPPING_H

/*!
 * \file MicromegasMapping.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <trackbase/TrkrDefs.h>

#include <map>
#include <string>
#include <vector>

/// micromegas mapping
class MicromegasMapping
{
  public:

  /// constructor
  MicromegasMapping();

  /// get list of fee ids
  std::vector<int> get_fee_id_list() const;
  
  /// get hitsetkey from fiber_id (fee_id)
  TrkrDefs::hitsetkey get_hitsetkey( int /*fee_id*/ ) const;

  /// get detector name (saclay) from fiber_id (fee_id)
  /** saclay detector name are of type MxxP and MxxZ, with xx the module number */
  std::string get_detname_saclay( int /*fee_id*/) const;

  /// get detector name (saclay) from fiber_id (fee_id)
  /** sphenix detector name are of type SWP, SWZ, etc. */
  std::string get_detname_sphenix( int /*fee_id*/) const;

  /// get physical strip number from channel_id
  /**
   * physical strip number correspond to a position in the
   * detector's geant implementation, with increasing number corresponding to strips further away from the detector's edge,
   * as defined in CylinderGeomMicromegas
   */
  int get_physical_strip( int /*fee_id*/, int /*channel_id*/) const;

  private:

  /// contains all relevant detector information
  /** this effectively implements mapping between fee_id as defined in EDBC,    * detector names (in both Saclay and sPHENIX convention),
   * and hitsetkey which is the detector unique identifier
   */
  class DetectorId
  {
    public:

    /// constructor
    DetectorId( int fee_id, TrkrDefs::hitsetkey hitsetkey, const std::string& detname_saclay, const std::string& detname_sphenix ):
      m_fee_id( fee_id ),
      m_hitsetkey( hitsetkey ),
      m_detname_saclay( detname_saclay ),
      m_detname_sphenix( detname_sphenix )
    {}

    /// fee_id
    int m_fee_id = 0;

    /// hitset key
    TrkrDefs::hitsetkey m_hitsetkey = 0;

    /// detector name
    std::string m_detname_saclay;

    /// detector name
    std::string m_detname_sphenix;
  };

  /// list of defined detectors
  std::vector<DetectorId> m_detectors;

  /// map detector_id to fee_id
  std::map<int, DetectorId> m_detector_map;

};

#endif
