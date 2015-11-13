#ifndef RAWTOWERDEFS_H
#define RAWTOWERDEFS_H

#include <iostream>

/*! Namespace with functions to encode / decode CaloTowerID. The highest 8 bits of the tower ID encode a unique ID
 * for the calorimeter the tower is in. The lower 24 bits uniquely identify the tower within a calorimeter.
 *
 */
namespace RawTowerDefs
{
  /*! Define data type of unique tower ID, i.e. for CaloTowerID
   */
  typedef unsigned int keytype;

  /*! Bit ranges for encoding calorimeter ID and tower indices in combined tower ID
   */
  static unsigned int calo_idbits = 8;
  static unsigned int tower_idbits = sizeof(keytype)*8 - calo_idbits;
  static unsigned int index1_idbits = tower_idbits/2;

  /*! Enum with all available calorimeter IDs. This enum can be extended up to 254 entries.
   * If adding new CalorimeterIDs, please also add them to the decode_caloname function below.
   */
  enum CalorimeterId {
    NONE,
    CEMC,
    HCALOUT,
    HCALIN,
    EEMC,
    FEMC,
    FHCAL,
  };

  /*! Returns CaloTowerID for given calorimeter ID, tower index 1, and tower index 2
   */
  inline RawTowerDefs::keytype encode_towerid( const CalorimeterId calo_id , const unsigned int tower_index_1 = 0 , const unsigned int tower_index_2 = 0)
  {
    RawTowerDefs::keytype calo_tower_id = 0;

    if ( calo_id < 0xFF && tower_index_1 < 0xFFF && tower_index_2 < 0xFFF )
      {
	calo_tower_id = ( calo_id << RawTowerDefs::tower_idbits )
	  + ( tower_index_1 << RawTowerDefs::index1_idbits ) + tower_index_2;
      }
    else
      {
	std::cout << "too large caloid, index1 and/or index2; caloid: " << calo_id << " (max val " << 0xFF << ")"
		  << ", index1: " << tower_index_1 << " (max val " << 0xFFF << ")"
		  << ", index2: " << tower_index_2 << " (max val " << 0xFFF << ")" << std::endl;
	exit(1);
      }

    return calo_tower_id;
  }

  /*! Extract calorimeter ID from CaloTowerID
   */
  inline unsigned int decode_caloid( const unsigned int calo_tower_id )
  {
    return ( calo_tower_id >> RawTowerDefs::tower_idbits ) &0xFFF;
  }

  /*! Extract tower index 1 of calorimeter tower from CaloTowerID
   */
  inline unsigned int decode_index1( const unsigned int calo_tower_id )
  {
    return ( calo_tower_id >> RawTowerDefs::index1_idbits ) &0xFFF;
  }

  /*! Extract tower index 2 of calorimeter tower from CaloTowerID
   */
  inline unsigned int decode_index2( const unsigned int calo_tower_id )
  {
    return calo_tower_id &0xFFF;
  }

}

#endif
