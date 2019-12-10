#ifndef CALOBASE_RAWTOWERDEFS_H
#define CALOBASE_RAWTOWERDEFS_H

#include <iostream>
#include <string>

#if !defined(__CINT__) || defined(__CLING__)
#include <cstdlib>
#else
#include <stdlib.h>
#endif

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
static unsigned int tower_idbits = sizeof(keytype) * 8 - calo_idbits;
static unsigned int index1_idbits = tower_idbits / 2;

/*! Enum with all available calorimeter IDs. This enum can be extended up to 254 entries.
   * If adding new CalorimeterIDs, please also add them to the decode_caloname function below.
   */
enum CalorimeterId
{
  NONE,
  CEMC,
  HCALOUT,
  HCALIN,
  EEMC,
  FEMC,
  FHCAL
};

/*! Returns CaloTowerID for given calorimeter ID, tower index 1, and tower index 2
   */
inline RawTowerDefs::keytype
encode_towerid(const CalorimeterId calo_id, const unsigned int tower_index_1,
               const unsigned int tower_index_2)
{
  RawTowerDefs::keytype calo_tower_id = 0;

  if (tower_index_1 < 0xFFF && tower_index_2 < 0xFFF)
  {
    calo_tower_id = (calo_id << RawTowerDefs::tower_idbits) + (tower_index_1 << RawTowerDefs::index1_idbits) + tower_index_2;
  }
  else
  {
    std::cout << "too large index1 and/or index2; index1: "
              << tower_index_1 << " (max val " << 0xFFF << ")"
              << ", index2: "
              << tower_index_2 << " (max val " << 0xFFF << ")" << std::endl;
    exit(1);
  }

  return calo_tower_id;
}

/*! Returns CaloTowerID for given calorimeter ID, tower index
   */
inline RawTowerDefs::keytype
encode_towerid(const CalorimeterId calo_id, const unsigned int tower_index)
{
  RawTowerDefs::keytype calo_tower_id = 0;

  if (tower_index < 0xFFFFFF)
  {
    calo_tower_id = (calo_id << RawTowerDefs::tower_idbits) + tower_index;
  }
  else
  {
    std::cout << "too large index; index: " << tower_index
              << " (max val " << 0xFFFFFF << ")" << std::endl;
    exit(1);
  }

  return calo_tower_id;
}

/*! Extract calorimeter ID from CaloTowerID
   */
inline CalorimeterId
decode_caloid(const unsigned int calo_tower_id)
{
  return static_cast<CalorimeterId>((calo_tower_id >> RawTowerDefs::tower_idbits) & 0xFFF);
}

/*! Extract tower index of calorimeter tower from CaloTowerID
   */
inline unsigned int
decode_index(const unsigned int calo_tower_id)
{
  return (calo_tower_id) &0xFFFFFF;
}

/*! Extract tower index 1 of calorimeter tower from CaloTowerID
   */
inline unsigned int
decode_index1(const unsigned int calo_tower_id)
{
  return (calo_tower_id >> RawTowerDefs::index1_idbits) & 0xFFF;
}

/*! Extract tower index 2 of calorimeter tower from CaloTowerID
   */
inline unsigned int
decode_index2(const unsigned int calo_tower_id)
{
  return calo_tower_id & 0xFFF;
}

/*! Convert calorimeter ID to name string
   */
inline std::string
convert_caloid_to_name(const RawTowerDefs::CalorimeterId calo_id)
{
  switch (calo_id)
  {
  case NONE:
    return "NONE";
    break;

  case CEMC:
    return "CEMC";
    break;

  case HCALIN:
    return "HCALIN";
    break;

  case HCALOUT:
    return "HCALOUT";
    break;

  case EEMC:
    return "EEMC";
    break;

  case FEMC:
    return "FEMC";
    break;

  case FHCAL:
    return "FHCAL";
    break;

  default:
    std::cout
        << "Invalid calorimeter ID passed to RawTowerDefs::convert_caloid_to_name"
        << std::endl;
    exit(1);
  }
}

/*! Convert name string to calorimeter ID
   */
inline RawTowerDefs::CalorimeterId
convert_name_to_caloid(const std::string &caloname)
{
  if (caloname == "NONE")
    return NONE;

  else if (caloname == "CEMC")
    return CEMC;

  else if (caloname == "HCALIN")
    return HCALIN;

  else if (caloname == "HCALOUT")
    return HCALOUT;

  else if (caloname == "EEMC")
    return EEMC;

  else if (caloname == "FEMC")
    return FEMC;

  else if (caloname == "FHCAL")
    return FHCAL;

  else
  {
    std::cout << "Invalid calorimeter name " << caloname
              << " passed to RawTowerDefs::convert_name_to_caloid" << std::endl;
    exit(1);
  }
}

}  // end namespace RawTowerDefs

#endif
