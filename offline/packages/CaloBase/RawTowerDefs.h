#ifndef CALOBASE_RAWTOWERDEFS_H
#define CALOBASE_RAWTOWERDEFS_H

#include <bitset>
#include <cstdlib>
#include <iostream>
#include <string>

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
    NONE = 0,
    CEMC = 1,
    HCALOUT = 2,
    HCALIN = 3,
    EEMC = 4,
    FEMC = 5,
    FHCAL = 6,
    DRCALO = 7,
    EHCAL = 8,
    EEMC_crystal = 9,
    EEMC_glass = 10,
    LFHCAL = 11,
    BECAL = 12,
    ZDC = 13,
    B0ECAL = 14,
    BWD_0 = 15,
    BWD_1 = 16,
    BWD_2 = 17,
    BWD_3 = 18,
    BWD_4 = 19
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

  /*! Extract tower index 1 of calorimeter tower from CaloTowerID with 3 indices
   */
  inline unsigned int
  decode_index1v2(const unsigned int calo_tower_id)
  {
    //     static unsigned int bitsIndex1 = 10;
    static unsigned int bitsIndex2 = 10;  // max 0x3FF (1023)
    static unsigned int bitsIndex3 = 4;   // max 0xF (15)

    //     std::cout << std::bitset<32>(calo_tower_id) << "\t index 1: " << ((calo_tower_id >> (bitsIndex2+bitsIndex3) & 0x3FF)) <<  "\t"<< std::bitset<32>((calo_tower_id >> (bitsIndex2+bitsIndex3))) << "\t"<< std::bitset<32>((calo_tower_id >> (bitsIndex2+bitsIndex3))& 0x3FF) <<std::endl;
    return (calo_tower_id >> (bitsIndex2 + bitsIndex3)) & 0x3FF;
  }

  /*! Extract tower index 2 of calorimeter tower from CaloTowerID with 3 indices
   */
  inline unsigned int
  decode_index2v2(const unsigned int calo_tower_id)
  {
    static unsigned int bitsIndex3 = 4;  // max 0xF (15)
                                         //     std::cout << std::bitset<32>(calo_tower_id) << "\t index 2: " << ((calo_tower_id >> (bitsIndex3) & 0x3FF)) <<  "\t"<< std::bitset<32>((calo_tower_id >> (bitsIndex3))) << "\t"<< std::bitset<32>((calo_tower_id >> (bitsIndex3))& 0x3FF) <<std::endl;
    return (calo_tower_id >> (bitsIndex3)) & 0x3FF;
  }

  /*! Extract tower index 3 of calorimeter tower from CaloTowerID with 3 indices
   */
  inline unsigned int
  decode_index3v2(const unsigned int calo_tower_id)
  {
    //     std::cout << std::bitset<32>(calo_tower_id) << "\t index 3: " << (calo_tower_id & 0xF) <<  "\t"<< std::bitset<32>((calo_tower_id & 0xF)) <<std::endl;
    return calo_tower_id & 0xF;
  }

  /*! Returns CaloTowerID for given calorimeter ID, tower index 1, tower index 2 and tower index 3
   */
  inline RawTowerDefs::keytype
  encode_towerid(const CalorimeterId calo_id, const unsigned int tower_index_1,
                 const unsigned int tower_index_2, const unsigned int tower_index_3)
  {
    RawTowerDefs::keytype calo_tower_id = 0;

    //     static unsigned int bitsIndex1 = 10; // max 0x3FF (1023)
    static unsigned int bitsIndex2 = 10;  // max 0x3FF (1023)
    static unsigned int bitsIndex3 = 4;   // max 0xF (15)

    if (tower_index_1 < 0x3FF && tower_index_2 < 0x3FF && tower_index_3 < 0xF)
    {
      calo_tower_id = (calo_id << RawTowerDefs::tower_idbits) + (tower_index_1 << (bitsIndex2 + bitsIndex3)) + (tower_index_2 << bitsIndex3) + tower_index_3;
    }
    else
    {
      std::cout << "too large index1 and/or index2; index1: "
                << tower_index_1 << " (max val " << 0x3FF << ")"
                << ", index2: "
                << tower_index_2 << " (max val " << 0x3FF << ")"
                << ", index3: "
                << tower_index_3 << " (max val " << 0xF << ")" << std::endl;
      exit(1);
    }
    //     std::cout << std::bitset<32>(calo_tower_id) << "\t" << std::bitset<8>(calo_id) << "\t"<< tower_index_1 << "\t"<<  std::bitset<10>(tower_index_1) << "\t"<< tower_index_2 << "\t"<<  std::bitset<10>(tower_index_2) << "\t"<< tower_index_3<< "\t"<<  std::bitset<4>(tower_index_3)  << std::endl;
    //     decode_index1v2(calo_tower_id);
    //     decode_index2v2(calo_tower_id);
    //     decode_index3v2(calo_tower_id);
    return calo_tower_id;
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

    case DRCALO:
      return "DRCALO";
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

    case EHCAL:
      return "EHCAL";
      break;

    case FEMC:
      return "FEMC";
      break;

    case FHCAL:
      return "FHCAL";
      break;

    case BECAL:
      return "BECAL";
      break;

    case EEMC_crystal:
      return "EEMC_crystal";
      break;

    case EEMC_glass:
      return "EEMC_glass";
      break;

    case LFHCAL:
      return "LFHCAL";
      break;

    case ZDC:
      return "ZDC";
      break;

    case B0ECAL:
      return "B0ECAL";
      break;

    case BWD_0:
      return "BWD_0";
      break;

    case BWD_1:
      return "BWD_1";
      break;

    case BWD_2:
      return "BWD_2";
      break;

    case BWD_3:
      return "BWD_3";
      break;

    case BWD_4:
      return "BWD_4";
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

    else if (caloname == "DRCALO")
      return DRCALO;

    else if (caloname == "HCALIN")
      return HCALIN;

    else if (caloname == "HCALOUT")
      return HCALOUT;

    else if (caloname == "EEMC")
      return EEMC;

    else if (caloname == "EHCAL")
      return EHCAL;

    else if (caloname == "FEMC")
      return FEMC;

    else if (caloname == "FHCAL")
      return FHCAL;

    else if (caloname == "EEMC_crystal")
      return EEMC_crystal;

    else if (caloname == "EEMC_glass")
      return EEMC_glass;

    else if (caloname == "LFHCAL")
      return LFHCAL;

    else if (caloname == "BECAL")
      return BECAL;

    else if (caloname == "ZDC")
      return ZDC;

    else if (caloname == "B0ECAL")
      return B0ECAL;

    else if (caloname == "BWD_0")
      return BWD_0;

    else if (caloname == "BWD_1")
      return BWD_1;

    else if (caloname == "BWD_2")
      return BWD_2;

    else if (caloname == "BWD_3")
      return BWD_3;

    else if (caloname == "BWD_4")
      return BWD_4;

    else
    {
      std::cout << "Invalid calorimeter name " << caloname
                << " passed to RawTowerDefs::convert_name_to_caloid" << std::endl;
      exit(1);
    }
  }

}  // end namespace RawTowerDefs

#endif
