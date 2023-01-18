/*
sEPD Geometry class

Tristan Protzman
tlprotzman@gmail.com
December 15, 2022

Converts from r/phi index to location in r/phi space

Key format:
100000000000000000000_2
211111111110000000000
098765432109876543210
*/

#ifndef EPD_GEOM_H
#define EPD_GEOM_H

#include <map>
#include <utility>


// sEPD geometry class
class EPD_GEOM {

public:
  const uint NUM_TOWERS = 768;
  const uint MAX_R = 32;
  const uint MAX_PHI = 24;

  EPD_GEOM();
  ~EPD_GEOM();

  uint side_r_phi_to_id(uint side, uint r_index, uint phi_index);
  uint side_sector_tile_to_id(uint side, uint sector, uint tile);
  std::tuple<uint, uint, uint> id_to_side_r_phi(uint id);
  std::tuple<uint, uint, uint> id_to_side_sector_tile(uint id);
  float r(uint id);
  float r_from_side_r_phi(uint side, uint r_index, uint phi_index);
  float r_from_side_sector_tile(uint side, uint sector, uint tile);
  float phi(uint id);
  float phi_from_side_r_phi(uint side, uint r_index, uint phi_index);
  float phi_from_side_sector_tile(uint side, uint sector, uint tile);
  float z(uint id);
  float z_from_side_r_phi(uint side, uint r_inxed, uint phi_index);
  float z_from_side_sector_tile(uint side, uint sector, uint tile);
  
  
  bool test_id_mapping();

private:
  std::map<int, float> r_map;
  std::map<int, float> phi_map;
  std::map<int, float> z_map;

  // const float NAIVE_R_INDEX[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  const float NAIVE_R_LOC[16] = {6.8, 11.2, 15.6, 20.565, 26.095, 31.625, 37.155, 42.685, 48.215, 53.745, 59.275, 64.805, 70.335, 75.865, 81.395, 86.925};
  
  const float NAIVE_PHI_LOC_RING_0[12] = {0.26179939, 0.78539816, 1.30899694, 1.83259571, 2.35619449,
                                          2.87979327, 3.40339204, 3.92699082, 4.45058959, 4.97418837,
                                          5.49778714, 6.02138592};

  const float NAIVE_PHI_LOC[24] = {0.13089969, 0.39269908, 0.65449847, 0.91629786, 1.17809725,
                                   1.43989663, 1.70169602, 1.96349541, 2.2252948 , 2.48709418,
                                   2.74889357, 3.01069296, 3.27249235, 3.53429174, 3.79609112,
                                   4.05789051, 4.3196899 , 4.58148929, 4.84328867, 5.10508806,
                                   5.36688745, 5.62868684, 5.89048623, 6.15228561};

  const float NAIVE_Z_MAP[2] = {-325, 325};

  void build_map();
  

};

#endif // EPD_GEOM_H