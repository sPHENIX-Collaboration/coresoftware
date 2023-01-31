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
#include <iostream>

#include <phool/PHObject.h>


// sEPD geometry class
class EPD_GEOM {

public:
  const unsigned int NUM_TOWERS = 768;
  const unsigned int MAX_R = 32;
  const unsigned int MAX_PHI = 24;

  EPD_GEOM();
  ~EPD_GEOM();

  unsigned int side_r_phi_to_id(unsigned int side, unsigned int r_index, unsigned int phi_index);
  unsigned int side_sector_tile_to_id(unsigned int side, unsigned int sector, unsigned int tile);
  std::tuple<unsigned int, unsigned int, unsigned int> id_to_side_r_phi(unsigned int id);
  std::tuple<unsigned int, unsigned int, unsigned int> id_to_side_sector_tile(unsigned int id);
  float r(unsigned int id);
  float r_from_side_r_phi(unsigned int side, unsigned int r_index, unsigned int phi_index);
  float r_from_side_sector_tile(unsigned int side, unsigned int sector, unsigned int tile);
  float phi(unsigned int id);
  float phi_from_side_r_phi(unsigned int side, unsigned int r_index, unsigned int phi_index);
  float phi_from_side_sector_tile(unsigned int side, unsigned int sector, unsigned int tile);
  float z(unsigned int id);
  float z_from_side_r_phi(unsigned int side, unsigned int r_inxed, unsigned int phi_index);
  float z_from_side_sector_tile(unsigned int side, unsigned int sector, unsigned int tile);
  
  
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

  const float NAIVE_Z_MAP[2] = {-316, 316};

  void build_map();
  

};

#endif // EPD_GEOM_H