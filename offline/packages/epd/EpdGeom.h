/*
sEPD Geometry class

Tristan Protzman
tlprotzman@gmail.com
December 15, 2022

Converts from r/phi index to location in r/phi space

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

  uint r_phi_to_id(uint side, uint r_index, uint phi_index);
  uint side_sector_tile_to_id(uint side, uint sector, uint tile);
  std::tuple<uint, uint, uint> id_to_r_phi(uint id);
  std::tuple<uint, uint, uint> id_to_side_sector_tile(uint id);
  float r(uint id);
  float r_from_side_r_phi(uint side, uint r_index, uint phi_index);
  float r_from_side_sector_tile(uint side, uint sector, uint tile);
  float phi(uint id);
  float phi_from_side_r_phi(uint side, uint r_index, uint phi_index);
  float phi_from_side_sector_tile(uint side, uint sector, uint tile);

private:
  std::map<int, float> r_map;
  std::map<int, float> phi_map;

  // const float NAIVE_R_INDEX[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  const float NAIVE_R_LOC[16] = {6.8, 11.2, 15.6, 20.565, 26.095, 31.625, 37.155, 42.685, 48.215, 53.745, 59.275, 64.805, 70.335, 75.865, 81.395, 86.925};
  
  const float NAIVE_PHI_LOC_RING_0[12] = {0.28559091, 0.85677273, 1.42795455, 1.99913636,
                                            2.57031818, 3.1415    , 3.71268182, 4.28386364,
                                            4.85504546, 5.42622727, 5.99740909, 6.56859091};

  const float NAIVE_PHI_LOC[24] = {0.13658696, 0.40976087, 0.68293478, 0.95610869, 1.22928261,
                                     1.50245652, 1.77563043, 2.04880435, 2.32197826, 2.59515217,
                                     2.86832609, 3.1415    , 3.41467391, 3.68784782, 3.96102174,
                                     4.23419565, 4.50736956, 4.78054348, 5.05371739, 5.3268913 ,
                                     5.60006522, 5.87323913, 6.14641304, 6.41958696};

  bool test_id_mapping();
  void build_map();
  

};

#endif // EPD_GEOM_H