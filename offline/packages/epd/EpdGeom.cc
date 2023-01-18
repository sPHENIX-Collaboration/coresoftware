#include "EpdGeom.h"

#include <map>
#include <utility>

// Initializes the EPD Geometry
EPD_GEOM::EPD_GEOM() {
  this->build_map();
}

// Cleanup
EPD_GEOM::~EPD_GEOM() {};


// Calculates the appropriate tower id from a given side/r/phi index
uint EPD_GEOM::side_r_phi_to_id(uint side, uint r_index, uint phi_index) {
  uint id = 0x0;
  id = id | (side << 20);
  id = id | (r_index << 10);
  id = id | phi_index;
  return id;
};

// Calculates the appropriate tower id from a given side/sector/tile index
uint EPD_GEOM::side_sector_tile_to_id(uint side, uint sector, uint tile) {
  int rmap[31] = {0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15};
  int phimap[31] = {0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
  uint phi = phimap[tile] + (sector * 2);
  return this->side_r_phi_to_id(side, rmap[tile], phi);    
}

// Calculates the appropriate side/r/phi index from a tower ID
std::tuple<uint, uint, uint> EPD_GEOM::id_to_side_r_phi(uint id) {
  uint side, r_index, phi_index;
  side = (id >> 20) & 0x1;
  r_index = (id >> 10) & 0x3ff;
  phi_index = id & 0x3ff;
  return std::tuple<uint, uint, uint>(side, r_index, phi_index);
};

// Calculates the appropriate side/sector/tile from a tower ID
std::tuple<uint, uint, uint> EPD_GEOM::id_to_side_sector_tile(uint id) {
  uint side, r_index, phi_index;
  std::tie(side, r_index, phi_index) = this->id_to_side_r_phi(id);
  uint sector = phi_index / 2; 
  uint tile = r_index * 2 - (phi_index % 2);
  return std::tuple<uint, uint, uint>(side, sector, tile);
}

// Returns the r location of the specified tile.
// Arguments: 
// int id: The tile's ID
// Returns:
// float: the tile's location in r/phi space
float EPD_GEOM::r(uint id) {
  return this->r_map.at(id);
};

// Returns the r location of the specified tile.
// Arguments: 
// int r_index: The r index of the desired tile
// int phi_index: The phi index of the desired tile
// Returns:
// float: the tile's location in r/phi space
float EPD_GEOM::r_from_side_r_phi(uint side, uint r_index, uint phi_index) {
  uint id;
  id = this->side_r_phi_to_id(side, r_index, phi_index);
  return this->r(id);
};


float EPD_GEOM::r_from_side_sector_tile(uint side, uint sector, uint tile) {
  uint id = this->side_sector_tile_to_id(side, sector, tile);
  return this->r(id);
}

// Returns the phi location of the specified tile.
// Arguments: 
// int id: The tile's ID
// Returns:
// float: the tile's location in r/phi space
float EPD_GEOM::phi(uint id) {
  return this->phi_map.at(id);
};

// Returns the phi location of the specified tile.
// Arguments: 
// int r_index: The r index of the desired tile
// int phi_index: The phi index of the desired tile
// Returns:
// float: the tile's location in r/phi space
float EPD_GEOM::phi_from_side_r_phi(uint side, uint r_index, uint phi_index) {
  uint id;
  id = this->side_r_phi_to_id(side, r_index, phi_index);
  return this->phi(id);
};

float EPD_GEOM::phi_from_side_sector_tile(uint side, uint sector, uint tile) {
  uint id = this->side_sector_tile_to_id(side, sector, tile);
  return this->phi(id);
}

// Generates the maps returning the r and phi for a particular tile
void EPD_GEOM::build_map() {
  for (uint id = 0; id < this->NUM_TOWERS; id++) {
    uint side, r_index, phi_index;
    std::tie(side, r_index, phi_index) = this->id_to_side_r_phi(id);
    this->r_map.at(id) = this->NAIVE_R_LOC[r_index];
    if (r_index == 0) {
      this->phi_map.at(id) = this->NAIVE_PHI_LOC_RING_0[phi_index]; // Inner ring only has 12 tiles
    } else {
      this->phi_map.at(id) = this->NAIVE_PHI_LOC[phi_index];        // The rest have 24
    }
  }
}

