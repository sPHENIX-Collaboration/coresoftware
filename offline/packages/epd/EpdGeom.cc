#include "EpdGeom.h"

// Initializes the EPD Geometry
EPD_GEOM::EPD_GEOM() {
  this->r_map = new std::map<int, float>;
  this->phi_map = new std::map<int, float>;
  this->build_map();
}

// Cleanup
EPD_GEOM::~EPD_GEOM() {};

// Calculates the appropriate r/phi index from a tower ID
void EPD_GEOM::id_to_r_phi(uint id, uint &side, uint &r_index, uint &phi_index) {
  side = id >> 20;
  r_index = (id - (side << 20)) >> 10; 
  phi_index = id - (side << 20) - (r_index << 10);
};

// Calculates the appropriate tower id from a given r/phi index
void EPD_GEOM::r_phi_to_id(uint side, uint r_index, uint phi_index, uint &id) {
  
};

void EPD_GEOM::id_to_side_sector_tile(uint id, uint &side, uint &sector, uint &tile) {
}

// Returns the r location of the specified tile.
// Arguments: 
// int id: The tile's ID
// Returns:
// float: the tile's location in r/phi space
float EPD_GEOM::r(uint id) {
  return this->r_map->at(id);
};

// Returns the r location of the specified tile.
// Arguments: 
// int r_index: The r index of the desired tile
// int phi_index: The phi index of the desired tile
// Returns:
// float: the tile's location in r/phi space
float EPD_GEOM::r(uint r_index, uint phi_index) {
  uint id;
  this->r_phi_to_id(r_index, phi_index, id);
  return this->r(id);
};

// Returns the phi location of the specified tile.
// Arguments: 
// int id: The tile's ID
// Returns:
// float: the tile's location in r/phi space
float EPD_GEOM::phi(uint id) {
  return this->phi_map->at(id);
};

// Returns the phi location of the specified tile.
// Arguments: 
// int r_index: The r index of the desired tile
// int phi_index: The phi index of the desired tile
// Returns:
// float: the tile's location in r/phi space
float EPD_GEOM::phi(uint r_index, uint phi_index) {
  uint id;
  this->r_phi_to_id(r_index, phi_index, id);
  return this->phi(id);
};

// Generates the maps returning the r and phi for a particular tile
void EPD_GEOM::build_map() {
  
}

