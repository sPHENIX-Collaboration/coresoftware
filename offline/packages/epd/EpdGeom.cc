#include "EpdGeom.h"

#include <map>
#include <utility>
#include <iostream>

// Initializes the EPD Geometry
EPD_GEOM::EPD_GEOM() {
  this->build_map();
}

// Cleanup
EPD_GEOM::~EPD_GEOM() {};


// Calculates the appropriate tower id from a given side/r/phi index
unsigned int EPD_GEOM::side_r_phi_to_id(unsigned int side, unsigned int r_index, unsigned int phi_index) {
  unsigned int id = 0x0;
  id = id | (side << 20);
  id = id | (r_index << 10);
  id = id | phi_index;
  return id;
};

// Calculates the appropriate tower id from a given side/sector/tile index
unsigned int EPD_GEOM::side_sector_tile_to_id(unsigned int side, unsigned int sector, unsigned int tile) {
  int rmap[31] = {0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15};
  int phimap[31] = {0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
  unsigned int phi = phimap[tile] + (sector * 2);
  return this->side_r_phi_to_id(side, rmap[tile], phi);    
}

// Calculates the appropriate side/r/phi index from a tower ID
std::tuple<unsigned int, unsigned int, unsigned int> EPD_GEOM::id_to_side_r_phi(unsigned int id) {
  unsigned int side, r_index, phi_index;
  side = (id >> 20) & 0x1;
  r_index = (id >> 10) & 0x3ff;
  phi_index = id & 0x3ff;
  return std::tuple<unsigned int, unsigned int, unsigned int>(side, r_index, phi_index);
};

// Calculates the appropriate side/sector/tile from a tower ID
std::tuple<unsigned int, unsigned int, unsigned int> EPD_GEOM::id_to_side_sector_tile(unsigned int id) {
  unsigned int side, r_index, phi_index;
  std::tie(side, r_index, phi_index) = this->id_to_side_r_phi(id);
  unsigned int sector = phi_index / 2; 
  unsigned int tile = r_index * 2;
  if (r_index && (phi_index % 2 == 0)) {
    tile --;
  }
  return std::tuple<unsigned int, unsigned int, unsigned int>(side, sector, tile);
}

// Returns the r location of the specified tile.
// Arguments: 
// int id: The tile's ID
// Returns:
// float: the tile's location in r/phi space
float EPD_GEOM::r(unsigned int id) {
  return this->r_map.at(id);
};

// Returns the r location of the specified tile.
// Arguments: 
// int r_index: The r index of the desired tile
// int phi_index: The phi index of the desired tile
// Returns:
// float: the tile's location in r/phi space
float EPD_GEOM::r_from_side_r_phi(unsigned int side, unsigned int r_index, unsigned int phi_index) {
  unsigned int id;
  id = this->side_r_phi_to_id(side, r_index, phi_index);
  return this->r(id);
};


float EPD_GEOM::r_from_side_sector_tile(unsigned int side, unsigned int sector, unsigned int tile) {
  unsigned int id = this->side_sector_tile_to_id(side, sector, tile);
  return this->r(id);
}

// Returns the phi location of the specified tile.
// Arguments: 
// int id: The tile's ID
// Returns:
// float: the tile's location in r/phi space
float EPD_GEOM::phi(unsigned int id) {
  return this->phi_map.at(id);
};

// Returns the phi location of the specified tile.
// Arguments: 
// int r_index: The r index of the desired tile
// int phi_index: The phi index of the desired tile
// Returns:
// float: the tile's location in r/phi space
float EPD_GEOM::phi_from_side_r_phi(unsigned int side, unsigned int r_index, unsigned int phi_index) {
  unsigned int id;
  id = this->side_r_phi_to_id(side, r_index, phi_index);
  return this->phi(id);
};

float EPD_GEOM::phi_from_side_sector_tile(unsigned int side, unsigned int sector, unsigned int tile) {
  unsigned int id = this->side_sector_tile_to_id(side, sector, tile);
  return this->phi(id);
}


float EPD_GEOM::z(unsigned int id) {
  return this->z_map.at(id);
}

float EPD_GEOM::z_from_side_r_phi(unsigned int side, unsigned int r_index, unsigned int phi_index) {
  unsigned int id = this->side_r_phi_to_id(side, r_index, phi_index);
  return this->z(id);
}

float EPD_GEOM::z_from_side_sector_tile(unsigned int side, unsigned int sector, unsigned int tile) {
  unsigned int id = this->side_sector_tile_to_id(side, sector, tile);
  return this->z(id);
}


// Generates the maps returning the r and phi for a particular tile
void EPD_GEOM::build_map() {
  for (unsigned int side = 0; side < 2; side++) {
    for (unsigned int r_index = 0; r_index < 16; r_index++) {
      for (unsigned int phi_index = 0; phi_index < 24; phi_index++) {
        unsigned int id = this->side_r_phi_to_id(side, r_index, phi_index);
        this->z_map[id] = this->NAIVE_Z_MAP[side];
        this->r_map[id] = this->NAIVE_R_LOC[r_index];
        if (r_index == 0) {
          this->phi_map[id] = this->NAIVE_PHI_LOC_RING_0[phi_index]; // Inner ring only has 12 tiles
        } else {
          this->phi_map[id] = this->NAIVE_PHI_LOC[phi_index];        // The rest have 24
        }
      }
    }
  }
}

bool EPD_GEOM::test_id_mapping() {
  bool pass = true;
  for (unsigned int side = 0; side < 2; side++) {
    for (unsigned int r_index = 0; r_index < 16; r_index++) {
      for (unsigned int phi_index = 0; phi_index < 24; phi_index++) {
        if (r_index == 0 && (phi_index & 0x1)) {
          continue;
        }
        std::cout << "side: " << side << "\tr: " << r_index << "\tphi: " << phi_index << std::endl;
        unsigned int id_from_r_phi = this->side_r_phi_to_id(side, r_index, phi_index);  // side r phi to id
        std::cout << "id 1: " << id_from_r_phi << std::endl;
        unsigned int temp_side, sector, tile;
        std::tie(temp_side, sector, tile) = this->id_to_side_sector_tile(id_from_r_phi);  // id to side sector tile
        std::cout << "side: " << side << "\tsector: " << sector << "\ttile: " << tile << std::endl;
        unsigned int id_from_sector_tile = this->side_sector_tile_to_id(temp_side, sector, tile); // side sector tile to id
        std::cout << "id 2: " << id_from_sector_tile << std::endl;
        unsigned int final_side, final_r_index, final_phi_index;
        std::tie(final_side, final_r_index, final_phi_index) = this->id_to_side_r_phi(id_from_sector_tile); /// id to side r phi
        std::cout << "side: " << final_side << "\tr: " << final_r_index << "\tphi: " << final_phi_index << std::endl;
        if (side != final_side || r_index != final_r_index || phi_index != final_phi_index) {
          pass = false;
          std::cout << "COORDINATE FAILED" << std::endl;
        } else{
          std::cout << "coordinate passed" << std::endl;
        }
        std::cout << "\n\n" << std::endl;
      }
    }
  }
  return pass;
}

