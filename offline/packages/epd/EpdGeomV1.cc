#include "EpdGeomV1.h"

#include <map>
#include <utility>
#include <tuple>
#include <iostream>
#include <fstream>
#include <sstream>

// Initializes the EPD Geometry
EpdGeomV1::EpdGeomV1() {
  build_lookup();
}

// Cleanup
EpdGeomV1::~EpdGeomV1() {};


// Calculates the appropriate tower id from a given side/r/phi index
unsigned int EpdGeomV1::side_r_phi_to_id(unsigned int side, unsigned int r_index, unsigned int phi_index) {
  unsigned int id = 0x0;
  id = id | (side << 20);
  id = id | (r_index << 10);
  id = id | phi_index;
  return id;
};

// Calculates the appropriate tower id from a given side/sector/tile index
unsigned int EpdGeomV1::side_sector_tile_to_id(unsigned int side, unsigned int sector, unsigned int tile) {
  int rmap[31] = {0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15};
  int phimap[31] = {0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
  unsigned int phi = phimap[tile] + (sector * 2);
  return side_r_phi_to_id(side, rmap[tile], phi);    
}

// Calculates the appropriate side/r/phi index from a tower ID
std::tuple<unsigned int, unsigned int, unsigned int> EpdGeomV1::id_to_side_r_phi(unsigned int id) {
  unsigned int side, r_index, phi_index;
  side = (id >> 20) & 0x1;
  r_index = (id >> 10) & 0x3ff;
  phi_index = id & 0x3ff;
  return std::tuple<unsigned int, unsigned int, unsigned int>(side, r_index, phi_index);
};

// Calculates the appropriate side/sector/tile from a tower ID
std::tuple<unsigned int, unsigned int, unsigned int> EpdGeomV1::id_to_side_sector_tile(unsigned int id) {
  unsigned int side, r_index, phi_index;
  std::tie(side, r_index, phi_index) = id_to_side_r_phi(id);
  unsigned int sector = phi_index / 2; 
  unsigned int tile = r_index * 2;
  if (r_index && (phi_index % 2 == 0)) {
    tile --;
  }
  return std::tuple<unsigned int, unsigned int, unsigned int>(side, sector, tile);
}

float EpdGeomV1::r(unsigned int id) {
  return r_lookup[id];
}

// Returns the r location of the specified tile.
// Arguments: 
// int id: The tile's ID
// Returns:
// float: the tile's location in r/phi space
float EpdGeomV1::r_from_key(unsigned int key) {
  return r(decode_epd(key));
};

// Returns the r location of the specified tile.
// Arguments: 
// int r_index: The r index of the desired tile
// int phi_index: The phi index of the desired tile
// Returns:
// float: the tile's location in r/phi space
float EpdGeomV1::r_from_side_r_phi(unsigned int side, unsigned int r_index, unsigned int phi_index) {
  unsigned int key;
  key = side_r_phi_to_id(side, r_index, phi_index);
  return r_from_key(key);
};


float EpdGeomV1::r_from_side_sector_tile(unsigned int side, unsigned int sector, unsigned int tile) {
  unsigned int key = side_sector_tile_to_id(side, sector, tile);
  return r_from_key(key);
}

float EpdGeomV1::phi(unsigned int id) {
  return phi_lookup[id];
}

// Returns the phi location of the specified tile.
// Arguments: 
// int id: The tile's ID
// Returns:
// float: the tile's location in r/phi space
float EpdGeomV1::phi_from_key(unsigned int key) {
  return phi(decode_epd(key));
};

// Returns the phi location of the specified tile.
// Arguments: 
// int r_index: The r index of the desired tile
// int phi_index: The phi index of the desired tile
// Returns:
// float: the tile's location in r/phi space
float EpdGeomV1::phi_from_side_r_phi(unsigned int side, unsigned int r_index, unsigned int phi_index) {
  unsigned int key;
  key = side_r_phi_to_id(side, r_index, phi_index);
  return phi_from_key(key);
};

float EpdGeomV1::phi_from_side_sector_tile(unsigned int side, unsigned int sector, unsigned int tile) {
  unsigned int key = side_sector_tile_to_id(side, sector, tile);
  return phi_from_key(key);
}

float EpdGeomV1::z(unsigned int id) {
  return z_lookup[id];
}

float EpdGeomV1::z_from_key(unsigned int key) {
  return z(decode_epd(key));
}

float EpdGeomV1::z_from_side_r_phi(unsigned int side, unsigned int r_index, unsigned int phi_index) {
  unsigned int key = side_r_phi_to_id(side, r_index, phi_index);
  return z_from_key(key);
}

float EpdGeomV1::z_from_side_sector_tile(unsigned int side, unsigned int sector, unsigned int tile) {
  unsigned int key = side_sector_tile_to_id(side, sector, tile);
  return z_from_key(key);
}

unsigned int EpdGeomV1::decode_epd(unsigned int tower_key) {
  int channels_per_sector = 31;
  int supersector = channels_per_sector * 12;
  unsigned int ns_sector = tower_key >> 20U;
  unsigned int rbin = (tower_key - (ns_sector << 20U)) >> 10U;
  unsigned int phibin = tower_key - (ns_sector << 20U) - (rbin << 10U);
  int epdchnlmap[16][2] = {{0, 0}, {1, 2}, {3, 4}, {5, 6}, {7, 8}, {9, 10}, {11, 12}, {13, 14}, {15, 16}, {17, 18}, {19, 20}, {21, 22}, {23, 24}, {25, 26}, {27, 28}, {29, 30}};
  int sector = phibin / 2;
  int channel = 0;
  if (rbin > 0)
    {
      channel = epdchnlmap[rbin][phibin - 2 * sector];
    }
  else
    {
      channel = 0;
    }
  unsigned int index = 0;
  index = ns_sector * supersector + sector * channels_per_sector + channel;
  return index;
}


// Generates the maps returning the r and phi for a particular tile
void EpdGeomV1::build_lookup() {
  r_lookup = std::vector<float>(NUM_TOWERS, -999);
  phi_lookup = std::vector<float>(NUM_TOWERS, -999);
  z_lookup = std::vector<float>(NUM_TOWERS, -999);
  std::string epd_geom_file_path("epd_geometry.tsv");
  std::ifstream epd_geom_file(epd_geom_file_path);
  if (!epd_geom_file.is_open()) {
    std::cerr << "Could not read EPD Geometry file: " << epd_geom_file_path << std::endl;
    return; // Should this throw an error?  Not sure what the sPHENIX way to do things like that is
  }
  std::string line;
  uint filled = 0;
  while(std::getline(epd_geom_file, line)) {
    std::stringstream tokens(line);
    int tower;
    tokens >> tower;  // Read the first column into the tower name;
    float r, phi, z;
    tokens >> r;
    tokens >> phi;
    tokens >> z;  // This could all use error checking
    r_lookup[tower] = r;
    phi_lookup[tower] = phi;
    z_lookup[tower] = z;
    filled++;
  }
  if (filled != NUM_TOWERS) {
    std::cerr << "Did not find geometry values for all EPD towers!" << std::endl;
  }
  // for (uint i = 0; i < NUM_TOWERS; i++) {
  //   printf("%d\t%f\t%f\t%f\n", i, r_lookup[i], phi_lookup[i], z_lookup[i]);
  // }
}

bool EpdGeomV1::test_id_mapping() {
  bool pass = true;
  for (unsigned int side = 0; side < 2; side++) {
    for (unsigned int r_index = 0; r_index < 16; r_index++) {
      for (unsigned int phi_index = 0; phi_index < 24; phi_index++) {
        if (r_index == 0 && (phi_index & 0x1)) {
          continue;
        }
        std::cout << "side: " << side << "\tr: " << r_index << "\tphi: " << phi_index << std::endl;
        unsigned int id_from_r_phi = side_r_phi_to_id(side, r_index, phi_index);  // side r phi to id
        std::cout << "id 1: " << id_from_r_phi << std::endl;
        unsigned int temp_side, sector, tile;
        std::tie(temp_side, sector, tile) = id_to_side_sector_tile(id_from_r_phi);  // id to side sector tile
        std::cout << "side: " << side << "\tsector: " << sector << "\ttile: " << tile << std::endl;
        unsigned int id_from_sector_tile = side_sector_tile_to_id(temp_side, sector, tile); // side sector tile to id
        std::cout << "id 2: " << id_from_sector_tile << std::endl;
        unsigned int final_side, final_r_index, final_phi_index;
        std::tie(final_side, final_r_index, final_phi_index) = id_to_side_r_phi(id_from_sector_tile); /// id to side r phi
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
