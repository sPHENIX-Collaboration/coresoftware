#ifndef PHG4MICROMEGASSURVEY_H
#define PHG4MICROMEGASSURVEY_H

/*!
 * \file PHG4MicromegasSurvey.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 * \brief implements survey data for TPOT definition
 */

/**
 * survey data are provided in the form of a G4Transform3D object for each TPOT detector
 * it is to be applied on top of the default GEANT tranformation as defined in PHG4MicromegasDetector
 * to move the strips from their "ideal" to their surveyed position
 */

#include <Geant4/G4Transform3D.hh>

#include <array>
#include <string>
#include <unordered_map>

class PHG4MicromegasSurvey
{
 public:
  /// constructor
  PHG4MicromegasSurvey();

  /// get module name from tile and layer
  std::string get_module_name(int layer, uint tile) const;

  /// get transformation from tile and layer
  G4Transform3D get_transformation(int layer, uint tile) const;

 private:
  /// internal detector definition (tile number, detector)
  struct tile_id_t
  {
    tile_id_t(int layer, uint tile)
      : m_layer(layer)
      , m_tile(tile)
    {
    }
    int m_layer = 0;
    uint m_tile = 0;

    bool operator==(const tile_id_t& other) const
    {
      return other.m_layer == m_layer && other.m_tile == m_tile;
    }
  };

  struct tile_id_hash_t
  {
    std::size_t operator()(const tile_id_t& id) const noexcept
    {
      return id.m_tile + (id.m_layer << 4);
    }
  };

  /// map tile_id to module name
  using tile_map_t = std::unordered_map<tile_id_t, std::string, tile_id_hash_t>;
  tile_map_t m_tile_map;

  /// internal G4Transformation format
  using rotation_t = std::array<double, 3>;
  using translation_t = std::array<double, 3>;
  struct transformation_t
  {
    transformation_t(const rotation_t& rotation, const translation_t& translation)
      : m_rotation(rotation)
      , m_translation(translation)
    {
    }

    /// x, y and z axis rotation in order, degrees
    rotation_t m_rotation = {{0}};

    /// x, y, z translation, cm
    translation_t m_translation = {{0}};
  };

  /// map module name to transformation
  using transformation_map_t = std::unordered_map<std::string, transformation_t>;
  transformation_map_t m_transformation_map;
};

#endif
