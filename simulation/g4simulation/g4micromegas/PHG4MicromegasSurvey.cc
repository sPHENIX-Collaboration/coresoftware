/*!
 * \file PHG4MicromegasSurvey.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 * \brief implements survey data for TPOT definition
 */

#include "PHG4MicromegasSurvey.h"

#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>

//____________________________________________________________________________________________________
PHG4MicromegasSurvey::PHG4MicromegasSurvey()
  :  // map layer and tile number to detector name
  m_tile_map{
      {{55, 0}, "M5P"}, {{56, 0}, "M5Z"}, {{55, 1}, "M8P"}, {{56, 1}, "M8Z"}, {{55, 2}, "M4P"}, {{56, 2}, "M4Z"}, {{55, 3}, "M10P"}, {{56, 3}, "M10Z"}, {{55, 4}, "M9P"}, {{56, 4}, "M9Z"}, {{55, 5}, "M2P"}, {{56, 5}, "M2Z"}, {{55, 6}, "M6P"}, {{56, 6}, "M6Z"}, {{55, 7}, "M7P"}, {{56, 7}, "M7Z"}}
  ,

  // map module name to transformation
  /*
   * these are constructing by mapping the first and last strip positions as given by GEANT4 default implementation of TPOT
   * to the ones measured by the survey crew
   */
  m_transformation_map{
      {"M10P",
       {{0.0238291, -0.0610202, 359.729},
        {0.32958, 0.104316, 0.491085}}},
      {"M10Z",
       {{0.0240982, -0.01546, -0.267002},
        {0.196309, -0.00618558, 0.126347}}},
      {"M4P",
       {{0.10541, 0.135221, -0.691686},
        {0.775549, 0.281716, 0.603817}}},
      {"M4Z",
       {{0.10526, 0.122168, -0.691598},
        {0.691176, 0.259581, 0.233525}}},
      {"M8P",
       {{-0.417662, 0.0569648, -0.0964277},
        {-0.0666588, 0.366775, -0.147532}}},
      {"M8Z",
       {{-0.417529, 0.120316, -0.0951897},
        {-0.086922, 0.333419, -0.512579}}},
      {"M5P",
       {{-0.202661, 0.0403903, -0.235734},
        {0.153994, 0.313934, 0.13618}}},
      {"M5Z",
       {{-0.203231, 0.0481545, -0.228361},
        {0.0652422, 0.30468, -0.241307}}},
      {"M2P",
       {{0.0243278, 0.0186967, 0.524609},
        {-1.1679, 0.746204, 0.704047}}},
      {"M2Z",
       {{-0.00159013, -0.0326792, 0.522416},
        {-1.21045, 0.7655, 0.315458}}},
      {"M9P",
       {{-0.172856, 0.156516, 0.630981},
        {-1.23611, 0.828348, 0.349842}}},
      {"M9Z",
       {{-0.174169, 0.153218, 0.632272},
        {-1.32091, 0.861492, -0.0240633}}},
      {"M7P",
       {{-0.0802472, 0.12362, 0.804231},
        {-1.15284, -0.140113, 0.102525}}},
      {"M7Z",
       {{-0.0822581, 0.127216, 0.803837},
        {-1.20763, -0.242328, -0.269243}}},
      {"M6P",
       {{-0.221227, -0.0357901, 0.597703},
        {-0.951645, 0.183264, -0.169542}}},
      {"M6Z",
       {{-0.169038, -0.126228, 0.59627},
        {-0.99097, -0.0452179, -0.574241}}}}
{
}

//____________________________________________________________________________________________________
std::string PHG4MicromegasSurvey::get_module_name(int layer, uint tile) const
{
  const auto iter = m_tile_map.find({layer, tile});
  if (iter == m_tile_map.end())
  {
    std::cout << " PHG4MicromegasSurvey::get_module_name - module not found. layer: " << layer << " tile: " << tile << std::endl;
    return std::string();
  }
  else
  {
    return iter->second;
  }
}

//____________________________________________________________________________________________________
G4Transform3D PHG4MicromegasSurvey::get_transformation(int layer, uint tile) const
{
  const auto tile_iter = m_tile_map.find({layer, tile});
  if (tile_iter == m_tile_map.end())
  {
    std::cout << " PHG4MicromegasSurvey::get_transformation - module not found. layer: " << layer << " tile: " << tile << std::endl;
    return G4Transform3D();
  }

  const auto& tile_name = tile_iter->second;
  const auto transform_iter = m_transformation_map.find(tile_name);
  if (transform_iter == m_transformation_map.end())
  {
    std::cout
        << " PHG4MicromegasSurvey::get_transformation - transformation not found."
        << " layer: " << layer
        << " tile: " << tile
        << " tile_name: " << tile_name
        << std::endl;
    return G4Transform3D();
  }

  // get transformation
  const auto& transformation = transform_iter->second;

  // translation
  const G4ThreeVector translation(
      transformation.m_translation[0] * cm,
      transformation.m_translation[1] * cm,
      transformation.m_translation[2] * cm);

  // rotation
  G4RotationMatrix rotation;
  rotation.rotateX(transformation.m_rotation[0] * deg);
  rotation.rotateY(transformation.m_rotation[1] * deg);
  rotation.rotateZ(transformation.m_rotation[2] * deg);

  return G4Transform3D(rotation, translation);
}
