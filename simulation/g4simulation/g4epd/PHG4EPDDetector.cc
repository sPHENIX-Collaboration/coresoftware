/* vim: set sw=2 ft=cpp: */

#include "PHG4EPDDetector.h"

#include "PHG4EPDDisplayAction.h"

#include <epd/EPDDefs.h>

#include <g4main/PHG4Detector.h>
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <phparameter/PHParameters.h>

#include <Geant4/G4ExtrudedSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4TwoVector.hh>        // for G4TwoVector
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume

#include <algorithm>  // for max
#include <cmath>
#include <vector>  // for vector

PHG4EPDDetector::PHG4EPDDetector(PHG4Subsystem* subsys,
                                 PHCompositeNode* node,
                                 PHParameters* parameters,
                                 std::string const& name)
  : PHG4Detector(subsys, node, name)
  , m_DisplayAction(dynamic_cast<PHG4EPDDisplayAction*>(subsys->GetDisplayAction()))
  , m_Params(parameters)
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_SupportActiveFlag(m_Params->get_int_param("supportactive"))
{
}

void PHG4EPDDetector::ConstructMe(G4LogicalVolume* world)
{
  G4Material* material = GetDetectorMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  G4ThreeVector positive(0., 0., m_Params->get_double_param("place_z") * cm);
  G4ThreeVector negative(0., 0., -m_Params->get_double_param("place_z") * cm);

  constexpr int32_t ntiles = 31;
  constexpr int32_t nsectors = 12;

  for (int32_t i = 0; i < ntiles; ++i)
  {
    std::string label = "EPD_tile_" + std::to_string(i);

    G4ExtrudedSolid* block = construct_block(i);
    G4LogicalVolume* volume = new G4LogicalVolume(block, material, label, nullptr, nullptr, nullptr);

    GetDisplayAction()->AddVolume(volume, volume->GetName());
    m_ActiveLogVolSet.insert(volume);

    for (int32_t k = 0; k < nsectors; ++k)
    {

      G4RotationMatrix* rotate = new G4RotationMatrix();
     
      double phi_shift = (k + 9) * 2 * M_PI / nsectors;
     
      if(phi_shift >= (2.0 * M_PI))
      {
         phi_shift -= (2.0 * M_PI);
      }     
      else if (phi_shift < 0.0)
      {
         phi_shift += (2.0 * M_PI);
      }

      rotate->rotateZ(-1*phi_shift);

      m_volumes.emplace(
          new G4PVPlacement( rotate, negative, volume, label, world, false, 2 * k + 0, OverlapCheck()),
          module_id_for(i, k, 0));

      m_volumes.emplace(
          new G4PVPlacement( rotate, positive, volume, label, world, false, 2 * k + 1, OverlapCheck()),
          module_id_for(i, k, 1));
    }
  }
}

int PHG4EPDDetector::IsInDetector(G4VPhysicalVolume* volume) const
{
  G4LogicalVolume* mylogvol = volume->GetLogicalVolume();
  if (m_ActiveFlag)
  {
    if (m_ActiveLogVolSet.find(mylogvol) != m_ActiveLogVolSet.end())
    {
      return 1;
    }
  }
  if (m_SupportActiveFlag)
  {
    if (m_SupportLogVolSet.find(mylogvol) != m_SupportLogVolSet.end())
    {
      return -2;
    }
  }
  return 0;
}

uint32_t PHG4EPDDetector::module_id_for(uint32_t tile_id, uint32_t sector, uint32_t arm)
{
  return EPDDefs::make_epd_key(arm,sector,tile_id);
}

uint32_t PHG4EPDDetector::module_id_for(G4VPhysicalVolume* volume)
{
  return m_volumes[volume];
}

static constexpr double dz = 6.;


static constexpr double coordinates[31][5][2] =
{
  {{ -22.43,  40.39},   { -44.41,  78.05},   { -23.29,  86.91},   {  -0.57,  89.80},   {-0.77,46.19}},
  {{ -44.96,  79.66},   { -66.53, 116.62},   { -35.41, 129.51},   { -24.54,  88.12},   { 0.00, 0.00}},
  {{ -22.81,  88.58},   { -34.08, 129.86},   {  -0.69, 134.26},   {  -0.90,  91.47},   { 0.00, 0.00}},
  {{ -67.15, 118.09},   { -88.72, 155.05},   { -46.90, 172.37},   { -36.03, 130.98},   { 0.00, 0.00}},
  {{ -34.29, 131.45},   { -45.57, 172.73},   {  -0.69, 178.64},   {  -0.90, 135.85},   { 0.00, 0.00}},
  {{ -89.34, 156.53},   {-116.61, 203.35},   { -61.34, 226.25},   { -47.51, 173.85},   { 0.00, 0.00}},
  {{ -45.78, 174.32},   { -60.01, 226.61},   {  -0.69, 234.42},   {  -0.90, 180.22},   { 0.00, 0.00}},
  {{-117.22, 204.83},   {-144.50, 251.66},   { -75.77, 280.13},   { -61.95, 227.73},   { 0.00, 0.00}},
  {{ -60.21, 228.19},   { -74.44, 280.48},   {  -0.69, 290.19},   {  -0.90, 236.00},   { 0.00, 0.00}},
  {{-145.11, 253.14},   {-172.39, 299.96},   { -90.21, 334.00},   { -76.38, 281.60},   { 0.00, 0.00}},
  {{ -74.65, 282.07},   { -88.88, 334.36},   {  -0.69, 345.97},   {  -0.90, 291.78},   { 0.00, 0.00}},
  {{-173.00, 301.44},   {-200.28, 348.27},   {-104.65, 387.88},   { -90.82, 335.48},   { 0.00, 0.00}},
  {{ -89.09, 335.95},   {-103.31, 388.24},   {  -0.69, 401.75},   {  -0.90, 347.56},   { 0.00, 0.00}},
  {{-200.89, 349.75},   {-228.17, 396.57},   {-119.08, 441.76},   {-105.26, 389.36},   { 0.00, 0.00}},
  {{-103.52, 389.82},   {-117.75, 442.11},   {  -0.69, 457.52},   {  -0.90, 403.33},   { 0.00, 0.00}},
  {{-228.78, 398.05},   {-256.05, 444.88},   {-133.52, 495.63},   {-119.69, 443.23},   { 0.00, 0.00}},
  {{-117.96, 443.70},   {-132.19, 495.99},   {  -0.69, 513.30},   {  -0.90, 459.11},   { 0.00, 0.00}},
  {{-256.67, 446.35},   {-283.94, 493.18},   {-147.95, 549.51},   {-134.13, 497.11},   { 0.00, 0.00}},
  {{-132.40, 497.58},   {-146.62, 549.87},   {  -0.69, 569.08},   {  -0.90, 514.89},   { 0.00, 0.00}},
  {{-284.56, 494.66},   {-311.83, 541.49},   {-162.39, 603.39},   {-148.57, 550.99},   { 0.00, 0.00}},
  {{-146.83, 551.45},   {-161.06, 603.74},   {  -0.69, 624.86},   {  -0.90, 570.66},   { 0.00, 0.00}},
  {{-312.44, 542.96},   {-339.72, 589.79},   {-176.83, 657.26},   {-163.00, 604.86},   { 0.00, 0.00}},
  {{-161.27, 605.33},   {-175.50, 657.62},   {  -0.69, 680.63},   {  -0.90, 626.44},   { 0.00, 0.00}},
  {{-340.33, 591.27},   {-367.61, 638.09},   {-191.26, 711.14},   {-177.44, 658.74},   { 0.00, 0.00}},
  {{-175.70, 659.21},   {-189.93, 711.50},   {  -0.69, 736.41},   {  -0.90, 682.22},   { 0.00, 0.00}},
  {{-368.22, 639.57},   {-395.50, 686.40},   {-205.70, 765.02},   {-191.87, 712.62},   { 0.00, 0.00}},
  {{-190.14, 713.08},   {-204.37, 765.37},   {  -0.69, 792.19},   {  -0.90, 738.00},   { 0.00, 0.00}},
  {{-396.11, 687.88},   {-423.39, 734.70},   {-220.13, 818.89},   {-206.31, 766.49},   { 0.00, 0.00}},
  {{-204.58, 766.96},   {-218.80, 819.25},   {  -0.69, 847.96},   {  -0.90, 793.77},   { 0.00, 0.00}},
  {{-424.00, 736.18},   {-451.58, 783.75},   {-234.88, 873.51},   {-220.75, 820.37},   { 0.00, 0.00}},
  {{-219.01, 820.84},   {-233.34, 873.92},   {  -0.79, 904.53},   {  -0.90, 849.55},   { 0.00, 0.00}},
};


G4ExtrudedSolid* PHG4EPDDetector::construct_block(int32_t index)
{
  std::string label("tile_" + std::to_string(index));

  const double(*coords)[5][2] = &coordinates[index];

  std::vector<G4TwoVector> vertices;

  for (int32_t i = 0; i < 5; ++i)
  {
    double x = (*coords)[i][0];
    double y = (*coords)[i][1];

    if (x == 0. && y == 0.)
    {
      continue;
    }
    vertices.emplace_back(x * mm, y * mm);
  }

  std::vector<G4ExtrudedSolid::ZSection> zsections = {
      G4ExtrudedSolid::ZSection(-dz, G4TwoVector(), 1.),
      G4ExtrudedSolid::ZSection(dz, G4TwoVector(), 1.),
  };

  return new G4ExtrudedSolid(label, vertices, zsections);
}
