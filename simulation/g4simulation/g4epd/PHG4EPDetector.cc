/* vim: set sw=2 ft=cpp: */

#include "PHG4EPDetector.h"

#include <g4main/PHG4Detector.h>

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <Geant4/G4Colour.hh>
#include <Geant4/G4ExtrudedSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4String.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4TwoVector.hh>  // for G4TwoVector
#include <Geant4/G4VisAttributes.hh>

#include <cmath>
#include <iterator>  // for end
#include <vector>    // for vector

PHG4EPDetector::PHG4EPDetector(PHG4Subsystem* subsys,
                               PHCompositeNode* node,
                               PHParametersContainer* params,
                               std::string const& name)
  : PHG4Detector(subsys, node, name)
{
  PHParameters const* pars = params->GetParameters(-1);
  m_z_position = pars->get_double_param("z_position");
}

void PHG4EPDetector::ConstructMe(G4LogicalVolume* world)
{
  G4Material* material = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  G4VisAttributes* attrs = new G4VisAttributes();

  attrs->SetVisibility(true);
  attrs->SetForceSolid(true);
  attrs->SetColour(G4Colour::Red());

  G4ThreeVector positive(0., 0., m_z_position * cm);
  G4ThreeVector negative(0., 0., -m_z_position * cm);

  constexpr int32_t ntiles = 31;
  constexpr int32_t nslices = 12;

  for (int32_t i = 0; i < ntiles; ++i)
  {
    std::string label = "EPD_tile_" + std::to_string(i);

    G4ExtrudedSolid* block = construct_block(i);
    G4LogicalVolume* volume = new G4LogicalVolume(
        block, material, label.data(), 0, 0, 0);

    volume->SetVisAttributes(attrs);

    for (int32_t k = 0; k < nslices; ++k)
    {
      G4RotationMatrix* rotate = new G4RotationMatrix();

      rotate->rotateZ(k * 2 * M_PI / nslices);

      m_volumes.emplace(
          new G4PVPlacement(
              rotate, positive, volume, label, world,
              false, 2 * k + 0, OverlapCheck()),
          module_id_for(i, k, 0));

      m_volumes.emplace(
          new G4PVPlacement(
              rotate, negative, volume, label, world,
              false, 2 * k + 1, OverlapCheck()),
          module_id_for(i, k, 1));
    }
  }
}

bool PHG4EPDetector::contains(G4VPhysicalVolume* volume) const
{
  return m_volumes.find(volume) != std::end(m_volumes);
}

uint32_t PHG4EPDetector::module_id_for(int32_t index, int32_t slice,
                                       int32_t side)
{
  return (side << 9 | slice << 5 | index) & 0x3FF;
}

uint32_t PHG4EPDetector::module_id_for(G4VPhysicalVolume* volume)
{
  return m_volumes[volume];
}

static constexpr double dz = 6.;

static constexpr double coordinates[31][6][2] = {
    {{-20.5, 82.4}, {0.0, 85.1}, {20.5, 82.4}, {9.7, 41.9}, {-9.7, 41.9}, {0.0, 0.0}},
    {{31.9, 124.9}, {1.9, 128.8}, {0.8, 129.0}, {0.8, 86.6}, {0.0, 0.0}, {20.9, 84.0}},
    {{-31.9, 124.9}, {-1.9, 128.8}, {-0.8, 129.0}, {-0.8, 86.6}, {0.0, 0.0}, {-20.9, 84.0}},
    {{43.3, 167.4}, {3.0, 172.7}, {0.8, 173.0}, {0.8, 130.6}, {1.9, 130.5}, {32.3, 126.5}},
    {{-43.3, 167.4}, {-3.0, 172.7}, {-0.8, 173.0}, {-0.8, 130.6}, {-1.9, 130.5}, {-32.3, 126.5}},
    {{57.6, 220.8}, {4.4, 227.8}, {0.8, 228.3}, {0.8, 174.6}, {3.1, 174.3}, {43.7, 169.0}},
    {{-57.6, 220.8}, {-4.4, 227.8}, {-0.8, 228.3}, {-0.8, 174.6}, {-3.1, 174.3}, {-43.7, 169.0}},
    {{71.9, 274.2}, {5.9, 282.9}, {0.8, 283.6}, {0.8, 229.9}, {4.5, 229.4}, {58.0, 222.4}},
    {{-71.9, 274.2}, {-5.9, 282.9}, {-0.8, 283.6}, {-0.8, 229.9}, {-4.5, 229.4}, {-58.0, 222.4}},
    {{86.2, 327.6}, {7.3, 338.0}, {0.8, 338.9}, {0.8, 285.2}, {5.9, 284.6}, {72.3, 275.8}},
    {{-86.2, 327.6}, {-7.3, 338.0}, {-0.8, 338.9}, {-0.8, 285.2}, {-5.9, 284.6}, {-72.3, 275.8}},
    {{100.5, 381.1}, {8.3, 394.2}, {0.8, 394.2}, {0.8, 340.5}, {7.3, 339.7}, {86.7, 329.2}},
    {{-100.5, 381.1}, {-8.3, 394.2}, {-0.8, 394.2}, {-0.8, 340.5}, {-7.3, 339.7}, {-86.7, 329.2}},
    {{114.9, 434.5}, {8.3, 448.5}, {0.8, 449.5}, {0.8, 395.9}, {8.3, 394.9}, {101.0, 382.7}},
    {{-114.9, 434.5}, {-8.3, 448.5}, {-0.8, 449.5}, {-0.8, 395.9}, {-8.3, 394.9}, {-101.0, 382.7}},
    {{129.2, 487.9}, {8.3, 503.8}, {0.8, 504.8}, {0.8, 451.2}, {8.3, 450.2}, {115.3, 436.1}},
    {{-129.2, 487.9}, {-8.3, 503.8}, {-0.8, 504.8}, {-0.8, 451.2}, {-8.3, 450.2}, {-115.3, 436.1}},
    {{143.5, 541.3}, {8.3, 559.1}, {0.8, 560.1}, {0.8, 506.5}, {8.3, 505.5}, {129.6, 489.5}},
    {{-143.5, 541.3}, {-8.3, 559.1}, {-0.8, 560.1}, {-0.8, 506.5}, {-8.3, 505.5}, {-129.6, 489.5}},
    {{157.8, 594.8}, {8.3, 614.4}, {0.8, 615.4}, {0.8, 561.8}, {8.3, 560.8}, {143.9, 542.9}},
    {{-157.8, 594.8}, {-8.3, 614.4}, {-0.8, 615.4}, {-0.8, 561.8}, {-8.3, 560.8}, {-143.9, 542.9}},
    {{172.1, 648.2}, {8.3, 669.7}, {0.8, 670.7}, {0.8, 617.1}, {8.3, 616.1}, {158.2, 596.4}},
    {{-172.1, 648.2}, {-8.3, 669.7}, {-0.8, 670.7}, {-0.8, 617.1}, {-8.3, 616.1}, {-158.2, 596.4}},
    {{186.4, 701.6}, {8.3, 725.0}, {0.8, 726.0}, {0.8, 672.4}, {8.3, 671.4}, {172.5, 649.8}},
    {{-186.4, 701.6}, {-8.3, 725.0}, {-0.8, 726.0}, {-0.8, 672.4}, {-8.3, 671.4}, {-172.5, 649.8}},
    {{200.7, 755.0}, {8.3, 780.4}, {0.8, 781.3}, {0.8, 727.7}, {8.3, 726.7}, {186.8, 703.2}},
    {{-200.7, 755.0}, {-8.3, 780.4}, {-0.8, 781.3}, {-0.8, 727.7}, {-8.3, 726.7}, {-186.8, 703.2}},
    {{215.0, 808.4}, {8.3, 835.7}, {0.8, 836.6}, {0.8, 783.0}, {8.3, 782.0}, {201.2, 756.6}},
    {{-215.0, 808.4}, {-8.3, 835.7}, {-0.8, 836.6}, {-0.8, 783.0}, {-8.3, 782.0}, {-201.2, 756.6}},
    {{229.4, 861.8}, {17.5, 890.6}, {0.8, 892.7}, {0.8, 838.3}, {8.3, 837.3}, {215.5, 810.0}},
    {{-229.4, 861.8}, {-17.5, 890.6}, {-0.8, 892.7}, {-0.8, 838.3}, {-8.3, 837.3}, {-215.5, 810.0}},
};

G4ExtrudedSolid* PHG4EPDetector::construct_block(int32_t index)
{
  G4String label("tile_" + std::to_string(index));

  const double(*coords)[6][2] = &coordinates[index];

  std::vector<G4TwoVector> vertices;

  for (int32_t i = 0; i < 6; ++i)
  {
    double x = (*coords)[i][0];
    double y = (*coords)[i][1];

    if (x == 0. && y == 0.)
      continue;

    vertices.emplace_back(x * mm, y * mm);
  }

  std::vector<G4ExtrudedSolid::ZSection> zsections = {
      G4ExtrudedSolid::ZSection(-dz, G4TwoVector(), 1.),
      G4ExtrudedSolid::ZSection(dz, G4TwoVector(), 1.),
  };

  return new G4ExtrudedSolid(label, vertices, zsections);
}
