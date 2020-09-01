/* vim: set sw=2 ft=cpp: */

#ifndef G4EPD_PHG4EPDETECTOR_H
#define G4EPD_PHG4EPDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <Geant4/G4ExtrudedSolid.hh>

#include <cstdint>
#include <map>
#include <string>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHParametersContainer;
class PHG4Subsystem;

class PHG4EPDetector : public PHG4Detector {
  public:
    PHG4EPDetector(PHG4Subsystem* subsys,
                   PHCompositeNode* node,
                   PHParametersContainer* params,
                   std::string const& name);

    void ConstructMe(G4LogicalVolume* world) override;

    bool contains(G4VPhysicalVolume*) const;

    uint32_t module_id_for(int32_t index, int32_t slice, int32_t side);
    uint32_t module_id_for(G4VPhysicalVolume* volume);

    void SuperDetector(std::string const& name) { superdetector = name; }
    const std::string SuperDetector() const { return superdetector; }

  private:
    G4ExtrudedSolid* construct_block(int32_t index);

    PHParametersContainer* m_params;

    std::map<G4VPhysicalVolume*, uint32_t> m_volumes;

    std::string superdetector;
};

#endif /* G4EPD_PHG4EPDETECTOR_H */
