#ifndef PHG4SiliconTrackerDetector_h
#define PHG4SiliconTrackerDetector_h

#include "PHG4SiliconTrackerDefs.h"

#include <g4main/PHG4Detector.h>

#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Types.hh>
#include <Geant4/globals.hh>

#include <map>
#include <set>
#include <utility>
#include <vector>

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSolid;
class PHParametersContainer;

class PHG4SiliconTrackerDetector : public PHG4Detector
{
 public:
  typedef std::vector<std::pair<int, int>> vpair;
  //! constructor
  PHG4SiliconTrackerDetector(PHCompositeNode *Node, PHParametersContainer *parameters, const std::string &dnam, const std::pair<std::vector<std::pair<int, int>>::const_iterator, std::vector<std::pair<int, int>>::const_iterator> &layer_b_e);

  //! destructor
  virtual ~PHG4SiliconTrackerDetector() {}
  //! construct
  virtual void Construct(G4LogicalVolume *world);

  //!@name volume accessors
  //@{
  int IsInSiliconTracker(G4VPhysicalVolume *) const;
  //@}

  void SuperDetector(const std::string &name)
  {
    superdetector = name;
  }
  const std::string SuperDetector() const
  {
    return superdetector;
  }
  void Detector(const std::string &name)
  {
    detector_type = name;
  }
  const std::string Detector() const
  {
    return detector_type;
  }

 private:
  void AddGeometryNode();
  int ConstructSiliconTracker(G4LogicalVolume *sandwich);
  int DisplayVolume(G4VSolid *volume, G4LogicalVolume *logvol, G4RotationMatrix *rotm = nullptr);

  PHParametersContainer *paramscontainer;
  //  vpair layerconfig_;
  /* unsigned int nlayer_; */
  /* int layermin_; */
  /* int layermax_; */

//  G4double overlap_fraction;

  std::string detector_type;
  std::string superdetector;

  G4double sensor_radius_inner[4];
  G4double sensor_radius_outer[4];
  G4double ladder_radius_inner[4];
  G4double ladder_radius_outer[4];
  G4double posz[4][2];
  G4double strip_x_offset[4];

  std::set<G4LogicalVolume *> absorberlogvols;
  std::set<G4LogicalVolume *> activelogvols;
  std::map<int, int> IsActive;
  std::map<int, int> IsAbsorberActive;
  std::pair<std::vector<std::pair<int, int>>::const_iterator, std::vector<std::pair<int, int>>::const_iterator> layer_begin_end;
};

#endif
