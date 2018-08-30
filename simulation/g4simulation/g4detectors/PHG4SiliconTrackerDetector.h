// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4DETECTORS_PHG4SILICONTRACKERDETECTOR_H
#define G4DETECTORS_PHG4SILICONTRACKERDETECTOR_H

#include "PHG4SiliconTrackerDefs.h"

#include <g4main/PHG4Detector.h>

#include <map>
#include <set>
#include <tuple>
#include <vector>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHParametersContainer;

class PHG4SiliconTrackerDetector : public PHG4Detector
{
 public:
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
    m_SuperDetector = name;
  }
  const std::string SuperDetector() const
  {
    return m_SuperDetector;
  }
  void Detector(const std::string &name)
  {
    m_DetectorType = name;
  }
  const std::string Detector() const
  {
    return m_DetectorType;
  }

  std::map<G4VPhysicalVolume *,std::tuple<int, int, int, int>>::const_iterator get_ActiveVolumeTuple(G4VPhysicalVolume *physvol) const;
  std::map<G4LogicalVolume *,std::tuple<int, int>>::const_iterator get_PassiveVolumeTuple(G4LogicalVolume *logvol) const;

 private:
  void AddGeometryNode();
  int ConstructSiliconTracker(G4LogicalVolume *sandwich);

  PHParametersContainer *m_ParamsContainer;

  std::string m_DetectorType;
  std::string m_SuperDetector;

  G4double m_SensorRadiusInner[4];
  G4double m_SensorRadiusOuter[4];
  G4double m_PosZ[4][2];
  G4double m_StripOffsetX[4];

  std::set<G4LogicalVolume *> m_ActiveLogVols;
  std::map<int, int> m_IsActiveMap;
  std::map<int, int> m_IsAbsorberActiveMap;
  std::pair<std::vector<std::pair<int, int>>::const_iterator, std::vector<std::pair<int, int>>::const_iterator> m_LayerBeginEndIteratorPair;
  std::map<G4VPhysicalVolume *,std::tuple<int, int, int, int>> m_ActiveVolumeTuple;
  std::map<G4LogicalVolume *,std::tuple<int, int>> m_PassiveVolumeTuple;
};

#endif
