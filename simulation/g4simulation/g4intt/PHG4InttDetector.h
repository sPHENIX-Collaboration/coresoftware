// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4INTT_PHG4INTTDETECTOR_H
#define G4INTT_PHG4INTTDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <map>
#include <set>
#include <string>  // for string
#include <tuple>
#include <utility>  // for pair
#include <vector>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4InttDisplayAction;
class PHG4Subsystem;
class PHParametersContainer;

class PHG4InttDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4InttDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParametersContainer *parameters, const std::string &dnam, const std::pair<std::vector<std::pair<int, int>>::const_iterator, std::vector<std::pair<int, int>>::const_iterator> &layer_b_e);

  //! destructor
  ~PHG4InttDetector() override {}
  //! construct
  void ConstructMe(G4LogicalVolume *world) override;

  //!@name volume accessors
  //@{
  int IsInIntt(G4VPhysicalVolume *) const;
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

  std::map<G4VPhysicalVolume *, std::tuple<int, int, int, int>>::const_iterator get_ActiveVolumeTuple(G4VPhysicalVolume *physvol) const;
  std::map<G4LogicalVolume *, std::tuple<int, int>>::const_iterator get_PassiveVolumeTuple(G4LogicalVolume *logvol) const;

 private:
  void AddGeometryNode();
  int ConstructIntt(G4LogicalVolume *sandwich);

  PHG4InttDisplayAction *m_DisplayAction = nullptr;
  PHParametersContainer *m_ParamsContainer = nullptr;

  std::string m_DetectorType;
  std::string m_SuperDetector;

  int m_IsSupportActive = 0;
  int m_IsEndcapActive = 0;

  double m_PosZ[8][2];
  double m_SensorRadius[8];
  double m_StripOffsetX[8];

  std::set<G4LogicalVolume *> m_ActiveLogVols;
  std::map<int, int> m_IsActiveMap;
  std::map<int, int> m_IsAbsorberActiveMap;
  std::pair<std::vector<std::pair<int, int>>::const_iterator, std::vector<std::pair<int, int>>::const_iterator> m_LayerBeginEndIteratorPair;
  std::map<G4VPhysicalVolume *, std::tuple<int, int, int, int>> m_ActiveVolumeTuple;
  std::map<G4LogicalVolume *, std::tuple<int, int>> m_PassiveVolumeTuple;
};

#endif
