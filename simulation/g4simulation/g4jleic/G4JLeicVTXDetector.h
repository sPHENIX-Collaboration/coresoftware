// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4JLEIC_G4JLEICVTXDETECTOR_H
#define G4JLEIC_G4JLEICVTXDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <map>
#include <string>                 // for string

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHParametersContainer;

class G4JLeicVTXDetector : public PHG4Detector
{
 public:
  //! constructor
  G4JLeicVTXDetector(PHCompositeNode *Node, PHParametersContainer *params_array, const std::string &dnam = "VTX");

  //! destructor
  virtual ~G4JLeicVTXDetector(){}

  //! construct
  virtual void Construct(G4LogicalVolume *world);

  virtual void Print(const std::string &what = "ALL") const;

  //!@name volume accessors
  //@{
  int IsInVTX(G4VPhysicalVolume *) const;
  //@}

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }
 
protected:
  int m_IsActiveFlag;
  int m_IsAbsorberActiveFlag;
  int m_Layers;
  PHParametersContainer *m_ParamsContainer;
  std::map<G4VPhysicalVolume *, int> m_PhysicalVolumesMap;

  std::string m_SuperDetector;
};

#endif // G4JLEIC_G4JLEICVTXDETECTOR_H
