// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4JLEIC_G4JLEICDIRCDETECTOR_H
#define G4JLEIC_G4JLEICDIRCDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <set>
#include <string>                 // for string

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHParameters;

class G4JLeicDIRCDetector : public PHG4Detector
{
 public:
  //! constructor
  G4JLeicDIRCDetector(PHCompositeNode *Node, PHParameters *params_array, const std::string &dnam = "DIRC");

  //! destructor
  virtual ~G4JLeicDIRCDetector(){}

  //! construct
  virtual void Construct(G4LogicalVolume *world);

  virtual void Print(const std::string &what = "ALL") const;

  //!@name volume accessors
  //@{
  int IsInDIRC(G4VPhysicalVolume *) const;
  //@}

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }
 
protected:
  int m_IsActiveFlag;
  int m_IsAbsorberActiveFlag;
  int m_Layers;
  PHParameters *m_Params;
  std::set<G4VPhysicalVolume *> m_PhysicalVolumesSet;

  std::string m_SuperDetector;
};

#endif // G4JLEIC_G4JLEICDIRCDETECTOR_H
