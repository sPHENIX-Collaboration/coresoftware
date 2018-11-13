// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef G4DETECTORS_PHG4BLOCKDETECTOR_H
#define G4DETECTORS_PHG4BLOCKDETECTOR_H

#include <g4main/PHG4Detector.h>

class G4LogicalVolume;
class PHParameters;
class G4VPhysicalVolume;

class PHG4BlockDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4BlockDetector(PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam = "BLOCK", const int lyr = 0);

  //! destructor
  virtual ~PHG4BlockDetector(void)
  {
  }

  //! construct
  virtual void Construct(G4LogicalVolume *world);

  //!@name volume accessors
  //@{
  bool IsInBlock(G4VPhysicalVolume *) const;
  //@}

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }
  int get_Layer() const { return m_Layer; }

 private:
  PHParameters *m_Params;

  G4VPhysicalVolume *m_BlockPhysi;

  int m_Layer;
  std::string m_SuperDetector;
};

#endif
