// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4PSTOFDETECTOR_H
#define G4DETECTORS_PHG4PSTOFDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <map>
#include <string>  // for string

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHParametersContainer;
class PHG4Subsystem;

class PHG4PSTOFDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4PSTOFDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParametersContainer *params_array, const std::string &dnam);

  //! destructor
  ~PHG4PSTOFDetector() override {}

  //! construct
  void ConstructMe(G4LogicalVolume *world) override;

  void Print(const std::string &what = "ALL") const override;

  //!@name volume accessors
  //@{
  int IsInPSTOF(G4VPhysicalVolume *) const;
  //@}

  void SuperDetector(const std::string &name) { superdetector = name; }
  const std::string SuperDetector() const { return superdetector; }

 protected:
  int IsActive;
  int IsAbsorberActive;
  int nmod;
  int nrows;
  PHParametersContainer *paramscontainer;
  std::map<G4VPhysicalVolume *, int> active_phys_vols;

  std::string superdetector;
};

#endif
