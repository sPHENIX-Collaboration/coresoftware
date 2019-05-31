#ifndef G4DETECTORS_PHG4PSTOFDETECTOR_H
#define G4DETECTORS_PHG4PSTOFDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <map>
#include <set>
#include <vector>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHParametersContainer;

class PHG4PSTOFDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4PSTOFDetector(PHCompositeNode *Node, PHParametersContainer *params_array, const std::string &dnam = "PSTOF");

  //! destructor
  virtual ~PHG4PSTOFDetector(){}

  //! construct
  virtual void Construct(G4LogicalVolume *world);

  virtual void Print(const std::string &what = "ALL") const;

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
