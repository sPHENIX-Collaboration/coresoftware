#ifndef PHG4PSTOFDetector_h
#define PHG4PSTOFDetector_h

#include <g4main/PHG4Detector.h>

#include <map>
#include <set>
#include <vector>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHG4ParametersContainer;

class PHG4PSTOFDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4PSTOFDetector(PHCompositeNode *Node, PHG4ParametersContainer *params_array, const std::string &dnam = "PSTOF");

  //! destructor
  virtual ~PHG4PSTOFDetector();

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
  std::set<G4VPhysicalVolume *> active_phys_vols;
  G4LogicalVolume *active_volume;
  PHG4ParametersContainer *paramscontainer;

  std::string superdetector;
};

#endif
