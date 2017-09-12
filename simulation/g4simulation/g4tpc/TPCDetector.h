#ifndef __TPC_PHG4DETECTOR_H__
#define __TPC_PHG4DETECTOR_H__

#include <g4main/PHG4Detector.h>

class G4UserLimits;
class G4LogicalVolume;
class G4VPhysicalVolume;

class TPCDetector: public PHG4Detector {
 public:
  TPCDetector(PHCompositeNode *Node);
  virtual ~TPCDetector() {}
  void BuildCage(G4LogicalVolume *wrld);
  void BuildFiducial(G4LogicalVolume *wrld);
  void Construct(G4LogicalVolume *wrld);
  void SetVerbosity(int v) {fVerbosity=v;}
  bool IsThisActive(G4VPhysicalVolume *test);
  void SetMaxStep(G4double maxStep);

 protected:
  int fVerbosity;
  G4VPhysicalVolume *fTPC;
  G4UserLimits *fLimits;
};

#endif
