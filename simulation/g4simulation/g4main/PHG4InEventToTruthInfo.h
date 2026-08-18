// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4INEVENTTOTRUTHINFO_H
#define G4MAIN_PHG4INEVENTTOTRUTHINFO_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class PHG4Particle;
class PHG4TruthInfoContainer;
class PHG4VtxPoint;

class PHG4InEventToTruthInfo : public SubsysReco
{
 public:
  explicit PHG4InEventToTruthInfo(const std::string &name = "PHG4INEVENTTOTRUTHINFO");
  ~PHG4InEventToTruthInfo() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  int CreateNodeTree(PHCompositeNode *topNode);

 private:
  static PHG4Particle *makeParticleCopy(PHG4Particle *particle);
  static PHG4VtxPoint *makePrimaryVertexCopy(const PHG4VtxPoint *vertex, int vtxid);

};

#endif  // G4MAIN_PHG4INEVENTTOTRUTHINFO_H
