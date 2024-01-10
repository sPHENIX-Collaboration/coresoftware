// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef REACTIONPLANEAFTERBURNER_H
#define REACTIONPLANEAFTERBURNER_H

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>

#include <string>

class PHCompositeNode;

class ReactionPlaneAfterburner : public SubsysReco
{
 public:

  ReactionPlaneAfterburner(const std::string &name = "ReactionPlaneAfterburner");

  ~ReactionPlaneAfterburner() override;

  int Init(PHCompositeNode *topNode) override;


  int process_event(PHCompositeNode *topNode) override;


 private:
  gsl_rng *RandomGenerator {nullptr};

};

#endif // REACTIONPLANEAFTERBURNER_H
