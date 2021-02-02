// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FERMIMOTIONAFTERBURNER_H
#define FERMIMOTIONAFTERBURNER_H

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>

#include <string>

class PHCompositeNode;

class FermimotionAfterburner : public SubsysReco
{
 public:
  FermimotionAfterburner(const std::string &name = "FermimotionAfterburner");

  virtual ~FermimotionAfterburner();

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

 private:
  void AddpF(PHCompositeNode *);

  gsl_rng *RandomGenerator;
};

#endif  // FERMIMOTIONAFTERBURNER_H
