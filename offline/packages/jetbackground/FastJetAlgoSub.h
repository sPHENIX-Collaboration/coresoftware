#ifndef JETBACKGROUND_FASTJETALGOSUB_H
#define JETBACKGROUND_FASTJETALGOSUB_H

#include <jetbase/Jet.h>
#include <jetbase/JetAlgo.h>

#include <iostream>
#include <vector>

class JetContainer;

class FastJetAlgoSub : public JetAlgo
{
 public:
  FastJetAlgoSub(Jet::ALGO algo, float par, int verbosity = 0);
  ~FastJetAlgoSub() override {}

  void identify(std::ostream& os = std::cout) override;
  Jet::ALGO get_algo() override { return _algo; }
  float get_par() override { return _par; }

  /* std::vector<Jet*> get_jets(std::vector<Jet*> particles) override; */
  void cluster_and_fill(std::vector<Jet*>& part_in, JetContainer* jets_out) override;

 private:
  int _verbosity;
  Jet::ALGO _algo;
  float _par;
};

#endif
