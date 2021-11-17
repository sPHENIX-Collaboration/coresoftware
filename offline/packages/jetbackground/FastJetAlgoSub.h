#ifndef JETBACKGROUND_FASTJETALGOSUB_H
#define JETBACKGROUND_FASTJETALGOSUB_H

#include <g4jets/JetAlgo.h>

#include <g4jets/Jet.h>

#include <iostream>
#include <vector>

class FastJetAlgoSub : public JetAlgo
{
 public:
  FastJetAlgoSub(Jet::ALGO algo, float par, float verbosity = 0);
  ~FastJetAlgoSub() override {}

  void identify(std::ostream& os = std::cout) override;
  Jet::ALGO get_algo() override { return _algo; }
  float get_par() override { return _par; }

  std::vector<Jet*> get_jets(std::vector<Jet*> particles) override;

 private:
  int _verbosity;
  Jet::ALGO _algo;
  float _par;
};

#endif
