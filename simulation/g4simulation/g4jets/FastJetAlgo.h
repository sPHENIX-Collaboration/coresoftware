#ifndef G4JET_FASTJETALGO_H
#define G4JET_FASTJETALGO_H

#include "Jet.h"
#include "JetAlgo.h"

#include <cmath>     // for NAN
#include <iostream>  // for cout, ostream
#include <vector>    // for vector

class FastJetAlgo : public JetAlgo
{
 public:
  FastJetAlgo(Jet::ALGO algo, float par, int verbosity = 0);
  ~FastJetAlgo() override {}

  void identify(std::ostream& os = std::cout) override;
  Jet::ALGO get_algo() override { return m_AlgoFlag; }
  float get_par() override { return m_Par; }

  void set_do_SoftDrop(bool do_SD) { m_SDFlag = do_SD; }

  void set_SoftDrop_beta(float beta) { m_SDBeta = beta; }

  void set_SoftDrop_zcut(float zcut) { m_SDZCut = zcut; }

  std::vector<Jet*> get_jets(std::vector<Jet*> particles) override;

 private:
  int m_Verbosity = 0;
  Jet::ALGO m_AlgoFlag = Jet::NONE;
  float m_Par = NAN;

  bool m_SDFlag = false;
  float m_SDBeta = NAN;
  float m_SDZCut = NAN;
};

#endif
