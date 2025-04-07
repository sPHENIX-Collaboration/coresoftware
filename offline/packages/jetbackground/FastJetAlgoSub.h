#ifndef JETBACKGROUND_FASTJETALGOSUB_H
#define JETBACKGROUND_FASTJETALGOSUB_H

#include <jetbase/FastJetOptions.h>
#include <jetbase/Jet.h>
#include <jetbase/JetAlgo.h>

#include <iostream>
#include <vector>

class JetContainer;

class FastJetAlgoSub : public JetAlgo
{
 public:
  FastJetAlgoSub(const FastJetOptions& options);
  ~FastJetAlgoSub() override = default;

  //----------------------------------------------------------------------
  //  Legacy code interface. It is better to use FastJetOptions, but
  //  there is no harm is using these, as well.
  //----------------------------------------------------------------------
  FastJetAlgoSub(Jet::ALGO algo, float par, int verbosity = 0)
    : FastJetAlgoSub({{algo, JET_R, par, VERBOSITY, static_cast<float>(verbosity)}})
  {
  }
  //--end-legacy-code-interface-------------------------------------------

  void identify(std::ostream& os = std::cout) override;
  Jet::ALGO get_algo() override { return m_opt.algo; }
  float get_par() override { return m_opt.jet_R; }

  /* std::vector<Jet*> get_jets(std::vector<Jet*> particles) override; */
  void cluster_and_fill(std::vector<Jet*>& part_in, JetContainer* jets_out) override;

 private:
  FastJetOptions m_opt{};
};

#endif
