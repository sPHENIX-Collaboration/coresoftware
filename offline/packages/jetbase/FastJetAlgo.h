#ifndef JETBASE_FASTJETALGO_H
#define JETBASE_FASTJETALGO_H

#include "Jet.h"
#include "JetAlgo.h"

#include <cmath>     // for NAN
#include <climits>     // for NAN
#include <iostream>  // for cout, ostream
#include <vector>    // for vector

namespace fastjet {
  class PseudoJet;
}
class JetContainer;

class FastJetAlgo : public JetAlgo
{
 public:
  FastJetAlgo(Jet::ALGO algo, float par, int verbosity=0, Jet::SORT sort=Jet::SORT::PT);
  ~FastJetAlgo() override {}


  void identify(std::ostream& os = std::cout) override;
  Jet::ALGO get_algo() override { return m_AlgoFlag; }
  /* float get_par() override { return m_Par; } */

  void set_do_SoftDrop(bool do_SD) { m_SDFlag = do_SD; }

  void set_SoftDrop_beta(float beta) { m_SDBeta = beta; }

  void set_SoftDrop_zcut(float zcut) { m_SDZCut = zcut; }

  std::vector<Jet*> get_jets(std::vector<Jet*> particles) override;

  void cluster_and_fill(std::vector<Jet*>& part_in, JetContainer* jets_out) override;

  // Fastjet calculate jet area
  void set_do_JetArea(bool do_JA) { m_JetAreaFlag = do_JA; }
  void set_do_JetArea(float ghostArea, float ghostMaxRap=0) { 
    m_GhostArea = ghostArea; m_GhostMaxRap=ghostMaxRap;
  }

  // Set if doing fastjet sort on jets
  // takes only fastjet sorting algorithms for E, ETA, PT;
  // defaults to PT, but can also select Jet:NO_SORT
  void set_jet_sorting (Jet::SORT sort) { m_whichsort = sort; }

  // Fastjet Rho Median
  void set_do_RhoMedian(bool do_RM) { m_RhoMedianFlag = do_RM; }

  float get_RhoMedian() { return m_RhoMedian; }

  void set_RapCutHardest (float _) { m_RapCutHardest = _; }
  void set_CutNHardest   (int   _) { m_CutNHardest   = _; }

 private:

  bool m_first_cluster_call { true };

  int       m_Verbosity = 0;
  Jet::ALGO m_AlgoFlag  = Jet::NONE;
  float     m_Par       = NAN; // Jet_R value

  bool      m_SDFlag    = false;
  float     m_SDBeta    = NAN;
  float     m_SDZCut    = NAN;
  Jet::SORT m_whichsort;

  // options to calculate jet areas
  bool  m_JetAreaFlag { false };
  float m_GhostArea   { 0.01  };
  float m_GhostMaxRap { 1.2   };

  // option to calculate background median density
  bool  m_RhoMedianFlag  { false };
  int   m_CutNHardest    { 2 };
  float m_RapCutHardest  { 0 }; // if 0 default to 1.1 - m_Par
  float m_RhoMedian      { NAN };
  // otherwise default ghose area to m_GhostArea and kt_algorithm

  // for convenience save indices of the zg, Rg, mu, and area for jets in the JetContainer
  unsigned int m_zg_index   { UINT_MAX };
  unsigned int m_Rg_index   { UINT_MAX };
  unsigned int m_mu_index   { UINT_MAX };
  unsigned int m_area_index { UINT_MAX };

  void fillJetContainer(std::vector<fastjet::PseudoJet>* jets, JetContainer* jetconn, std::vector<Jet*>& particles);
};

#endif
