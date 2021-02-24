#ifndef __SHEPGEN_H__
#define __SHEPGEN_H__

#include <fun4all/SubsysReco.h>
#include <phhepmc/PHHepMCGenHelper.h>

#if !defined(__CINT__) || defined(__CLING__)
#include <gsl/gsl_rng.h>
#endif

#include <iostream>
#include <string>

class PHCompositeNode;
class PHHepMCGenEvent;

class HGenManager;
class HLorentzVector;

namespace HepMC
{
  class GenEvent;
};

class sHEPGen : public SubsysReco, public PHHepMCGenHelper
{
 public:
  sHEPGen(const std::string &name = "sHEPGen");
  virtual ~sHEPGen();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void set_datacard_file(const char *cfg_file)
  {
    if (cfg_file) _datacardFile = cfg_file;
  }

  void set_momentum_electron(double emom)
  {
    _p_electron_lab = emom;
  }

  void set_momentum_hadron(double hmom)
  {
    _p_hadron_lab = hmom;
  }

  void beam_vertex_parameters(double beamX,
                              double beamY,
                              double beamZ,
                              double beamXsigma,
                              double beamYsigma,
                              double beamZsigma)
  {
    set_vertex_distribution_mean(beamX, beamY, beamZ, 0);
    set_vertex_distribution_width(beamXsigma, beamYsigma, beamZsigma, 0);
  }

 private:
  /** Print HEPGen++ logo to screen
   */
  void printlogo();


  int _eventcount;

  double _p_electron_lab;
  double _p_hadron_lab;

  HLorentzVector *_p4_electron_lab;
  HLorentzVector *_p4_hadron_lab;
  HLorentzVector *_p4_hadron_lab_invert;
  HLorentzVector *_p4_electron_prest;
  HLorentzVector *_p4_hadron_prest;

  // HEPGen++
  HGenManager *_hgenManager;

  std::string _datacardFile;
};

#endif /* __SHEPGEN_H__ */
