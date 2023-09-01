#ifndef JETBASE_FASTJETOPTIONS
#define JETBASE_FASTJETOPTIONS

#include "Jet.h"

#include <ostream>
#include <vector>

// This is a structure that the FastJetAlgo class uses for its optoins. The users
// can initialize it as a double-brace enclosed list which is a mix of
// Something like this:
//
// FastJetOptions fj_opt {{ Jet::ALGO::ANTIKT, DO_SOFTDROP, SD_BETA, 0.0, SD_ZCUT 0.1, JET_R, 0.4 }};
// jet_reco_obj->add_algo(new FastJetAlgo(fj_opt), "AntiKt_Tower_r04");
// /* can also update it */
// jet_reco_obj->add_algo(new FastJetAlgo(fj_opt({{JET_R,0.5}}), "AntiKt_Tower_r05");

enum FastJetOptEnum
{
  JET_R,
  JET_MIN_PT,
  JET_MAX_ETA  // defaults to 1.1
  ,
  CUT_EDGE_ETA  // set jet_max_eta to 1.1-jet_R

  ,
  CONSTITUENT_MIN_PT

  ,
  DO_SOFTDROP,
  SD_BETA,
  SD_ZCUT,
  SD_JET_MIN_PT  // defaults to 5.

  ,
  CALC_AREA,
  GHOST_AREA,
  GHOST_MAX_RAP

  ,
  CALC_RhoMedDens,
  CUT_RhoMedNHardest,
  NONE

  ,
  SAVE_JET_COMPONENTS,
  DONT_SAVE_JET_COMPONENTS

  ,
  VERBOSITY
};

struct FastJetOptItem
{  // All the things you can feed into FastJetAlgo
  FastJetOptEnum opt{FastJetOptEnum::NONE};
  bool is_opt{false};

  float val{0.};
  bool is_val{false};

  Jet::ALGO algo{};
  bool is_algo{false};

  FastJetOptItem(float _val)
    : val{_val}
    , is_val{true} {};
  FastJetOptItem(Jet::ALGO _algo)
    : algo{_algo}
    , is_algo{true} {};
  FastJetOptItem(FastJetOptEnum _opt)
    : opt{_opt}
    , is_opt{true} {};
};

struct FastJetOptions
{
  FastJetOptions(){};
  FastJetOptions(const std::vector<FastJetOptItem>& _vitem) { update(_vitem); };
  FastJetOptions& update(std::vector<FastJetOptItem>);
  FastJetOptions& operator()(const std::vector<FastJetOptItem>& _vitem) { return update(_vitem); };

  void print(std::ostream& os = std::cout);

  float next_val(int& i, std::vector<FastJetOptItem>&);
  // default options
  float jet_R = 0.4;
  Jet::ALGO algo = Jet::ALGO::ANTIKT;

  float jet_max_eta = 1.1;
  bool handset_maxeta = false;  // If user uses jet_max_eta, then cut_edge_eta won't reset eta at all
  bool cut_edge_eta = false;

  float jet_min_pt = 5;

  float constituent_min_pt = 0.;

  bool save_jet_components = true;

  bool doSoftDrop = false;
  float SD_beta = 0;
  float SD_zcut = 0;
  float SD_jet_min_pt = 5.;

  bool calc_area = false;
  float ghost_area = 0.001;
  float ghost_max_rap = 1.2;

  bool calc_jetmedbkgdens = false;
  float nhardestcut_jetmedbkgdens = 2;
  float etahardestcut_jetmedbkgdens = 0.;
  float rhoMedEtaNHardCut = 0.;

  int verbosity = 0;

  // for convenience when running FastJetAlgo
  bool has_jet_selection = false;  // set when initialized
  bool jet_eta_iscut = false;
  bool min_const_pt_iscut = false;
  bool jet_min_pt_iscut = false;

  void initialize();  // updates run with the first call
};

#endif
