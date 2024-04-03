#include "FastJetOptions.h"

#include <phool/phool.h>  // for PHWHERE
                         
#include <cassert>
#include <iostream>

/* Need simpler input logic -- it is too clever by half in this first iteration:
  float jet_max_eta  = 10; // essentially no cut, unless the user sets it




*/

FastJetOptions& FastJetOptions::update(std::vector<FastJetOptItem> input)
{
  int size = input.size();
  int i = 0;
  while (i < size)
  {
    auto& item = input[i];
    if (item.is_algo)
    {
      algo = item.algo;
    }
    else if (item.is_opt)
    {
      if (item.opt == JET_R)
      {
        jet_R = next_val(i, input);
      }
      else if (item.opt == JET_MIN_PT)
      {
        use_jet_min_pt = true;
        jet_min_pt = next_val(i, input);
      }
      else if (item.opt == JET_MAX_ETA)
      {
        use_jet_max_eta = true;
        jet_max_eta = next_val(i, input);
      }
      else if (item.opt == CONSTITUENT_MIN_PT)
      {
        use_constituent_min_pt = true;
        constituent_min_pt = next_val(i, input);
      }
      else if (item.opt == DO_SOFTDROP)
      {
        doSoftDrop = true;
      }
      else if (item.opt == SD_BETA)
      {
        SD_beta = next_val(i, input);
      }
      else if (item.opt == SD_ZCUT)
      {
        SD_zcut = next_val(i, input);
      }
      else if (item.opt == SD_JET_MIN_PT)
      {
        SD_jet_min_pt = next_val(i, input);
      }
      else if (item.opt == CALC_AREA)
      {
        calc_area = true;
      }
      else if (item.opt == GHOST_AREA)
      {
        ghost_area = next_val(i, input);
      }
      else if (item.opt == GHOST_MAX_RAP)
      {
        ghost_max_rap = next_val(i, input);
      }
      else if (item.opt == CALC_RhoMedDens)
      {
        calc_jetmedbkgdens = true;
      }
      else if (item.opt == CUT_RhoMedNHardest)
      {
        nhardestcut_jetmedbkgdens = static_cast<int>(next_val(i, input));
      }
      else if (item.opt == FJCS_doConstSub)
      {
        cs_calc_constsub = true;
      }
      else if (item.opt == FJCS_max_eta)
      {
        cs_max_eta = static_cast<float>(next_val(i, input));
      }
      else if (item.opt == FJCS_GridMedBkgEst_Size)
      {
        cs_gridmedestsize = static_cast<float>(next_val(i, input));
      }
      else if (item.opt == FJCS_max_dist)
      {
        cs_max_dist = static_cast<float>(next_val(i, input));
      }
      else if (item.opt == FJCS_alpha)
      {
        cs_alpha = static_cast<float>(next_val(i, input));
      }
      else if (item.opt == FJCS_max_pt)
      {
        cs_max_pt = static_cast<float>(next_val(i, input));
      }
      else if (item.opt == FJCS_ghost_area)
      {
        cs_ghost_area = static_cast<float>(next_val(i, input));
      }
      else if (item.opt == VERBOSITY)
      {
        verbosity = static_cast<int>(next_val(i, input));
      }
      else if (item.opt == DONT_SAVE_JET_COMPONENTS)
      {
        save_jet_components = false;
      }
      else if (item.opt == SAVE_JET_COMPONENTS)
      {
        save_jet_components = true;
      }
    }
    ++i;
  }
  return *this;
}

float FastJetOptions::next_val(int& i, std::vector<FastJetOptItem>& inputs)
{
  if (inputs.size() > static_cast<size_t>(i + 1) && inputs[i + 1].is_val)
  {
    ++i;
    return inputs[i].val;
  }
  else
  {
    std::cout << PHWHERE << std::endl;
    std::cout << "Error in FastJetOptions, option required to have a value which isn't provided." << std::endl;
    assert(false);
  }
}

void FastJetOptions::print(std::ostream& os)
{
  initialize();

  os << "FastJetOptions (input options for fastjet in FastJetAlgp)" << std::endl;
  os << " FastJet input options: " << std::endl
     << " - R: " << jet_R << std::endl
     << " - algorithm: " << (algo == Jet::ALGO::ANTIKT ? "ANTIKT" : algo == Jet::ALGO::KT      ? "KT"
                                                                : algo == Jet::ALGO::CAMBRIDGE ? "CAMBRIDGE"
                                                                                               : "none")
     << std::endl
     << " - save jet components ids: " << (save_jet_components ? "yes" : "no") << std::endl
     << " - verbosity: " << verbosity << std::endl;
  if (use_constituent_min_pt)
  {
    os << " - minimum constituent pT cut: " << constituent_min_pt << std::endl;
  }
  if (use_jet_min_pt) {
    os << " - minimum jet pt: " << jet_min_pt << std::endl;
  }
  if (use_jet_max_eta) {
    os << " - maximum |eta_jet|: " << jet_max_eta << std::endl;
  }
  if (doSoftDrop)
  {
    os << " - do softdrop with Beta(" << SD_beta << ") and Zcut(" << SD_zcut
       << ") for jets w/pT>" << SD_jet_min_pt << std::endl;
  }
  if (calc_area)
  {
    os << " - calculate jet areas (using KT jets) with ghost_area("
       << ghost_area << ") and ghost_max_rap(" << ghost_max_rap << ")" << std::endl;
  }
  if (calc_jetmedbkgdens)
  {
    os << " - calculate jet median background estimator density with cutting "
       << nhardestcut_jetmedbkgdens << " hardest jets within |eta|<" << etahardestcut_jetmedbkgdens << std::endl;
  }
  
  if (cs_calc_constsub) {
    os << " - calculate jet background constituent subtractor " << std::endl;
  }
}

void FastJetOptions::initialize()
{
  // set some required derived options when first running FastJetAlgo
  if (calc_jetmedbkgdens) {
    calc_area = true;
  }

  if (calc_area && ghost_max_rap == 0) {
    if (use_jet_max_eta) ghost_max_rap = (jet_max_eta + jet_R);
    else ghost_max_rap = 5.;
  }

  if (calc_jetmedbkgdens) {
    if (nhardestcut_jetmedbkgdens > 0 && etahardestcut_jetmedbkgdens == 0) {
      if (use_jet_max_eta) etahardestcut_jetmedbkgdens = jet_max_eta;
      else etahardestcut_jetmedbkgdens = ghost_max_rap;
    }
  }
  
  use_jet_selection = (use_jet_max_eta || use_jet_min_pt);
}
