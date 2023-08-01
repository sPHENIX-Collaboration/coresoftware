#include "FastJetOptions.h"
#include <phool/phool.h>  // for PHWHERE
#include <iostream>
#include <cassert>

FastJetOptions& FastJetOptions::update(std::vector<FastJetOptItem> input) {
  int size = input.size();
  int i =0;
  while (i<size) {
    auto& item = input[i];
    if (item.is_algo) {
      algo = item.algo;
    } else if (item.is_opt) {
      if (item.opt == JET_R) {
        jet_R = next_val(i, input);
      } else if (item.opt == JET_MIN_PT) {
        jet_min_pt = next_val(i, input);
      } else if (item.opt == JET_MAX_ETA) {
        jet_max_eta = next_val(i, input);
        handset_maxeta = true;
      } else if (item.opt == CUT_EDGE_ETA) {
        cut_edge_eta = true;
      } else if (item.opt == CONSTITUENT_MIN_PT) {
        constituent_min_pt = next_val(i, input);
        min_const_pt_iscut = (constituent_min_pt != 0.);
      } else if (item.opt == DO_SOFTDROP) {
        doSoftDrop = true;
      } else if (item.opt == SD_BETA) {
        SD_beta = next_val(i, input);
      } else if (item.opt == SD_ZCUT) {
        SD_zcut = next_val(i, input);
      } else if (item.opt == SD_JET_MIN_PT) {
        SD_jet_min_pt = next_val(i, input);
      } else if (item.opt == CALC_AREA) {
        calc_area = true;
      } else if (item.opt == GHOST_AREA) {
        ghost_area = next_val(i, input);
      } else if (item.opt == GHOST_MAX_RAP) {
        ghost_max_rap = next_val(i, input);
      } else if (item.opt == CALC_RhoMedDens) {
        calc_jetmedbkgdens = true;
      } else if (item.opt == CUT_RhoMedNHardest) {
        nhardestcut_jetmedbkgdens = static_cast<int>(next_val(i, input));
      } else if (item.opt == VERBOSITY) {
        verbosity = static_cast<int>(next_val(i, input));
      } else if (item.opt == DONT_SAVE_JET_COMPONENTS) {
        save_jet_components = false;
      } else if (item.opt == SAVE_JET_COMPONENTS) {
        save_jet_components = true;
      }
    }
    ++i;
  }
  return *this;
}

float FastJetOptions::next_val(int& i, std::vector<FastJetOptItem>& inputs) {
  if ( inputs.size() > static_cast<size_t>(i+1) && inputs[i+1].is_val ) {
    ++i;
    return inputs[i].val;
  } else {
    std::cout << PHWHERE << std::endl;
    std::cout << "Error in FastJetOptions, option required to have a value which isn't provided." << std::endl;
    assert(false);    
  }
}

void FastJetOptions::print(std::ostream& os) {
  os << " FastJet input options: " << std::endl
     << " - R: " << jet_R << std::endl
     << " - algorithm: " <<
        (algo==Jet::ALGO::ANTIKT ? "ANTIKT"
       : algo==Jet::ALGO::KT     ? "KT"
       : algo==Jet::ALGO::CAMBRIDGE ? "CAMBRIDGE"
       : "none") << std::endl
      << " - save jet components ids: " << (save_jet_components ? "yes" : "no") << std::endl
      << " - verbosity: " << verbosity << std::endl;
  if (min_const_pt_iscut) {
    os << " - minimum constituent pT cut: " << constituent_min_pt << std::endl;
  }
  os << " - minimum jet pt: " << jet_min_pt << std::endl;
  if (jet_max_eta != 1.1) os << " - maximum |eta_jet|: " << jet_max_eta << std::endl;
  if (doSoftDrop) {
    os << " - do softdrop with Beta(" << SD_beta <<") and Zcut(" << SD_zcut 
      <<") for jets w/pT>" << SD_jet_min_pt << std::endl;
  }
  if (calc_area) {
    os << " - calculate jet areas (using KT jets) with ghost_area("
      <<ghost_area<<") and ghost_max_rap("<<ghost_max_rap<<")"<<std::endl;
  }
  if (calc_jetmedbkgdens) {
    os << " - calculate jet median background estimator density with cutting " 
       << nhardestcut_jetmedbkgdens << " hardest jets" << std::endl;
  }
}

void FastJetOptions::initialize() {
  // set some required derivefd options when first running FastJetAlgo

  // set values if calculating jet areas and rapidities
  if (calc_area && ghost_max_rap == 0) {
    ghost_max_rap = 1.1;
  }

  if (calc_jetmedbkgdens) {
    if (jet_eta_iscut) etahardestcut_jetmedbkgdens = jet_max_eta;
    else etahardestcut_jetmedbkgdens = 1.1 - jet_R;
  }

  jet_eta_iscut = (jet_max_eta != 1.1);
  jet_min_pt_iscut = (jet_min_pt != 0.);
  has_jet_selection = (jet_eta_iscut || jet_min_pt_iscut);
}
