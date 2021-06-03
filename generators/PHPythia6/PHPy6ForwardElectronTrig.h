#ifndef PHPYTHIA6_PHPY6FORWARDELECTRONTRIG_H
#define PHPYTHIA6_PHPY6FORWARDELECTRONTRIG_H

#include "PHPy6GenTrigger.h"

#include <string>             // for string

namespace HepMC
{
  class GenEvent;
}

class PHPy6ForwardElectronTrig : public PHPy6GenTrigger
{
 public:
  //! constructor
  PHPy6ForwardElectronTrig(const std::string& name = "PHPy6ForwardElectronTrigger");

  //! destructor
  ~PHPy6ForwardElectronTrig(void) override {}

  bool Apply(const HepMC::GenEvent* evt) override;

  void set_electrons_required(int n)
  {
    n_ep_required = n;
  }
  void set_positrons_required(int n) { n_em_required = n; }
  void set_combined_required(int n) { n_comb_required = n; }
  void set_pt_required(float set_pt) { pt_required = set_pt; }
  void set_eta_range(float set_eta_low, float set_eta_high)
  {
    eta_low = set_eta_low;
    eta_high = set_eta_high;
  }

  void PrintConfig();

  void SetRequireElectron()
  {
    RequireElectron = true;
    RequirePositron = false;
    RequireOR = false;
    RequireAND = false;
    RequireCOMBO = false;
  }
  void SetRequirePositron()
  {
    RequireElectron = false;
    RequirePositron = true;
    RequireOR = false;
    RequireAND = false;
    RequireCOMBO = false;
  }
  void SetRequireOR()
  {
    RequireElectron = false;
    RequirePositron = false;
    RequireOR = true;
    RequireAND = false;
    RequireCOMBO = false;
  }
  void SetRequireAND()
  {
    RequireElectron = false;
    RequirePositron = false;
    RequireOR = false;
    RequireAND = true;
    RequireCOMBO = false;
  }
  void SetRequireCOMBINED()
  {
    RequireElectron = false;
    RequirePositron = false;
    RequireOR = false;
    RequireAND = false;
    RequireCOMBO = true;
  }

 protected:
  int ntriggered_forward_electron;
  int nconsidered_forward_electron;

  // trigger variables
  unsigned int n_ep_required;
  unsigned int n_em_required;
  unsigned int n_comb_required;
  float pt_required;
  float eta_low;
  float eta_high;

  bool RequireElectron;
  bool RequirePositron;
  bool RequireOR;
  bool RequireAND;
  bool RequireCOMBO;
};

#endif
