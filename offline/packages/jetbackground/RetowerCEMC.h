#ifndef JETBACKGROUND_RETOWERCEMC_H
#define JETBACKGROUND_RETOWERCEMC_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

class RetowerCEMC : public SubsysReco
{
 public:
  RetowerCEMC(const std::string &name = "RetowerCEMC");
  ~RetowerCEMC() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  void SetEnergyDistribution(int val) { _weighted_energy_distribution = val; }
  void set_frac_cut(double frac_cut) { _frac_cut = frac_cut; }
  void set_towerinfo(bool use_towerinfo) { m_use_towerinfo = use_towerinfo; }
  void set_towerNodePrefix(const std::string &prefix)
  {
    m_towerNodePrefix = prefix;
    return;
  }

 private:
  int CreateNode(PHCompositeNode *topNode);
  void get_first_phi_index(PHCompositeNode *topNode);
  void get_fraction(PHCompositeNode *topNode);
  void get_weighted_fraction(PHCompositeNode *topNode);

  int _weighted_energy_distribution{1};
  double _frac_cut{0.5};
  bool m_use_towerinfo{false};
  std::string m_towerNodePrefix{"TOWERINFO_CALIB"};

  static const int neta_ihcal = 24;
  static const int neta_emcal = 96;
  static const int nphi_ihcal = 64;
  static const int nphi_emcal = 256;

  int retower_lowerbound_originaltower_ieta[neta_ihcal] = {0};
  int retower_upperbound_originaltower_ieta[neta_ihcal] = {0};
  double retower_lowerbound_originaltower_fraction[neta_ihcal] = {0.0};
  double retower_upperbound_originaltower_fraction[neta_ihcal] = {0.0};
  double retower_totalarea[neta_ihcal] = {0.0};
  int retower_first_lowerbound_originaltower_iphi{-1};

  double rawtower_e[neta_emcal][nphi_emcal] = {{0.0}};
  double rawtower_time[neta_emcal][nphi_emcal] = {{0.0}};
  int rawtower_status[neta_emcal][nphi_emcal] = {{0}};

  std::string EMTowerName;
  std::string IHTowerName;
  std::string EMRetowerName;
};

#endif
