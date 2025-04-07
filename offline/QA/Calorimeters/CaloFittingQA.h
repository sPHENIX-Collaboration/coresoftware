#ifndef CALOFITTINGQA_CALOFITTINGQA_H
#define CALOFITTINGQA_CALOFITTINGQA_H

#include <fun4all/SubsysReco.h>

#include <caloreco/CaloTowerDefs.h>

#include <cdbobjects/CDBTTree.h>  // for CDBTTree

#include <string>
#include <vector>

// Forward declarations
class CDBTTree;
class PHCompositeNode;
class TH1;
class TH2;
class TProfile2D;
class TProfile;

class CaloFittingQA : public SubsysReco
{
 public:
  //! constructor
  CaloFittingQA(const std::string& name = "CaloFittingQA");  // const std::string &filename = "testQA.root"); //int nevents = 100);

  //! destructor
  ~CaloFittingQA() override;

  //! full initialization
  int Init(PHCompositeNode*) override;
  int InitRun(PHCompositeNode *) override;

  //! event processing method
  int process_event(PHCompositeNode*) override;

  //! end of run method
  int End(PHCompositeNode*) override;

  int process_towers(PHCompositeNode*);
  int process_data(PHCompositeNode *topNode, CaloTowerDefs::DetectorSystem dettype, std::vector<std::vector<float>> &waveforms);
  bool skipChannel(int ich, int pid, CaloTowerDefs::DetectorSystem dettype);

  void set_offlineflag(const bool f = true)
  {
    m_UseOfflinePacketFlag = f;
  }
  void set_simflag(const bool f = false)
  {
    m_SimFlag = f;
  }
  void set_cemc_lowadcthreshold(int lath)
  {
    m_cemc_adc_threshold = lath;
  }
  void set_cemc_highadcthreshold(int hath)
  {
    m_cemc_high_adc_threshold = hath;
  }
  void set_hcal_lowadcthreshold(int lath)
  {
    m_hcal_adc_threshold = lath;
  }
  void set_hcal_highadcthreshold(int hath)
  {
    m_hcal_high_adc_threshold = hath;
  }
  void set_cemc_hit_threshold(int hthres)
  {
    m_cemc_hit_threshold = hthres;
  }
  void set_ihcal_hit_threshold(int hthres)
  {
    m_ihcal_hit_threshold = hthres;
  }
  void set_ohcal_hit_threshold(int hthres)
  {
    m_ohcal_hit_threshold = hthres;
  }

 private:
  void createHistos();
  std::string getHistoPrefix() const;
  TProfile2D* h_cemc_etaphi_ZScrosscalib{nullptr};
  TProfile2D* h_ihcal_etaphi_ZScrosscalib{nullptr};
  TProfile2D* h_ohcal_etaphi_ZScrosscalib{nullptr};
  TProfile2D* h_cemc_etaphi_pedestal{nullptr};
  TProfile2D* h_ihcal_etaphi_pedestal{nullptr};
  TProfile2D* h_ohcal_etaphi_pedestal{nullptr};
  TH2* h_cemc_zs_frac_vs_multiplicity{nullptr};
  TH2* h_ihcal_zs_frac_vs_multiplicity{nullptr};
  TH2* h_ohcal_zs_frac_vs_multiplicity{nullptr};
  TH1* h_packet_events{nullptr};

  CDBTTree *cdbttree = nullptr;

  int _eventcounter{0};

  bool m_UseOfflinePacketFlag{true};
  bool m_SimFlag{false};

  float m_cemc_adc_threshold = 200.;
  float m_cemc_high_adc_threshold = 2000.;
  float m_hcal_adc_threshold = 100.;
  float m_hcal_high_adc_threshold = 2000.;

  float m_cemc_hit_threshold = 200; // ~ 300 MeV
  float m_ihcal_hit_threshold = 600; // ~ 300 MeV
  float m_ohcal_hit_threshold = 100; // ~ 300 MeV

  std::string m_calibName;
  std::string m_fieldname;
  std::string m_outputFileName;
  std::string OutputFileName;
};

#endif
