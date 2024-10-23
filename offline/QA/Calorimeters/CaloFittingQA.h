#ifndef CALOFITTINGQA_CALOFITTINGQA_H
#define CALOFITTINGQA_CALOFITTINGQA_H

#include <fun4all/SubsysReco.h>

#include <caloreco/CaloTowerDefs.h>

#include <string>
#include <vector>

// Forward declarations
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
  ~CaloFittingQA() override = default;

  //! full initialization
  int Init(PHCompositeNode*) override;

  //! event processing method
  int process_event(PHCompositeNode*) override;

  //! end of run method
  int End(PHCompositeNode*) override;

  int process_towers(PHCompositeNode*);
  int process_data(PHCompositeNode *topNode, CaloTowerDefs::DetectorSystem dettype, std::vector<std::vector<float>> &waveforms);
  bool skipChannel(int ich, int pid, CaloTowerDefs::DetectorSystem dettype);

  void set_debug(bool debug) { m_debug = debug; }

  void set_offlineflag(const bool f = true)
  {
    m_UseOfflinePacketFlag = f;
  }
  void set_simflag(const bool f = false)
  {
    m_SimFlag = f;
  }

 private:
  void createHistos();
  std::string getHistoPrefix() const;
  TProfile2D* h_cemc_etaphi_ZScrosscalib{nullptr};
  TProfile2D* h_ihcal_etaphi_ZScrosscalib{nullptr};
  TProfile2D* h_ohcal_etaphi_ZScrosscalib{nullptr};

  int _eventcounter{0};

  bool m_debug{false};
  bool m_UseOfflinePacketFlag{true};
  bool m_SimFlag{false};

  std::string m_outputFileName;
  std::string OutputFileName;
};

#endif
