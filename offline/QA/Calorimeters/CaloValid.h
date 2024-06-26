#ifndef CALOVALID_CALOVALID_H
#define CALOVALID_CALOVALID_H

#include <fun4all/SubsysReco.h>

#include <string>

// Forward declarations
class PHCompositeNode;
class TH1;
class TH2;

class CaloValid : public SubsysReco
{
 public:
  //! constructor
  CaloValid(const std::string& name = "CaloValid");  // const std::string &filename = "testQA.root"); //int nevents = 100);

  //! destructor
  ~CaloValid() override;

  //! full initialization
  int Init(PHCompositeNode*) override;

  //! event processing method
  int process_event(PHCompositeNode*) override;

  //! end of run method
  int End(PHCompositeNode*) override;

  int process_g4hits(PHCompositeNode*);
  int process_g4cells(PHCompositeNode*);
  int process_towers(PHCompositeNode*);
  int process_clusters(PHCompositeNode*);

  void Detector(const std::string& name) { detector = name; }
  void set_timing_cut_width(const int& t) { _range = t; }

  void set_debug(bool debug) { m_debug = debug; }
  TH2* LogYHist2D(const std::string& name, const std::string& title, int, double, double, int, double, double);

 private:
  int Getpeaktime(TH1* h);
  void createHistos();
  void MirrorHistogram(TH1* histogram);
  std::string getHistoPrefix() const;

  TH1* h_cemc_channel_pedestal[128 * 192]{nullptr};
  TH1* h_ihcal_channel_pedestal[32 * 48]{nullptr};
  TH1* h_ohcal_channel_pedestal[32 * 48]{nullptr};

  TH1* h_cemc_channel_energy[128 * 192]{nullptr};
  TH1* h_ihcal_channel_energy[32 * 48]{nullptr};
  TH1* h_ohcal_channel_energy[32 * 48]{nullptr};

  int _eventcounter{0};
  int _range{1};

  bool m_debug{false};

  std::string detector;
  std::string m_outputFileName;
  std::string OutputFileName;
};

#endif
