#ifndef GLOBALQA_GLOBALQA_H
#define GLOBALQA_GLOBALQA_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector> 

// Forward declarations
class PHCompositeNode;
class TH1;
class TH2;
class TProfile;

class GlobalQA : public SubsysReco {
  public:
    //! constructor
    GlobalQA(const std::string &name =
        "GlobalQA"); // const std::string &filename = "testQA.root");
    // //int nevents = 100);

    //! destructor
    ~GlobalQA() override;

    //! full initialization
    int Init(PHCompositeNode *) override;

    //! event processing method
    int process_event(PHCompositeNode *) override;

    //! end of run method
    int End(PHCompositeNode *) override;

    int process_g4hits(PHCompositeNode *);
    int process_g4cells(PHCompositeNode *);
    int process_towers(PHCompositeNode *);
    int process_clusters(PHCompositeNode *);

    void Detector(const std::string &name) { detector = name; }
    void set_timing_cut_width(const int &t) { _range = t; }

    void set_debug(bool debug) { m_debug = debug; }
    TH2 *LogYHist2D(const std::string &name, const std::string &title, int,
        double, double, int, double, double);

  private:
    int evtcount = 0;
    int Getpeaktime(TH1 *h);
    void createHistos();
    float TSAMPLE = 1.0 / (9.4e+6 * 6);

    // MBD histos
    TH1 *h_GlobalQA_mbd_zvtx = nullptr;
    TH1 *h_GlobalQA_mbd_zvtx_wide = nullptr;
    TH1 *h_GlobalQA_calc_zvtx = nullptr;
    TH1 *h_GlobalQA_calc_zvtx_wide = nullptr;
    TH1 *h_GlobalQA_mbd_charge_s = nullptr;
    TH1 *h_GlobalQA_mbd_charge_n = nullptr;

    TH1 *h_GlobalQA_mbd_charge_sum = nullptr;     //MBD north + south charge sum
    TH2 *h2_GlobalQA_mbd_charge_NS_correlation = nullptr;    // North charge vs south charge

    TH1 *h_GlobalQA_mbd_nhit_s = nullptr;
    TH1 *h_GlobalQA_mbd_nhit_n = nullptr;

    TH2 *h2_GlobalQA_mbd_nhits_NS_correlation = nullptr; //Number of event that have hits on South & North   

    TH1 *h_GlobalQA_mbd_zvtxq = nullptr;




    // ZDC histos
    TH1 *h_GlobalQA_zdc_zvtx = nullptr;
    TH1 *h_GlobalQA_zdc_zvtx_wide = nullptr;
    TH1 *h_GlobalQA_zdc_energy_s = nullptr;
    TH1 *h_GlobalQA_zdc_energy_n = nullptr;
    TH1* h_GlobalQA_triggerVec{nullptr};

    int _eventcounter{0}; // num events processed

    int _range{1};

    bool m_debug{false};

    std::string detector = "";
    std::string m_outputFileName = "";
    std::string OutputFileName = "";
};

#endif

