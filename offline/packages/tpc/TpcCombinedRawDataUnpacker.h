// Tell emacs that this is a C++ source
// //  -*- C++ -*-.
#ifndef TPC_COMBINEDRAWDATAUNPACKER_H
#define TPC_COMBINEDRAWDATAUNPACKER_H

#include <fun4all/SubsysReco.h>

#include <limits>
#include <map>
#include <string>
#include <vector>

class PHCompositeNode;
class CDBTTree;
class CDBInterface;
class TH2I;
class TFile;
class TNtuple;

class TpcCombinedRawDataUnpacker : public SubsysReco
{
 public:
  TpcCombinedRawDataUnpacker(std::string const &name = "TpcCombinedRawDataUnpacker", std::string const &outF = "TpcCombinedRawDataUnpackerOutput.root");

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *) override;
  int process_event(PHCompositeNode *) override;
  int End(PHCompositeNode *topNode) override;
  void writeTree() { m_writeTree = true; }
  void do_zero_suppression(bool b) { m_do_zerosup = b; }
  void set_pedestalSigmaCut(float b) { m_ped_sig_cut = b; }
  void do_noise_rejection(bool b) { m_do_noise_rejection = b; }
  void doBaselineCorr(bool val) { m_do_baseline_corr = val; }
  void doZSEmulation(bool val) { m_do_zs_emulation = val; }
  void ReadZeroSuppressedData() { 
    m_do_zs_emulation = true;
    m_zs_threshold = 10;
  }
  void set_presampleShift(int b) { m_presampleShift = b; }
  void set_zs_threshold(int b) { m_zs_threshold = b; }
  void skipNevent(int b) { startevt = b; }
  void useRawHitNodeName(const std::string &name) { m_TpcRawNodeName = name; }

  void event_range(int a, int b)
  {
    startevt = a;
    endevt = b;
  }
  struct chan_info
  {
    unsigned int fee = std::numeric_limits<unsigned int>::max();
    float ped = -1;
    float width = -1;
  };
  unsigned int get_rx(unsigned int layer)
  {
    return (layer - 7) / 16;
  }

  unsigned int create_fee_key(unsigned int side, unsigned int sector, unsigned int rx, unsigned int fee)
  {
    unsigned int key = rx;
    key += side * 10;
    key += sector * 100;
    key += fee * 10000;
    return key;
  }
  void unpack_fee_key(unsigned int &side, unsigned int &sector, unsigned int &rx, unsigned int &fee, unsigned int fee_key)
  {
    rx = fee_key % 10;
    side = (fee_key % 100) / 10;
    sector = ((fee_key - side * 10 - rx) / 100) / 100;
    fee = fee_key / 10000;
    return;
  }
  unsigned int create_pad_key(unsigned int side, unsigned int layer, unsigned int padnum)
  {
    unsigned int key = side;
    key += layer * 10;
    key += padnum * 1000;
    return key;
  }
  void unpack_pad_key(unsigned int &side, unsigned int &layer, unsigned int &pad_num, unsigned int pad_key)
  {
    side = pad_key % 10;
    layer = (pad_key % 1000) / 10;
    pad_num = (pad_key - layer * 10 - side) / 1000;
    return;
  }

 private:
  TNtuple *m_ntup{nullptr};
  TNtuple *m_ntup_hits = nullptr;
  TNtuple *m_ntup_hits_corr = nullptr;
  TFile *m_file{nullptr};
  CDBTTree *m_cdbttree{nullptr};
  CDBInterface *m_cdb{nullptr};

  int m_presampleShift = 40;  // number of presamples shifted to line up t0
  int _ievent{0};
  int startevt{-1};
  int endevt{9999999};
  int mc_sectors[12]{5, 4, 3, 2, 1, 0, 11, 10, 9, 8, 7, 6};
  int FEE_map[26]{4, 5, 0, 2, 1, 11, 9, 10, 8, 7, 6, 0, 1, 3, 7, 6, 5, 4, 3, 2, 0, 2, 1, 3, 5, 4};
  int FEE_R[26]{2, 2, 1, 1, 1, 3, 3, 3, 3, 3, 3, 2, 2, 1, 2, 2, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3};

  float m_ped_sig_cut{4.0};

  bool m_writeTree{false};
  bool m_do_zerosup{true};
  bool m_do_noise_rejection{true};
  bool m_do_baseline_corr{false};
  bool m_do_zs_emulation{false};
  int pedestal_offset{30};
  int m_zs_threshold{30};
  std::string m_TpcRawNodeName{"TPCRAWHIT"};
  std::string outfile_name;
  std::map<unsigned int, chan_info> chan_map;                  // stays in place
  std::map<unsigned int, TH2I *> feeadc_map;                   // histos reset after each event
  std::map<unsigned int, std::vector<float>> feebaseline_map;  // cleared after each event
};

#endif  // TPC_COMBINEDRAWDATAUNPACKER_H
