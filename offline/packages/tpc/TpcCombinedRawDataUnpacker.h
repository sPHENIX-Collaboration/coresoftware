// Tell emacs that this is a C++ source
// //  -*- C++ -*-.
#ifndef TPC_COMBINEDRAWDATAUNPACKER_H
#define TPC_COMBINEDRAWDATAUNPACKER_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class CDBTTree;
class CDBInterface;
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
  void skipNevent(int b) { startevt = b; }
  void event_range(int a, int b)
  {
    startevt = a;
    endevt = b;
  }

 private:
  TNtuple *m_ntup{nullptr};
  TFile *m_file{nullptr};
  CDBTTree *m_cdbttree{nullptr};
  CDBInterface *m_cdb{nullptr};

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

  std::string m_TpcRawNodeName{"TPCRAWHIT"};
  std::string outfile_name;
};

#endif  // TPC_COMBINEDRAWDATAUNPACKER_H
