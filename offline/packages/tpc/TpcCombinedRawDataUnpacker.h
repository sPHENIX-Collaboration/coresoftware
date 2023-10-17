// Tell emacs that this is a C++ source
// //  -*- C++ -*-.
#ifndef TPC_COMBINEDRAWDATAUNPACKER_H
#define TPC_COMBINEDRAWDATAUNPACKER_H

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSet.h>

#include <fun4all/SubsysReco.h>

#include <ffarawobjects/TpcRawHit.h>
#include <ffarawobjects/TpcRawHitContainer.h>

#include <TFile.h>
#include <memory>
#include <string>

class PHCompositeNode;
class CDBTTree;
class CDBInterface;
class TNtuple;
//class TpcRawHit;
//class TpcRawHitContainer;

class TpcCombinedRawDataUnpacker : public SubsysReco
{
 public:
  TpcCombinedRawDataUnpacker(std::string const& name = "TpcCombinedRawDataUnpacker", std::string const& outF = "TpcCombinedRawDataUnpackerOutput.root");

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode *topNode) override;  

 protected:
  std::string outfile_name;
  
  CDBTTree *m_cdbttree = nullptr;
  CDBInterface *m_cdb = nullptr;

  int _ievent = 0;

 private:
  std::string m_TpcRawNodeName = "TPCRAWHIT";

  TNtuple *m_ntup = nullptr;
  TFile *m_file = nullptr;

  int mc_sectors[12] = { 5, 4, 3, 2, 1, 0, 11, 10, 9, 8, 7, 6};
  int FEE_map[26] = {4, 5, 0, 2, 1, 11, 9, 10, 8, 7, 6, 0, 1, 3, 7, 6, 5, 4, 3, 2, 0, 2, 1, 3, 5, 4};
  int FEE_R[26] = {2, 2, 1, 1, 1, 3, 3, 3, 3, 3, 3, 2, 2, 1, 2, 2, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3};
};

#endif  // TPC_COMBINEDRAWDATAUNPACKER_H
