// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef DIJETQA_H
#define DIJETQA_H

#include <fun4all/Fun4AllBase.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>
#include <jetqa/JetQADefs.h>
#include <qautils/QAHistManagerDef.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <jetbase/JetContainer.h>
#include <jetbase/JetMap.h>
#include <jetbase/Jetv1.h>
#include <jetbase/Jetv2.h>

#include <centrality/CentralityInfo.h>
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <math.h>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
#include <boost/format.hpp>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#define PI 3.1415926535
class PHCompositeNode;

class DijetQA : public SubsysReco
{
 public:
  DijetQA(const std::string &name = "DijetQA", const std::string &recojetname="AntiKt_Tower_r04");

  ~DijetQA() override;

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;
  void FindPairs(JetContainer *);

  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode *topNode) override;

  /// Called at the end of each run.
  int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;
  //////////////////////////////////////////////////////////////
  //							    //
  //        X_j = (p_(T, 1))/(p_(T,2))			    //
  //        A_j = (p_(T, 1) -p_(T,2))/P_T			    //
  //							    //
  //////////////////////////////////////////////////////////////

 private:
  Fun4AllHistoManager *m_manager{nullptr};
  std::string m_moduleName = "";
  std::pair<float, float> m_etaRange, m_ptRange;
  float DeltaPhiOne = 3.141529694 / 32.;  // cut on the opening angle of phi for the identified jets
                                          // Should set to integer multilple of hcal phi tower size ->Pi/32
  int ntowers_opening = 2;
  float DeltaPhi = ntowers_opening * DeltaPhiOne;
  int m_nJet = 0, m_nJetPair = 0;
  float /* m_centrality = 0.,*/ m_zvtx = 0.,/* m_impactparam = 0., */m_Ajj = 0., m_xj = 0., m_ptl = 0., m_ptsl = 0.;
  float m_phil = 0., m_phisl = 0., m_dphil = 0., m_dphi = 0., m_etal = 0., m_etasl = 0., m_deltaeta = 0.;
  TH1F *h_Ajj = nullptr, *h_xj = nullptr, *h_pt = nullptr, *h_dphi = nullptr;
  TH2F *h_Ajj_pt = nullptr, *h_xj_pt = nullptr, *h_dphi_pt = nullptr, *h_dphi_Ajj = nullptr;
  TH1F *h_Ajj_l = nullptr, *h_xj_l = nullptr, *h_pt_l = nullptr, *h_dphi_l = nullptr;
  TH2F *h_Ajj_pt_l = nullptr, *h_xj_pt_l = nullptr, *h_dphi_pt_l = nullptr, *h_dphi_Ajj_l = nullptr;
  bool m_doTrgSelect;
  uint32_t m_trgToSelect;
  std::string m_recoJetName = "AntiKT_Truth_r04";
};

#endif  // DIJETQA_H
