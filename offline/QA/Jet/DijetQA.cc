//////////////////////////////////////////////////////////////////
//							    	//
//		Dijet QA module for JetQA		    	//
//		QA module that looks at dijet Ajj, xj		//
//								//
//								//
//								//
//		Author:Skadi 				  	//
//		First commit:  		13 Nov 24		//
//		Most recent update:	02 Dec 24		//
//		version:		v1.0			//
//////////////////////////////////////////////////////////////////
#include "DijetQA.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

//____________________________________________________________________________..
DijetQA::DijetQA(const std::string& name, const std::string& recojetname)
  : SubsysReco(name)
  , m_moduleName(name)
  , m_etaRange(-1.1, 1.1)
  , m_ptLeadRange(1, 100)
  , m_ptSubRange(1, 100)
  , m_nJet(-1)
  , m_nJetPair(-1)
  /*, m_centrality(-1)*/
  , m_zvtx(-1)
  /*, m_impactparam(-1)*/
  , m_Ajj(-1)
  , m_xj(-1)
  , m_ptl(-1)
  , m_ptsl(-1)
  , m_phil(-1)
  , m_phisl(-1)
  , m_dphil(-1)
  , m_dphi(-1)
  , m_etal(-1)
  , m_etasl(-1)
  , m_deltaeta(-1)
  , m_doTrgSelect(false)
  , m_trgToSelect(JetQADefs::GL1::MBDNSJet1)
  , m_recoJetName(recojetname)
{
  if(Verbosity() > 1 )
  {
    std::cout << "DijetQA::DijetQA(const std::string &name) Calling ctor" << std::endl;
  }
}

//____________________________________________________________________________..
DijetQA::~DijetQA()
{
  if (Verbosity() > 1)
  {
    std::cout << "DijetQA::~DijetQA() Calling dtor" << std::endl;
  }
}

//____________________________________________________________________________..
int DijetQA::Init(PHCompositeNode* /*topNode*/)
{
  //  std::cout << "DijetQA::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  delete m_analyzer;  // make cppcheck happy
  m_analyzer = new TriggerAnalyzer();
  m_manager = QAHistManagerDef::getHistoManager();  // get the histogram anager

	if(!m_manager){
		std::cerr<<PHWHERE <<": PANIC: couldn't grab histogram manager!" <<std::endl;
		assert(m_manager);
	}
	std::string smallModuleName = m_moduleName; //make sure name is lowercase
	std::transform(
	      smallModuleName.begin(),
	      smallModuleName.end(),
	      smallModuleName.begin(),
	      ::tolower);
	h_Ajj=new TH1F(boost::str(boost::format("h_%s_Ajj") % smallModuleName).c_str(), boost::str(boost::format("A_{jj} for identified jet pairs for %s; A_{jj}; N_{pairs}")% m_recoJetName).c_str(), 100, -0.005, 0.995);
	h_xj=new TH1F(boost::str(boost::format("h_%s_xj") % smallModuleName).c_str(), boost::str(boost::format("x_{j} for identified jet pairs for %s; x_{j}; N_{pairs}")% m_recoJetName).c_str(), 100, -0.005, 0.995);
	h_pt=new TH1F(boost::str(boost::format("h_%s_pt") % smallModuleName).c_str(), boost::str(boost::format("p_{T} for leading jets in identified pairs for %s; p_{T} [GeV/c]; N_{jet}")% m_recoJetName).c_str(), 70, -0.5, 69.5);
	h_dphi=new TH1F(boost::str(boost::format("h_%s_dphi")% smallModuleName).c_str(), boost::str(boost::format("|#Delta #varphi| for identified jet pairs for %s; |#Delta #phi|; N_{pairs}")% m_recoJetName).c_str(), 64, 0, 6.2831);
	h_Ajj_pt=new TH2F(boost::str(boost::format("h_%s_Ajj_pt")% smallModuleName).c_str(), boost::str(boost::format("A_{jj} as a function of leading jet $p_{T}$ for %s; p_{T}^{leading} [GeV/c]; A_{jj}; N_{pairs}")% m_recoJetName).c_str(), 70, -0.5, 69.5, 100, -0.005, 0.995);
	h_xj_pt=new TH2F(boost::str(boost::format("h_%s_xj_pt")% smallModuleName).c_str(), boost::str(boost::format("x_{j} as a function of leading jet $p_{T}$ for %s; p_{T}^{leading} [GeV]; x_{j}; N_{pairs}")% m_recoJetName).c_str(), 70, -0.5, 69.5, 100, -0.005, 0.995);
	h_dphi_pt=new TH2F(boost::str(boost::format("h_%s_dphi_pt")% smallModuleName).c_str(), boost::str(boost::format("|#Delta #varphi| of dijet pair as a function of leading jet p_{T} for %s; p_{T}^{leading} [GeV/c]; |#Delta #varphi|; N_{pairs}")% m_recoJetName).c_str(), 70, -0.5, 69.5, 64, 0, 6.2832);
	h_dphi_Ajj=new TH2F(boost::str(boost::format("h_%s_dphi_Ajj")% smallModuleName).c_str(), boost::str(boost::format("A_{jj} of dijet pair as a function of |#Delta #varphi| for %s; |#Delta #varphi|; A_{jj}; N_{pairs}")% m_recoJetName).c_str(), 64, 0, 6.2831, 100, -0.005, 0.995);
	h_Ajj_l=new TH1F(boost::str(boost::format("h_%s_Ajj_l")% smallModuleName).c_str(), boost::str(boost::format("A_{jj} for event leading jet pairs for %s; A_{jj}; N_{pairs}")% m_recoJetName).c_str(), 100, -0.005, 0.995);
	h_xj_l=new TH1F(boost::str(boost::format("h_%s_xj_l")% smallModuleName).c_str(), boost::str(boost::format("x_{j} for event leading jet pairs for %s; x_{j}; N_{pairs}")% m_recoJetName).c_str(), 100, -0.005, 0.995);
	h_pt_l=new TH1F(boost::str(boost::format("h_%s_pt_l")% smallModuleName).c_str(), boost::str(boost::format("p_{T} for leading jets in event leading pair for %s; p_{T} [GeV/c]; N_{jet}")% m_recoJetName).c_str(), 70, -0.5, 69.5);
	h_dphi_l=new TH1F(boost::str(boost::format("h_%s_dphi_l")% smallModuleName).c_str(), boost::str(boost::format("|#Delta #varphi| for leading jet pairs for %s; |#Delta #varphi|; N_{pairs}")% m_recoJetName).c_str(), 64, 0, 6.2831);
	h_Ajj_pt_l=new TH2F(boost::str(boost::format("h_%s_Ajj_pt_l")% smallModuleName).c_str(), boost::str(boost::format("A_{jj} of event leading dijet pair as a function of leading jet p_{T} for %s; p_{T}^{leading} [GeV/c]; A_{jj}; N_{pairs}")% m_recoJetName).c_str(), 70, -0.5, 69.5, 100, -0.005, 0.995);
	h_xj_pt_l=new TH2F(boost::str(boost::format("h_%s_xj_pt_l")% smallModuleName).c_str(), boost::str(boost::format("x_{j} of event leading dijet pair as a function of leading jet p_{T} for %s; p_{T}^{leading} [GeV/c]; x_{j}; N_{pairs}")% m_recoJetName).c_str(), 70, -0.5, 69.5, 100, -0.005, 0.995);
	h_dphi_pt_l=new TH2F(boost::str(boost::format("h_%s_dphi_pt_l")% smallModuleName).c_str(), boost::str(boost::format("|#Delta #varphi| of event leading dijet pair as a function of leading jet p_{T} for %s; p_{T}^{leading} [GeV/c]; |#Delta #varphi|; N_{pairs}")% m_recoJetName).c_str(), 70, -0.5, 69.5, 64, 0, 6.2831);
	h_dphi_Ajj_l=new TH2F(boost::str(boost::format("h_%s_dphi_Ajj_l")% smallModuleName).c_str(), boost::str(boost::format("A_{jj} of event leading dijet pair as a function of |#Delta #varphi| for %s; |#Delta #varphi|^{leading}; A_{jj}; N_{pairs}")% m_recoJetName).c_str(), 64, 0, 6.2831, 100, -0.005, 0.995);
	
	return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int DijetQA::InitRun(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 1)
  {
    std::cout << "DijetQA::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int DijetQA::process_event(PHCompositeNode* topNode)
{
  if(Verbosity() > 1) {
	std::cout << "DijetQA::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  }
  if (m_doTrgSelect)
  {
    m_analyzer->decodeTriggers(topNode);
    bool hasTrigger = JetQADefs::DidTriggerFire(m_trgToSelect, m_analyzer);
    if (!hasTrigger)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }
/* //removing the centrality for now, will add back in with a conditional flag for pp
  CentralityInfo* cent_node = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
  if (!cent_node)
  {
    if (Verbosity() > 1)
    {
      std::cout << "No centrality info found" << std::endl;
    }
    m_centrality = 1.;
    m_impactparam = 0.;
  }
  else
  {
    m_centrality = cent_node->get_centile(CentralityInfo::PROP::bimp);
    m_impactparam = cent_node->get_quantity(CentralityInfo::PROP::bimp);
  }*/
  GlobalVertex* vtx = nullptr;
  GlobalVertexMap* vtxmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vtxmap || vtxmap->empty())
  {
    if( Verbosity() > 1 ){
	 std::cout << "No vertex map found, assuming the vertex has z=0" << std::endl;
    }
    m_zvtx = 0;
  }
  else
  {
    vtx = vtxmap->begin()->second;
    m_zvtx = vtx->get_z();
  }
  JetContainer* jets = findNode::getClass<JetContainer>(topNode, m_recoJetName);
  if (!jets)
  {
    if (Verbosity() > 1)
    {
      std::cerr << "No Jet container found" << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }
  else
  {
    FindPairs(jets);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
//____________________________________________________________________________..
void DijetQA::FindPairs(JetContainer* jets)
{
  // find all pairs that are potenital dijets and measure the kinematics
  if (Verbosity() > 1)
  {
    std::cout << "JetKinematicCheck::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  }
  Jet* jet_leading = nullptr;
  float pt_leading = 0, pt1 = 0, pt2 = 0;
  m_nJet = jets->size();
  if( Verbosity() > 2 ){
	std::cout << "number of jets is" << m_nJet << std::endl;
  }
  std::vector<std::pair<Jet*, Jet*> > jet_pairs;
  bool set_leading = false;
  for (auto j1 : *jets)
  {
    // assert(j1);
    Jet *jet_pair1 = nullptr, *jet_pair2 = nullptr;
    if (j1->get_pt() < m_ptLeadRange.first || j1->get_eta() < m_etaRange.first || j1->get_pt() > m_ptLeadRange.second || j1->get_eta() > m_etaRange.second)
    {
      continue;  // cut on 1 GeV jets
    }
    if (j1->get_pt() > pt_leading)
    {
      set_leading = true;
      pt_leading = j1->get_pt();
      jet_leading = j1;
    }
    for (auto j2 : (*jets))
    {
      if (j2 < j1)
      {
        continue;
      }
      if (/*j2 == j1 ||*/ j2->get_pt() < m_ptSubRange.first || j2->get_eta() < m_etaRange.first || j2->get_pt() > m_ptSubRange.second || j2->get_eta() > m_etaRange.second)
      {
        continue;
      }
      if (std::abs(j2->get_phi() - j1->get_phi()) > PI - DeltaPhi)
      {
        if (j2->get_pt() > j1->get_pt())
        {
          jet_pair1 = j2;
          jet_pair2 = j1;
        }
        else
        {
          jet_pair1 = j1;
          jet_pair2 = j2;
        }
        jet_pairs.emplace_back(jet_pair1, jet_pair2);
        if (set_leading && jet_pair1 && jet_pair2)
        {
          pt1 = jet_pair1->get_pt();
          pt2 = jet_pair2->get_pt();
          m_Ajj = (pt1 - pt2) / (pt1 + pt2);
          m_xj = pt2 / pt1;
          m_ptl = pt1;
          m_ptsl = pt2;
          m_phil = jet_pair1->get_phi();
          m_phisl = jet_pair2->get_phi();
          m_dphil = m_phil - m_phisl;
          m_etal = jet_pair1->get_eta();
          m_etasl = jet_pair2->get_eta();
          m_deltaeta = m_etal - m_etasl;
          h_Ajj_l->Fill(m_Ajj);
          h_xj_l->Fill(m_xj);
          h_pt_l->Fill(m_ptl);
          h_dphi_l->Fill(m_dphil);
          h_Ajj_pt_l->Fill(m_ptl, m_Ajj);
          h_xj_pt_l->Fill(m_ptl, m_xj);
          h_dphi_pt_l->Fill(m_ptl, std::abs(m_dphi));
          h_dphi_Ajj_l->Fill(std::abs(m_dphi), m_Ajj);
          //	m_T->Fill();
        }
        set_leading = false;
      }
    }
  }
  if(Verbosity() > 2) {
	std::cout << "Finished search for pairs" << std::endl;
  }
  m_nJetPair = jet_pairs.size();
  float Ajj = 0., xj = 0.;
  if (jet_pairs.size() > 0)
  {
    for (auto js : jet_pairs)
    {
      Jet *jet_pair1 = js.first, *jet_pair2 = js.second;
      if (!jet_pair1 || !jet_pair2)
      {
        continue;
      }
      if (jet_pair1 && Verbosity() > 2)
      {
        std::cout << "jetpair 1 object has a pt of " << jet_pair1->get_pt() << std::endl;
      }
      pt1 = jet_pair1->get_pt();
      pt2 = jet_pair2->get_pt();
      if (pt1 < pt2)
      {
        auto j = jet_pair1;
        jet_pair1 = jet_pair2;
        jet_pair2 = j;
        pt1 = jet_pair1->get_pt();
        pt2 = jet_pair2->get_pt();
      }
      float dphi = jet_pair1->get_phi() - jet_pair2->get_phi();
      Ajj = (pt1 - pt2) / (pt1 + pt2);
      xj = pt2 / pt1;
      h_Ajj->Fill(Ajj);
      h_xj->Fill(xj);
      h_pt->Fill(pt1);
      h_dphi->Fill(std::abs(dphi));
      h_Ajj_pt->Fill(pt1, Ajj);
      h_xj_pt->Fill(pt1, xj);
      h_dphi_pt->Fill(pt1, std::abs(dphi));
      h_dphi_Ajj->Fill(std::abs(dphi), Ajj);
      if (Verbosity() > 2)
      {
        std::cout << "highest pt jet is " << jet_leading->get_pt() << " and highest pt in a pair is " << jet_pair1->get_pt() << std::endl;
      }
    }
  }
  else
  {
    if (Verbosity() > 2)
    {
      std::cout << "Did not find a pair of jets" << std::endl;
    }
  }
  return;
}

//____________________________________________________________________________..
int DijetQA::ResetEvent(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 1)
  {
    std::cout << "DijetQA::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int DijetQA::EndRun(const int runnumber)
{
  if (Verbosity() > 1)
  {
    std::cout << "DijetQA::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int DijetQA::End(PHCompositeNode* /*topNode*/)
{
  //  std::cout << "DijetQA::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  /*h_Ajj->SetStats(0);
  h_xj->SetStats(0);
  h_dphi->SetStats(0);
  h_Ajj_pt->SetStats(0);
  h_xj_pt->SetStats(0);
  h_dphi_pt->SetStats(0);
  h_dphi_Ajj->SetStats(0);
  h_Ajj_l->SetStats(0);
  h_xj_l->SetStats(0);
  h_dphi_l->SetStats(0);
  h_Ajj_pt_l->SetStats(0);
  h_xj_pt_l->SetStats(0);
  h_dphi_pt_l->SetStats(0);
  h_dphi_Ajj_l->SetStats(0);
  TLegend* l1=new TLegend(0.7, 0.9, 0.9, 1.0); //TLegend for the Ajj plots
  l1->SetFillStyle(0);
  l1->SetBorderSize(0);
  l1->SetTextSize(0.06f);
  l1->AddEntry((TObject*) nullptr, boost::str(boost::format("A_{jj} dijet pairs with pt_{l} #geq %d", m_ptLeadRange.first),"");
  h_Ajj->GetListOfFunctions()->Add(l1);*/

  m_manager->registerHisto(h_Ajj);
  m_manager->registerHisto(h_xj);
  // m_manager->registerHisto(h_pt); //this is turned off but can be turned on to diagnose issues
  m_manager->registerHisto(h_dphi);
  m_manager->registerHisto(h_Ajj_pt);
  m_manager->registerHisto(h_xj_pt);
  m_manager->registerHisto(h_dphi_pt);
  m_manager->registerHisto(h_dphi_Ajj);
  m_manager->registerHisto(h_Ajj_l);
  m_manager->registerHisto(h_xj_l);
  // m_manager->registerHisto(h_pt_l); //this is turned off but can be turned on to diagnose issues
  m_manager->registerHisto(h_dphi_l);
  m_manager->registerHisto(h_Ajj_pt_l);
  m_manager->registerHisto(h_xj_pt_l);
  m_manager->registerHisto(h_dphi_pt_l);
  m_manager->registerHisto(h_dphi_Ajj_l);
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int DijetQA::Reset(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 1)
  {
    std::cout << "DijetQA::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void DijetQA::Print(const std::string& what) const
{
 if(Verbosity() > 1){
	 std::cout << "DijetQA::Print(const std::string &what) const Printing info for " << what << std::endl;
	}
}
