#include "TpcLaminationFitting.h"

#include <trackbase/CMFlashDifferenceContainerv1.h>
#include <trackbase/CMFlashDifferencev1.h>
#include <trackbase/LaserClusterContainerv1.h>
#include <trackbase/LaserClusterv1.h>
#include <trackbase/TpcDefs.h>

#include <ffaobjects/EventHeader.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2.h>
#include <TH3.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include <TVector3.h>

#include <boost/format.hpp>

#include <cmath>
#include <iomanip>
#include <set>
#include <string>
#include <vector>

//________________________________
TpcLaminationFitting::TpcLaminationFitting(const std::string &name)
  : SubsysReco(name)
{
}

//___________________________________________________________
void TpcLaminationFitting::set_grid_dimensions(int phibins, int rbins)
{
  m_phibins = phibins;
  m_rbins = rbins;
}

//____________________________________________
int TpcLaminationFitting::InitRun(PHCompositeNode *topNode)
{
  for (int s = 0; s < 2; s++)
  {
    for (int l = 0; l < 18; l++)
    {
      double shift = l * M_PI / 9;
      if (s == 0)
      {
        shift += M_PI / 18;
      }

      m_hLamination[l][s] = new TH2D((boost::format("hLamination%d_%s") %l %(s == 1 ? "North" : "South")).str().c_str(), (boost::format("Lamination %d %s, #phi_{expected}=%.2f;R [cm];#phi") %l %(s == 1 ? "North" : "South") %shift).str().c_str(), 200, 30, 80, 200, shift - 0.2, shift + 0.2);
      //m_fLamination[l][s] = new TF1((boost::format("fLamination%d_%s") %l %(s == 1 ? "North" : "South")).str().c_str(), "[0]+[1]*exp(-[2]*x)", 30, 80);
      //m_fLamination[l][s] = new TF1((boost::format("fLamination%d_%s") %l %(s == 1 ? "North" : "South")).str().c_str(), "[0]*(1+exp(-[2]*(x-[1])))", 30, 80);
      m_fLamination[l][s] = new TF1((boost::format("fLamination%d_%s") %l %(s == 1 ? "North" : "South")).str().c_str(), "[3]+[0]*(1-exp(-[2]*(x-[1])))", 30, 80);
      //m_fLamination[l][s]->SetParameters(-0.022 + shift, log(3.0/(-0.22 + shift)), 0.12);
      m_fLamination[l][s]->SetParameters(-0.011, 30, 0.16, 0.0);
      m_fLamination[l][s]->SetParLimits(0, -0.22, 0.0);
      m_fLamination[l][s]->SetParLimits(1, 0, 80);
      m_fLamination[l][s]->SetParLimits(2, 0.0, 3);
      m_fLamination[l][s]->FixParameter(3, shift);
      m_laminationCenter[l][s] = shift;
    }
  }

  //Make map for run and ZDC rate for pp mode
  m_run_ZDC_map_pp.insert(std::pair<int, float>(49709, 555.0));
  m_run_ZDC_map_pp.insert(std::pair<int, float>(52077, 0.0));
  m_run_ZDC_map_pp.insert(std::pair<int, float>(52078, 0.0));
  m_run_ZDC_map_pp.insert(std::pair<int, float>(53534, 3013.5));
  m_run_ZDC_map_pp.insert(std::pair<int, float>(53630, 6849.3));
  m_run_ZDC_map_pp.insert(std::pair<int, float>(53631, 5577.8));
  m_run_ZDC_map_pp.insert(std::pair<int, float>(53632, 5151.2));
  m_run_ZDC_map_pp.insert(std::pair<int, float>(53652, 4600.0));
  m_run_ZDC_map_pp.insert(std::pair<int, float>(53687, 3967.2));
  m_run_ZDC_map_pp.insert(std::pair<int, float>(53716, 3070.1));
  m_run_ZDC_map_pp.insert(std::pair<int, float>(53738, 4510.7));
  m_run_ZDC_map_pp.insert(std::pair<int, float>(53739, 4165.0));
  m_run_ZDC_map_pp.insert(std::pair<int, float>(53741, 3738.1));
  m_run_ZDC_map_pp.insert(std::pair<int, float>(53742, 3721.4));
  m_run_ZDC_map_pp.insert(std::pair<int, float>(53743, 3693.4));
  m_run_ZDC_map_pp.insert(std::pair<int, float>(53744, 3581.9));
  m_run_ZDC_map_pp.insert(std::pair<int, float>(53756, 4471.4));
  m_run_ZDC_map_pp.insert(std::pair<int, float>(53783, 4825.7));
  m_run_ZDC_map_pp.insert(std::pair<int, float>(53871, 6871.5));
  m_run_ZDC_map_pp.insert(std::pair<int, float>(53876, 5082.3));
  m_run_ZDC_map_pp.insert(std::pair<int, float>(53877, 4758.5));
  m_run_ZDC_map_pp.insert(std::pair<int, float>(53879, 4315.0));

  //beam off go into pp
  m_run_ZDC_map_pp.insert(std::pair<int, float>(53098, 0.0));
  m_run_ZDC_map_pp.insert(std::pair<int, float>(53271, 0.0));

  m_run_ZDC_map_auau.insert(std::pair<int, float>(54966, 12400.));
  m_run_ZDC_map_auau.insert(std::pair<int, float>(54967, 11600.));
  m_run_ZDC_map_auau.insert(std::pair<int, float>(54968, 10500.));
  m_run_ZDC_map_auau.insert(std::pair<int, float>(54969, 9680.));
  
  int ret = GetNodes(topNode);
  return ret;
}

//______________________________________
int TpcLaminationFitting::GetNodes(PHCompositeNode *topNode)
{
  m_correctedCMcluster_map = findNode::getClass<LaserClusterContainer>(topNode, "LAMINATION_CLUSTER");
  if (!m_correctedCMcluster_map)
  {
    std::cout << PHWHERE << "CORRECTED_CM_CLISTER Node missing, abort." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_dcc_in_module_edge = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerModuleEdge");
  if (m_dcc_in_module_edge)
  {
    std::cout << "TpcLaminationFitting::GetNodes - found TPC distortion correction container module edge" << std::endl;
  }

  m_dcc_in_static = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerStatic");
  if (m_dcc_in_static)
  {
    std::cout << "TpcLaminationFitting::GetNodes - found TPC distortion correction container static" << std::endl;
  }

  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  auto flashDiffContainer = findNode::getClass<CMFlashDifferenceContainerv1>(topNode, "CM_FLASH_DIFFERENCES");
  if (!flashDiffContainer)
  {
    PHNodeIterator dstIter(dstNode);
    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    flashDiffContainer = new CMFlashDifferenceContainerv1;
    PHIODataNode<PHObject> *CMFlashDifferenceNode = new PHIODataNode<PHObject>(flashDiffContainer, "CM_FLASH_DIFFERENCES", "PHObject");
    DetNode->addNode(CMFlashDifferenceNode);
  }

  m_cm_flash_diffs = findNode::getClass<CMFlashDifferenceContainerv1>(topNode, "CM_FLASH_DIFFERENCES");
  if (!m_cm_flash_diffs)
  {
    std::cout << PHWHERE << " ERROR: Can't find CM_FLASH_DIFFERENCES." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  const std::string dcc_out_node_name = "TpcDistortionCorrectionContainerAverage";
  m_dcc_out = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, dcc_out_node_name);
  if (!m_dcc_out)
  {
    auto runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
    if (!runNode)
    {
      std::cout << "TpcLaminationFitting::InitRun = RUN Node missing, quitting" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    std::cout << "TpcLaminationFitting::GetNodes - creating TpcDistortionCorrectionContainer in node " << dcc_out_node_name << std::endl;
    m_dcc_out = new TpcDistortionCorrectionContainer;
    auto node = new PHDataNode<TpcDistortionCorrectionContainer>(m_dcc_out, dcc_out_node_name);
    runNode->addNode(node);
  }

  const float phiMin = m_phiMin - (m_phiMax - m_phiMin) / m_phibins;
  const float phiMax = m_phiMax + (m_phiMax - m_phiMin) / m_phibins;

  const float rMin = m_rMin - (m_rMax - m_rMin) / m_rbins;
  const float rMax = m_rMax + (m_rMax - m_rMin) / m_rbins;

  const std::array<const std::string, 2> extension = {{"_negz", "_posz"}};

  m_dcc_out->m_dimensions = 2;

  for (int i = 0; i < 2; ++i)
  {
    delete m_dcc_out->m_hDPint[i];
    m_dcc_out->m_hDPint[i] = new TH2F((boost::format("hIntDistortionP%s") % extension[i]).str().c_str(), (boost::format("hIntDistortionP%s") % extension[i]).str().c_str(), m_phibins + 2, phiMin, phiMax, m_rbins + 2, rMin, rMax);
    delete m_dcc_out->m_hDRint[i];
    m_dcc_out->m_hDRint[i] = new TH2F((boost::format("hIntDistortionR%s") % extension[i]).str().c_str(), (boost::format("hIntDistortionR%s") % extension[i]).str().c_str(), m_phibins + 2, phiMin, phiMax, m_rbins + 2, rMin, rMax);
    delete m_dcc_out->m_hDZint[i];
    m_dcc_out->m_hDZint[i] = new TH2F((boost::format("hIntDistortionZ%s") % extension[i]).str().c_str(), (boost::format("hIntDistortionZ%s") % extension[i]).str().c_str(), m_phibins + 2, phiMin, phiMax, m_rbins + 2, rMin, rMax);
    delete m_dcc_out->m_hentries[i];
    m_dcc_out->m_hentries[i] = new TH2I((boost::format("hEntries%s") % extension[i]).str().c_str(), (boost::format("hEntries%s") % extension[i]).str().c_str(), m_phibins + 2, phiMin, phiMax, m_rbins + 2, rMin, rMax);
  }

  m_laminationTree = new TTree("laminationTree","laminationTree");
  m_laminationTree->Branch("side",&m_side);
  m_laminationTree->Branch("lamIndex",&m_lamIndex);
  m_laminationTree->Branch("lamPhi",&m_lamPhi);
  m_laminationTree->Branch("goodFit",&m_goodFit);
  m_laminationTree->Branch("A",&m_A);
  m_laminationTree->Branch("B",&m_B);
  m_laminationTree->Branch("C",&m_C);
  m_laminationTree->Branch("A_err",&m_A_err);
  m_laminationTree->Branch("B_err",&m_B_err);
  m_laminationTree->Branch("C_err",&m_C_err);
  m_laminationTree->Branch("distanceToFit",&m_dist);
  m_laminationTree->Branch("nBinsFit",&m_nBins);


  return Fun4AllReturnCodes::EVENT_OK;
}

//__________________________________________
int TpcLaminationFitting::process_event(PHCompositeNode *topNode)
{
  eventHeader = findNode::getClass<EventHeader>(topNode, "EventHeader");
  if (!eventHeader)
  {
    std::cout << PHWHERE << " EventHeader Node missing, abort" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_runnumber = eventHeader->get_RunNumber();
  
  if (m_useHeader && eventHeader->get_EvtSequence() == 0)
  {
    m_useHeader = false;
  }

  if (m_useHeader)
  {
    m_event_index = eventHeader->get_EvtSequence();
  }

  if (!m_correctedCMcluster_map || m_correctedCMcluster_map->size() < 1000)
  {
    if (!m_useHeader)
    {
      m_event_index++;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_nClusters += m_correctedCMcluster_map->size();
  m_nEvents++;

  auto clusrange = m_correctedCMcluster_map->getClusters();
  for (auto cmitr = clusrange.first; cmitr != clusrange.second; ++cmitr)
  {
    const auto &[cmkey, cmclus_orig] = *cmitr;
    LaserCluster *cmclus = cmclus_orig;
    // const unsigned int adc = cmclus->getAdc();
    bool side = (bool) TpcDefs::getSide(cmkey);
    if (cmclus->getNLayers() <= m_nLayerCut)
    {
      continue;
    }

    //Acts::Vector3 pos(cmclus->getX(), cmclus->getY(), cmclus->getZ());
    Acts::Vector3 pos(cmclus->getX(), cmclus->getY(), (side ? 1.0 : -1.0));
    if (m_dcc_in_module_edge)
    {
      pos = m_distortionCorrection.get_corrected_position(pos, m_dcc_in_module_edge);
    }
    if (m_dcc_in_static)
    {
      pos = m_distortionCorrection.get_corrected_position(pos, m_dcc_in_static);
    }

    TVector3 tmp_pos(pos[0], pos[1], pos[2]);

    for (int l = 0; l < 18; l++)
    {
      double shift = m_laminationCenter[l][side];

      double phi2pi = tmp_pos.Phi();
      if (side && phi2pi < -0.2)
      {
        phi2pi += 2 * M_PI;
      }
      if (!side && phi2pi < M_PI / 18 - 0.2)
      {
        phi2pi += 2 * M_PI;
      }

      if (phi2pi > shift - 0.2 && phi2pi < shift + 0.2)
      {
        m_hLamination[l][side]->Fill(tmp_pos.Perp(), phi2pi);
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________
int TpcLaminationFitting::fitLaminations()
{
  //int nBinAvg = 6;
  int nBinAvg = 4;
  // double contentCut = 0.0;

  /*
  for(int s=0; s<2; s++)
  {
          for(int l=0; l<18; l++)
          {
                  if(m_hLamination[l][s]->GetMaximum()/5.0 > contentCut)
                  {
                          contentCut = m_hLamination[l][s]->GetMaximum()/5.0;
                  }
          }
  }
  */
  //double seedScale = (m_nClusters / m_nEvents) / 3718.8030;

  float ZDC = 4500.0;
  TF1 *Af[2] = {new TF1("AN","pol1",0,100000), new TF1("AS","pol1",0,100000)};
  TF1 *Bf[2] = {new TF1("BN","pol1",0,100000), new TF1("BS","pol1",0,100000)};
  double Cseed[2] = {0.16, 0.125};

  if(ppMode)
  {
    auto it = m_run_ZDC_map_pp.find(m_runnumber);
    if( it != m_run_ZDC_map_pp.end() )
    {
      std::cout << "pp runnumber " << m_runnumber << " found. It has ZDC NS rate of " << it->second << std::endl;
      ZDC = it->second;
    }
    else
    {
      std::cout << "pp runnumber " << m_runnumber << " not found. Using default value of " << ZDC << std::endl;
    }

    Af[0]->SetParameters(-0.007999,-1.783e-6);
    Af[1]->SetParameters(-0.003288,-2.297e-6);
    
    Bf[0]->SetParameters(31.55,0.0006141);
    Bf[1]->SetParameters(34.7,0.0005226);
  }
  else
  {
    ZDC = 10000.0;
    
    auto it = m_run_ZDC_map_auau.find(m_runnumber);
    if( it != m_run_ZDC_map_auau.end() )
    {
      std::cout << "AuAu runnumber " << m_runnumber << " found. It has ZDC NS rate of " << it->second << std::endl;
      ZDC = it->second;
    }
    else
    {
      std::cout << "AuAu runnumber " << m_runnumber << " not found. Using default value of " << ZDC << std::endl;
    }

    Af[0]->SetParameters(-0.003836,-1.025e-6);
    Af[1]->SetParameters(-0.003283,-8.176e-7);
    
    Bf[0]->SetParameters(32.96,0.0002997);
    Bf[1]->SetParameters(31.19,0.0005622);

    Cseed[0] = 0.125;
    Cseed[1] = 0.122;
  }
  

  
  for (int s = 0; s < 2; s++)
  {
    for (int l = 0; l < 18; l++)
    {
      if (m_hLamination[l][s]->GetEntries() < 100)
      {
        m_laminationGoodFit[l][s] = false;
        continue;
      }

      // double histMax = m_hLamination[l][s]->GetMaximum();

      TGraph *gr = new TGraph();
      TGraph *proj = new TGraph();

      //m_fLamination[l][s]->SetParameters(-0.022 + m_laminationCenter[l][s], 4.595 * seedScale, 0.138);
      //m_fLamination[l][s]->SetParameters(-0.022 + m_laminationCenter[l][s], log(4.595 * seedScale/(-0.022 + m_laminationCenter[l][s])), 0.138);
      //m_fLamination[l][s]->SetParameters(-0.011 + m_laminationCenter[l][s], 0.025, 0.16);
      //m_fLamination[l][s]->SetParameters(-0.011, 30, 0.16, m_laminationCenter[l][s]);
      m_fLamination[l][s]->SetParameters(Af[s]->Eval(ZDC), Bf[s]->Eval(ZDC), Cseed[s], m_laminationCenter[l][s]);
      m_fLamination[l][s]->FixParameter(3, m_laminationCenter[l][s]);
      
      TF1 *fitSeed = (TF1 *) m_fLamination[l][s]->Clone();
      fitSeed->SetName((boost::format("fitSeed%d_%s") %l %(s == 1 ? "North" : "South")).str().c_str());

      TF1 *localFit = (TF1 *) m_fLamination[l][s]->Clone();
      localFit->SetName((boost::format("localFit%d_%s") %l %(s == 1 ? "North" : "South")).str().c_str());

      for (int i = m_hLamination[l][s]->GetNbinsX(); i >= 1; i--)
      {
        double R = m_hLamination[l][s]->GetXaxis()->GetBinCenter(i);
        double sum = 0.0;
        double c = 0.0;
        double noFitPred = fitSeed->Eval(R);
        int noFitPredBin = m_hLamination[l][s]->GetYaxis()->FindBin(noFitPred);

        double pred2 = -999;
        int pred2Bin = -999;

        if (gr->GetN() >= 5)
        {
          gr->Fit(localFit);
          pred2 = localFit->Eval(R);
          pred2Bin = m_hLamination[l][s]->GetYaxis()->FindBin(pred2);
        }

        if (pred2 == -999)
        {
          proj->AddPoint(R, noFitPred);
        }
        else
        {
          proj->AddPoint(R, pred2);
        }

        for (int j = 1; j <= m_hLamination[l][s]->GetNbinsY(); j++)
        {
          double phi = m_hLamination[l][s]->GetYaxis()->GetBinCenter(j);

          if (fabs(phi - m_laminationCenter[l][s]) > 0.05)
          {
            continue;
          }
          if (pred2 == -999 && std::abs(j - noFitPredBin) > nBinAvg)
          {
            continue;
          }
          if (pred2 != -999 && std::abs(j - noFitPredBin) > nBinAvg && std::abs(j - pred2Bin) > nBinAvg)
          {
            continue;
          }

          double content = m_hLamination[l][s]->GetBinContent(i, j);
          // if(content < histMax/20.0 || content > histMax/10.0)
          if (content <= 0)
          {
            continue;
          }

          sum += content * phi;
          c += content;
        }

        if (c > 0)
        {
          gr->AddPoint(R, sum / c);
        }
      }

      gr->Fit(m_fLamination[l][s]);

      // m_hLamination[l][s]->GetXaxis()->SetRangeUser(contentCut, 2*m_hLamination[l][s]->GetMaximum());
      // m_hLamination[l][s]->GetXaxis()->SetRangeUser(histMax/20.0,histMax/10.0);

      double distToFunc = 0.0;
      int nBinsUsed = 0;

      for (int i = 1; i <= m_hLamination[l][s]->GetNbinsX(); i++)
      {
        double R = m_hLamination[l][s]->GetXaxis()->GetBinCenter(i);
        double funcVal = m_fLamination[l][s]->Eval(R);
        int funcBin = m_hLamination[l][s]->GetYaxis()->FindBin(funcVal);
        for (int j = 0; j <= nBinAvg; j++)
        {
          double contentp = -999;
          double contentm = -999;

          if (funcBin + j < m_hLamination[l][s]->GetNbinsY())
          {
            contentp = m_hLamination[l][s]->GetBinContent(i, funcBin + j);
          }
          if (funcBin - j >= 1)
          {
            contentm = m_hLamination[l][s]->GetBinContent(i, funcBin - j);
          }
          // if(contentp >= contentCut || contentm >= contentCut)
          // if(contentp >= histMax/20.0 || contentm >= histMax/20.0)
          if (contentp >= 1 || contentm >= 1)
          {
            distToFunc += j;
            nBinsUsed++;
            break;
          }
        }
      }

      m_distanceToFit[l][s] = distToFunc / nBinsUsed;
      m_nBinsFit[l][s] = nBinsUsed;
      if (nBinsUsed < 10 || distToFunc / nBinsUsed > 1.0)
      {
        m_laminationGoodFit[l][s] = false;
      }
      else
      {
        m_laminationGoodFit[l][s] = true;
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcLaminationFitting::InterpolatePhiDistortions(TH2 *simPhiDistortion[2])
{
  phiDistortionLamination[0] = (TH2 *) simPhiDistortion[0]->Clone();
  phiDistortionLamination[0]->Reset();
  phiDistortionLamination[0]->SetName("phiDistortionLamination0");
  phiDistortionLamination[1] = (TH2 *) simPhiDistortion[1]->Clone();
  phiDistortionLamination[1]->Reset();
  phiDistortionLamination[1]->SetName("phiDistortionLamination1");

  for (int s = 0; s < 2; s++)
  {
    for (int l = 0; l < 18; l++)
    {
      if (!m_laminationGoodFit[l][s])
      {
        continue;
      }

      for (int i = 1; i <= phiDistortionLamination[s]->GetNbinsY(); i++)
      {
        double R = phiDistortionLamination[s]->GetYaxis()->GetBinCenter(i);
        if (R < 30 || R > 80)
        {
          continue;
        }
        double phi = m_fLamination[l][s]->Eval(R);
        if (phi < 0.0)
        {
          phi += 2 * M_PI;
        }
        if (phi > 2 * M_PI)
        {
          phi -= 2 * M_PI;
        }
        int phiBin = phiDistortionLamination[s]->GetXaxis()->FindBin(phi);
        //m_fLamination[l][s]->SetParameter(0, m_fLamination[l][s]->GetParameter(0) - m_laminationCenter[l][s]);
	m_fLamination[l][s]->SetParameter(3, m_laminationOffset[l][s]);
	/*
	if(s==0)
	{
	  m_fLamination[l][s]->SetParameter(3, 0.0);
	}
        else
	{
	  m_fLamination[l][s]->SetParameter(3, 0.0);
	}
	*/
	//m_fLamination[l][s]->SetParameter(3, 0.0);
        double phiDistortion = R * m_fLamination[l][s]->Integral(phiDistortionLamination[s]->GetYaxis()->GetBinLowEdge(i), phiDistortionLamination[s]->GetYaxis()->GetBinLowEdge(i + 1)) / (phiDistortionLamination[s]->GetYaxis()->GetBinLowEdge(i + 1) - phiDistortionLamination[s]->GetYaxis()->GetBinLowEdge(i));
        //m_fLamination[l][s]->SetParameter(0, m_fLamination[l][s]->GetParameter(0) + m_laminationCenter[l][s]);
        m_fLamination[l][s]->SetParameter(3, m_laminationCenter[l][s]);
        phiDistortionLamination[s]->SetBinContent(phiBin, i, phiDistortion);
      }
    }
  }

  for (int s = 0; s < 2; s++)
  {
    m_dcc_out->m_hDPint[s] = (TH2 *) phiDistortionLamination[s]->Clone();
    m_dcc_out->m_hDPint[s]->SetName((boost::format("hIntDistortionP%s") %(s == 0 ? "_negz" : "_posz")).str().c_str());
  }

  /*
  for(int s=0; s<2; s++)
  {
          for(int i=1; i<=m_dcc_out->m_hDPint[s]->GetNbinsX(); i++)
          {
                  for(int j=1; j<=m_dcc_out->m_hDPint[s]->GetNbinsY(); j++)
                  {
                          if(phiDistortionLamination[s]->GetBinContent(i,j) != 0.0)
                          {
                                  m_dcc_out->m_hDPint[s]->SetBinContent(i,j, phiDistortionLamination[s]->GetBinContent(i,j));
                          }
                  }
          }
  }
  */

  // m_dcc_out->m_hDPint[0] = (TH2*)phiDistortionLamination[0]->Clone();
  // m_dcc_out->m_hDPint[1] = (TH2*)phiDistortionLamination[1]->Clone();

  for (int s = 0; s < 2; s++)
  {
    for (int i = 1; i <= m_dcc_out->m_hDPint[s]->GetNbinsY(); i++)
    {
      double R = m_dcc_out->m_hDPint[s]->GetYaxis()->GetBinCenter(i);
      if (R < 30 || R > 80)
      {
        continue;
      }
      std::vector<int> laminationPhiBins;
      for (int j = 1; j <= m_dcc_out->m_hDPint[s]->GetNbinsX(); j++)
      {
        if (m_dcc_out->m_hDPint[s]->GetBinContent(j, i) != 0.0)
        {
          laminationPhiBins.push_back(j);
        }
      }

      if (laminationPhiBins.size() > 1)
      {
        laminationPhiBins.push_back(laminationPhiBins[0]);
      }

      for (int lamPair = 0; lamPair < (int) laminationPhiBins.size() - 1; lamPair++)
      {
        double dist0 = m_dcc_out->m_hDPint[s]->GetBinContent(laminationPhiBins[lamPair], i);
        double dist1 = m_dcc_out->m_hDPint[s]->GetBinContent(laminationPhiBins[lamPair + 1], i);

        int nEmptyBins = laminationPhiBins[lamPair + 1] - laminationPhiBins[lamPair] - 1;
        if (laminationPhiBins[lamPair] > laminationPhiBins[lamPair + 1])
        {
          double binW = m_dcc_out->m_hDPint[s]->GetXaxis()->GetBinWidth(1);
          nEmptyBins = (int) round((m_dcc_out->m_hDPint[s]->GetXaxis()->GetBinCenter(laminationPhiBins[lamPair + 1]) + 2 * M_PI - m_dcc_out->m_hDPint[s]->GetXaxis()->GetBinCenter(laminationPhiBins[lamPair])) / binW) - 1;
        }

        bool wrap = false;
        int wrapBin = -1;
        for (int j = 1; j <= nEmptyBins; j++)
        // for(int j=+ 1; j<laminationPhiBins[lamPair+1]; j++)
        {
          if (!wrap && laminationPhiBins[lamPair] + j > m_dcc_out->m_hDPint[s]->GetNbinsX() - 1)
          {
            wrap = true;
            wrapBin = j;
          }
          double dist = dist0 + ((dist1 - dist0) * (1.0 * j / (nEmptyBins + 1)));
          if (wrap)
          {
            m_dcc_out->m_hDPint[s]->SetBinContent(2 + (j - wrapBin), i, dist);
          }
          else
          {
            m_dcc_out->m_hDPint[s]->SetBinContent(laminationPhiBins[lamPair] + j, i, dist);
          }
        }
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcLaminationFitting::End(PHCompositeNode * /*topNode*/)
{
  int fitSuccess = fitLaminations();
  if (fitSuccess != Fun4AllReturnCodes::EVENT_OK)
  {
    std::cout << PHWHERE << " Return code for lamination fitting was " << fitSuccess << " and not successful" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if(m_fitFileName != "")
  {
    TCanvas *c1 = new TCanvas();
    gStyle->SetPalette(56);
    c1->SetLogz();
    gStyle->SetOptStat(0);
    c1->SaveAs((boost::format("%s[") % m_fitFileName).str().c_str());
    for (int s = 0; s < 2; s++)
    {
      for (int l = 0; l < 18; l++)
      {
	c1->Clear();
	m_hLamination[l][s]->Draw("COLZ");
	m_fLamination[l][s]->SetLineColor(kRed);
	if (!m_laminationGoodFit[l][s])
	{
	  m_fLamination[l][s]->SetLineColor(kMagenta);
	}
	m_fLamination[l][s]->Draw("same");
	
	TLine *line = new TLine(30,m_laminationCenter[l][s],80,m_laminationCenter[l][s]);
	line->SetLineColor(kBlue);
	line->SetLineStyle(2);
	line->Draw("same");
	
	TPaveText *pars = new TPaveText(0.6, 0.55, 0.85, 0.85, "NDC");
	pars->AddText("#phi = #phi_{ideal} + A#times (1 - e^{-C#times (R - B)})");
	pars->AddText((boost::format("A=%.3f#pm %.3f") %m_fLamination[l][s]->GetParameter(0) %m_fLamination[l][s]->GetParError(0)).str().c_str());
	pars->AddText((boost::format("#phi_{ideal}=%.3f#pm %.3f") %m_fLamination[l][s]->GetParameter(3) %m_fLamination[l][s]->GetParError(3)).str().c_str());
	pars->AddText((boost::format("B=%.3f#pm %.3f") %m_fLamination[l][s]->GetParameter(1) %m_fLamination[l][s]->GetParError(1)).str().c_str());
	pars->AddText((boost::format("C=%.3f#pm %.3f") %m_fLamination[l][s]->GetParameter(2) %m_fLamination[l][s]->GetParError(2)).str().c_str());
	pars->AddText((boost::format("Distance to line=%.2f") %m_distanceToFit[l][s]).str().c_str());
	pars->AddText((boost::format("Number of Bins used=%d") %m_nBinsFit[l][s]).str().c_str());
	pars->Draw("same");
	c1->SaveAs(m_fitFileName.c_str());
      }
    }
    c1->SaveAs((boost::format("%s]") %m_fitFileName).str().c_str());
  }
  
  TFile *simDistortion = new TFile("/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/release/release_new/new.10/share/calibrations/distortion_maps/average_minus_static_distortion_inverted_10-new.root", "READ");
  TH3 *hIntDistortionP_posz = (TH3 *) simDistortion->Get("hIntDistortionP_posz");
  hIntDistortionP_posz->GetZaxis()->SetRange(2, 2);
  TH2 *simPhiDistortion[2];
  simPhiDistortion[1] = (TH2 *) hIntDistortionP_posz->Project3D("yx");
  TH3 *hIntDistortionP_negz = (TH3 *) simDistortion->Get("hIntDistortionP_negz");
  hIntDistortionP_negz->GetZaxis()->SetRange(hIntDistortionP_negz->GetNbinsZ() - 1, hIntDistortionP_negz->GetNbinsZ() - 1);
  simPhiDistortion[0] = (TH2 *) hIntDistortionP_negz->Project3D("yx");

  int interpolateSuccess = InterpolatePhiDistortions(simPhiDistortion);
  if (interpolateSuccess != Fun4AllReturnCodes::EVENT_OK)
  {
    std::cout << PHWHERE << " Return code for lamination interpolation was " << interpolateSuccess << " and not successful" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  for (int s = 0; s < 2; s++)
  {
    scaleFactorMap[s] = (TH2 *) m_dcc_out->m_hDPint[s]->Clone();
    scaleFactorMap[s]->SetName((boost::format("scaleFactorMap%d") %s).str().c_str());
    scaleFactorMap[s]->Divide(simPhiDistortion[s]);
  }

  TH3 *hIntDistortionR_posz = (TH3 *) simDistortion->Get("hIntDistortionR_posz");
  hIntDistortionR_posz->GetZaxis()->SetRange(2, 2);
  TH2 *simRDistortion[2];
  simRDistortion[1] = (TH2 *) hIntDistortionR_posz->Project3D("yx");
  TH3 *hIntDistortionR_negz = (TH3 *) simDistortion->Get("hIntDistortionR_negz");
  hIntDistortionR_negz->GetZaxis()->SetRange(hIntDistortionR_negz->GetNbinsZ() - 1, hIntDistortionR_negz->GetNbinsZ() - 1);
  simRDistortion[0] = (TH2 *) hIntDistortionR_negz->Project3D("yx");

  for (int s = 0; s < 2; s++)
  {
    /*
    for(int i=1; i<=m_dcc_out->m_hDRint[s]->GetNbinsX(); i++)
    {
            for(int j=1; j<=m_dcc_out->m_hDRint[s]->GetNbinsY(); j++)
            {
                    if(simRDistortion[s]->GetBinContent(i,j) != 0.0)
                    {
                            m_dcc_out->m_hDRint[s]->SetBinContent(i,j, simRDistortion[s]->GetBinContent(i,j));
                    }
            }
    }
    */
    m_dcc_out->m_hDRint[s] = (TH2 *) simRDistortion[s]->Clone();
    m_dcc_out->m_hDRint[s]->SetName((boost::format("hIntDistortionR%s") %(s == 0 ? "_negz" : "_posz")).str().c_str());
    m_dcc_out->m_hDRint[s]->Multiply(scaleFactorMap[s]);
  }


  fill_guarding_bins(m_dcc_out);


  for(int s=0; s<2; s++)
  {
    for(int l=0; l<18; l++)
    {
      m_side = s;
      m_lamIndex = s*18 + l;
      m_lamPhi = m_laminationCenter[l][s];
      m_goodFit = m_laminationGoodFit[l][s];
      m_A = m_fLamination[l][s]->GetParameter(0);
      m_B = m_fLamination[l][s]->GetParameter(1);
      m_C = m_fLamination[l][s]->GetParameter(2);
      m_A_err = m_fLamination[l][s]->GetParError(0);
      m_B_err = m_fLamination[l][s]->GetParError(1);
      m_C_err = m_fLamination[l][s]->GetParError(2);
      m_dist = m_distanceToFit[l][s];
      m_nBins = m_nBinsFit[l][s];
      m_laminationTree->Fill();
    }
  }

  TFile *outputfile = new TFile(m_outputfile.c_str(), "RECREATE");
  outputfile->cd();
  for (int s = 0; s < 2; s++)
  {
    for (const auto &h : {m_dcc_out->m_hDRint[s], m_dcc_out->m_hDPint[s], m_dcc_out->m_hDZint[s], m_dcc_out->m_hentries[s]})
    {
      if (h)
      {
        h->Write();
      }
    }
    phiDistortionLamination[s]->Write();
    scaleFactorMap[s]->Write();
  }
  m_laminationTree->Write();

  m_hLamination[13][0]->Write();
  m_hLamination[13][1]->Write();
  m_hLamination[14][1]->Write();

  
  outputfile->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

//________________________________
void TpcLaminationFitting::fill_guarding_bins(TpcDistortionCorrectionContainer *dcc)
{
  for (int s = 0; s < 2; s++)
  {
    for (const auto &h : {dcc->m_hDRint[s], dcc->m_hDPint[s]})
    {
      const auto phibins = h->GetNbinsX();
      const auto rbins = h->GetNbinsY();

      int r30Bin = h->GetYaxis()->FindBin(30.0) + 1;

      for (int ir = 0; ir < rbins; ir++)
      {
        h->SetBinContent(1, ir + 1, h->GetBinContent(phibins - 1, ir + 1));
        h->SetBinError(1, ir + 1, h->GetBinError(phibins - 1, ir + 1));

        h->SetBinContent(phibins, ir + 1, h->GetBinContent(2, ir + 1));
        h->SetBinError(phibins, ir + 1, h->GetBinError(2, ir + 1));
      }

      for (int iphi = 0; iphi < phibins; iphi++)
      {
        h->SetBinContent(iphi + 1, r30Bin-1, h->GetBinContent(iphi + 1, r30Bin));
        h->SetBinError(iphi + 1, r30Bin-1, h->GetBinError(iphi + 1, r30Bin));

        h->SetBinContent(iphi + 1, rbins, h->GetBinContent(iphi + 1, rbins - 1));
        h->SetBinError(iphi + 1, rbins, h->GetBinError(iphi + 1, rbins - 1));
      }
    }
  }
}


//_____________________________
void TpcLaminationFitting::setOutputfile(const std::string &outputfile)
{
  m_outputfile = outputfile;
}
