#include "PHTpcCentralMembraneMatcher.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <string>
#include <TVector3.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <cmath>

/// Tracking includes
#include <trackbase/CMFlashClusterv1.h>
#include <trackbase/CMFlashClusterContainerv1.h>
#include <trackbase/CMFlashDifferencev1.h>
#include <trackbase/CMFlashDifferenceContainerv1.h>

using namespace std;

//____________________________________________________________________________..
PHTpcCentralMembraneMatcher::PHTpcCentralMembraneMatcher(const std::string &name):
  SubsysReco(name)
{

  // Make some histograms on an output file for diagnostics
  char temp[500];
  sprintf(temp,
	  "./eval_output/Matcher_Histograms_%i.root", _process);
  fout = new TFile(temp,"RECREATE");
  
  hxy_reco = new TH2F("hxy_reco","reco cluster x:y",800,-100,100,800,-80,80);
  hxy_truth = new TH2F("hxy_truth","truth cluster x:y",800,-100,100,800,-80,80);
  hdrdphi = new TH2F("hdrdphi","dr vs dphi",800,-0.5,0.5,800,-0.001,0.001);
  hdrdphi->GetXaxis()->SetTitle("dr");  
  hdrdphi->GetYaxis()->SetTitle("dphi");  
  hrdr = new TH2F("hrdr","dr vs r",800,0.0,80.0,800,-0.5,0.5);
  hrdr->GetXaxis()->SetTitle("r");  
  hrdr->GetYaxis()->SetTitle("dr");  
  hrdphi = new TH2F("hrdphi","dphi vs r",800,0.0,80.0,800,-0.001,0.001);
  hrdphi->GetXaxis()->SetTitle("r");  
  hrdphi->GetYaxis()->SetTitle("dphi");  
  
  // Get truth cluster positions
  //=====================
  
  CalculateCenters(nPads_R1, R1_e, nGoodStripes_R1_e, keepUntil_R1_e, nStripesIn_R1_e, nStripesBefore_R1_e, cx1_e, cy1_e);
  CalculateCenters(nPads_R1, R1, nGoodStripes_R1, keepUntil_R1, nStripesIn_R1, nStripesBefore_R1, cx1, cy1);
  CalculateCenters(nPads_R2, R2, nGoodStripes_R2, keepUntil_R2, nStripesIn_R2, nStripesBefore_R2, cx2, cy2);
  CalculateCenters(nPads_R3, R3, nGoodStripes_R3, keepUntil_R3, nStripesIn_R3, nStripesBefore_R3, cx3, cy3);
  
  TVector3 dummyPos;
  const double phi_petal = M_PI / 9.0;  // angle span of one petal
  
  // k is the petal ID, rotate the stripe center to this petal	      
  
  // inner region extended is the 8 layers inside 30 cm    
  for(int j = 0; j < nRadii; ++j)
    for(int i=0; i < nGoodStripes_R1_e[j]; ++i)
      for(int k =0; k<18; ++k)
	{
	  dummyPos.SetXYZ(cx1_e[i][j], cy1_e[i][j], 0.0);
	  dummyPos.RotateZ(k * phi_petal);
	  
	  truth_pos.push_back(dummyPos);
	  
	  std::cout << " i " << i << " j " << j << " k " << k << " x1 " << dummyPos.X() << " y1 " << dummyPos.Y()
		    <<  " theta " << atan2(dummyPos.Y(), dummyPos.X())
		    << " radius " << sqrt(pow(dummyPos.X(),2)+pow(dummyPos.Y(),2)) << std::endl; 
	  if(_histos) hxy_truth->Fill(dummyPos.X(),dummyPos.Y());      	  
	}  
  
  // inner region is the 8 layers outside 30 cm  
  for(int j = 0; j < nRadii; ++j)
    for(int i=0; i < nGoodStripes_R1[j]; ++i)
      for(int k =0; k<18; ++k)
	{
	  dummyPos.SetXYZ(cx1[i][j], cy1[i][j], 0.0);
	  dummyPos.RotateZ(k * phi_petal);
	  
	  truth_pos.push_back(dummyPos);
	  
	  std::cout << " i " << i << " j " << j << " k " << k << " x1 " << dummyPos.X() << " y1 " << dummyPos.Y()
		    <<  " theta " << atan2(dummyPos.Y(), dummyPos.X())
		    << " radius " << sqrt(pow(dummyPos.X(),2)+pow(dummyPos.Y(),2)) << std::endl; 
	  if(_histos) hxy_truth->Fill(dummyPos.X(),dummyPos.Y());      	  
	}  
  
  for(int j = 0; j < nRadii; ++j)
    for(int i=0; i < nGoodStripes_R2[j]; ++i)
      for(int k =0; k<18; ++k)
	{
	  dummyPos.SetXYZ(cx2[i][j], cy2[i][j], 0.0);
	  dummyPos.RotateZ(k * phi_petal);
	  
	  truth_pos.push_back(dummyPos);
	  
	  std::cout << " i " << i << " j " << j << " k " << k << " x1 " << dummyPos.X() << " y1 " << dummyPos.Y()
		    <<  " theta " << atan2(dummyPos.Y(), dummyPos.X())
		    << " radius " << sqrt(pow(dummyPos.X(),2)+pow(dummyPos.Y(),2)) << std::endl; 
	  if(_histos) hxy_truth->Fill(dummyPos.X(),dummyPos.Y());      	  
	}      	  
  
  for(int j = 0; j < nRadii; ++j)
    for(int i=0; i < nGoodStripes_R3[j]; ++i)
      for(int k =0; k<18; ++k)
	{
	  dummyPos.SetXYZ(cx3[i][j], cy3[i][j], 0.0);
	  dummyPos.RotateZ(k * phi_petal);

	  truth_pos.push_back(dummyPos);
	  
	  std::cout << " i " << i << " j " << j << " k " << k << " x1 " << dummyPos.X() << " y1 " << dummyPos.Y()
		    <<  " theta " << atan2(dummyPos.Y(), dummyPos.X())
		    << " radius " << sqrt(pow(dummyPos.X(),2)+pow(dummyPos.Y(),2)) << std::endl; 
	  if(_histos) hxy_truth->Fill(dummyPos.X(),dummyPos.Y());      	  
	}      	  
}

//____________________________________________________________________________..
PHTpcCentralMembraneMatcher::~PHTpcCentralMembraneMatcher()
{

}

//____________________________________________________________________________..
int PHTpcCentralMembraneMatcher::InitRun(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);
  return ret;
}

//____________________________________________________________________________..
int PHTpcCentralMembraneMatcher::process_event(PHCompositeNode * /*topNode*/)
{
 
  // read the reconstructed CM clusters

  auto clusrange = _corrected_CMcluster_map->getClusters();
  for (auto cmitr = clusrange.first;
       cmitr !=clusrange.second;
       ++cmitr)
    {
      auto cmkey = cmitr->first;
      auto cmclus = cmitr->second;
      TVector3 tmp_pos(cmclus->getX(), cmclus->getY(), cmclus->getZ());

      // Do we want to do the static + average distortion corrections here?
      // That makes the most sense, since the pattern matching becomes far easier if the differences are small

      reco_pos.push_back(tmp_pos);

      if(Verbosity() > 0)
	std::cout << "found cluster " << cmkey << " with adc " << cmclus->getAdc() 
		  << " x " << cmclus->getX() << " y " << cmclus->getY() << cmclus->getZ()
		  << std::endl; 

      if(_histos)
	{      
	  hxy_reco->Fill(cmclus->getX(), cmclus->getY());
	}
    }

  // Match reco and truth positions
  std::map<unsigned int, unsigned int> matched_pair;
  // loop over truth positions
  for(unsigned int i=0; i<truth_pos.size(); ++i)
    {
      double rad1=truth_pos[i].Mag();
      double phi1 = truth_pos[i].Phi();
      for(unsigned int j = 0; j < reco_pos.size(); ++j)
	{
	  double rad2=reco_pos[j].Mag();
	  double phi2 = reco_pos[j].Phi();

	  if(fabs(rad1-rad2) < _rad_cut && fabs(phi1-phi2) < _phi_cut)
	    {
	      matched_pair.insert(std::make_pair(i,j));
	      break;
	    }
	}      
    }
  
 
  for(const auto &p : matched_pair)
    {
      if(_histos)
	{
	  double dr = truth_pos[p.first].Mag() - reco_pos[p.second].Mag();
	  double dphi = truth_pos[p.first].Phi() - reco_pos[p.second].Phi();
	  double r =  truth_pos[p.first].Mag();
	  
	  hdrdphi->Fill(dr, dphi);
	  hrdr->Fill(r,dr);
	  hrdphi->Fill(r,dphi);
	}

      unsigned int key = p.first;
      std::pair<float, float> reco_xy(reco_pos[p.second].X(), reco_pos[p.second].Y());

      std::cout << " key " << key << " reco_xy x value " << reco_xy.first << " reco_xy y value " << reco_xy.second << std::endl;

      // add to node tree
     auto cmdiff = new CMFlashDifferencev1();

     cmdiff->setTruthX(truth_pos[p.first].X());
     cmdiff->setTruthY(truth_pos[p.first].Y());      
     cmdiff->setRecoX(reco_pos[p.second].X());
     cmdiff->setRecoY(reco_pos[p.second].Y());
      
      _cm_flash_diffs->addDifferenceSpecifyKey(key, cmdiff);

    } 
  
  // read back differences from node tree as a check
  auto diffrange = _cm_flash_diffs->getDifferences();
  for (auto cmitr = diffrange.first;
       cmitr !=diffrange.second;
       ++cmitr)
    {
      auto key = cmitr->first;
      auto cmreco = cmitr->second;

      std::cout << " key " << key << " truth X " << cmreco->getTruthX() << " reco X " << cmreco->getRecoX()
		<< " truth Y " << cmreco->getTruthY() << " reco Y " << cmreco->getRecoY()
		<< std::endl;

    }


 
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTpcCentralMembraneMatcher::End(PHCompositeNode * /*topNode*/ )
{

  if(_histos)
    {
      fout->cd();
      
      hxy_reco->Write();
      hxy_truth->Write();
      hdrdphi->Write();
      hrdr->Write();
      hrdphi->Write();
            
      fout->Close();
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..

int  PHTpcCentralMembraneMatcher::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  _corrected_CMcluster_map  = findNode::getClass<CMFlashClusterContainer>(topNode, "CORRECTED_CM_CLUSTER");
  if(!_corrected_CMcluster_map)
    {
      std::cout << PHWHERE << "CORRECTED_CM_CLUSTER Node missing, abort." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }      


  std::cout << "Creating node CM_FLASH_DIFFERENCES" << std::endl;  
  PHNodeIterator iter(topNode);
  
  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }      
  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode =
    dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
  if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }
  
  _cm_flash_diffs = new CMFlashDifferenceContainerv1;
  PHIODataNode<PHObject> *CMFlashDifferenceNode =
    new PHIODataNode<PHObject>(_cm_flash_diffs, "CM_FLASH_DIFFERENCES", "PHObject");
  DetNode->addNode(CMFlashDifferenceNode);
  
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________
void PHTpcCentralMembraneMatcher::CalculateCenters(
    int nPads,
    const std::array<double, nRadii>& R,
    std::array<int, nRadii>& nGoodStripes,
    const std::array<int, nRadii>& keepUntil,
    std::array<int, nRadii>& nStripesIn,
    std::array<int, nRadii>& nStripesBefore,
    double cx[][nRadii], double cy[][nRadii])
{
  const double phi_module = M_PI / 6.0;  // angle span of a module
  const int pr_mult = 3;                 // multiples of intrinsic resolution of pads
  const int dw_mult = 8;                 // multiples of diffusion width
  const double diffwidth = 0.6 * mm;     // diffusion width
  const double adjust = 0.015;           //arbitrary angle to center the pattern in a petal

  double theta = 0.0;

  //center coords

  //calculate spacing first:
  std::array<double, nRadii> spacing;
  for (int i = 0; i < nRadii; i++)
  {
    spacing[i] = 2.0 * ((dw_mult * diffwidth / R[i]) + (pr_mult * phi_module / nPads));
  }

  //center calculation
  for (int j = 0; j < nRadii; j++)
  {
    int i_out = 0;
    for (int i = keepThisAndAfter[j]; i < keepUntil[j]; i++)
    {
      if (j % 2 == 0)
      {
        theta = i * spacing[j] + (spacing[j] / 2) - adjust;
        cx[i_out][j] = R[j] * cos(theta) / cm;
        cy[i_out][j] = R[j] * sin(theta) / cm;
      }
      else
      {
        theta = (i + 1) * spacing[j] - adjust;
        cx[i_out][j] = R[j] * cos(theta) / cm;
        cy[i_out][j] = R[j] * sin(theta) / cm;
      }

      std::cout << " j " << j << " i " << i << " i_out " << i_out << " theta " << theta << " cx " << cx[i_out][j] << " cy " << cy[i_out][j] 
		<< " radius " << sqrt(pow(cx[i_out][j],2)+pow(cy[i_out][j],2)) << std::endl; 

      i_out++;

      nStripesBefore_R1_e[0] = 0;

      nStripesIn[j] = keepUntil[j] - keepThisAndAfter[j];
      if (j==0)
      {
	nStripesBefore[j] = 0;
      }
      else
      {
        nStripesBefore[j] = nStripesIn[j - 1] + nStripesBefore[j - 1];
      }
      nStripesBefore_R1_e[0] = 0;
    }
    nGoodStripes[j] = i_out;
  }
}


