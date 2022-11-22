/**
 * \file PHTpcCentralMembraneMatcher.cc
 * \brief match reconstructed CM clusters to CM pads, calculate differences, store on the node tree and compute distortion reconstruction maps
 * \author Tony Frawley <frawley@fsunuc.physics.fsu.edu>, Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "PHTpcCentralMembraneMatcher.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <trackbase/CMFlashClusterv1.h>
#include <trackbase/CMFlashClusterContainerv1.h>
#include <trackbase/CMFlashDifferencev1.h>
#include <trackbase/CMFlashDifferenceContainerv1.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TNtuple.h>
#include <TString.h>
#include <TVector3.h>

#include <cmath>
#include <set>
#include <string>

namespace 
{
  template<class T> inline constexpr T delta_phi(const T& phi)
  {
    if (phi > M_PI) return phi - 2. * M_PI;
    else if (phi <= -M_PI) return phi + 2.* M_PI;
    else return phi;
  }
  
  template<class T> inline constexpr T square( const T& x ) { return x*x; }
  
  template<class T> inline T get_r( const T& x, const T& y ) { return std::sqrt( square(x) + square(y) ); }
  
  // stream acts vector3
  [[maybe_unused]] std::ostream& operator << (std::ostream& out, const Acts::Vector3& v )
  { 
    out << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
    return out;
  }

  /// normalize distortions based on the number of entries in each cell, as recorded in the m_hentries histogram
  [[maybe_unused]] void normalize_distortions( TpcDistortionCorrectionContainer* dcc )
  {
    // loop over side
    for( unsigned int i = 0; i<2; ++i )
    {

      // loop over bins in entries
      for( int ip = 0; ip < dcc->m_hentries[i]->GetNbinsX(); ++ip )
        for( int ir = 0; ir < dcc->m_hentries[i]->GetNbinsY(); ++ir )
      {
        // count number of times a given cell was filled
        const auto entries = dcc->m_hentries[i]->GetBinContent( ip+1, ir+1 );
        if( entries <= 1 ) continue;

        // normalize histograms
        for( const auto& h:{dcc->m_hDRint[i], dcc->m_hDPint[i], dcc->m_hDZint[i]} )
        {
          h->SetBinContent( ip+1, ir+1, h->GetBinContent( ip+1, ir+1 )/entries );
          h->SetBinError( ip+1, ir+1, h->GetBinError( ip+1, ir+1 )/entries );
        }
      }
    }
  }

  /// fill distortion correction histograms' guarding bins, to allow ::Interpolate to work over the full acceptance
  [[maybe_unused]] void fill_guarding_bins( TpcDistortionCorrectionContainer* dcc )
  {

    // loop over side
    for( unsigned int i = 0; i<2; ++i )
    {
      for( const auto& h:{dcc->m_hDRint[i], dcc->m_hDPint[i], dcc->m_hDZint[i], dcc->m_hentries[i]} )
      {

        // fill guarding phi bins
        /*
        * we use 2pi periodicity to do that:
        * - last valid bin is copied to first guarding bin;
        * - first valid bin is copied to last guarding bin
        */
        const auto phibins = h->GetNbinsX();
        const auto rbins = h->GetNbinsY();
        for( int ir = 0; ir < rbins; ++ir )
        {
          // copy last valid bin to first guarding bin
          h->SetBinContent( 1, ir+1, h->GetBinContent( phibins-1, ir+1 ) );
          h->SetBinError( 1, ir+1, h->GetBinError( phibins-1, ir+1 ) );

          // copy first valid bin to last guarding bin
          h->SetBinContent( phibins, ir+1, h->GetBinContent( 2, ir+1 ) );
          h->SetBinError( phibins, ir+1, h->GetBinError( 2, ir+1 ) );
        }

        // fill guarding r bins
        for( int iphi = 0; iphi < phibins; ++iphi )
        {
          // copy first valid bin to first guarding bin
          h->SetBinContent( iphi+1, 1, h->GetBinContent( iphi+1, 2 ) );
          h->SetBinError( iphi+1, 1, h->GetBinError( iphi+1, 2 ) );

          // copy last valid bin to last guarding bin
          h->SetBinContent( iphi+1, rbins, h->GetBinContent( iphi+1, rbins-1 ) );
          h->SetBinError( iphi+1, rbins, h->GetBinError( iphi+1, rbins-1 ) );
        }
      }
    }
  }

}

//____________________________________________________________________________..
PHTpcCentralMembraneMatcher::PHTpcCentralMembraneMatcher(const std::string &name):
  SubsysReco(name)
{
  // calculate stripes center positions
  CalculateCenters(nPads_R1, R1_e, nGoodStripes_R1_e, keepUntil_R1_e, nStripesIn_R1_e, nStripesBefore_R1_e, cx1_e, cy1_e);
  CalculateCenters(nPads_R1, R1, nGoodStripes_R1, keepUntil_R1, nStripesIn_R1, nStripesBefore_R1, cx1, cy1);
  CalculateCenters(nPads_R2, R2, nGoodStripes_R2, keepUntil_R2, nStripesIn_R2, nStripesBefore_R2, cx2, cy2);
  CalculateCenters(nPads_R3, R3, nGoodStripes_R3, keepUntil_R3, nStripesIn_R3, nStripesBefore_R3, cx3, cy3);
}


//___________________________________________________________
void PHTpcCentralMembraneMatcher::set_grid_dimensions( int phibins, int rbins )
{
  m_phibins = phibins;
  m_rbins = rbins;
}

//____________________________________________________________________________..
int PHTpcCentralMembraneMatcher::InitRun(PHCompositeNode *topNode)
{
  if( m_savehistograms )
  { 

    static constexpr float max_dr = 1.0;
    static constexpr float max_dphi = 0.05;

    fout.reset( new TFile(m_histogramfilename.c_str(),"RECREATE") ); 
    hxy_reco = new TH2F("hxy_reco","reco cluster x:y",800,-100,100,800,-80,80);
    hxy_truth = new TH2F("hxy_truth","truth cluster x:y",800,-100,100,800,-80,80);
    
    hdrdphi = new TH2F("hdrdphi","dr vs dphi",800,-max_dr,max_dr,800,-max_dphi,max_dphi);
    hdrdphi->GetXaxis()->SetTitle("dr");  
    hdrdphi->GetYaxis()->SetTitle("dphi");  
    
    hrdr = new TH2F("hrdr","dr vs r",800,0.0,80.0,800,-max_dr, max_dr);
    hrdr->GetXaxis()->SetTitle("r");  
    hrdr->GetYaxis()->SetTitle("dr");  
    
    hrdphi = new TH2F("hrdphi","dphi vs r",800,0.0,80.0,800,-max_dphi,max_dphi);
    hrdphi->GetXaxis()->SetTitle("r");  
    hrdphi->GetYaxis()->SetTitle("dphi");
    
    hdphi = new TH1F("hdphi","dph",800,-max_dphi,max_dphi);
    hdphi->GetXaxis()->SetTitle("dphi");

    hdr1_single = new TH1F("hdr1_single", "innner dr single", 200,-max_dr, max_dr);
    hdr2_single = new TH1F("hdr2_single", "mid dr single", 200,-max_dr, max_dr);
    hdr3_single = new TH1F("hdr3_single", "outer dr single", 200,-max_dr, max_dr);
    hdr1_double = new TH1F("hdr1_double", "innner dr double", 200,-max_dr, max_dr);
    hdr2_double = new TH1F("hdr2_double", "mid dr double", 200,-max_dr, max_dr);
    hdr3_double = new TH1F("hdr3_double", "outer dr double", 200,-max_dr, max_dr);
    hdrphi = new TH1F("hdrphi","r * dphi", 200, -0.05, 0.05);
    hnclus = new TH1F("hnclus", " nclusters ", 3, 0., 3.);
  }
  
  // Get truth cluster positions
  //=====================
  
  const double phi_petal = M_PI / 9.0;  // angle span of one petal
   
  /*
   * utility function to
   * - duplicate generated truth position to cover both sides of the central membrane
   * - assign proper z,
   * - insert in container
   */
  auto save_truth_position = [&](TVector3 source) 
  {
    source.SetZ( +1 );
    m_truth_pos.push_back( source );
    
    source.SetZ( -1 );
    m_truth_pos.push_back( source );
  };
  
  // inner region extended is the 8 layers inside 30 cm    
  for(int j = 0; j < nRadii; ++j)
    for(int i=0; i < nGoodStripes_R1_e[j]; ++i)
      for(int k =0; k<18; ++k)
	{
	  TVector3 dummyPos(cx1_e[i][j], cy1_e[i][j], 0.0);
	  dummyPos.RotateZ(k * phi_petal);
	  save_truth_position(dummyPos);

	  if(Verbosity() > 2)	  
	    std::cout << " i " << i << " j " << j << " k " << k << " x1 " << dummyPos.X() << " y1 " << dummyPos.Y()
		      <<  " theta " << std::atan2(dummyPos.Y(), dummyPos.X())
		      << " radius " << get_r( dummyPos.X(), dummyPos.y()) << std::endl; 
	  if(m_savehistograms) hxy_truth->Fill(dummyPos.X(),dummyPos.Y());      	  
	}  
  
  // inner region is the 8 layers outside 30 cm  
  for(int j = 0; j < nRadii; ++j)
    for(int i=0; i < nGoodStripes_R1[j]; ++i)
      for(int k =0; k<18; ++k)
	{
	  TVector3 dummyPos(cx1[i][j], cy1[i][j], 0.0);
	  dummyPos.RotateZ(k * phi_petal);
	  save_truth_position(dummyPos);

	  if(Verbosity() > 2)	  
	    std::cout << " i " << i << " j " << j << " k " << k << " x1 " << dummyPos.X() << " y1 " << dummyPos.Y()
		      <<  " theta " << std::atan2(dummyPos.Y(), dummyPos.X())
		      << " radius " << get_r( dummyPos.X(), dummyPos.y()) << std::endl; 
	  if(m_savehistograms) hxy_truth->Fill(dummyPos.X(),dummyPos.Y());      	  
	}  

  for(int j = 0; j < nRadii; ++j)
    for(int i=0; i < nGoodStripes_R2[j]; ++i)
      for(int k =0; k<18; ++k)
	{
	  TVector3 dummyPos(cx2[i][j], cy2[i][j], 0.0);
	  dummyPos.RotateZ(k * phi_petal);
	  save_truth_position(dummyPos);

	  if(Verbosity() > 2)	  
	    std::cout << " i " << i << " j " << j << " k " << k << " x1 " << dummyPos.X() << " y1 " << dummyPos.Y()
		      <<  " theta " << std::atan2(dummyPos.Y(), dummyPos.X())
		      << " radius " << get_r( dummyPos.X(), dummyPos.y()) << std::endl; 
	  if(m_savehistograms) hxy_truth->Fill(dummyPos.X(),dummyPos.Y());      	  
	}      	  
  
  for(int j = 0; j < nRadii; ++j)
    for(int i=0; i < nGoodStripes_R3[j]; ++i)
      for(int k =0; k<18; ++k)
	{
	  TVector3 dummyPos(cx3[i][j], cy3[i][j], 0.0);
	  dummyPos.RotateZ(k * phi_petal);
	  save_truth_position(dummyPos);
	  
	  if(Verbosity() > 2)
	    std::cout << " i " << i << " j " << j << " k " << k << " x1 " << dummyPos.X() << " y1 " << dummyPos.Y()
		      <<  " theta " << std::atan2(dummyPos.Y(), dummyPos.X())
		      << " radius " << get_r( dummyPos.X(), dummyPos.y()) << std::endl; 
	  if(m_savehistograms) hxy_truth->Fill(dummyPos.X(),dummyPos.Y());      	  
	}   
  
  int ret = GetNodes(topNode);
  return ret;
}

//____________________________________________________________________________..
int PHTpcCentralMembraneMatcher::process_event(PHCompositeNode * /*topNode*/)
{
  std::vector<TVector3> reco_pos;
  std::vector<unsigned int> reco_nclusters;
 
  // reset output distortion correction container histograms
  for( const auto& harray:{m_dcc_out->m_hDRint, m_dcc_out->m_hDPint, m_dcc_out->m_hDZint, m_dcc_out->m_hentries} )
  { for( const auto& h:harray ) { h->Reset(); } }
  
  // read the reconstructed CM clusters
  auto clusrange = m_corrected_CMcluster_map->getClusters();
  for (auto cmitr = clusrange.first;
       cmitr !=clusrange.second;
       ++cmitr)
    {
      const auto& [cmkey, cmclus] = *cmitr;
      const unsigned int nclus = cmclus->getNclusters();
      
      // Do the static + average distortion corrections if the container was found
      Acts::Vector3 pos(cmclus->getX(), cmclus->getY(), cmclus->getZ());
      if( m_dcc_in) pos = m_distortionCorrection.get_corrected_position( pos, m_dcc_in ); 
      
      TVector3 tmp_pos(pos[0], pos[1], pos[2]);
      reco_pos.push_back(tmp_pos);      
      reco_nclusters.push_back(nclus);

      if(Verbosity())
	{
	  double raw_rad = sqrt( cmclus->getX()*cmclus->getX() + cmclus->getY()*cmclus->getY() );
	  double corr_rad = sqrt( tmp_pos.X()*tmp_pos.X() + tmp_pos.Y()*tmp_pos.Y() );
	  std::cout << "found raw cluster " << cmkey << " with x " << cmclus->getX() << " y " << cmclus->getY() << " z " << cmclus->getZ()   << " radius " << raw_rad << std::endl; 
	  std::cout << "                --- corrected positions: " << tmp_pos.X() << "  " << tmp_pos.Y() << "  " << tmp_pos.Z() << " radius " << corr_rad << std::endl; 
	}

      if(m_savehistograms)
	{      
	  hxy_reco->Fill(tmp_pos.X(), tmp_pos.Y());
	}
    }

  // Match reco and truth positions
  //std::map<unsigned int, unsigned int> matched_pair;
  std::vector<std::pair<unsigned int, unsigned int>> matched_pair;
  std::vector<unsigned int> matched_nclus;
      
  // loop over truth positions
  for(unsigned int i=0; i<m_truth_pos.size(); ++i)
  {
    const double z1 = m_truth_pos[i].Z(); 
    const double rad1= get_r( m_truth_pos[i].X(),m_truth_pos[i].Y());
    const double phi1 = m_truth_pos[i].Phi();
    
    // loop over cluster positions
    for(unsigned int j = 0; j < reco_pos.size(); ++j)
    {
      
      const auto& nclus = reco_nclusters[j];
      const double z2 = reco_pos[j].Z(); 
      const double rad2=get_r(reco_pos[j].X(), reco_pos[j].Y());
      const double phi2 = reco_pos[j].Phi();
    
      // only match pairs that are on the same side of the TPC
      const bool accepted_z = ((z1>0)==(z2>0));
      if( !accepted_z ) continue;
      
      const auto dr = rad1-rad2;
      const bool accepted_r = std::abs(dr) < m_rad_cut;

      const auto dphi = delta_phi(phi1-phi2);
      const bool accepted_phi = std::abs(dphi) < m_phi_cut;
      
      if(m_savehistograms)
      {
       
        hnclus->Fill( (float) nclus);

        double r =  rad2;

        if( accepted_r )
        {
          hdrphi->Fill(r * dphi);
          hdphi->Fill(dphi);
          hrdphi->Fill(r,dphi);
        }

        if( accepted_r && accepted_phi)
        { hdrdphi->Fill(dr, dphi); }

        if( accepted_phi )
        {
          hrdr->Fill(r,dr);
          if(nclus==1)
          {
            if(r < 40.0) hdr1_single->Fill(dr); 
            if(r >= 40.0 && r < 58.0) hdr2_single->Fill(dr); 
            if(r >= 58.0) hdr3_single->Fill(dr); 	  
          }
          else
          {
            if(r < 40.0) hdr1_double->Fill(dr); 
            if(r >= 40.0 && r < 58.0) hdr2_double->Fill(dr); 
            if(r >= 58.0) hdr3_double->Fill(dr); 	  
          }
        }
      }
      
      if( accepted_r && accepted_phi ) 
      {
        matched_pair.emplace_back(i,j);
        matched_nclus.push_back(nclus);
        break;
      }
    }      
  }
  
  // print some statistics: 
  if( Verbosity() )
  {
    const auto n_valid_truth = std::count_if( m_truth_pos.begin(), m_truth_pos.end(), []( const TVector3& pos ) { return get_r( pos.x(), pos.y() ) >  30; } );
    const auto n_reco_size1 = std::count_if( reco_nclusters.begin(), reco_nclusters.end(), []( const unsigned int& value ) { return value==1; } );
    const auto n_reco_size2 = std::count_if( reco_nclusters.begin(), reco_nclusters.end(), []( const unsigned int& value ) { return value==2; } );
    std::cout << "PHTpcCentralMembraneMatcher::process_event - m_truth_pos size: " << m_truth_pos.size() << std::endl;
    std::cout << "PHTpcCentralMembraneMatcher::process_event - m_truth_pos size, r>30cm: " << n_valid_truth << std::endl;
    std::cout << "PHTpcCentralMembraneMatcher::process_event - reco_pos size: " << reco_pos.size() << std::endl;
    std::cout << "PHTpcCentralMembraneMatcher::process_event - reco_pos size (nclus==1): " << n_reco_size1 << std::endl;
    std::cout << "PHTpcCentralMembraneMatcher::process_event - reco_pos size (nclus==2): " << n_reco_size2 << std::endl;
    std::cout << "PHTpcCentralMembraneMatcher::process_event - matched_pair size: " << matched_pair.size() << std::endl;
  }
  
  for(unsigned int ip = 0; ip < matched_pair.size(); ++ip)
  {
    const std::pair<unsigned int, unsigned int>& p = matched_pair[ip];
    const unsigned int& nclus = matched_nclus[ip];

    // add to node tree
    unsigned int key = p.first;
    auto cmdiff = new CMFlashDifferencev1();
    cmdiff->setTruthPhi(m_truth_pos[p.first].Phi());
    cmdiff->setTruthR(m_truth_pos[p.first].Perp());
    cmdiff->setTruthZ(m_truth_pos[p.first].Z());
    
    cmdiff->setRecoPhi(reco_pos[p.second].Phi());
    cmdiff->setRecoR(reco_pos[p.second].Perp());
    cmdiff->setRecoZ(reco_pos[p.second].Z());
    
    cmdiff->setNclusters(nclus);
    
    m_cm_flash_diffs->addDifferenceSpecifyKey(key, cmdiff);
    
    // store cluster position
    const double clus_r = reco_pos[p.second].Perp();
    double clus_phi = reco_pos[p.second].Phi();
    if ( clus_phi < 0 ) clus_phi += 2*M_PI;

    const double clus_z = reco_pos[p.second].z();
    const unsigned int side = (clus_z<0) ? 0:1;
    
    // calculate residuals (cluster - truth)
    const double dr = reco_pos[p.second].Perp() - m_truth_pos[p.first].Perp();
    const double dphi = delta_phi( reco_pos[p.second].Phi() - m_truth_pos[p.first].Phi() );
    const double rdphi = reco_pos[p.second].Perp() * dphi;
    const double dz = reco_pos[p.second].z() - m_truth_pos[p.first].z();

    // fill distortion correction histograms
    /* 
     * TODO: 
     * - we might need to only fill the histograms for cm clusters that have 2 clusters only
     * - we might need a smoothing procedure to fill the bins that have no entries using neighbors
     */
    for( const auto& dcc:{m_dcc_out, m_dcc_out_aggregated.get()} )
    {
      static_cast<TH2*>(dcc->m_hDRint[side])->Fill( clus_phi, clus_r, dr );
      static_cast<TH2*>(dcc->m_hDPint[side])->Fill( clus_phi, clus_r, rdphi );
      static_cast<TH2*>(dcc->m_hDZint[side])->Fill( clus_phi, clus_r, dz );
      static_cast<TH2*>(dcc->m_hentries[side])->Fill( clus_phi, clus_r );
    }
    
  }
  
  if(Verbosity())
  {
    std::cout << "PHTpcCentralMembraneMatcher::process_events - cmclusters: " << m_corrected_CMcluster_map->size() << std::endl;
    std::cout << "PHTpcCentralMembraneMatcher::process_events - matched pairs: " << matched_pair.size() << std::endl;
    std::cout << "PHTpcCentralMembraneMatcher::process_events - differences: " << m_cm_flash_diffs->size() << std::endl;
    std::cout << "PHTpcCentralMembraneMatcher::process_events - entries: " << m_dcc_out->m_hentries[0]->GetEntries() << ", " << m_dcc_out->m_hentries[1]->GetEntries() << std::endl;  
  }

  // normalize per-event distortion correction histograms and fill guarding bins
  normalize_distortions( m_dcc_out );
  fill_guarding_bins( m_dcc_out );
    
  if(Verbosity())
    {	
      // read back differences from node tree as a check
      auto diffrange = m_cm_flash_diffs->getDifferences();
      for (auto cmitr = diffrange.first;
	   cmitr !=diffrange.second;
	   ++cmitr)
	{
	  auto key = cmitr->first;
	  auto cmreco = cmitr->second;
	  
	  std::cout << " key " << key 
		    << " nclus " << cmreco->getNclusters() 
		    << " truth Phi " << cmreco->getTruthPhi() << " reco Phi " << cmreco->getRecoPhi()
		    << " truth R " << cmreco->getTruthR() << " reco R " << cmreco->getRecoR() 
		    << " truth Z " << cmreco->getTruthZ() << " reco Z " << cmreco->getRecoZ() 
		    << std::endl;
	}
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTpcCentralMembraneMatcher::End(PHCompositeNode * /*topNode*/ )
{

  // write distortion corrections
  if( m_dcc_out_aggregated )
  {
 
    // normalize aggregated distortion correction histograms and fill guarding bins
    normalize_distortions( m_dcc_out_aggregated.get() );
    fill_guarding_bins( m_dcc_out_aggregated.get() );

    // create TFile and write all histograms
    std::unique_ptr<TFile> outputfile( TFile::Open( m_outputfile.c_str(), "RECREATE" ) );
    outputfile->cd();

    // loop over side
    for( unsigned int i = 0; i<2; ++i )
    {
      for( const auto& h:{m_dcc_out_aggregated->m_hDRint[i], m_dcc_out_aggregated->m_hDPint[i], m_dcc_out_aggregated->m_hDZint[i], m_dcc_out_aggregated->m_hentries[i]} )
      { if( h ) h->Write(); }
    }
  }
  
  // write evaluation histograms
  if(m_savehistograms && fout)
  {
    fout->cd();
    
    hxy_reco->Write();
    hxy_truth->Write();
    hdrdphi->Write();
    hrdr->Write();
    hrdphi->Write();
    hdphi->Write();
    hdrphi->Write();
    hdr1_single->Write();
    hdr2_single->Write();
    hdr3_single->Write();
    hdr1_double->Write();
    hdr2_double->Write();
    hdr3_double->Write();
    hnclus->Write();
    
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

  m_corrected_CMcluster_map  = findNode::getClass<CMFlashClusterContainer>(topNode, "CORRECTED_CM_CLUSTER");
  if(!m_corrected_CMcluster_map)
    {
      std::cout << PHWHERE << "CORRECTED_CM_CLUSTER Node missing, abort." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }      

  // input tpc distortion correction
  m_dcc_in = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainer");
  if( m_dcc_in )
    { 
      std::cout << "PHTpcCentralMembraneMatcher:   found TPC distortion correction container" << std::endl; 
    }

  // create node for results of matching
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
  
  m_cm_flash_diffs = new CMFlashDifferenceContainerv1;
  PHIODataNode<PHObject> *CMFlashDifferenceNode =
    new PHIODataNode<PHObject>(m_cm_flash_diffs, "CM_FLASH_DIFFERENCES", "PHObject");
  DetNode->addNode(CMFlashDifferenceNode);
  
//   // output tpc fluctuation distortion container
//   // this one is filled on the fly on a per-CM-event basis, and applied in the tracking chain
//   const std::string dcc_out_node_name = "TpcDistortionCorrectionContainerFluctuation";
//   m_dcc_out = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,dcc_out_node_name);
//   if( !m_dcc_out )
//   { 
//   
//     /// Get the DST node and check
//     auto dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
//     if (!dstNode)
//     {
//       std::cout << "PHTpcCentralMembraneMatcher::InitRun - DST Node missing, quitting" << std::endl;
//       return Fun4AllReturnCodes::ABORTRUN;
//     }
//     
//     // Get the tracking subnode and create if not found
//     auto svtxNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "SVTX"));
//     if (!svtxNode)
//     {
//       svtxNode = new PHCompositeNode("SVTX");
//       dstNode->addNode(svtxNode);
//     }
// 
//     std::cout << "PHTpcCentralMembraneMatcher::GetNodes - creating TpcDistortionCorrectionContainer in node " << dcc_out_node_name << std::endl;
//     m_dcc_out = new TpcDistortionCorrectionContainer;
//     auto node = new PHDataNode<TpcDistortionCorrectionContainer>(m_dcc_out, dcc_out_node_name);
//     svtxNode->addNode(node);
//   }

  // create per event distortions. Do not put on the node tree
  m_dcc_out = new TpcDistortionCorrectionContainer;

  // also prepare the local distortion container, used to aggregate multple events 
  m_dcc_out_aggregated.reset( new TpcDistortionCorrectionContainer );

  // compute axis limits to include guarding bins, needed for TH2::Interpolate to work
  const float phiMin = m_phiMin - (m_phiMax-m_phiMin)/m_phibins;
  const float phiMax = m_phiMax + (m_phiMax-m_phiMin)/m_phibins;
  
  const float rMin = m_rMin - (m_rMax-m_rMin)/m_rbins;
  const float rMax = m_rMax + (m_rMax-m_rMin)/m_rbins;

  // reset all output distortion container so that they match the requested grid size
  const std::array<const std::string,2> extension = {{ "_negz", "_posz" }};
  for( const auto& dcc:{m_dcc_out, m_dcc_out_aggregated.get()} )
  {
    // set dimensions to 2, since central membrane flashes only provide distortions at z = 0
    dcc->dimensions = 2;
    
    // create all histograms
    for( int i =0; i < 2; ++i )
    {
      delete dcc->m_hDPint[i]; dcc->m_hDPint[i] = new TH2F( Form("hIntDistortionP%s", extension[i].c_str()), Form("hIntDistortionP%s", extension[i].c_str()), m_phibins+2, phiMin, phiMax, m_rbins+2, rMin, rMax );
      delete dcc->m_hDRint[i]; dcc->m_hDRint[i] = new TH2F( Form("hIntDistortionR%s", extension[i].c_str()), Form("hIntDistortionR%s", extension[i].c_str()), m_phibins+2, phiMin, phiMax, m_rbins+2, rMin, rMax );
      delete dcc->m_hDZint[i]; dcc->m_hDZint[i] = new TH2F( Form("hIntDistortionZ%s", extension[i].c_str()), Form("hIntDistortionZ%s", extension[i].c_str()), m_phibins+2, phiMin, phiMax, m_rbins+2, rMin, rMax );
      delete dcc->m_hentries[i]; dcc->m_hentries[i] = new TH2I( Form("hEntries%s", extension[i].c_str()), Form("hEntries%s", extension[i].c_str()), m_phibins+2, phiMin, phiMax, m_rbins+2, rMin, rMax );
    }
  }
  
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

      if( Verbosity() > 2 )
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


