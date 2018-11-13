#include "PHG4TPCClusterizer.h"
#include "SvtxCluster.h"
#include "SvtxClusterMap.h"
#include "SvtxClusterMap_v1.h"
#include "SvtxCluster_v1.h"
#include "SvtxHit.h"
#include "SvtxHitMap.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/getClass.h>

#include <TMath.h>

#include <TF1.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TSpectrum2.h>
#include <TProfile2D.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TStopwatch.h>
#include <TH1D.h>

#include <TMatrixF.h>

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;

PHG4TPCClusterizer::PHG4TPCClusterizer(const char *name) : 
  SubsysReco(name),
  fNPhiBins(1),
  fNZBins(1),
  fGeoLayer(NULL),
  fFitW(1.0),
  fFitSumP(0.0),
  fFitSumZ(0.0),
  fFitSumP2(0.0),
  fFitSumZ2(0.0),
  fFitSumPZ(0.0),
  fFitP0(0.0),
  fFitZ0(0.0),
  fFitRangeP(1),
  fFitRangeZ(1),
  fFitRangeMP(1),
  fFitRangeMZ(1),
  fClusterCut(20), 
  fPedestal(74.4),
  fFitSizeP(0.0),
  fFitSizeZ(0),
  fShapingLead(32.0*8.0/1000.0),
  fShapingTail(48.0*8.0/1000.0),
  fMinLayer(0),
  fMaxLayer(0),
  fEnergyCut(0.0),
  fClusterWindow(2),
  fClusterZSplit(true),
  fDeconMode(false),
  fDCT(0.006),
  fDCL(0.012),
  _inv_sqrt12( 1.0/TMath::Sqrt(12) ),
  _twopi( TMath::TwoPi() ),
//zz_shaping_correction(0.0557),  // correction for 80 ns SAMPA
  zz_shaping_correction(0.0754),  // correction for 80 ns SAMPA
  fHClusterEnergy(NULL),
  fHClusterSizePP(NULL),
  fHClusterSizeZZ(NULL),
  fHClusterErrorPP(NULL),
  fHClusterErrorZZ(NULL),
  fHClusterDensity(NULL),
  fHClusterSizePP2(NULL),
  fHClusterSizeZZ2(NULL),
  fHClusterErrorPP2(NULL),
  fHClusterErrorZZ2(NULL),
  fHClusterDensity2(NULL),
  fHClusterWindowP(NULL),
  fHClusterWindowZ(NULL),
  fSW(NULL),
  fHTime(NULL),
  fSource(NULL),
  fResponse(NULL)
{
}
//===================
PHG4TPCClusterizer::~PHG4TPCClusterizer() {
  if(fHClusterEnergy) delete fHClusterEnergy;
  if(fHClusterSizePP) delete fHClusterSizePP;
  if(fHClusterSizeZZ) delete fHClusterSizeZZ;
  if(fHClusterErrorPP) delete fHClusterErrorPP;
  if(fHClusterErrorZZ) delete fHClusterErrorZZ;
  if(fHClusterDensity) delete fHClusterDensity;
  if(fHClusterSizePP2) delete fHClusterSizePP2;
  if(fHClusterSizeZZ2) delete fHClusterSizeZZ2;
  if(fHClusterErrorPP2) delete fHClusterErrorPP2;
  if(fHClusterErrorZZ2) delete fHClusterErrorZZ2;
  if(fHClusterDensity2) delete fHClusterDensity2;
}
//===================
int PHG4TPCClusterizer::wrap_phibin(int bin) {
  if(bin < 0) bin += fNPhiBins;
  if(bin >= fNPhiBins) bin -= fNPhiBins;
  return bin;
}
//===================
bool PHG4TPCClusterizer::is_local_maximum(int phi, int z) {
  bool is_max = true;

  if(fDeconMode){
    if(fSource[z][phi]<=0) return false;

    float cent_val = fSource[z][phi];
    
    for(int iz=-1; iz<=1; ++iz) {
      int cz = z + iz;
      if(cz < 0) continue; // skip edge
      if(cz >= fNZBins) continue; // skip edge
      for(int ip=-1; ip<=1; ++ip) {
	if((iz == 0) && (ip == 0)) continue; // skip center
	int cp = wrap_phibin(phi + ip);
	if(fSource[cz][cp] > cent_val) {
	  is_max = false;
	  break;
	}
      }
      if(!is_max) break;
    }
  }
  else{
    if(fAmps[z * fNPhiBins + phi] <= 0.) return false;

    float cent_val = fAmps[z*fNPhiBins + phi];
    int FitRangeZ = fFitRangeZ;
    int FitRangeP = fFitRangeP;
    if(FitRangeZ>1)FitRangeZ-=1;
    if(FitRangeP>1)FitRangeP-=1;

    for(int iz=-FitRangeZ; iz!=FitRangeZ+1; ++iz) {
      int cz = z + iz;
      if(cz < 0) continue; // skip edge
      if(cz >= fNZBins) continue; // skip edge
      for(int ip=-FitRangeP; ip!=FitRangeP+1; ++ip) {
	if((iz == 0) && (ip == 0)) continue; // skip center
	int cp = wrap_phibin(phi + ip);
	if(fAmps[cz*fNPhiBins + cp] > cent_val) {
	  is_max = false;
	  break;
	}
      }
      if(!is_max) break;
    }
  }
  return is_max;
}

void PHG4TPCClusterizer::prepare_layer(float radius){
  double m_sqrt2 = sqrt(2.0);
  TH2F *hdum = new TH2F("hdum","",fNZBins,-105.5,105.5,fNPhiBins,-3.14159,3.14159);

  double zpos = 5.0;
  //  printf("zpos: %f \n",zpos);
  {

    double phi = 0;
    double z = zpos;
    int phibin = 500;
    phi = hdum->GetYaxis()->GetBinCenter(phi);
    int zbin = hdum->GetXaxis()->FindBin(zpos);
    double zbinwidth  = hdum->GetXaxis()->GetBinWidth(2);

    double layer_depth = 1.2;
    
    double zrange = layer_depth*sqrt(zpos*zpos+radius*radius)/radius;
    
    // converting Edep to Total Number Of Electrons
    Int_t nelec = 1000;
    
    double sigmaL[2];
    double sigmaT = 0.04;
    sigmaL[0] = 32*3.0/1000;
    sigmaL[1] = 48*3.0/1000;
    double fDiffusionT = 0.0130;
    double fDiffusionL = 0.0130;
    
    double cloud_sig_rp = sqrt( fDiffusionT*fDiffusionT*(105.5 - fabs(zpos)) + sigmaT*sigmaT );
    double cloud_sig_zz[2];
    cloud_sig_zz[0] = sqrt( fDiffusionL*fDiffusionL*(105.5 - fabs(zpos)) + sigmaL[0]*sigmaL[0] );
    cloud_sig_zz[1] = sqrt( fDiffusionL*fDiffusionL*(105.5 - fabs(zpos)) + sigmaL[1]*sigmaL[1] );
    
    int nseg = zrange/cloud_sig_zz[1];// must be odd number
    if(nseg%2==0)nseg+=1;

    for(int izr =0;izr<nseg;izr++)
      {

	double zsegoffset = (izr - nseg/2) * zrange/nseg;
	int zbinseg = zbin = hdum->GetXaxis()->FindBin( z + zsegoffset);
	if(zbinseg < 0 || zbinseg >= fNZBins){continue;}
	//define the window
	
	double zdispseg = z + zsegoffset - hdum->GetXaxis()->GetBinCenter(zbin);
	int n_rp = int(3*cloud_sig_rp/(radius*hdum->GetXaxis()->GetBinWidth(2))+1);
	int n_zz = int(3*(cloud_sig_zz[0]+cloud_sig_zz[1])/(2.0*zbinwidth)+1);
	
	double cloud_sig_rp_inv = 1./cloud_sig_rp;
	double  cloud_sig_zz_inv[2]; 
	cloud_sig_zz_inv[0] = 1./cloud_sig_zz[0];
	cloud_sig_zz_inv[1] = 1./cloud_sig_zz[1];
	
	//loop over bins in window
	
	for( int iphi = -n_rp; iphi != n_rp+1; ++iphi ) {

	  double phiLim1 = 0.5*m_sqrt2*( (iphi+0.5)*radius*hdum->GetXaxis()->GetBinWidth(2))*cloud_sig_rp_inv;
	  double phiLim2 = 0.5*m_sqrt2*( (iphi-0.5)*radius*hdum->GetXaxis()->GetBinWidth(2))*cloud_sig_rp_inv;
	  double phi_integral = 0.5*( TMath::Erf(phiLim1) - TMath::Erf(phiLim2) );
	  
	  for( int iz = -n_zz; iz != n_zz+1; ++iz ) {
	    double zLim1 = 0.5*m_sqrt2*( (iz+0.5)*zbinwidth - zdispseg)*cloud_sig_zz_inv[0];
	    double zLim2 = 0.5*m_sqrt2*( (iz-0.5)*zbinwidth - zdispseg)*cloud_sig_zz_inv[0];

	    if(zLim1 > 0)
	      zLim1 =  0.5*m_sqrt2*( (iz+0.5)*zbinwidth - zdispseg )*cloud_sig_zz_inv[1];
	    if(zLim2 > 0)
	      zLim2 = 0.5*m_sqrt2*( (iz-0.5)*zbinwidth - zdispseg )*cloud_sig_zz_inv[1];

	    double z_integral = 0.5*( TMath::Erf(zLim1) - TMath::Erf(zLim2) );
	    
	    //float neffelectrons = (2000/nseg)*nelec*( phi_integral * z_integral ); 
	    float neffelectrons = (1/nseg)*nelec*( phi_integral * z_integral ); 
	    
	    if(neffelectrons < 0) continue; // skip no signals
	    
	    //Fill at cur_z_bin, cur_phi_bin neffelectrons
	    double thisPhi = hdum->GetYaxis()->GetBinCenter(phibin + iphi);
	    double thisZ = hdum->GetXaxis()->GetBinCenter(zbinseg + iz);
	    hdum->Fill(thisZ,thisPhi,neffelectrons);
	  } //iz
	} //iphi
      } // izr
  }

  if(fResponse!=NULL){
    for (Int_t i=0;i<fNZBins;i++){
      delete[] fSource[i];
      delete[] fResponse[i];
    }
    delete[] fSource;
    delete[] fResponse;
  }

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 10, 4)
  fSource = new double *[fNZBins];  
  fResponse = new double *[fNZBins];  
#else
  fSource = new float *[fNZBins];  
  fResponse = new float *[fNZBins];  
#endif

  for (Int_t i=0;i<fNZBins;i++){

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 10, 4)
    fSource[i]=new double[fNPhiBins];
    fResponse[i]=new double[fNPhiBins];
#else
    fSource[i]=new float[fNPhiBins];
    fResponse[i]=new float[fNPhiBins];
#endif
    for(Int_t j = 0;j<fNPhiBins;j++){
      fResponse[i][j] = 0;
      fSource[i][j] = 0;
    }
  }
  TH2F *htest2D = new TH2F("htest2D","",fNZBins,-105.5,105.5,fNPhiBins,-3.14159,3.14159);
  //prep response matrix
  
  //Get x range:
  TH1D *hist1d = hdum->ProjectionX("hist1d",1,fNPhiBins);
  Int_t xmin = 0; 
  Int_t xmax = fNPhiBins;

  for(Int_t bin = 1;bin<=fNPhiBins;bin++){
    float val = hist1d->GetBinContent(bin);
    if(xmin==0&&val>0)
      xmin = bin;
    if(xmax==fNPhiBins&&xmin>0&&val==0.0){
      xmax = bin-1;
    }
  }

  hist1d = hdum->ProjectionY("hist1d",1,fNZBins);
  Int_t ymin = 0; 
  Int_t ymax = fNZBins;

  for(Int_t bin = 1;bin<=fNZBins;bin++){
    float val = hist1d->GetBinContent(bin);
    if(ymin==0&&val>0)
      ymin = bin;
    if(ymax==fNZBins&&ymin>0&&val==0.0){
      ymax = bin-1;
    }
  }

  Int_t iy = 0;
  for(Int_t ybin = ymin;ybin<=ymax;ybin++){
    TH1D *hist1d = hdum->ProjectionX("hist1d",ybin,ybin);
    Int_t ix = 0;
    for(Int_t xbin = xmin;xbin<=xmax;xbin++){

      fResponse[ix][iy] = hist1d->GetBinContent(xbin);
      htest2D->SetBinContent(ix,iy,hist1d->GetBinContent(xbin));
      ix++;
    }
    iy++;
    //    delete hist1d;
  }

  iy = 0;
  for(Int_t ybin = ymin;ybin<=ymax;ybin++){
    TH1D *hist1d = hdum->ProjectionX("hist1d",ybin,ybin);
    Int_t ix = 0;
    for(Int_t xbin = xmin;xbin<=xmax;xbin++){
      fResponse[ix][iy] += hist1d->GetBinContent(xbin);
      htest2D->SetBinContent(ix,iy,fResponse[ix][iy]);
      ix++;
    }
    iy++;
    //    delete hist1d;
  }
  delete hdum;
  delete htest2D;
}

//===================
void PHG4TPCClusterizer::deconvolution() {
  TSpectrum2 pfinder;
  pfinder.Deconvolution(fSource,fResponse,fNZBins,fNPhiBins,10,10,5);
  return;
}

//===================
void PHG4TPCClusterizer::find_phi_range(int zbin, int phibin, int phimax, float peak, int& phiup, int& phidown){

  
  for(int ip=0; ip<=fFitRangeP; ++ip) {
    int cp = wrap_phibin(phibin + ip);
    int bin = zbin * fNPhiBins + cp;
    if(fAmps[bin] == 0){
      phiup = ip;
      break; // skip small (include empty)
    }
    int bin1 = zbin * fNPhiBins + wrap_phibin(phibin + ip + 1);
    int bin2 = zbin * fNPhiBins + wrap_phibin(phibin + ip + 2);
    int bin3 = zbin * fNPhiBins + wrap_phibin(phibin + ip + 3);
    if(ip<fFitRangeP){
      if(fAmps[bin]+fAmps[bin1] < fAmps[bin2]+fAmps[bin3]){//rising again
	phiup = ip+1;
	break;
      }
    }
  }
  for(int ip=0; ip<=fFitRangeP; ++ip) {
    int cp = wrap_phibin(phibin - ip);
    int bin = zbin * fNPhiBins + cp;
    if(fAmps[bin] == 0){
      phidown = ip;
      break; // skip small (include empty)
    }
    int bin1 = zbin * fNPhiBins + wrap_phibin(phibin - ip - 1);
    int bin2 = zbin * fNPhiBins + wrap_phibin(phibin - ip - 2);
    int bin3 = zbin * fNPhiBins + wrap_phibin(phibin - ip - 3);
    if(ip<fFitRangeP){
      if(fAmps[bin]+fAmps[bin1] < fAmps[bin2]+fAmps[bin3]){//rising again
	phidown = ip+1;
	break;
      }
    }
  }
  /*
    if(phiup<1)phiup = 1;
    if(phidown<1)phidown = 1;
  */
}

//===================
void PHG4TPCClusterizer::find_z_range(int zbin, int phibin, int zmax, float peak, int& zup, int& zdown){

  zup   = fFitRangeZ;
  zdown = fFitRangeZ;
    
  for(int iz=0; iz< zmax; iz++)
    {
      int cz = zbin + iz;
      if(cz < 0) continue; // truncate edge
      if(cz >= fNZBins-3) continue; // truncate edge
      
      // consider only the peak bin in phi when searching for Z limit     
      int cp = wrap_phibin(phibin);
      int bin = cz * fNPhiBins + cp;
      int bin1 = (cz+1) * fNPhiBins + cp;
      int bin2 = (cz+2) * fNPhiBins + cp;
      int bin3 = (cz+3) * fNPhiBins + cp;
      if(fAmps[bin] == 0) {
	zup = iz;
	if(Verbosity() > 1000) cout << " failed threshold cut, set izup to " << zup << endl;
	break;
      }
      if(fClusterZSplit){
	//check local minima and break at minimum.
	if(iz<zmax-4){//make sure we stay clear from the edge
	  if(fAmps[bin]+fAmps[bin1] < fAmps[bin2]+fAmps[bin3]){//rising again
	    zup = iz+1;
	    break;
	  }
	}
      }
      
    }

  for(int iz=0; iz< zmax; iz++)
    {
      int cz = zbin - iz;
      if(cz <= 3) continue; // truncate edge
      if(cz >= fNZBins) continue; // truncate edge
      
      int cp = wrap_phibin(phibin);
      int bin = cz * fNPhiBins + cp;
      int bin1 = (cz-1) * fNPhiBins + cp;
      int bin2 = (cz-2) * fNPhiBins + cp;
      int bin3 = (cz-3) * fNPhiBins + cp;
      if(fAmps[bin] == 0) {
	zdown = iz;
	if(Verbosity() > 1000) cout << " failed threshold cut, set izdown to " << zdown << endl;
	break;
      }
      if(fClusterZSplit){
	if(iz<zmax-4){//make sure we stay clear from the edge
	  if(fAmps[bin]+fAmps[bin1] < fAmps[bin2]+fAmps[bin3]){//rising again
	    zdown = iz+1;
	    break;
	  }
	}
      }
    }
  /*  if(zup<2)zup = 2;
      if(zdown<2)zdown = 2;
  */
}

//===================
void PHG4TPCClusterizer::fit(int pbin, int zbin, int& nhits_tot) {
  float peak = fAmps[zbin * fNPhiBins + pbin];
  fFitW = 0.0;
  fFitP0 = fGeoLayer->get_phicenter( pbin );
  fFitZ0 = fGeoLayer->get_zcenter( zbin );
  fFitSumP = 0;
  fFitSumZ = 0;
  fFitSumP2 = 0;
  fFitSumZ2 = 0;
  fFitSumPZ = 0;
  fFitSizeP = 0;
  fFitSizeZ = 0;
  // Here we have to allow for highly angled tracks, where the Z range can be much larger than the Z diffusion
  // Check that we do not have to increase this range
  // Start at the peak and go out in both directions until the yield falls below threshold

  int izmost = 30; // don't look more than izmost Z bins in each direction
  int ipmost   = fFitRangeP;
  int ipup   = fFitRangeP;
  int ipdown = fFitRangeP;

  //Get Phi dimension of cluster at maximum
  find_phi_range(zbin,pbin,ipmost,peak,ipup,ipdown);

  //Get Z Dimension per phi slice
  int izup[2*fFitRangeP+1];
  int izdown[2*fFitRangeP+1];
  int izupmax = -1;
  int izdownmax = -1;

  for(int iz= 0;iz<2*fFitRangeP+1;iz++){
    izup[iz] = -1;
    izdown[iz] = -1;
  }

  for(int ip  = -1*ipdown;ip<=ipup;ip++){
    find_z_range(zbin,pbin+ip,izmost,peak,izup[ip+fFitRangeP],izdown[ip+fFitRangeP]);
    if(izup[ip+fFitRangeP]>izupmax)
      izupmax = izup[ip+fFitRangeP];
    if(izdown[ip+fFitRangeP]>izdownmax)
      izdownmax = izdown[ip+fFitRangeP];
  }

  if(Verbosity()>1000) 
    std::cout << "max " << fAmps[zbin*fNPhiBins+pbin] << std::endl;
  if(Verbosity()>1000) 
    std::cout << "izdown " << izdown << " izup " << izup << std::endl;

  for(int iz=-izdownmax; iz<izupmax; ++iz) {
    int cz = zbin + iz;
    if(cz < 0) continue; // truncate edge
    if(cz >= fNZBins) continue; // truncate edge
    bool used = false;
    int nphis = 0;

    for(int ip=-ipdown; ip<=ipup; ++ip) {
      if(iz<-1*izdown[ip+fFitRangeP])continue;
      if(iz>izup[ip+fFitRangeP])continue;
      int cp = wrap_phibin(pbin + ip);
      int bin = cz * fNPhiBins + cp;
      if(Verbosity()>1000) {
	std::cout << Form("%.2f | ",fAmps[bin]);
	if(ip==fFitRangeP) std::cout << std::endl;
      }
      if(fAmps[bin] == 0) continue; // skip small (include empty)
      used = true;
      nphis++;
      float ee = fAmps[bin];
      float dz = fGeoLayer->get_zcenter(cz) - fFitZ0;
      float dp = fGeoLayer->get_phicenter(cp) - fFitP0;
      if(dp<-3.0) dp = dp + _twopi; // binphi has been reduced;
      if(dp>+3.0) dp = dp - _twopi; // binphi has been increased;
      fFitW += ee;
      fFitSumP += ee*dp;
      fFitSumZ += ee*dz;
      fFitSumP2 += ee*dp*dp;
      fFitSumZ2 += ee*dz*dz;
      fFitSumPZ += ee*dp*dz;
      // add this cell to the list of contributing cells for this cluster
      fCellz.push_back(cz);
      fCellphi.push_back(cp);
      nhits_tot -= 1; //taken
      fNHitsPerZ[cz] -= 1; //taken
      fAmps[bin] = 0.; //removed
    }
    if(used) fFitSizeZ++;
    fFitSizeP = TMath::Max(fFitSizeP,float(nphis));
  }
  if(Verbosity()>1000) {
    std::cout << " FIT | phi " << fit_p_mean() << " from " << fFitP0;
    std::cout << " | z " << fit_z_mean() << " from " << fFitZ0 << std::endl;
  }
}
//===================
int PHG4TPCClusterizer::InitRun(PHCompositeNode* topNode) {
  if(Verbosity()>1) {
    fHClusterEnergy = new TH1F("CLUSTER_Energy","CLUSTER_Energy",1000,0,1000);
    fHClusterDensity = new TProfile2D("CLUSTER_Density","CLUSTER_Density;LayerNo;ZZ;<E>",50,-0.5,49.5,220,-110,+110);
    fHClusterSizePP = new TProfile2D("CLUSTER_SizePP","CLUSTER_SizePP;LayerNo;ZZ;<rphisize>",50,-0.5,49.5,220,-110,+110);
    fHClusterSizeZZ = new TProfile2D("CLUSTER_SizeZZ","CLUSTER_SizeZZ;LayerNo;ZZ;<zsize>",50,-0.5,49.5,220,-110,+110);
    fHClusterErrorPP = new TProfile2D("CLUSTER_ErrorPP","CLUSTER_ErrorPP;LayerNo;ZZ;<rphierror>",50,-0.5,49.5,220,-110,+110);
    fHClusterErrorZZ = new TProfile2D("CLUSTER_ErrorZZ","CLUSTER_ErrorZZ;LayerN0;ZZ;<zerror>",50,-0.5,49.5,220,-110,+110);
    fHClusterDensity2 = new TProfile2D("CLUSTER_Density2","CLUSTER_Density2;Phi;ZZ;<E>",100,-TMath::Pi(),TMath::Pi(),220,-110,+110);
    fHClusterSizePP2 = new TProfile2D("CLUSTER_SizePP2","CLUSTER_SizePP2;Phi;ZZ;<rphisize>",100,-TMath::Pi(),TMath::Pi(),220,-110,+110);
    fHClusterSizeZZ2 = new TProfile2D("CLUSTER_SizeZZ2","CLUSTER_SizeZZ2;Phi;ZZ;<zsize>",100,-TMath::Pi(),TMath::Pi(),220,-110,+110);
    fHClusterErrorPP2 = new TProfile2D("CLUSTER_ErrorP2P","CLUSTER_ErrorPP2;Phi;ZZ;<rphierror>",100,-TMath::Pi(),TMath::Pi(),220,-110,+110);
    fHClusterErrorZZ2 = new TProfile2D("CLUSTER_ErrorZZ2","CLUSTER_ErrorZZ2;Phi;ZZ;<zerror>",100,-TMath::Pi(),TMath::Pi(),220,-110,+110);
    fHClusterWindowP = new TProfile2D("CLUSTER_WindowP","CLUSTER_WindowP;LayerNo;ZZ",50,-0.5,49.5,220,-110,+110);
    fHClusterWindowZ = new TProfile2D("CLUSTER_WindowZ","CLUSTER_WindowZ;LayerNo;ZZ",50,-0.5,49.5,220,-110,+110);
    fSW = new TStopwatch();
    fHTime = new TH1F("TIME_CLUSTER","CLUSTER_TIME;sec per event",1000,0,50);
    Fun4AllServer *se = Fun4AllServer::instance();
    se->registerHisto( fHClusterEnergy );
    se->registerHisto( fHClusterDensity );
    se->registerHisto( fHClusterSizePP );
    se->registerHisto( fHClusterSizeZZ );
    se->registerHisto( fHClusterErrorPP );
    se->registerHisto( fHClusterErrorZZ );
    se->registerHisto( fHClusterDensity2 );
    se->registerHisto( fHClusterSizePP2 );
    se->registerHisto( fHClusterSizeZZ2 );
    se->registerHisto( fHClusterErrorPP2 );
    se->registerHisto( fHClusterErrorZZ2 );
    se->registerHisto( fHClusterWindowP );
    se->registerHisto( fHClusterWindowZ );
    se->registerHisto( fHTime );
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
//===================
void PHG4TPCClusterizer::reset() {}
//===================
int PHG4TPCClusterizer::process_event(PHCompositeNode* topNode) {
  if(Verbosity()>1000) std::cout << "PHG4TPCClusterizer::Process_Event" << std::endl;
  if(Verbosity()>1) {
    fSW->Reset();
    fSW->Start();
  }
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode) {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator iter_dst(dstNode);
  SvtxHitMap* hits = findNode::getClass<SvtxHitMap>(dstNode, "SvtxHitMap");
  if (!hits) {
    cout << PHWHERE << "ERROR: Can't find node SvtxHitMap" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHCompositeNode* svxNode = dynamic_cast<PHCompositeNode*>(iter_dst.findFirst("PHCompositeNode", "SVTX"));
  if (!svxNode) {
    svxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svxNode);
  }
  SvtxClusterMap* svxclusters = findNode::getClass<SvtxClusterMap>(dstNode, "SvtxClusterMap");
  if (!svxclusters) {
    svxclusters = new SvtxClusterMap_v1();
    PHIODataNode<PHObject>* SvtxClusterMapNode =
        new PHIODataNode<PHObject>(svxclusters, "SvtxClusterMap", "PHObject");
    svxNode->addNode(SvtxClusterMapNode);
  }
  PHG4CylinderCellGeomContainer* geom_container =
    findNode::getClass<PHG4CylinderCellGeomContainer>(topNode,"CYLINDERCELLGEOM_SVTX");
  if (!geom_container) {
    std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHG4CellContainer* cells =  findNode::getClass<PHG4CellContainer>(dstNode,"G4CELL_SVTX");
  if (!cells) {
    std::cout << PHWHERE << "ERROR: Can't find node G4CELL_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // ==> making layer_sorted, initially an empty vector of hits for each layer 
  std::vector<std::vector<const SvtxHit*> > layer_sorted;
  PHG4CylinderCellGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for(PHG4CylinderCellGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter) {
    //cout << "Layer " << layeriter->second->get_layer() << endl;
    if( (unsigned int) layeriter->second->get_layer() < fMinLayer) {
      if(Verbosity()>1000) std::cout << "Skipping layer " << layeriter->second->get_layer() << std::endl;
      continue;
    }
    layer_sorted.push_back(std::vector<const SvtxHit*>());
  }
  for(SvtxHitMap::Iter iter = hits->begin(); iter != hits->end(); ++iter) {
    SvtxHit* hit = iter->second;
    //cout << " hit_get_layer = " << hit->get_layer() << " fMinLayer " << fMinLayer << endl;
    if( (unsigned int) hit->get_layer() < fMinLayer) continue;
    layer_sorted[hit->get_layer() - fMinLayer].push_back(hit);
  }
  // <==

  for(PHG4CylinderCellGeomContainer::ConstIterator layeriter = layerrange.first;
      layeriter != layerrange.second; ++layeriter) {
    unsigned int layer = (unsigned int)layeriter->second->get_layer();
    if (layer < fMinLayer) continue;
    if (layer > fMaxLayer) continue;
    fGeoLayer = geom_container->GetLayerCellGeom(layer);
    fNPhiBins = layeriter->second->get_phibins();
    fNZBins = layeriter->second->get_zbins();
    if(Verbosity()>0) {
      std::cout << "************ Layer " << layer;
      std::cout << " fNPhiBins " << fNPhiBins;
      std::cout << " fNZBins " << fNZBins;
      std::cout << std::endl;
    }
    fNHitsPerZ.clear();
    fNHitsPerZ.assign(fNZBins, 0);
    fAmps.clear();
    fAmps.assign(fNPhiBins * fNZBins, 0.);
    fCellIDs.clear();
    fCellIDs.assign(fNPhiBins * fNZBins, 0);
    // ==>unpacking information
    for(unsigned int i = 0; i < layer_sorted[layer - fMinLayer].size(); ++i) {
      const SvtxHit* hit = layer_sorted[layer - fMinLayer][i];
      if(hit->get_e() <= 0.) continue;
      if(Verbosity()>2000) 
	std::cout << hit->get_cellid();
      PHG4Cell* cell = cells->findCell(hit->get_cellid()); //not needed once geofixed
      int phibin = PHG4CellDefs::SizeBinning::get_phibin(cell->get_cellid());//cell->get_binphi();
      int zbin = PHG4CellDefs::SizeBinning::get_zbin(cell->get_cellid());//cell->get_binz();
      if(Verbosity()>0) std::cout << " phibin " << phibin << " zbin " << zbin << " z " << fGeoLayer->get_zcenter( zbin ) << " adc " << hit->get_adc() - fPedestal << std::endl;
      fNHitsPerZ[zbin] += 1;
      fAmps[zbin * fNPhiBins + phibin] += hit->get_adc() - fPedestal;  // subtract pedestal in ADC counts, determined elsewhere
      if(fAmps[zbin * fNPhiBins + phibin] < 0)  fAmps[zbin * fNPhiBins + phibin]  = 0;  // our simple clustering algorithm does not handle negative bins well
      fCellIDs[zbin * fNPhiBins + phibin] = hit->get_id();
      //if(Verbosity() > 100)
      if(Verbosity() > 100)
	//if(layer == 47) 
	  {
	    cout << "Clusterizer: adding input SvtxHit " <<  hit->get_id() << endl;;
	    //hit->identify();
	    cout << "      layer " << layer << " zbin " << zbin << " phibin " << phibin << " cellid " << hit->get_cellid() << " adc " << fAmps[zbin * fNPhiBins + phibin] << endl;
	  }
    }
    if(fDeconMode){
      cout << "deconvoluting layer: " << layer << endl;
      deconvolution();
      cout << "done deconvoluting" << endl;
    }
    int nhits_tot = 0;
    for(int zbin = 0; zbin!=fNZBins; ++zbin)
      nhits_tot += fNHitsPerZ[zbin];
    if(Verbosity()>0) 
      std::cout << " nhits_tot " << nhits_tot << std::endl;

    float stepz = fGeoLayer->get_zstep();
    float stepp = fGeoLayer->get_phistep();

    //while(nhits_tot > 2) {
    if(Verbosity()>0) 
      std::cout << " => nhits_tot " << nhits_tot << std::endl;
      for(int zbin = 0; zbin!=fNZBins; ++zbin) {
        if(fNHitsPerZ[zbin]<=0) continue;
	float abszbincenter = TMath::Abs(fGeoLayer->get_zcenter( zbin ));
	float sigmaZ = TMath::Sqrt(pow((fShapingTail), 2) + fDCL*fDCL*(105.5-abszbincenter));  // shaping time + drift diffusion, used only to calculate FitRangeZ
	float TPC_padgeo_sigma = 0.04;  // 0.4 mm (from Tom)
	float sigmaP = TMath::Sqrt(pow(TPC_padgeo_sigma, 2) + fDCT*fDCT*(105.5-abszbincenter)); // readout geometry + drift diffusion, used only to calculate FitRangeP
	fFitRangeZ = int( fClusterWindow*sigmaZ/stepz + 1);
	if(fFitRangeZ<1) fFitRangeZ = 1; // should never happen
	if(Verbosity() > 2000) 
	  cout << " sigmaZ " << sigmaZ << " fFitRangeZ " << fFitRangeZ << " sigmaP " << sigmaP << endl;

        for(int phibin = 0; phibin!=fNPhiBins; ++phibin) {
	  float radius = fGeoLayer->get_radius(); // returns center of layer
	  fFitRangeP = int( fClusterWindow*sigmaP/(radius*stepp) + 1);

	  if(fFitRangeP<1) fFitRangeP = 1;
	  if(fFitRangeP>fFitRangeMP) fFitRangeP = fFitRangeMP;
          if(!is_local_maximum(phibin, zbin)) continue;
	  if(Verbosity()>2000){
	    std::cout << " maxima found in window rphi z " << fFitRangeP << " " << fFitRangeZ << std::endl;
	    cout << " phibin: "  << phibin << " zbin: " << zbin 
		 << endl;
	  }
	  fCellz.clear();
	  fCellphi.clear();
          fit(phibin,zbin,nhits_tot);
          if(fFitW < fClusterCut) continue; // ignore this cluster
          SvtxCluster_v1 clus;
          clus.set_layer(layer);
 	  clus.set_e( fFitW);
 	  clus.set_adc( fFitW);
	  float phi = fit_p_mean();
	  float pp = radius*phi;
	  float zz = fit_z_mean();
	  // Correction for bias in electron z position due to asymmeteric SAMPA shaping
	  float zz_raw = zz;
	  if(zz < 0)
	    zz -= zz_shaping_correction;
	  else
	    zz += zz_shaping_correction;
	  float pp_err = radius * fGeoLayer->get_phistep() * _inv_sqrt12;
	  float zz_err = fGeoLayer->get_zstep() * _inv_sqrt12;

	  // fit_p_cov is essentially the weighted mean of dphi^2. The error is then:
	  // e_phi = sigma_dphi/sqrt(N) = sqrt( sigma_dphi^2 / N )  -- where N is the number of samples of the distribution with standard deviation sigma_dphi
	  //    - N is the number of electrons that drift to the readout
	  // We have to convert fFitW (sum of adc units for all bins in the cluster) to number of ionization electrons N
	  // Conversion gain is 20 mV/fC - relates total charge collected on pad to PEAK voltage out of ADC
	  // GEM gain is assumed to be 2000
	  // To get equivalent charge per Z bin, so that summing ADC input voltage over all Z bins returns total input charge, divide voltages by 2.4 for 80 ns SAMPA
	  // Equivalent charge per Z bin is then  (ADU x 2200 mV / 1024) / 2.4 x (1/20) fC/mV x (1/1.6e-04) electrons/fC x (1/2000) = ADU x 0.14
	  if(fFitSizeP>1) pp_err = radius * TMath::Sqrt( fit_p_cov()/(fFitW*0.14) );
	  if(fFitSizeZ>1) zz_err = TMath::Sqrt( fit_z_cov()/(fFitW*0.14) );
	  float pp_size = radius*fFitSizeP*fGeoLayer->get_phistep();
	  float zz_size = fFitSizeZ*fGeoLayer->get_zstep();

	  if(Verbosity() > 100)
	    //if(layer == 47) 
	      {
		cout << endl << " clusterizer: layer " << layer << " fFitW " << fFitW << " number of primary electrons (adc) = " << fFitW * 0.14 
		     << " zz_raw " << zz_raw << " zz " << zz
		     << " zz_size " << zz_size << " fFitsizeZ " << fFitSizeZ << " phi " << phi << endl;
		cout << "       zz_err " << zz_err << " fit_z_cov " << fit_z_cov() << endl;

	      }

	  if(Verbosity()>1) {
	    fHClusterEnergy->Fill(fFitW);
	    fHClusterDensity->Fill(layer,zz,fFitW);
	    fHClusterSizePP->Fill(layer,zz,pp_size);
	    fHClusterSizeZZ->Fill(layer,zz,zz_size);
	    fHClusterErrorPP->Fill(layer,zz,pp_err);
	    fHClusterErrorZZ->Fill(layer,zz,zz_err);
	    fHClusterDensity2->Fill(phi,zz,fFitW);
	    fHClusterSizePP2->Fill(phi,zz,pp_size);
	    fHClusterSizeZZ2->Fill(phi,zz,zz_size);
	    fHClusterErrorPP2->Fill(phi,zz,pp_err);
	    fHClusterErrorZZ2->Fill(phi,zz,zz_err);
	    fHClusterWindowP->Fill(layer,zz,fFitRangeP);
	    fHClusterWindowZ->Fill(layer,zz,fFitRangeZ);
	  }
	  if(Verbosity()>2000)
	    std::cout << " cluster fitted " << std::endl;
	    if(Verbosity()>1000)
	    {
	    std::cout << " | rad " << radius;
	    std::cout << " | size rp z " << pp_size << " " << zz_size;
	    std::cout << " | error_rphi error_z " << pp_err << " " << zz_err << std::endl;
	    std::cout << " | sgn " << fFitW;
	    std::cout << " | rphi z " << pp << " " << zz;
	    std::cout << " | phi_sig z_sig phiz_cov " << TMath::Sqrt(fit_p_cov()) << " " << TMath::Sqrt(fit_z_cov()) << " " << fit_pz_cov();
	    std::cout << std::endl;
	    }
          clus.set_position(0, radius*TMath::Cos( phi ) );
          clus.set_position(1, radius*TMath::Sin( phi ) );
          clus.set_position(2, zz);
	  for(unsigned int i=0;i<fCellz.size();i++)
	    {
	      if(Verbosity() > 10)
		//if(layer == 47)
		cout  << "   Fitted cluster contains SvtxHit " << fCellIDs[ fCellz[i] * fNPhiBins + fCellphi[i] ] << endl;
	      clus.insert_hit( fCellIDs[ fCellz[i] * fNPhiBins + fCellphi[i] ]);
	    }

	  TMatrixF DIM(3,3);
	  DIM[0][0] = 0.0;
	  DIM[0][1] = 0.0;
	  DIM[0][2] = 0.0;
	  DIM[1][0] = 0.0;
	  DIM[1][1] = pp_size/radius*pp_size/radius; //cluster_v1 expects polar coordinates covariance
	  DIM[1][2] = 0.0;
	  DIM[2][0] = 0.0;
	  DIM[2][1] = 0.0;
	  DIM[2][2] = zz_size*zz_size;

	  TMatrixF ERR(3,3);
	  ERR[0][0] = 0.0;
	  ERR[0][1] = 0.0;
	  ERR[0][2] = 0.0;
	  ERR[1][0] = 0.0;
	  ERR[1][1] = pp_err*pp_err; //cluster_v1 expects rad, arc, z as elementsof covariance
	  ERR[1][2] = 0.0;
	  ERR[2][0] = 0.0;
	  ERR[2][1] = 0.0;
	  ERR[2][2] = zz_err*zz_err;

	  TMatrixF ROT(3,3);
	  ROT[0][0] = cos(phi);
	  ROT[0][1] = -sin(phi);
	  ROT[0][2] = 0.0;
	  ROT[1][0] = sin(phi);
	  ROT[1][1] = cos(phi);
	  ROT[1][2] = 0.0;
	  ROT[2][0] = 0.0;
	  ROT[2][1] = 0.0;
	  ROT[2][2] = 1.0;

	  TMatrixF ROT_T(3,3);
	  ROT_T.Transpose(ROT);
      
	  TMatrixF COVAR_DIM(3,3);
	  COVAR_DIM = ROT * DIM * ROT_T;

	  clus.set_size( 0 , 0 , COVAR_DIM[0][0] );
	  clus.set_size( 0 , 1 , COVAR_DIM[0][1] );
	  clus.set_size( 0 , 2 , COVAR_DIM[0][2] );
	  clus.set_size( 1 , 0 , COVAR_DIM[1][0] );
	  clus.set_size( 1 , 1 , COVAR_DIM[1][1] );
	  clus.set_size( 1 , 2 , COVAR_DIM[1][2] );
	  clus.set_size( 2 , 0 , COVAR_DIM[2][0] );
	  clus.set_size( 2 , 1 , COVAR_DIM[2][1] );
	  clus.set_size( 2 , 2 , COVAR_DIM[2][2] );
	  //cout << " covar_dim[2][2] = " <<  COVAR_DIM[2][2] << endl;

	  TMatrixF COVAR_ERR(3,3);
	  COVAR_ERR = ROT * ERR * ROT_T;
	    
	  clus.set_error( 0 , 0 , COVAR_ERR[0][0] );
	  clus.set_error( 0 , 1 , COVAR_ERR[0][1] );
	  clus.set_error( 0 , 2 , COVAR_ERR[0][2] );
	  clus.set_error( 1 , 0 , COVAR_ERR[1][0] );
	  clus.set_error( 1 , 1 , COVAR_ERR[1][1] );
	  clus.set_error( 1 , 2 , COVAR_ERR[1][2] );
	  clus.set_error( 2 , 0 , COVAR_ERR[2][0] );
	  clus.set_error( 2 , 1 , COVAR_ERR[2][1] );
	  clus.set_error( 2 , 2 , COVAR_ERR[2][2] );
	  //cout << " covar_err[2][2] = " <<  COVAR_ERR[2][2] << endl;

	  /*
	  not do this yet, first clear cluster class
	  clus.set_size( 0 , 0 , xx_size );
	  clus.set_size( 0 , 1 , 0.0 );
	  clus.set_size( 0 , 2 , 0.0 );
	  clus.set_size( 1 , 0 , 0.0 );
	  clus.set_size( 1 , 1 , yy_size );
	  clus.set_size( 1 , 2 , 0.0 );
	  clus.set_size( 2 , 0 , 0.0 );
	  clus.set_size( 2 , 1 , 0.0 );
	  clus.set_size( 2 , 2 , zz_size );
	  clus.set_error( 0 , 0 , xx_err );
	  clus.set_error( 0 , 1 , 0.0 );
	  clus.set_error( 0 , 2 , 0.0 );
	  clus.set_error( 1 , 0 , 0.0 );
	  clus.set_error( 1 , 1 , yy_err );
	  clus.set_error( 1 , 2 , 0.0 );
	  clus.set_error( 2 , 0 , 0.0 );
	  clus.set_error( 2 , 1 , 0.0 );
	  clus.set_error( 2 , 2 , zz_err );
	  */
	  svxclusters->insert(&clus);
	  //}
      }
    }
  }
  reset();
  if(Verbosity()>1000) std::cout << "PHG4TPCClusterizer::Process_Event DONE" << std::endl;
  if(Verbosity()>1) fHTime->Fill( fSW->RealTime() );
  return Fun4AllReturnCodes::EVENT_OK;
}
