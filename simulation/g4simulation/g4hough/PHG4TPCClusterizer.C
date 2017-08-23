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
  fFitEnergyThreshold(0.05),
  fFitSizeP(0.0),
  fFitSizeZ(0),
  fShapingLead(32.0*6.0/1000.0),
  fShapingTail(48.0*6.0/1000.0),
  fMinLayer(0),
  fMaxLayer(0),
  fEnergyCut(0.1),
  fDCT(0.006),
  fDCL(0.012),
  _inv_sqrt12( 1.0/TMath::Sqrt(12) ),
  _twopi( TMath::TwoPi() ),
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
  fHTime(NULL)
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
  if(fAmps[z * fNPhiBins + phi] <= 0.) return false;
  float cent_val = fAmps[z*fNPhiBins + phi];
  bool is_max = true;
  for(int iz=-fFitRangeZ; iz!=fFitRangeZ+1; ++iz) {
    int cz = z + iz;
    if(cz < 0) continue; // skip edge
    if(cz >= fNZBins) continue; // skip edge
    for(int ip=-fFitRangeP; ip!=fFitRangeP+1; ++ip) {
      if((iz == 0) && (ip == 0)) continue; // skip center
      int cp = wrap_phibin(phi + ip);
      if(fAmps[cz*fNPhiBins + cp] > cent_val) {
        is_max = false;
        break;
      }
    }
    if(!is_max) break;
  }
  return is_max;
}
//===================
void PHG4TPCClusterizer::fit(int pbin, int zbin, int& nhits_tot) {
  float peak = fAmps[zbin * fNPhiBins + pbin];
  fFitW = peak;
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

  int izup = fFitRangeZ;  
  for(int iz=0; iz< izmost; iz++)
    {
      int cz = zbin + iz;
      if(cz < 0) continue; // truncate edge
      if(cz >= fNZBins) continue; // truncate edge
      
      // consider only the peak bin in phi when searching for Z limit     
      int cp = wrap_phibin(pbin);
      int bin = cz * fNPhiBins + cp;
      if(fAmps[bin] < fFitEnergyThreshold*peak) 
	{
	  izup = iz;
	  if(verbosity > 1000) cout << " failed threshold cut, set izup to " << izup << endl;
	  break;
	}
    }

  int izdown = fFitRangeZ;
  for(int iz=0; iz< izmost; iz++)
    {
      int cz = zbin - iz;
      if(cz < 0) continue; // truncate edge
      if(cz >= fNZBins) continue; // truncate edge
      
      int cp = wrap_phibin(pbin);
      int bin = cz * fNPhiBins + cp;
      if(fAmps[bin] < fFitEnergyThreshold*peak) 
	{
	  izdown = iz;
	  if(verbosity > 1000) cout << " failed threshold cut, set izdown to " << izdown << endl;
	  break;
	}
    }
  
  if(verbosity>1000) std::cout << "max " << fAmps[zbin*fNPhiBins+pbin] << std::endl;
  if(verbosity>1000) std::cout << "izdown " << izdown << " izup " << izup << std::endl;

  for(int iz=-izdown; iz!=izup; ++iz) {
    int cz = zbin + iz;
    if(cz < 0) continue; // truncate edge
    if(cz >= fNZBins) continue; // truncate edge
    bool used = false;
    int nphis = 0;
    for(int ip=-fFitRangeP; ip!=fFitRangeP+1; ++ip) {
      int cp = wrap_phibin(pbin + ip);
      int bin = cz * fNPhiBins + cp;
      if(verbosity>1000) {
	std::cout << Form("%.2f | ",fAmps[bin]);
	if(ip==fFitRangeP) std::cout << std::endl;
      }
      if(fAmps[bin] < fFitEnergyThreshold*peak) continue; // skip small (include empty)
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
      nhits_tot -= 1; //taken
      fNHitsPerZ[cz] -= 1; //taken
      fAmps[bin] = 0.; //removed
    }
    if(used) fFitSizeZ++;
    fFitSizeP = TMath::Max(fFitSizeP,float(nphis));
  }
  if(verbosity>1000) {
    std::cout << " FIT | phi " << fit_p_mean() << " from " << fFitP0;
    std::cout << " | z " << fit_z_mean() << " from " << fFitZ0 << std::endl;
  }
}
//===================
int PHG4TPCClusterizer::InitRun(PHCompositeNode* topNode) {
  if(verbosity>1) {
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
  if(verbosity>1000) std::cout << "PHG4TPCClusterizer::Process_Event" << std::endl;
  if(verbosity>1) {
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

  // ==> making layer_sorted
  std::vector<std::vector<const SvtxHit*> > layer_sorted;
  PHG4CylinderCellGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for(PHG4CylinderCellGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter) {
    if( (unsigned int) layeriter->second->get_layer() < fMinLayer) {
      if(verbosity>1000) std::cout << "Skipping layer " << layeriter->second->get_layer() << std::endl;
      continue;
    }
    layer_sorted.push_back(std::vector<const SvtxHit*>());
  }
  for(SvtxHitMap::Iter iter = hits->begin(); iter != hits->end(); ++iter) {
    SvtxHit* hit = iter->second;
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
    if(verbosity>2000) {
      std::cout << "Layer " << layer;
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
      if(verbosity>2000) std::cout << hit->get_cellid();
      PHG4Cell* cell = cells->findCell(hit->get_cellid()); //not needed once geofixed
      int phibin = PHG4CellDefs::SizeBinning::get_phibin(cell->get_cellid());//cell->get_binphi();
      int zbin = PHG4CellDefs::SizeBinning::get_zbin(cell->get_cellid());//cell->get_binz();
      if(verbosity>2000) std::cout << " phibin " << phibin << " zbin " << zbin << " energy " << hit->get_e() << std::endl;
      fNHitsPerZ[zbin] += 1;
      fAmps[zbin * fNPhiBins + phibin] += hit->get_e();
      fCellIDs[zbin * fNPhiBins + phibin] = hit->get_id();
    }
    int nhits_tot = 0;
    for(int zbin = 0; zbin!=fNZBins; ++zbin)
      nhits_tot += fNHitsPerZ[zbin];
    if(verbosity>2000) std::cout << " nhits_tot " << nhits_tot << std::endl;
    // <==
    float stepz = fGeoLayer->get_zstep();
    float stepp = fGeoLayer->get_phistep();
    //std::cout << "STEPP " << stepp << " STEPZ " << stepz << std::endl;
    while(nhits_tot > 0) {
      if(verbosity>2000) std::cout << " => nhits_tot " << nhits_tot << std::endl;
      for(int zbin = 0; zbin!=fNZBins; ++zbin) {
        if(fNHitsPerZ[zbin]<=0) continue;
	float abszbincenter = TMath::Abs(fGeoLayer->get_zcenter( zbin ));
	float sigmaZ = TMath::Sqrt(pow((fShapingTail), 2) + fDCL*fDCL*(105.5-abszbincenter));  // shaping time + drift diffusion, used only to calculate FitRangZ
	float TPC_padgeo_sigma = 0.04;  // 0.4 mm (from Tom)
	float sigmaP = TMath::Sqrt(pow(TPC_padgeo_sigma, 2) + fDCT*fDCT*(105.5-abszbincenter)); // readout geometry + drift diffusion, used only to calculate FitRangeP
	fFitRangeZ = int( 3.0*sigmaZ/stepz + 1);
	if(fFitRangeZ<1) fFitRangeZ = 1; // should never happen
	if(verbosity > 2000) cout << " sigmaZ " << sigmaZ << " fFitRangeZ " << fFitRangeZ << " sigmaP " << sigmaP << endl;
	//if(fFitRangeZ>fFitRangeMZ) fFitRangeZ = fFitRangeMZ;  // does not allow for high angle tracks, see mod to fit() method
        for(int phibin = 0; phibin!=fNPhiBins; ++phibin) {
          float radius = fGeoLayer->get_radius() + 0.5*fGeoLayer->get_thickness();
	  fFitRangeP = int( 3.0*sigmaP/(radius*stepp) + 1);
	  if(fFitRangeP<1) fFitRangeP = 1;
	  if(fFitRangeP>fFitRangeMP) fFitRangeP = fFitRangeMP;
          if(!is_local_maximum(phibin, zbin)) continue;
	  if(verbosity>2000) std::cout << " maxima found in window rphi z " << fFitRangeP << " " << fFitRangeZ << std::endl;
          fit(phibin,zbin,nhits_tot);
          if(fFitW/2000 < fEnergyCut) continue; // ignore this cluster
          SvtxCluster_v1 clus;
          clus.set_layer(layer);
          clus.set_e( fFitW/2000 );
	  float phi = fit_p_mean();
	  float pp = radius*phi;
	  float zz = fit_z_mean();
	  float pp_err = radius * fGeoLayer->get_phistep() * _inv_sqrt12;
	  float zz_err = fGeoLayer->get_zstep() * _inv_sqrt12;
	  if(fFitSizeP>1) pp_err = radius * TMath::Sqrt( fit_p_cov()/(fFitW/2000) );
	  if(fFitSizeZ>1) zz_err = TMath::Sqrt( fit_z_cov()/(fFitW/2000) );
	  //float rr_err = fGeoLayer->get_thickness() * _inv_sqrt12;
	  //float sinphi = TMath::Sin(phi);
	  //float cosphi = TMath::Cos(phi);
	  //float abscosphi = TMath::Abs(cosphi);
	  //float xx_err = rr_err*abscosphi;
	  //float yy_err = pp_err*abscosphi;
	  //float xx_err = TMath::Sqrt(pp_err*sinphi*pp_err*sinphi + rr_err*cosphi*rr_err*cosphi); // linearization
	  //float yy_err = TMath::Sqrt(pp_err*cosphi*pp_err*cosphi + rr_err*sinphi*rr_err*sinphi); // linearization
	  float pp_size = radius*fFitSizeP*fGeoLayer->get_phistep();
	  float zz_size = fFitSizeZ*fGeoLayer->get_zstep();

	  /*	  
	  cout << "layer " << layer 
	       << " get_zstep() returns "  << fGeoLayer->get_zstep() << " fFitSizeZ = " << fFitSizeZ << " fFitSizeZ*zstep = " << zz_size << endl
	       << " get_phistep() returns " << fGeoLayer->get_phistep() << " fFitSizeP = " << fFitSizeP << " radius*fFitSizeP*phistep = " << pp_size <<  endl;
	  */
	  //float xx_size = TMath::Sqrt(pp_size*sinphi*pp_size*sinphi + fGeoLayer->get_thickness()*cosphi*fGeoLayer->get_thickness()*cosphi); // linearization
	  //float yy_size = TMath::Sqrt(pp_size*cosphi*pp_size*cosphi + fGeoLayer->get_thickness()*sinphi*fGeoLayer->get_thickness()*sinphi); // linearization
	  if(verbosity>1) {
	    fHClusterEnergy->Fill(fFitW/2000);
	    fHClusterDensity->Fill(layer,zz,fFitW/2000);
	    fHClusterSizePP->Fill(layer,zz,pp_size);
	    fHClusterSizeZZ->Fill(layer,zz,zz_size);
	    fHClusterErrorPP->Fill(layer,zz,pp_err);
	    fHClusterErrorZZ->Fill(layer,zz,zz_err);
	    fHClusterDensity2->Fill(phi,zz,fFitW/2000);
	    fHClusterSizePP2->Fill(phi,zz,pp_size);
	    fHClusterSizeZZ2->Fill(phi,zz,zz_size);
	    fHClusterErrorPP2->Fill(phi,zz,pp_err);
	    fHClusterErrorZZ2->Fill(phi,zz,zz_err);
	    fHClusterWindowP->Fill(layer,zz,fFitRangeP);
	    fHClusterWindowZ->Fill(layer,zz,fFitRangeZ);
	  }
	  if(verbosity>2000)
	    std::cout << " cluster fitted " << std::endl;
	  if(verbosity>1000) {
	    std::cout << " | rad " << radius;
	    std::cout << " | size rp z " << pp_size << " " << zz_size;
	    std::cout << " | error_rphi error_z " << pp_err << " " << zz_err << std::endl;
	    std::cout << " | sgn " << fFitW/2000;
	    std::cout << " | rphi z " << pp << " " << zz;
	    std::cout << " | phi_sig z_sig phiz_cov " << TMath::Sqrt(fit_p_cov()) << " " << TMath::Sqrt(fit_z_cov()) << " " << fit_pz_cov();
	    std::cout << std::endl;
	  }
          clus.set_position(0, radius*TMath::Cos( phi ) );
          clus.set_position(1, radius*TMath::Sin( phi ) );
          clus.set_position(2, zz);
	  clus.insert_hit( fCellIDs[zbin * fNPhiBins + phibin] );


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
	  if(verbosity>1000) std::cout << " inserted " << std::endl;
        }
      }
    }
  }
  reset();
  if(verbosity>1000) std::cout << "PHG4TPCClusterizer::Process_Event DONE" << std::endl;
  if(verbosity>1) fHTime->Fill( fSW->RealTime() );
  return Fun4AllReturnCodes::EVENT_OK;
}
