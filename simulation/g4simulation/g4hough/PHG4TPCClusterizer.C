#include "PHG4TPCClusterizer.h"
#include "SvtxCluster.h"
#include "SvtxClusterMap.h"
#include "SvtxClusterMap_v1.h"
#include "SvtxCluster_v1.h"
#include "SvtxHit.h"
#include "SvtxHitMap.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <g4detectors/PHG4CylinderCell.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
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
  fMinLayer(0),
  fMaxLayer(0),
  fEnergyCut(0.1),
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
  fHClusterDensity2(NULL)
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
  for(int iz=-fFitRangeZ; iz!=fFitRangeZ; ++iz) {
    int cz = z + iz;
    if(cz < 0) continue; // skip edge
    if(cz >= fNZBins) continue; // skip edge
    for(int ip=-fFitRangeP; ip!=fFitRangeP; ++ip) {
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
  if(verbosity>1000) std::cout << "max " << fAmps[zbin*fNPhiBins+pbin] << std::endl;
  for(int iz=-fFitRangeZ; iz!=fFitRangeZ+1; ++iz) {
    int cz = zbin + iz;
    if(cz < 0) continue; // truncate edge
    if(cz >= fNZBins) continue; // truncate edge
    bool used = false;
    for(int ip=-fFitRangeP; ip!=fFitRangeP+1; ++ip) {
      int cp = wrap_phibin(pbin + ip);
      int bin = cz * fNPhiBins + cp;
      if(verbosity>1000) {
	std::cout << Form("%.2f | ",fAmps[bin]);
	if(ip==fFitRangeP) std::cout << std::endl;
      }
      if(fAmps[bin] < fFitEnergyThreshold*peak) continue; // skip small (include empty)
      used = true;
      fFitSizeP += 1.0;
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
  }
  if(fFitSizeZ>0) fFitSizeP /= fFitSizeZ; // avg size
  if(verbosity>1000) std::cout << " FIT " << fit_p_mean() << " from " << fFitP0 << std::endl;
}
//===================
int PHG4TPCClusterizer::InitRun(PHCompositeNode* topNode) {
  if(verbosity>1) {
    fHClusterEnergy = new TH1F("fHClusterEnergy","fHClusterEnergy",1000,0,1000);
    fHClusterDensity = new TProfile2D("fHClusterDensity","fHClusterDensity;LayerNo;ZZ;<E>",47,-0.5,46.5,210,-105,+105);
    fHClusterSizePP = new TProfile2D("fHClusterSizePP","fHClusterSizePP;LayerNo;ZZ;<rphisize>",47,-0.5,46.5,210,-105,+105);
    fHClusterSizeZZ = new TProfile2D("fHClusterSizeZZ","fHClusterSizeZZ;LayerNo;ZZ;<zsize>",47,-0.5,46.5,210,-105,+105);
    fHClusterErrorPP = new TProfile2D("fHClusterErrorPP","fHClusterErrorPP;LayerNo;ZZ;<rphierror>",47,-0.5,46.5,210,-105,+105);
    fHClusterErrorZZ = new TProfile2D("fHClusterErrorZZ","fHClusterErrorZZ;LayerN0;ZZ;<zerror>",47,-0.5,46.5,210,-105,+105);
    fHClusterDensity2 = new TProfile2D("fHClusterDensity2","fHClusterDensity2;Phi;ZZ;<E>",100,-TMath::Pi(),TMath::Pi(),210,-105,+105);
    fHClusterSizePP2 = new TProfile2D("fHClusterSizePP2","fHClusterSizePP2;Phi;ZZ;<rphisize>",100,-TMath::Pi(),TMath::Pi(),210,-105,+105);
    fHClusterSizeZZ2 = new TProfile2D("fHClusterSizeZZ2","fHClusterSizeZZ2;Phi;ZZ;<zsize>",100,-TMath::Pi(),TMath::Pi(),210,-105,+105);
    fHClusterErrorPP2 = new TProfile2D("fHClusterErrorP2P","fHClusterErrorPP2;Phi;ZZ;<rphierror>",100,-TMath::Pi(),TMath::Pi(),210,-105,+105);
    fHClusterErrorZZ2 = new TProfile2D("fHClusterErrorZZ2","fHClusterErrorZZ2;Phi;ZZ;<zerror>",100,-TMath::Pi(),TMath::Pi(),210,-105,+105);
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
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
//===================
void PHG4TPCClusterizer::reset() {}
//===================
int PHG4TPCClusterizer::process_event(PHCompositeNode* topNode) {
  if(verbosity>1000) std::cout << "PHG4TPCClusterizer::Process_Event" << std::endl;
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
  if (!geom_container) return Fun4AllReturnCodes::ABORTRUN;
  PHG4CylinderCellContainer* cells =  findNode::getClass<PHG4CylinderCellContainer>(dstNode,"G4CELL_SVTX");
  if (!cells) return Fun4AllReturnCodes::ABORTRUN;

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
      PHG4CylinderCell* cell = cells->findCylinderCell(hit->get_cellid()); //not needed once geofixed
      int phibin = cell->get_binphi();
      int zbin = cell->get_binz();
      if(verbosity>2000) std::cout << " phibin " << phibin << " zbin " << zbin << std::endl;
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
    while(nhits_tot > 0) {
      if(verbosity>2000) std::cout << " => nhits_tot " << nhits_tot << std::endl;
      for(int zbin = 0; zbin!=fNZBins; ++zbin) {
        if(fNHitsPerZ[zbin]<=0) continue;
	float abszbincenter = TMath::Abs(fGeoLayer->get_zcenter( zbin ));
	fFitRangeZ = 1+3*0.02*TMath::Sqrt(106-abszbincenter)/stepz;
	if(fFitRangeZ>fFitRangeMZ) fFitRangeZ = fFitRangeMZ;
        for(int phibin = 0; phibin!=fNPhiBins; ++phibin) {
	  fFitRangeP = 1+3*0.02*TMath::Sqrt(106-abszbincenter)/stepp;
	  if(fFitRangeP>fFitRangeMP) fFitRangeP = fFitRangeMP;
          if(!is_local_maximum(phibin, zbin)) continue;
	  if(verbosity>2000) std::cout << " maxima found " << std::endl;
          fit(phibin,zbin,nhits_tot);
          if(fFitW < fEnergyCut) continue; // ignore this cluster
          SvtxCluster_v1 clus;
          clus.set_layer(layer);
          clus.set_e( fFitW );
          float radius = fGeoLayer->get_radius() + 0.5*fGeoLayer->get_thickness();
	  float phi = fit_p_mean();
	  float pp = radius*phi;
	  float zz = fit_z_mean();
	  float pp_err = radius * fGeoLayer->get_phistep() * _inv_sqrt12;
	  float zz_err = fGeoLayer->get_zstep() * _inv_sqrt12;
	  if(fFitSizeP>1) pp_err = radius * TMath::Sqrt( fit_p_cov()/fFitW );
	  if(fFitSizeZ>1) zz_err = TMath::Sqrt( fit_z_cov()/fFitW );
	  float rr_err = fGeoLayer->get_thickness() * _inv_sqrt12;
	  float sinphi = TMath::Sin(phi);
	  float cosphi = TMath::Cos(phi);
	  float xx_err = TMath::Sqrt(pp_err*sinphi*pp_err*sinphi + rr_err*cosphi*rr_err*cosphi); // linearization
	  float yy_err = TMath::Sqrt(pp_err*cosphi*pp_err*cosphi + rr_err*sinphi*rr_err*sinphi); // linearization
	  float pp_size = radius*fFitSizeP*fGeoLayer->get_phistep();
	  float zz_size = fFitSizeZ*fGeoLayer->get_zstep();
	  float xx_size = TMath::Sqrt(pp_size*sinphi*pp_size*sinphi + fGeoLayer->get_thickness()*cosphi*fGeoLayer->get_thickness()*cosphi); // linearization
	  float yy_size = TMath::Sqrt(pp_size*cosphi*pp_size*cosphi + fGeoLayer->get_thickness()*sinphi*fGeoLayer->get_thickness()*sinphi); // linearization
	  if(verbosity>1) {
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
	  }
	  if(verbosity>2000)
	    std::cout << " cluster fitted " << std::endl;
	  if(verbosity>1000) {
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
	  clus.insert_hit( fCellIDs[zbin * fNPhiBins + phibin] );
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
	  svxclusters->insert(&clus);
	  if(verbosity>1000) std::cout << " inserted " << std::endl;
        }
      }
    }
  }
  reset();
  if(verbosity>1000) std::cout << "PHG4TPCClusterizer::Process_Event DONE" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
