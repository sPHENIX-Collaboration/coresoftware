#include "PHG4TPCClusterizer.h"
#include "SvtxCluster.h"
#include "SvtxClusterMap.h"
#include "SvtxClusterMap_v1.h"
#include "SvtxCluster_v1.h"
#include "SvtxHit.h"
#include "SvtxHitMap.h"

#include <fun4all/Fun4AllReturnCodes.h>
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
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1D.h>

#include <TMatrixF.h>

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;

static int phi_span = 10;
static int z_span = 10;

static int wrap_bin(int bin, int nbins) {
  if (bin < 0) {
    bin += nbins;
  }
  if (bin >= nbins) {
    bin -= nbins;
  }
  return bin;
}

static bool is_local_maximum(const std::vector<float>& amps, int nphibins,
                             int nzbins, int phi, int z) {
  int max_width = 32;
  if (amps[z * nphibins + phi] <= 0.) {
    return false;
  }
  float cent_val = amps[z * nphibins + phi];
  bool is_max = true;
  for (int iz = -max_width; iz <= max_width; ++iz) {
    int cz = z + iz;
    if (cz < 0) {
      continue;
    }
    if (cz >= (int)(nzbins)) {
      continue;
    }
    for (int ip = -max_width; ip <= max_width; ++ip) {
      if ((iz == 0) && (ip == 0)) {
        continue;
      }
      int cp = wrap_bin(phi + ip, nphibins);
      assert(cp >= 0);
      if (amps[cz * nphibins + cp] > cent_val) {
        is_max = false;
        break;
      }
    }
    if (is_max == false) {
      break;
    }
  }
  return is_max;
}

static void fit_cluster(std::vector<float>& amps, int nphibins, int nzbins,
                        int& nhits_tot, std::vector<int>& nhits, int phibin,
                        int zbin, PHG4CylinderCellGeom* geo,
			float& phi, float& z, float& e) {
  e = 0.;
  phi = 0.;
  z = 0.;
  
  float prop_cut = 0.05;
  float peak = amps[zbin * nphibins + phibin];

  for (int iz = -z_span; iz <= z_span; ++iz) {
    int cz = zbin + iz;
    if (cz < 0) {
      continue;
    }
    if (cz >= (int)(nzbins)) {
      continue;
    }
    for (int ip = -phi_span; ip <= phi_span; ++ip) {
      int cp = wrap_bin(phibin + ip, nphibins);
      assert(cp >= 0);
      if (amps[cz * nphibins + cp] <= 0.) {
        continue;
      }
      if (amps[cz * nphibins + cp] < prop_cut * peak) {
        continue;
      }
      e += amps[cz * nphibins + cp];
      phi += amps[cz * nphibins + cp] * geo->get_phicenter(cp);
      z += amps[cz * nphibins + cp] * geo->get_zcenter(cz);
      nhits_tot -= 1;
      nhits[cz] -= 1;
      amps[cz * nphibins + cp] = 0.;
    }
  }

  phi /= e;
  z /= e;
}

int PHG4TPCClusterizer::InitRun(PHCompositeNode* topNode) {
  phi_span = _phi_span;
  z_span = _z_span;
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4TPCClusterizer::reset() {}

int PHG4TPCClusterizer::process_event(PHCompositeNode* topNode) {
  if(verbosity>1000) std::cout << "PHG4TPCClusterizer::Process_Event" << std::endl;

  PHNodeIterator iter(topNode);

  PHCompositeNode* dstNode =
      static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
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

  std::vector<std::vector<const SvtxHit*> > layer_sorted;
  PHG4CylinderCellGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for (PHG4CylinderCellGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter) {
    // We only need TPC layers here, so skip the layers below _min_layer
    // This if statement is needed because although the maps ladder layers are not included in the cylinder cell geom container, 
    // the cylinder Svx layers are, so they have to be dropped here if they are present
    if( (unsigned int) layeriter->second->get_layer() < _min_layer) {
      if(verbosity>1000) std::cout << "Skipping layer " << layeriter->second->get_layer() << std::endl;
      continue;
    }
    layer_sorted.push_back(std::vector<const SvtxHit*>());
  }
  for (SvtxHitMap::Iter iter = hits->begin(); iter != hits->end(); ++iter) {
    SvtxHit* hit = iter->second;
    if( (unsigned int) hit->get_layer() < _min_layer) continue;
    layer_sorted[hit->get_layer() - _min_layer].push_back(hit);
  }
  
  for (PHG4CylinderCellGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second; ++layeriter) {
    unsigned int layer = (unsigned int)layeriter->second->get_layer();
    // exit on the MAPS layers...
    // needed in case cylinder svtx layers are present      
    if (layer < _min_layer) continue;
    if (layer > _max_layer) continue;
    
    PHG4CylinderCellGeom* geo = geom_container->GetLayerCellGeom(layer);
    const int nphibins = layeriter->second->get_phibins();
    const int nzbins = layeriter->second->get_zbins();
    if(verbosity>1000) {
      std::cout << "Layer " << layer;
      std::cout << " nphibins " << nphibins;
      std::cout << " nzbins " << nzbins;
      std::cout << std::endl;
    }
    nhits.clear();
    nhits.assign(nzbins, 0);
    amps.clear();
    amps.assign(nphibins * nzbins, 0.);
    cellids.clear();
    cellids.assign(nphibins * nzbins, 0);

    for (unsigned int i = 0; i < layer_sorted[layer - _min_layer].size(); ++i) {
      const SvtxHit* hit = layer_sorted[layer - _min_layer][i];
      if (hit->get_e() <= 0.) continue;
      if(verbosity>1000) std::cout << hit->get_cellid();
      PHG4CylinderCell* cell = cells->findCylinderCell(hit->get_cellid());
      int phibin = cell->get_binphi();
      int zbin = cell->get_binz();
      if(verbosity>1000) std::cout << " phibin " << phibin << " zbin " << zbin << std::endl;
      nhits[zbin] += 1;
      amps[zbin * nphibins + phibin] += hit->get_e();
      cellids[zbin * nphibins + phibin] = hit->get_id();
    }

    int nhits_tot = 0;
    for (int zbin = 0; zbin < nzbins; ++zbin) {
      nhits_tot += nhits[zbin];
    }
    if(verbosity>1000) std::cout << " nhits_tot " << nhits_tot << std::endl;
    while (nhits_tot > 0) {
      if(verbosity>1000) std::cout << " => nhits_tot " << nhits_tot << std::endl;
      for (int zbin = 0; zbin!=nzbins; ++zbin) {
        if (nhits[zbin] <= 0) continue;
        for (int phibin = 0; phibin!=nphibins; ++phibin) {
          if (is_local_maximum(amps, nphibins, nzbins, phibin, zbin) == false) continue;
          float phi = 0.;
          float z = 0.;
          float e = 0.;
	  if(verbosity>1000) std::cout << " maxima found " << std::endl;
          fit_cluster(amps, nphibins, nzbins, nhits_tot, nhits, phibin, zbin, geo,
		      phi, z, e);
          if ((layer > 2) && (e < energy_cut)) {
            continue;
          }
	  if(verbosity>1000) std::cout << " cluster fitted " << std::endl;
          SvtxCluster_v1 clus;
          clus.set_layer(layer);
          clus.set_e(e);
          double radius = geo->get_radius() + 0.5*geo->get_thickness();
          clus.set_position(0, radius * cos(phi));
          clus.set_position(1, radius * sin(phi));
          clus.set_position(2, z);
	  
          clus.insert_hit(cellids[zbin * nphibins + phibin]);

	  //float invsqrt12 = 1.0/sqrt(12.);
      
	  TMatrixF DIM(3,3);
	  DIM[0][0] = 0.0;//pow(0.0*0.5*thickness,2);
	  DIM[0][1] = 0.0;
	  DIM[0][2] = 0.0;
	  DIM[1][0] = 0.0;
	  DIM[1][1] = pow(0.5*0.011,2);
	  DIM[1][2] = 0.0;
	  DIM[2][0] = 0.0;
	  DIM[2][1] = 0.0;
	  DIM[2][2] = pow(0.5*0.03,2);

	  TMatrixF ERR(3,3);
	  ERR[0][0] = 0.0;//pow(0.0*0.5*thickness*invsqrt12,2);
	  ERR[0][1] = 0.0;
	  ERR[0][2] = 0.0;
	  ERR[1][0] = 0.0;
	  ERR[1][1] = pow(0.012,2);
	  ERR[1][2] = 0.0;
	  ERR[2][0] = 0.0;
	  ERR[2][1] = 0.0;
	  ERR[2][2] = pow(0.026,2);

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
