#include "PHTruthClustering.h"

#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex.h>     // for SvtxVertex
#include <trackbase_historic/SvtxVertex_v1.h>


#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrClusterv2.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>  // for hitkey, getLayer
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TpcDefs.h>

#include <micromegas/MicromegasDefs.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>               // for PHG4CylinderGeom
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <mvtx/CylinderGeom_Mvtx.h>
#include <intt/CylinderGeomIntt.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>                         // for PHIODataNode
#include <phool/PHNode.h>                               // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>                             // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <fun4all/Fun4AllReturnCodes.h>

// gsl
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <TVector3.h>
#include <TMatrixFfwd.h>    // for TMatrixF
//#include <TMatrixT.h>       // for TMatrixT, ope...
//#include <TMatrixTUtils.h>  // for TMatrixTRow

#include <iostream>                            // for operator<<, basic_ostream
#include <set>                                 // for _Rb_tree_iterator, set
#include <utility>                             // for pair

class PHCompositeNode;

using namespace std;

PHTruthClustering::PHTruthClustering(const std::string& name)
  : SubsysReco(name)
{
}

PHTruthClustering::~PHTruthClustering() 
{

}

int PHTruthClustering::InitRun(PHCompositeNode* topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  for(unsigned int layer = 0; layer < _nlayers_maps; ++layer)
    {
      clus_err_rphi[layer] = mvtx_clus_err_rphi;
      clus_err_z[layer] = mvtx_clus_err_z;
    }
  for(unsigned int layer = _nlayers_maps; layer < _nlayers_maps + _nlayers_intt; ++layer) 
    {
      clus_err_rphi[layer] = intt_clus_err_rphi;
      clus_err_z[layer] = intt_clus_err_z;
    }
  for(unsigned int layer = _nlayers_maps + _nlayers_intt; layer < _nlayers_maps + _nlayers_intt + 16; ++layer) 
    {
      clus_err_rphi[layer] = tpc_inner_clus_err_rphi;
      clus_err_z[layer] = tpc_inner_clus_err_z;
    }
  for(unsigned int layer = _nlayers_maps + _nlayers_intt + 16; layer < _nlayers_maps + _nlayers_intt +_nlayers_tpc; ++layer) 
    {
      clus_err_rphi[layer] = tpc_outer_clus_err_rphi;
      clus_err_z[layer] = tpc_outer_clus_err_z;
    }
  for(unsigned int layer = _nlayers_maps + _nlayers_intt +_nlayers_tpc; layer <  _nlayers_maps + _nlayers_intt +_nlayers_tpc + 1; ++layer) 
    {
      clus_err_rphi[layer] = mms_layer55_clus_err_rphi;
      clus_err_z[layer] = mms_layer55_clus_err_z;
    }

for(unsigned int layer = _nlayers_maps + _nlayers_intt +_nlayers_tpc + 1; layer <  _nlayers_maps + _nlayers_intt +_nlayers_tpc + 2; ++layer) 
  {
    clus_err_rphi[layer] = mms_layer56_clus_err_rphi;
    clus_err_z[layer] = mms_layer56_clus_err_z;
  }
 
 if(Verbosity() > 3)
   {
     for(unsigned int layer = 0; layer <  _nlayers_maps + _nlayers_intt +_nlayers_tpc + 2; ++layer)
       std::cout << " layer " << layer << " clus_err _rphi " << clus_err_rphi[layer] << " clus_err_z " << clus_err_z[layer] << std::endl;
   }
 
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTruthClustering::process_event(PHCompositeNode* topNode)
{
  if (Verbosity() > 0) 
    cout << "Filling truth cluster node " << endl;

  // get node for writing truth clusters
  auto m_clusterlist = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER_TRUTH");
  if (!m_clusterlist)
  {
    cout << PHWHERE << " ERROR: Can't find TRKR_CLUSTER_TRUTH" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  
  PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");      
  PHG4TruthInfoContainer::ConstRange range = truthinfo->GetParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
       iter != range.second;
       ++iter)
    {
      PHG4Particle* g4particle = iter->second;
      
      float gtrackID = g4particle->get_track_id();

      // switch to discard secondary clusters
      if(_primary_clusters_only && gtrackID < 0) continue;

      float gflavor = g4particle->get_pid();	  

      int gembed = 0;
      bool is_primary = false;
      if (g4particle->get_parent_id() == 0)
	{
	  // primary particle
	  is_primary = true;
	  gembed = truthinfo->isEmbeded(g4particle->get_track_id());
	}
      else
	{
	  PHG4Particle* primary = truthinfo->GetPrimaryParticle(g4particle->get_primary_id());
	  gembed = truthinfo->isEmbeded(primary->get_track_id());
	}
      
      if(Verbosity() > 0)
	cout << PHWHERE << " PHG4Particle ID " << gtrackID << " gflavor " << gflavor << " is_primary " << is_primary 
	     << " gembed " << gembed << endl;

      // Get the truth clusters from this particle
      auto truth_clusters =  all_truth_clusters(g4particle);

      // output the results
      if(Verbosity() > 0) std::cout << " Truth cluster summary for g4particle " << g4particle->get_track_id() << " by layer: " << std::endl;
      for ( const auto& [ckey, gclus]:truth_clusters )
	{	  
	  float gx = gclus->getX();
	  float gy = gclus->getY();
	  float gz = gclus->getZ();
	  float ng4hits = gclus->getAdc();
	  
	  TVector3 gpos(gx, gy, gz);
	  float gr = sqrt(gx*gx+gy*gy);
	  float gphi = gpos.Phi();
	  float geta = gpos.Eta();
	  
	  float gphisize = gclus->getSize(1,1);
	  float gzsize = gclus->getSize(2,2);

	  const unsigned int trkrId = TrkrDefs::getTrkrId(ckey);

	  if(Verbosity() > 0)
	    {
        const unsigned int layer = TrkrDefs::getLayer( ckey );
	      std::cout << PHWHERE << "  ****   truth: layer " << layer << "  truth cluster key " << ckey << " ng4hits " << ng4hits << std::endl;
	      std::cout << " gr " << gr << " gx " << gx << " gy " << gy << " gz " << gz 
			<< " gphi " << gphi << " geta " << geta << " gphisize " << gphisize << " gzsize " << gzsize << endl;
	    }

	  if(trkrId == TrkrDefs::tpcId)
	    {
        // add the filled out cluster to the truth cluster node for the TPC (and MM's)
        m_clusterlist->addClusterSpecifyKey(ckey, gclus);
	    }
	}
    }

  // For the other subsystems, we just copy over all of the the clusters from the reco map
  for(const auto& hitsetkey:_reco_cluster_map->getHitSetKeys())
  {
    auto range = _reco_cluster_map->getClusters(hitsetkey);
    unsigned int trkrid = TrkrDefs::getTrkrId(hitsetkey);

    // skip TPC
    if(trkrid == TrkrDefs::tpcId)  continue;
    
    for( auto clusIter = range.first; clusIter != range.second; ++clusIter ){
      TrkrDefs::cluskey cluskey = clusIter->first;

      // we have to make a copy of the cluster, to avoid problems later
      TrkrCluster* cluster = (TrkrCluster*) clusIter->second->CloneMe();
      
      unsigned int layer = TrkrDefs::getLayer(cluskey);
      if (Verbosity() >= 3)
	{
	  std::cout << PHWHERE <<" copying cluster in layer " << layer << " from reco clusters to truth clusters " << std::endl;;
	  cluster->identify();
	}
      
       m_clusterlist->addClusterSpecifyKey(cluskey, cluster);    
    }
  }

  if(Verbosity() >=3)
    {
      std::cout << "Final TRKR_CLUSTER_TRUTH clusters:";
      m_clusterlist->identify();
    }

  return Fun4AllReturnCodes::EVENT_OK;
}


std::map<TrkrDefs::cluskey, TrkrCluster* > PHTruthClustering::all_truth_clusters(PHG4Particle* particle)
{
  // get all g4hits for this particle
  std::set<PHG4Hit*> g4hits = all_truth_hits(particle);

  if(Verbosity() > 3)
    std::cout << PHWHERE << " Truth clustering for particle " << particle->get_track_id() << " with ng4hits " << g4hits.size() << std::endl;;

  // container for storing truth clusters
  //std::map<unsigned int, TrkrCluster*> truth_clusters;
  std::map<TrkrDefs::cluskey, TrkrCluster*> truth_clusters;
  if(g4hits.empty()) return truth_clusters;

  // convert truth hits for this particle to truth clusters in each layer
  // loop over layers
  unsigned int layer;
  for(layer = 0; layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc +_nlayers_mms; ++layer)
    {
      float gx = NAN;
      float gy = NAN;
      float gz = NAN;
      float gt = NAN;
      float gedep = NAN;
      
      std::vector<PHG4Hit*> contributing_hits;
      std::vector<double> contributing_hits_energy;
      std::vector<std::vector<double>> contributing_hits_entry;
      std::vector<std::vector<double>> contributing_hits_exit;
      LayerClusterG4Hits(g4hits, contributing_hits, contributing_hits_energy, contributing_hits_entry, contributing_hits_exit, layer, gx, gy, gz, gt, gedep);
      if(!(gedep > 0)) continue;
  
      // we have the cluster in this layer from this truth particle
      // add the cluster to a TrkrCluster object
      TrkrDefs::cluskey ckey;
      if(layer >= _nlayers_maps + _nlayers_intt && layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc)  // in TPC
	{      
	  unsigned int side = 0;
	  if(gz > 0) side = 1;	  
	  // need accurate sector for the TPC so the hitsetkey will be correct
	  unsigned int sector = getTpcSector(gx, gy);
	  ckey = TpcDefs::genClusKey(layer, sector, side, iclus);
	}
      else if(layer < _nlayers_maps)  // in MVTX
	{
	  unsigned int stave = 0;
	  unsigned int chip = 0;
	  unsigned int strobe = 0;
	  ckey = MvtxDefs::genClusKey(layer, stave, chip, strobe, iclus);
	}
      else if(layer >= _nlayers_maps && layer < _nlayers_maps  + _nlayers_intt)  // in INTT
	{
	  // dummy ladder and phi ID
	  unsigned int ladderzid = 0;
	  unsigned int ladderphiid = 0;
	  int crossing = 0;
	  ckey = InttDefs::genClusKey(layer, ladderzid, ladderphiid,crossing,iclus); 
	}
      else if(layer >= _nlayers_maps + _nlayers_intt + _nlayers_tpc)    // in MICROMEGAS
	{
	  unsigned int tile = 0;
	  MicromegasDefs::SegmentationType segtype;
	  segtype  =  MicromegasDefs::SegmentationType::SEGMENTATION_PHI;
	  TrkrDefs::hitsetkey hkey = MicromegasDefs::genHitSetKey(layer, segtype, tile);
    ckey = TrkrDefs::genClusKey(hkey, iclus);
	}
      else
	{
	  std::cout << PHWHERE << "Bad layer number: " << layer << std::endl;
	  continue;
	}
      
      auto clus = new TrkrClusterv2;
      iclus++;

      // need to convert gedep to ADC value
      unsigned int adc_output = getAdcValue(gedep);

      clus->setAdc(adc_output);
      clus->setPosition(0, gx);
      clus->setPosition(1, gy);
      clus->setPosition(2, gz);
      clus->setGlobal();

      /*
      // record which g4hits contribute to this truth cluster
      for(unsigned int i=0; i< contributing_hits.size(); ++i)
	{
	  _truth_cluster_truth_hit_map.insert(make_pair(ckey,contributing_hits[i]));      
	}
      */

      // Estimate the size of the truth cluster
      float g4phisize = NAN;
      float g4zsize = NAN;
      G4ClusterSize(ckey, layer, contributing_hits_entry, contributing_hits_exit, g4phisize, g4zsize);

      /*
      std::cout << PHWHERE << " g4trackID " << particle->get_track_id() << " gedep " << gedep << " adc value " << adc_output 
		<< " g4phisize " << g4phisize << " g4zsize " << g4zsize << std::endl;
      */

      // make an estimate of the errors
      // We expect roughly 150 microns in r-phi and 700 microns in z
      // we have to rotate the errors into (x,y,z) coords 

      double clusphi = atan2(gy, gx);

      TMatrixF DIM(3, 3);
      DIM[0][0] = 0.0;
      DIM[0][1] = 0.0;
      DIM[0][2] = 0.0;
      DIM[1][0] = 0.0;
      DIM[1][1] = pow(0.5 * g4phisize, 2);  //cluster_v1 expects 1/2 of actual size
      DIM[1][2] = 0.0;
      DIM[2][0] = 0.0;
      DIM[2][1] = 0.0;
      DIM[2][2] = pow(0.5 * g4zsize,2);
      
      TMatrixF ERR(3, 3);
      ERR[0][0] = 0.0;
      ERR[0][1] = 0.0;
      ERR[0][2] = 0.0;
      ERR[1][0] = 0.0;
      ERR[1][1] = pow(clus_err_rphi[layer], 2);  //cluster_v1 expects rad, arc, z as elementsof covariance
      ERR[1][2] = 0.0;
      ERR[2][0] = 0.0;
      ERR[2][1] = 0.0;
      ERR[2][2] = pow(clus_err_z[layer], 2);
  
      TMatrixF ROT(3, 3);
      ROT[0][0] = cos(clusphi);
      ROT[0][1] = -sin(clusphi);
      ROT[0][2] = 0.0;
      ROT[1][0] = sin(clusphi);
      ROT[1][1] = cos(clusphi);
      ROT[1][2] = 0.0;
      ROT[2][0] = 0.0;
      ROT[2][1] = 0.0;
      ROT[2][2] = 1.0;
      
      TMatrixF ROT_T(3, 3);
      ROT_T.Transpose(ROT);
      
      TMatrixF COVAR_DIM(3, 3);
      COVAR_DIM = ROT * DIM * ROT_T;
      
      clus->setSize(0, 0, COVAR_DIM[0][0]);
      clus->setSize(0, 1, COVAR_DIM[0][1]);
      clus->setSize(0, 2, COVAR_DIM[0][2]);
      clus->setSize(1, 0, COVAR_DIM[1][0]);
      clus->setSize(1, 1, COVAR_DIM[1][1]);
      clus->setSize(1, 2, COVAR_DIM[1][2]);
      clus->setSize(2, 0, COVAR_DIM[2][0]);
      clus->setSize(2, 1, COVAR_DIM[2][1]);
      clus->setSize(2, 2, COVAR_DIM[2][2]);
      //cout << " covar_dim[2][2] = " <<  COVAR_DIM[2][2] << endl;
      
      TMatrixF COVAR_ERR(3, 3);
      COVAR_ERR = ROT * ERR * ROT_T;
      
      clus->setError(0, 0, COVAR_ERR[0][0]);
      clus->setError(0, 1, COVAR_ERR[0][1]);
      clus->setError(0, 2, COVAR_ERR[0][2]);
      clus->setError(1, 0, COVAR_ERR[1][0]);
      clus->setError(1, 1, COVAR_ERR[1][1]);
      clus->setError(1, 2, COVAR_ERR[1][2]);
      clus->setError(2, 0, COVAR_ERR[2][0]);
      clus->setError(2, 1, COVAR_ERR[2][1]);
      clus->setError(2, 2, COVAR_ERR[2][2]);

      if(Verbosity() > 0)
	{
	  std::cout << "    layer " << layer << " cluskey " << ckey << " cluster phi " << clusphi << " local cluster error rphi  " << clus_err_rphi[layer] 
		    << " z " << clus_err_z[layer] << std::endl;
	  if(Verbosity() > 10)
	    {
	      std::cout << " global covariance matrix:" << std::endl;
	      for(int i1=0;i1<3;++i1)
		for(int i2=0;i2<3;++i2)
		  std::cout << "  " << i1 << "  " << i2 << " cov " << clus->getError(i1,i2)  << std::endl;      
	    }
	}
		
      truth_clusters.insert(std::make_pair(ckey, clus));

    }  // end loop over layers for this particle

  return truth_clusters;
}

void PHTruthClustering::LayerClusterG4Hits(std::set<PHG4Hit*> truth_hits, std::vector<PHG4Hit*> &contributing_hits, std::vector<double> &contributing_hits_energy, std::vector<std::vector<double>> &contributing_hits_entry, std::vector<std::vector<double>> &contributing_hits_exit, float layer, float &x, float &y, float &z,  float &t, float &e)
{
  bool use_geo = true;

  // Given a set of g4hits, cluster them within a given layer of the TPC

  float gx = 0.0;
  float gy = 0.0;
  float gz = 0.0;
  float gr = 0.0;
  float gt = 0.0;
  float gwt = 0.0;

  if (layer >= _nlayers_maps + _nlayers_intt && layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc)   // in TPC
    {
      //cout << "layer = " << layer << " _nlayers_maps " << _nlayers_maps << " _nlayers_intt " << _nlayers_intt << endl;

      // This calculates the truth cluster position for the TPC from all of the contributing g4hits from a g4particle, typically 2-4 for the TPC
      // Complicated, since only the part of the energy that is collected within a layer contributes to the position
      //===============================================================================
      
      PHG4TpcCylinderGeom* GeoLayer = _tpc_geom_container->GetLayerCellGeom(layer);
      // get layer boundaries here for later use
      // radii of layer boundaries
      float rbin = GeoLayer->get_radius() - GeoLayer->get_thickness() / 2.0;
      float rbout = GeoLayer->get_radius() + GeoLayer->get_thickness() / 2.0;

      if(Verbosity() > 3)
	cout << "      PHTruthClustering::LayerCluster hits for layer  " << layer << " with rbin " << rbin << " rbout " << rbout << endl;  

      // we do not assume that the truth hits know what layer they are in            
      for (std::set<PHG4Hit*>::iterator iter = truth_hits.begin();
	   iter != truth_hits.end();
	   ++iter)
	{
	  
	  PHG4Hit* this_g4hit = *iter;
	  float rbegin = sqrt(this_g4hit->get_x(0) * this_g4hit->get_x(0) + this_g4hit->get_y(0) * this_g4hit->get_y(0));
	  float rend = sqrt(this_g4hit->get_x(1) * this_g4hit->get_x(1) + this_g4hit->get_y(1) * this_g4hit->get_y(1));
	  //cout << " Eval: g4hit " << this_g4hit->get_hit_id() <<  " layer " << layer << " rbegin " << rbegin << " rend " << rend << endl;
	  
	  // make sure the entry point is at lower radius
	  float xl[2];
	  float yl[2];
	  float zl[2];
	  
	  if (rbegin < rend)
	    {
	      xl[0] = this_g4hit->get_x(0);
	      yl[0] = this_g4hit->get_y(0);
	      zl[0] = this_g4hit->get_z(0);
	      xl[1] = this_g4hit->get_x(1);
	      yl[1] = this_g4hit->get_y(1);
	      zl[1] = this_g4hit->get_z(1);
	    }
	  else
	    {
	      xl[0] = this_g4hit->get_x(1);
	      yl[0] = this_g4hit->get_y(1);
	      zl[0] = this_g4hit->get_z(1);
	      xl[1] = this_g4hit->get_x(0);
	      yl[1] = this_g4hit->get_y(0);
	      zl[1] = this_g4hit->get_z(0);
	      swap(rbegin, rend);
	      //cout << "swapped in and out " << endl;
	    }
	  
	  // check that the g4hit is not completely outside the cluster layer. Just skip this g4hit if it is
	  if (rbegin < rbin && rend < rbin)
	    continue;
	  if (rbegin > rbout && rend > rbout)
	    continue;


	  if(Verbosity() > 3)
	    {
	      cout << "     keep g4hit with rbegin " << rbegin << " rend " << rend  
		   << "         xbegin " <<  xl[0] << " xend " << xl[1]
		   << " ybegin " << yl[0] << " yend " << yl[1]
		   << " zbegin " << zl[0] << " zend " << zl[1]
		   << endl;
	    }


	  float xin = xl[0];
	  float yin = yl[0];
	  float zin = zl[0];
	  float xout = xl[1];
	  float yout = yl[1];
	  float zout = zl[1];
	  
	  float t = NAN;
	  
	  if (rbegin < rbin)
	    {
	      // line segment begins before boundary, find where it crosses
	      t = line_circle_intersection(xl, yl, zl, rbin);
	      if (t > 0)
		{
		  xin = xl[0] + t * (xl[1] - xl[0]);
		  yin = yl[0] + t * (yl[1] - yl[0]);
		  zin = zl[0] + t * (zl[1] - zl[0]);
		}
	    }
	  
	  if (rend > rbout)
	    {
	      // line segment ends after boundary, find where it crosses
	      t = line_circle_intersection(xl, yl, zl, rbout);
	      if (t > 0)
		{
		  xout = xl[0] + t * (xl[1] - xl[0]);
		  yout = yl[0] + t * (yl[1] - yl[0]);
		  zout = zl[0] + t * (zl[1] - zl[0]);
		}
	    }

	  double rin = sqrt(xin*xin + yin*yin);
	  double rout = sqrt(xout*xout + yout*yout);

	  // we want only the fraction of edep inside the layer
	  double efrac =  this_g4hit->get_edep() * (rout - rin) / (rend - rbegin);
	  gx += (xin + xout) * 0.5 * efrac;
	  gy += (yin + yout) * 0.5 * efrac;
	  gz += (zin + zout) * 0.5 * efrac;
	  gt += this_g4hit->get_avg_t() * efrac;
	  gr += (rin + rout) * 0.5 * efrac;
	  gwt += efrac;

	  if(Verbosity() > 3)
	    {
	      cout << "      rin  " << rin << " rout " << rout 
		   << " xin " << xin << " xout " << xout << " yin " << yin << " yout " << yout << " zin " << zin << " zout " << zout 
		   << " edep " << this_g4hit->get_edep() 
		   << " this_edep " <<  efrac << endl;
	      cout << "              xavge " << (xin+xout) * 0.5 << " yavge " << (yin+yout) * 0.5 << " zavge " << (zin+zout) * 0.5 << " ravge " << (rin+rout) * 0.5
		   << endl;
	    }

	  // Capture entry and exit points
	  std::vector<double> entry_loc;
	  entry_loc.push_back(xin);
	  entry_loc.push_back(yin);
	  entry_loc.push_back(zin);
	  std::vector<double> exit_loc;
	  exit_loc.push_back(xout);
	  exit_loc.push_back(yout);
	  exit_loc.push_back(zout);

	  // this_g4hit is inside the layer, add it to the vectors
	  contributing_hits.push_back(this_g4hit);
	  contributing_hits_energy.push_back( this_g4hit->get_edep() * (zout - zin) / (zl[1] - zl[0]) );
	  contributing_hits_entry.push_back(entry_loc);
	  contributing_hits_exit.push_back(exit_loc);

	}  // loop over contributing hits

      if(gwt == 0)
	{
	  e = gwt;	  
	  return;  // will be discarded 
	}

      gx /= gwt;
      gy /= gwt;
      gz /= gwt;
      gr = (rbin + rbout) * 0.5;
      gt /= gwt;

      if(Verbosity() > 3)
	{
	  cout << " weighted means:   gx " << gx << " gy " << gy << " gz " << gz << " gr " << gr << " e " << gwt << endl;
	}

      if(use_geo)
	{
	  // The energy weighted values above have significant scatter due to fluctuations in the energy deposit from Geant
	  // Calculate the geometric mean positions instead
	  float rentry = 999.0;
	  float xentry = 999.0;
	  float yentry = 999.0;
	  float zentry = 999.0;
	  float rexit = - 999.0;
	  float xexit = -999.0;
	  float yexit = -999.0;
	  float zexit = -999.0;
	  
	  for(unsigned int ientry = 0; ientry < contributing_hits_entry.size(); ++ientry)
	    {
	      float tmpx = contributing_hits_entry[ientry][0];
	      float tmpy = contributing_hits_entry[ientry][1];
	      float tmpr = sqrt(tmpx*tmpx + tmpy*tmpy);
	      
	      if(tmpr < rentry)
		{
		  rentry =  tmpr;
		  xentry = contributing_hits_entry[ientry][0];
		  yentry = contributing_hits_entry[ientry][1];
		  zentry = contributing_hits_entry[ientry][2];
		}
	      
	      tmpx = contributing_hits_exit[ientry][0];
	      tmpy = contributing_hits_exit[ientry][1];
	      tmpr = sqrt(tmpx*tmpx + tmpy*tmpy);
	      
	      if(tmpr > rexit)
		{
		  rexit =  tmpr;
		  xexit = contributing_hits_exit[ientry][0];
		  yexit = contributing_hits_exit[ientry][1];
		  zexit = contributing_hits_exit[ientry][2];
		}
	    }
	  
	  float geo_r = (rentry+rexit)*0.5;
	  float geo_x = (xentry+xexit)*0.5;
	  float geo_y = (yentry+yexit)*0.5;
	  float geo_z = (zentry+zexit)*0.5;
	  
	  if(rexit > 0)
	    {
	      gx = geo_x;
	      gy = geo_y;
	      gz = geo_z;
	      gr = geo_r;
	    }


	  if(Verbosity() > 3)
	    {
	      cout << "      rentry  " << rentry << " rexit " << rexit 
		   << " xentry " << xentry << " xexit " << xexit << " yentry " << yentry << " yexit " << yexit << " zentry " << zentry << " zexit " << zexit << endl;
	      
	      cout  << " geometric means: geo_x " << geo_x << " geo_y " << geo_y << " geo_z " << geo_z  << " geo r " << geo_r <<  " e " << gwt << endl << endl;
	    }
	}

    }  // if TPC
  else
    {
      // not TPC, one g4hit per cluster
      for (std::set<PHG4Hit*>::iterator iter = truth_hits.begin();
	   iter != truth_hits.end();
	   ++iter)
	{
	  
	  PHG4Hit* this_g4hit = *iter;

	  if(this_g4hit->get_layer() != (unsigned int) layer) continue;
	  
	  gx = this_g4hit->get_avg_x();
	  gy = this_g4hit->get_avg_y();
	  gz = this_g4hit->get_avg_z();
	  gt = this_g4hit->get_avg_t();
	  gwt += this_g4hit->get_edep();

	  // Capture entry and exit points
	  std::vector<double> entry_loc;
	  entry_loc.push_back(this_g4hit->get_x(0));
	  entry_loc.push_back(this_g4hit->get_y(0));
	  entry_loc.push_back(this_g4hit->get_z(0));
	  std::vector<double> exit_loc;
	  exit_loc.push_back(this_g4hit->get_x(1));
	  exit_loc.push_back(this_g4hit->get_y(1));
	  exit_loc.push_back(this_g4hit->get_z(1));

	  // this_g4hit is inside the layer, add it to the vectors
	  contributing_hits.push_back(this_g4hit);
	  contributing_hits_energy.push_back( this_g4hit->get_edep() );
	  contributing_hits_entry.push_back(entry_loc);
	  contributing_hits_exit.push_back(exit_loc);
	}
    }  // not TPC

  x = gx;
  y = gy;
  z = gz;
  t = gt;
  e = gwt;

  return;
}

void PHTruthClustering::G4ClusterSize(TrkrDefs::cluskey& ckey,unsigned int layer, std::vector<std::vector<double>> contributing_hits_entry,std::vector<std::vector<double>> contributing_hits_exit, float &g4phisize, float &g4zsize)
{

  // sort the contributing g4hits in radius
  double inner_radius = 100.;
  double inner_x = NAN;
  double inner_y = NAN;
  double inner_z = NAN;;

  double outer_radius = 0.;
  double outer_x = NAN;
  double outer_y = NAN;
  double outer_z = NAN;

  for(unsigned int ihit=0;ihit<contributing_hits_entry.size(); ++ihit)
    {
      double rad1 = sqrt(pow(contributing_hits_entry[ihit][0], 2) + pow(contributing_hits_entry[ihit][1], 2));      
      if(rad1 < inner_radius)
	{
	  inner_radius = rad1;
	  inner_x = contributing_hits_entry[ihit][0];
	  inner_y = contributing_hits_entry[ihit][1];
	  inner_z = contributing_hits_entry[ihit][2];    
	}

      double rad2 = sqrt(pow(contributing_hits_exit[ihit][0], 2) + pow(contributing_hits_exit[ihit][1], 2));
      if(rad2 > outer_radius)
	{
	  outer_radius = rad2;
	  outer_x = contributing_hits_exit[ihit][0];
	  outer_y = contributing_hits_exit[ihit][1];
	  outer_z = contributing_hits_exit[ihit][2];    
	}
    }

  double inner_phi =  atan2(inner_y, inner_x);
  double outer_phi =  atan2(outer_y, outer_x);
  double avge_z = (outer_z + inner_z) / 2.0;

  // Now fold these with the expected diffusion and shaping widths
  // assume spread is +/- equals this many sigmas times diffusion and shaping when extending the size
  double sigmas = 2.0;

  double radius = (inner_radius + outer_radius)/2.;
  if(radius > 28 && radius < 80)  // TPC
    {
      PHG4TpcCylinderGeom*layergeom = _tpc_geom_container->GetLayerCellGeom(layer);

      double tpc_length = 211.0;  // cm
      double drift_velocity = 8.0 / 1000.0;  // cm/ns

      // Phi size
      //======
      double diffusion_trans =  0.006;  // cm/SQRT(cm)
      double phidiffusion = diffusion_trans * sqrt(tpc_length / 2. - fabs(avge_z));

      double added_smear_trans = 0.085; // cm
      double gem_spread = 0.04;  // 400 microns

      if(outer_phi < inner_phi) swap(outer_phi, inner_phi);

      // convert diffusion from cm to radians
      double g4max_phi =  outer_phi + sigmas * sqrt(  pow(phidiffusion, 2) + pow(added_smear_trans, 2) + pow(gem_spread, 2) ) / radius;
      double g4min_phi =  inner_phi - sigmas * sqrt(  pow(phidiffusion, 2) + pow(added_smear_trans, 2) + pow(gem_spread, 2) ) / radius;

      // find the bins containing these max and min z edges
      unsigned int phibinmin = layergeom->get_phibin(g4min_phi);
      unsigned int phibinmax = layergeom->get_phibin(g4max_phi);
      unsigned int phibinwidth = phibinmax - phibinmin + 1;
      g4phisize = (double) phibinwidth * layergeom->get_phistep() * layergeom->get_radius();

      // Z size
      //=====
      double g4max_z = 0;
      double g4min_z = 0;
 
      outer_z = fabs(outer_z);
      inner_z = fabs(inner_z);

      double diffusion_long = 0.015;  // cm/SQRT(cm)
      double zdiffusion = diffusion_long * sqrt(tpc_length / 2. - fabs(avge_z)) ;
      double zshaping_lead = 32.0 * drift_velocity;  // ns * cm/ns = cm
      double zshaping_tail = 48.0 * drift_velocity;
      double added_smear_long = 0.105;  // cm

      // largest z reaches gems first, make that the outer z
      if(outer_z < inner_z) swap(outer_z, inner_z);
      g4max_z = outer_z  + sigmas*sqrt(pow(zdiffusion,2) + pow(added_smear_long,2) + pow(zshaping_lead, 2));
      g4min_z = inner_z  -  sigmas*sqrt(pow(zdiffusion,2) + pow(added_smear_long,2) + pow(zshaping_tail, 2));

      // find the bins containing these max and min z edges
      unsigned int binmin = layergeom->get_zbin(g4min_z);
      unsigned int binmax = layergeom->get_zbin(g4max_z);
      if(binmax < binmin) swap(binmax, binmin);
      unsigned int binwidth = binmax - binmin + 1;

      // multiply total number of bins that include the edges by the bin size
      g4zsize = (double) binwidth * layergeom->get_zstep();
    }
  else if(radius > 5 && radius < 20)  // INTT
    {
      // All we have is the position and layer number

      CylinderGeomIntt *layergeom = dynamic_cast<CylinderGeomIntt *>(_intt_geom_container->GetLayerGeom(layer));

      // inner location
      double world_inner[3] = {inner_x, inner_y, inner_z};
      TVector3 world_inner_vec = {inner_x, inner_y, inner_z};

      int segment_z_bin, segment_phi_bin;
      layergeom->find_indices_from_world_location(segment_z_bin, segment_phi_bin, world_inner);
      auto hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(ckey);
      auto surf = _tgeometry->maps().getSiliconSurface(hitsetkey);
      TVector3 local_inner_vec =  layergeom->get_local_from_world_coords(surf,_tgeometry, world_inner_vec);
      double yin = local_inner_vec[1];
      double zin = local_inner_vec[2];
      int strip_y_index, strip_z_index;
      layergeom->find_strip_index_values(segment_z_bin, yin, zin, strip_y_index, strip_z_index);

	// outer location
      double world_outer[3] = {outer_x, outer_y, outer_z};
      TVector3 world_outer_vec = {outer_x, outer_y, outer_z};

      layergeom->find_indices_from_world_location(segment_z_bin, segment_phi_bin, world_outer);
      auto outerhitsetkey = TrkrDefs::getHitSetKeyFromClusKey(ckey);
      auto outersurf = _tgeometry->maps().getSiliconSurface(outerhitsetkey);
      TVector3 local_outer_vec =  layergeom->get_local_from_world_coords(outersurf,_tgeometry, world_outer_vec);
      double yout = local_outer_vec[1];
      double zout = local_outer_vec[2];
      int strip_y_index_out, strip_z_index_out;
      layergeom->find_strip_index_values(segment_z_bin, yout, zout, strip_y_index_out, strip_z_index_out);
 
      int strips = abs(strip_y_index_out - strip_y_index) + 1;
      int cols = abs(strip_z_index_out - strip_z_index) + 1;


      double strip_width = (double) strips * layergeom->get_strip_y_spacing(); // cm
      double strip_length = (double) cols * layergeom->get_strip_z_spacing(); // cm

      g4phisize = strip_width;
      g4zsize = strip_length;

      /*
      if(Verbosity() > 1)
	cout << " INTT: layer " << layer << " strips " << strips << " strip pitch " <<  layergeom->get_strip_y_spacing() << " g4phisize "<< g4phisize 
	     << " columns " << cols << " strip_z_spacing " <<  layergeom->get_strip_z_spacing() << " g4zsize " << g4zsize << endl;
      */
    }
  else if(radius > 80)  // MICROMEGAS
    {
      // made up for now
      g4phisize = 300e-04;
      g4zsize = 300e-04;
    }
  else  // MVTX
    {
      unsigned int stave, stave_outer;
      unsigned int chip, chip_outer;
      int row, row_outer;
      int column, column_outer;

      // add diffusion to entry and exit locations
      double max_diffusion_radius = 25.0e-4;  // 25 microns
      double min_diffusion_radius = 8.0e-4;  // 8 microns

      CylinderGeom_Mvtx *layergeom = dynamic_cast<CylinderGeom_Mvtx *>(_mvtx_geom_container->GetLayerGeom(layer));

      TVector3 world_inner = {inner_x, inner_y, inner_z};
      std::vector<double> world_inner_vec = { world_inner[0], world_inner[1], world_inner[2] };
      layergeom->get_sensor_indices_from_world_coords(world_inner_vec, stave, chip);
      auto ihitsetkey = TrkrDefs::getHitSetKeyFromClusKey(ckey);
      auto isurf = _tgeometry->maps().getSiliconSurface(ihitsetkey);
      TVector3 local_inner = layergeom->get_local_from_world_coords(isurf,_tgeometry, world_inner);

      TVector3 world_outer = {outer_x, outer_y, outer_z};
      std::vector<double> world_outer_vec = { world_outer[0], world_outer[1], world_outer[2] };
      layergeom->get_sensor_indices_from_world_coords(world_outer_vec, stave_outer, chip_outer);
      auto ohitsetkey = TrkrDefs::getHitSetKeyFromClusKey(ckey);
      auto osurf = _tgeometry->maps().getSiliconSurface(ohitsetkey);
      TVector3 local_outer = layergeom->get_local_from_world_coords(osurf,_tgeometry, world_outer);

      double diff =  max_diffusion_radius * 0.6;  // factor of 0.6 gives decent agreement with low occupancy reco clusters
      if(local_outer[0] < local_inner[0]) 
	diff = -diff;
      local_outer[0] += diff;
      local_inner[0] -= diff;

      double diff_outer = min_diffusion_radius * 0.6;
      if(local_outer[2] < local_inner[2]) 
	diff_outer = -diff_outer;
      local_outer[2] += diff_outer;
      local_inner[2] -= diff_outer;

      layergeom->get_pixel_from_local_coords(local_inner, row, column);
      layergeom->get_pixel_from_local_coords(local_outer, row_outer, column_outer);

      if(row_outer < row) swap(row_outer, row);
      unsigned int rows = row_outer - row + 1;
      g4phisize = (double) rows * layergeom->get_pixel_x();

      if(column_outer < column) swap(column_outer, column);
      unsigned int columns = column_outer - column + 1;
      g4zsize = (double) columns * layergeom->get_pixel_z();

      /*
      if(Verbosity() > 1)
	cout << " MVTX: layer " << layer << " rows " << rows << " pixel x " <<  layergeom->get_pixel_x() << " g4phisize "<< g4phisize 
	     << " columns " << columns << " pixel_z " <<  layergeom->get_pixel_z() << " g4zsize " << g4zsize << endl;
      */

    }
}

std::set<PHG4Hit*> PHTruthClustering::all_truth_hits(PHG4Particle* particle)
{
  std::set<PHG4Hit*> truth_hits;

  // loop over all the g4hits in the cylinder layers
  if (_g4hits_svtx)
    {
      for (PHG4HitContainer::ConstIterator g4iter = _g4hits_svtx->getHits().first;
	   g4iter != _g4hits_svtx->getHits().second;
	   ++g4iter)
	{
	  PHG4Hit* g4hit = g4iter->second;
	  if (g4hit->get_trkid() == particle->get_track_id())
	    {
	      truth_hits.insert(g4hit);
	    }
	}
    }
  
  // loop over all the g4hits in the ladder layers
  if (_g4hits_tracker)
  {
    for (PHG4HitContainer::ConstIterator g4iter = _g4hits_tracker->getHits().first;
         g4iter != _g4hits_tracker->getHits().second;
         ++g4iter)
    {
      PHG4Hit* g4hit = g4iter->second;
      if (g4hit->get_trkid() == particle->get_track_id())
	{
	  truth_hits.insert(g4hit);
	}
    }
  }

  // loop over all the g4hits in the maps ladder layers
  if (_g4hits_maps)
  {
    for (PHG4HitContainer::ConstIterator g4iter = _g4hits_maps->getHits().first;
         g4iter != _g4hits_maps->getHits().second;
         ++g4iter)
    {
      PHG4Hit* g4hit = g4iter->second;
      if (g4hit->get_trkid() == particle->get_track_id())
	{
	  truth_hits.insert(g4hit);
	}
    }
  }

  // loop over all the g4hits in the micromegas layers
  if (_g4hits_mms)
  {
    for (PHG4HitContainer::ConstIterator g4iter = _g4hits_mms->getHits().first;
         g4iter != _g4hits_mms->getHits().second;
         ++g4iter)
    {
      PHG4Hit* g4hit = g4iter->second;
      if (g4hit->get_trkid() == particle->get_track_id())
	{
	  truth_hits.insert(g4hit);
	}
    }
  }

  return truth_hits;
}

float PHTruthClustering::line_circle_intersection(float x[], float y[], float z[], float radius)
{
  // parameterize the line in terms of t (distance along the line segment, from 0-1) as
  // x = x0 + t * (x1-x0); y=y0 + t * (y1-y0); z = z0 + t * (z1-z0)
  // parameterize the cylinder (centered at x,y = 0,0) as  x^2 + y^2 = radius^2,   then
  // (x0 + t*(x1-z0))^2 + (y0+t*(y1-y0))^2 = radius^2
  // (x0^2 + y0^2 - radius^2) + (2x0*(x1-x0) + 2y0*(y1-y0))*t +  ((x1-x0)^2 + (y1-y0)^2)*t^2 = 0 = C + B*t + A*t^2
  // quadratic with:  A = (x1-x0)^2+(y1-y0)^2 ;  B = 2x0*(x1-x0) + 2y0*(y1-y0);  C = x0^2 + y0^2 - radius^2
  // solution: t = (-B +/- sqrt(B^2 - 4*A*C)) / (2*A)

  float A = (x[1] - x[0]) * (x[1] - x[0]) + (y[1] - y[0]) * (y[1] - y[0]);
  float B = 2.0 * x[0] * (x[1] - x[0]) + 2.0 * y[0] * (y[1] - y[0]);
  float C = x[0] * x[0] + y[0] * y[0] - radius * radius;
  float tup = (-B + sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
  float tdn = (-B - sqrt(B * B - 4.0 * A * C)) / (2.0 * A);

  // The limits are 0 and 1, but we allow a little for floating point precision
  float t;
  if (tdn >= -0.0e-4 && tdn <= 1.0004)
    t = tdn;
  else if (tup >= -0.0e-4 && tup <= 1.0004)
    t = tup;
  else
  {
    cout << PHWHERE << "   **** Oops! No valid solution for tup or tdn, tdn = " << tdn << " tup = " << tup << endl;
    cout << "   radius " << radius << " rbegin " << sqrt(x[0] * x[0] + y[0] * y[0]) << " rend " << sqrt(x[1] * x[1] + y[1] * y[1]) << endl;
    cout << "   x0 " << x[0] << " x1 " << x[1] << endl;
    cout << "   y0 " << y[0] << " y1 " << y[1] << endl;
    cout << "   z0 " << z[0] << " z1 " << z[1] << endl;
    cout << "   A " << A << " B " << B << " C " << C << endl;

    t = -1;
  }

  return t;
}
unsigned int PHTruthClustering::getTpcSector(double x, double y)
{
  double phi = atan2(y, x);
  unsigned int sector = (int) (12.0 * (phi + M_PI) / (2.0 * M_PI) );
  return sector;
}

unsigned int PHTruthClustering::getAdcValue(double gedep)
{
  // see TPC digitizer for algorithm
  
  // drift electrons per GeV of energy deposited in the TPC
  double Ne_dEdx = 1.56;   // keV/cm
  double CF4_dEdx = 7.00;  // keV/cm
  double Ne_NTotal = 43;    // Number/cm
  double CF4_NTotal = 100;  // Number/cm
  double Tpc_NTot = 0.5*Ne_NTotal + 0.5*CF4_NTotal;
  double Tpc_dEdx = 0.5*Ne_dEdx + 0.5*CF4_dEdx;
  double Tpc_ElectronsPerKeV = Tpc_NTot / Tpc_dEdx;
  double electrons_per_gev = Tpc_ElectronsPerKeV * 1e6;
  
  double gem_amplification = 1400; // GEM output electrons per drifted electron
  double input_electrons = gedep * electrons_per_gev * gem_amplification;
  
  // convert electrons after GEM to ADC output
  double ChargeToPeakVolts = 20;
  double ADCSignalConversionGain = ChargeToPeakVolts * 1.60e-04 * 2.4;  // 20 (or 30) mV/fC * fC/electron * scaleup factor 
  double adc_input_voltage = input_electrons * ADCSignalConversionGain;  // mV, see comments above
  unsigned int adc_output = (unsigned int) (adc_input_voltage * 1024.0 / 2200.0);  // input voltage x 1024 channels over 2200 mV max range
  if (adc_output > 1023) adc_output = 1023;
    
  return adc_output;
}

int PHTruthClustering::GetNodes(PHCompositeNode* topNode)
{
  _tgeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  _g4truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_g4truth_container)
  {
    cerr << PHWHERE << " ERROR: Can't find node G4TruthInfo" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _g4hits_mms = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MICROMEGAS");
  _g4hits_svtx = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  _g4hits_tracker = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_INTT");
  _g4hits_maps = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MVTX");

  _mms_geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MICROMEGAS_FULL");
  _tpc_geom_container = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  _intt_geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  _mvtx_geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");

  // Create the truth cluster node if required
  PHNodeIterator iter(topNode);
  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  auto trkrclusters = findNode::getClass<TrkrClusterContainer>(dstNode, "TRKR_CLUSTER_TRUTH");
  if (!trkrclusters)
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode =
        dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    trkrclusters = new TrkrClusterContainerv4;
    PHIODataNode<PHObject> *TrkrClusterContainerNode =
        new PHIODataNode<PHObject>(trkrclusters, "TRKR_CLUSTER_TRUTH", "PHObject");
    DetNode->addNode(TrkrClusterContainerNode);
  }

  _reco_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_reco_cluster_map)
  {
    cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTruthClustering::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
