#include "PHTpcCentralMembraneClusterizer.h"

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

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

/// Tracking includes
#include <trackbase/TrkrDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/CMFlashClusterv1.h>
#include <trackbase/CMFlashClusterContainerv1.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase_historic/ActsTransformations.h>


namespace 
{
  
  // stream acts vector3
  [[maybe_unused]] std::ostream& operator << (std::ostream& out, const TVector3& v )
  { 
    out << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
    return out;
  }
  
}

//____________________________________________________________________________..
PHTpcCentralMembraneClusterizer::PHTpcCentralMembraneClusterizer(const std::string &name):
 SubsysReco(name)
{}

//____________________________________________________________________________..
int PHTpcCentralMembraneClusterizer::InitRun(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);
  
  if( _histos )
  {
    m_histogramfile.reset( new TFile(m_histogramfilename.c_str(), "RECREATE"));
    m_histogramfile->cd();   
    henergy = new TH1F("henergy", "cluster energy", 200, 0, 2000);
    hxy = new TH2F("hxy","cluster x:y",800,-100,100,800,-80,80);
    hz = new TH1F("hz","cluster z", 220, -2,2);
    
    hClustE[0]= new TH1F("hRawClusterEnergy","Cluster Energy Before Merging;E[?]",200,0,2000);
    hClustE[1] = new TH1F("hMatchedClusterEnergy","Pair Cluster Energy After Merging;E[?]",200,0,2000);
    hClustE[2] = new TH1F("hSoloClusterEnergy","Lone Cluster Energy After Merging;E[?]",200,0,2000);
    
    hDist=new TH1F("hDist","3D distance to nearby clusters on same padrow;dist[cm]",100,-1,10);
    hDistRow=new TH2F("hDistRow","phi distance to nearby clusters vs (lower)row;dist[rad];padrow",100,-0.001,0.01,60,-0.5,59.5);
    hDist2=new TH1F("hDist2","phi distance to nearby clusters on same padrow;dist[rad]",100,-0.001,0.01);
    hDistRowAdj=new TH2F("hDistRowAdj","phi distance to nearby clusters vs (lower)row;dist[rad];(lower) padrow",100,-0.001,0.01,60,-0.5,59.5);
    hDist2Adj=new TH1F("hDist2Adj","phi distance to nearby clusters on adjacent padrow;dist[rad]",100,-0.001,0.01);
  }
    
  return ret;
}

//____________________________________________________________________________..
int PHTpcCentralMembraneClusterizer::process_event(PHCompositeNode *topNode)
{
  
  //local coord conversion below
  auto tgeometry = findNode::getClass<ActsGeometry>(topNode,"ActsGeometry");
  
  if(!tgeometry)
  {
    std::cout << PHWHERE << "No Acts geometry on node tree. Can't  continue." << std::endl;
  } 
  
  if(Verbosity() > 0) std::cout << std::endl << "original size of cluster map: " << _cluster_map->size() << std::endl;  
  
  std::vector<TVector3>pos; //position vector in cartesian
  std::vector<int>layer; //cluster layer number
  std::vector<unsigned int> side; // cluster side
  std::vector<int>i_pair; //vector for pair matching
  std::vector<float>energy;//vector for energy values of clusters
  int nTpcClust = 0;

  for(const auto& hitsetkey:_cluster_map->getHitSetKeys(TrkrDefs::TrkrId::tpcId))
  {
    auto clusRange = _cluster_map->getClusters(hitsetkey);
    for (auto clusiter = clusRange.first; 
    clusiter != clusRange.second; ++clusiter)
    {
      
      ++m_total_clusters;
      
      const auto& [cluskey, cluster] = *clusiter;
      auto glob = tgeometry->getGlobalPosition(cluskey, cluster);
      // std::cout << "PHTpcCentralMembraneClusterizer::process_event - key: " << cluskey << "z: " << glob.z() << std::endl;
      
      float x = glob(0);
      float y = glob(1);
      float z = glob(2);
      
      if(Verbosity() > 0)
      {
        unsigned int lyr = TrkrDefs::getLayer(cluskey);
        unsigned short  side = TpcDefs::getSide(cluskey);
        std::cout << " z " << z << " side " << side << " layer " << lyr << " Adc " << cluster->getAdc() << " x " << x << " y " << y << std::endl;
      }
      
      if(cluster->getAdc() < _min_adc_value) continue;
      if( std::abs(z) < _min_z_value) continue;
      
      ++m_accepted_clusters;
      
      i_pair.push_back(-1);
      energy.push_back(cluster->getAdc());
      nTpcClust++;
      pos.push_back(TVector3(x,y,z));
      layer.push_back((int)(TrkrDefs::getLayer(cluskey)));
      side.push_back(TpcDefs::getSide(cluskey));
      
      if(Verbosity() > 0) std::cout << ":\t" << x << "\t" << y << "\t" << z <<std::endl;
    }
  }

  if(_histos)
  {
    // fill evaluation histograms
    for (int i=0; i<nTpcClust;  ++i)
      for (int j=i+1;  j<nTpcClust;  j++)
    { 
      // must match clusters that are on the same side
      if( side[i] != side[j] ) continue;
      
      if(layer[i]==layer[j])
      {
        const TVector3 delta=pos[i]-pos[j];
        const float dphi=std::abs(pos[i].DeltaPhi(pos[j]));
        hDist->Fill(delta.Mag());
        hDist2->Fill(dphi);
        hDistRow->Fill(dphi,layer[i]);
      }
      
      if (std::abs(layer[i]-layer[j])==1)
      {  
        //match those centers to the known/expected stripe positions
        const TVector3 delta=pos[i]-pos[j];
        const float dphi=std::abs(pos[i].DeltaPhi(pos[j]));
        hDist2Adj->Fill(dphi);
        hDistRowAdj->Fill(dphi,layer[i]);
      }
    }
  }
  
  //now for each cluster, find its nearest partner on an adjacent row
  
  /* TODO: check with Tony. Should be enough to loop for j=i+1 */
  
  const float maxphidist=0.003;//as read off the plots.
  for (int i=0;i<nTpcClust;i++)
  {
    float bestphidist=maxphidist;
    for (int j=0;j<nTpcClust;j++)
    {
      // redundant to the 'adjacent row' check:  if (i==j) continue; //can't match yourself.
      // must match to an adjacent row.
      if (std::abs(layer[i]-layer[j])!=1) continue;
      
      // must match clusters that are on the same side
      if( side[i] != side[j] ) continue;
      
      const float newphidist=std::abs(pos[i].DeltaPhi(pos[j]));
      if (newphidist<bestphidist)
      {
        i_pair[i]=j;
        bestphidist=newphidist;
      }
    }  	
  }

  // check to see if the cluster pairs each match each other
  std::vector<bool>goodPair;
  bool allGood=true;
  int  nGood=0;
  
  for (int i=0;i<nTpcClust;i++)
  {
    int myPair=i_pair[i];
    int itsPair=myPair<0?-1:i_pair[myPair];
    if (i!=itsPair )
    {
      
      if( Verbosity() )
      {
        std::cout << "PHTpcCentralMembraneClusterizer::process_event -"
          << " i: " << i 
          << " myPair: " << myPair 
          << " itsPair: " << itsPair 
          << std::endl;
      }
      
      goodPair.push_back(false);
      allGood=false;
      
    } else {
      if (i<myPair) ++nGood;
      goodPair.push_back(true);
    }
  }
  
  if(Verbosity())
  {
    if (allGood) std::cout << "PHTpcCentralMembraneClusterizer::process_event - all pairs are good" << std::endl;
    else std::cout << "PHTpcCentralMembraneClusterizer::process_event - nGood: " << nGood << " out of " << nTpcClust/2 << std::endl;
  }
  
  //build the weighted cluster centers
  //==========================
  std::vector<float>aveenergy;
  std::vector<TVector3> avepos;
  std::vector<unsigned int> nclusters;
  for (int i=0;i<nTpcClust;++i)
  {
    if(_histos)  hClustE[0]->Fill(energy[i]);
    
    if (goodPair[i])
    {
      if (i_pair[i]>i)
      {
        if(_histos)  hClustE[1]->Fill(energy[i]+energy[i_pair[i]]);
        
        aveenergy.push_back(energy[i]+energy[i_pair[i]]);
        
        // The pads measure phi and z accurately
        // They do not measure R! It is taken as the center of the padrow
        // The x and y values are derived from phi and R. Not normally a problem, since tracks cross entire padrow
        // CM flash clusters have limited radial extent, do not necessarily cross padrows completely - then their nominal R is wrong.
        // So:
        //     Get phi from the energy weighted cluster phi values.
        //     Get R from the procedure below
        
        // Get phi and z centroid
        double avePhi = (pos[i].Phi() * energy[i] + pos[i_pair[i]].Phi() * energy[i_pair[i]]) * (1./(energy[i]+energy[i_pair[i]]));	      
        double aveZ = (pos[i].Z() * energy[i] + pos[i_pair[i]].Z() * energy[i_pair[i]]) * (1./(energy[i]+energy[i_pair[i]])); 	      
        
        // Single padrow CM flash clusters, R position is not well defined because cluster is smaller than padrow
        // 2-cluster case: Weighting by padrow center radius is not correct because distribution does not fill padrow (needs to be approximately linearly)
        //      Use ratio of component cluster energies to estimate number of sigmas at row boundary
        float efrac = energy[i] / (energy[i] + energy[i_pair[i]]);
        
        PHG4TpcCylinderGeom *layergeom1 = _geom_container->GetLayerCellGeom(layer[i]);
        double rad1 = layergeom1->get_radius();
        PHG4TpcCylinderGeom *layergeom2 = _geom_container->GetLayerCellGeom(layer[i_pair[i]]);
        double rad2 = layergeom2->get_radius();
        PHG4TpcCylinderGeom *layergeom0;
        double layer_dr;
        if(layer[i] != 7 && layer[i] != 23 && layer[i] != 39)
        {
          layergeom0 = _geom_container->GetLayerCellGeom(layer[i]-1);
          layer_dr = rad1 - layergeom0->get_radius();
        }
        else
        {
          layergeom0 = _geom_container->GetLayerCellGeom(layer[i]+1);
          layer_dr = layergeom0->get_radius() - rad1; 
        }
        double rad_lyr_boundary = rad1 + layer_dr / 2.0;	 
        
        
        // We have to (temporarily) use distortion corrected cluster positions to determine which stripe this came from
        Acts::Vector3 dist_pos(pos[i].X(), pos[i].Y(), pos[i].Z());
        if( _dcc)  dist_pos = _distortionCorrection.get_corrected_position( dist_pos, _dcc ); 
        double dist_r = sqrt(dist_pos[0]*dist_pos[0] + dist_pos[1] * dist_pos[1]);
        double cmclus_dr = _cmclus_dr_outer; 
        if(dist_r < 41.0)
          cmclus_dr = _cmclus_dr_inner;
        else if(dist_r >= 41.0 && rad2 < 58.0)
          cmclus_dr = _cmclus_dr_mid; 
        // Use radial width of stripe and efrac to determine where radius at center of distribution must be
        double aveR = rad_lyr_boundary - efrac * cmclus_dr + cmclus_dr/2.0;
        
        if(Verbosity() > 0)
          std::cout << " efrac " << efrac << " _cmclus_dr "<< cmclus_dr << " rad_lyr_boundary " << rad_lyr_boundary << " aveR " << aveR 
          << " layer i " << layer[i] << " R i " << rad1 << " layer i_pair " << layer[i_pair[i]] << " R i_pair " << rad2 << " layer_dr " << layer_dr << std::endl;	   
        
        TVector3 temppos(aveR*cos(avePhi), aveR*sin(avePhi), aveZ);
        avepos.push_back(temppos);
        nclusters.push_back(2);
        
        if(Verbosity() > 0)
          std::cout << " layer i " << layer[i] << " energy " << energy[i] << " pos i " << pos[i].X() << "  " << pos[i].Y() << "  " << pos[i].Z()
          << " layer i_pair " << layer[i_pair[i]] << " energy i_pair " << energy[i_pair[i]] 
          << " pos i_pair " << pos[i_pair[i]].X() << "  " <<  pos[i_pair[i]].Y() << "  " <<  pos[i_pair[i]].Z()
          << " reco pos " << temppos.x() << "  " << temppos.Y() << "  " << temppos.Z() 
          << std::endl;
      }
    } else {
      if(_histos)  hClustE[2]->Fill(energy[i]);
      // These single cluster cases have good phi, but do not have a good radius centroid estimate - may want to skip them, record nclusters
      aveenergy.push_back(energy[i]);
      avepos.push_back(pos[i]);
      nclusters.push_back(1);
    }
  }      
      
  // Loop over the vectors and put the clusters on the node tree
  //==============================================	
  if(Verbosity() > 1)  std::cout << " vector size is " << avepos.size() << std::endl; 
  
  for(unsigned int iv = 0; iv <avepos.size(); ++iv)
  {
    auto cmfc = new CMFlashClusterv1();
    
    cmfc->setX(avepos[iv].X());
    cmfc->setY(avepos[iv].Y());
    cmfc->setZ(avepos[iv].Z());
    cmfc->setAdc(aveenergy[iv]);
    cmfc->setNclusters(nclusters[iv]);
    
    _corrected_CMcluster_map->addClusterSpecifyKey(iv, cmfc);
    
    ++m_cm_clusters;
    
    if( nclusters[iv]==1 ) ++m_cm_clusters_size1;
    else if( nclusters[iv]==2 ) ++m_cm_clusters_size2;
    
  }

  // read back the clusters and make some histograms

  auto clusrange = _corrected_CMcluster_map->getClusters();
  for (auto cmitr = clusrange.first;
    cmitr !=clusrange.second;
    ++cmitr)
  {
    auto cmkey = cmitr->first;
    auto cmclus = cmitr->second;
    
    if(Verbosity() > 0)
      std::cout << "found CM cluster " << cmkey << " with adc " << cmclus->getAdc() 
      << " x " << cmclus->getX() << " y " << cmclus->getY() << " z " << cmclus->getZ() 
      << " nclusters " << cmclus->getNclusters()
      << std::endl; 
    
    if(_histos)
    {      
      henergy->Fill(cmclus->getAdc());
      hxy->Fill(cmclus->getX(), cmclus->getY());
      hz->Fill(cmclus->getZ());
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTpcCentralMembraneClusterizer::End(PHCompositeNode * /*topNode*/ )
{

  if(_histos)
  {
    m_histogramfile->cd();
    
    henergy->Write();
    hxy->Write();
    hz->Write();
    
    hClustE[0]->Write();
    hClustE[1]->Write();
    hClustE[2]->Write();
    hDist->Write();
    hDistRow->Write();
    hDist2->Write();
    hDistRowAdj->Write();
    hDist2Adj->Write();
    
    m_histogramfile->Close();
  }
  
  // print statistics
  if( Verbosity() )
  {
    std::cout
      << "PHTpcCentralMembraneClusterizer::End -"
      << " cluster statistics total: " << m_total_clusters
      << " accepted: " << m_accepted_clusters << " fraction: "
      << 100.*m_accepted_clusters/m_total_clusters << "%"
      << std::endl;
    
    std::cout
      << "PHTpcCentralMembraneClusterizer::End -"
      << " cm clusters: " << m_cm_clusters
      << std::endl;
    
    std::cout
      << "PHTpcCentralMembraneClusterizer::End -"
      << " cm clusters size 1: " << m_cm_clusters_size1
      << std::endl;

    std::cout
      << "PHTpcCentralMembraneClusterizer::End -"
      << " cm clusters size 2: " << m_cm_clusters_size2
      << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int  PHTpcCentralMembraneClusterizer::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _geom_container = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!_geom_container)
  {
    std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // tpc distortion correction
  _dcc = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainer");
  if( _dcc )
  { 
    std::cout << "PHTpcCentralMembraneMatcher:   found TPC distortion correction container" << std::endl; 
  }

  _corrected_CMcluster_map  = findNode::getClass<CMFlashClusterContainer>(topNode, "CORRECTED_CM_CLUSTER");
  if(!_corrected_CMcluster_map)
  {
    std::cout << "Creating node CORRECTED_CM_CLUSTER" << std::endl;
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
    
    _corrected_CMcluster_map = new CMFlashClusterContainerv1;
    PHIODataNode<PHObject> *TrkrClusterContainerNode = new PHIODataNode<PHObject>(_corrected_CMcluster_map, "CORRECTED_CM_CLUSTER", "PHObject");
    DetNode->addNode(TrkrClusterContainerNode);
  }    
  
  return Fun4AllReturnCodes::EVENT_OK;

}
