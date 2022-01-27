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

/// Tracking includes
#include <trackbase/TrkrDefs.h>
#include <tpc/TpcDefs.h>
#include <trackbase/CMFlashClusterv1.h>
#include <trackbase/CMFlashClusterContainerv1.h>
#include <trackbase/TrkrClusterv3.h>
#include <trackbase/TrkrClusterContainerv3.h>
#include <trackbase/TrkrClusterHitAssocv2.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase_historic/SvtxTrack_v2.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitv2.h>
#include <trackbase/TrkrHitSetv1.h>

#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase/ActsSurfaceMaps.h>
#include <trackbase_historic/ActsTransformations.h>


using namespace std;

//____________________________________________________________________________..
PHTpcCentralMembraneClusterizer::PHTpcCentralMembraneClusterizer(const std::string &name):
 SubsysReco(name)
{


  // Make some histograms on an output file for diagnostics
  char temp[500];
  sprintf(temp,
	  "./eval_output/Energy_Histograms_%i.root", _process);
  fout = new TFile(temp,"RECREATE");
  
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

//____________________________________________________________________________..
PHTpcCentralMembraneClusterizer::~PHTpcCentralMembraneClusterizer()
{

}

//____________________________________________________________________________..
int PHTpcCentralMembraneClusterizer::InitRun(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);
  return ret;
}

//____________________________________________________________________________..
int PHTpcCentralMembraneClusterizer::process_event(PHCompositeNode *topNode)
{
 
  //local coord conversion below
  ActsTrackingGeometry *tgeometry = findNode::getClass<ActsTrackingGeometry>(topNode,"ActsTrackingGeometry");
  std::cout << "got tgeom" << std::endl;
  ActsSurfaceMaps *surfmaps = findNode::getClass<ActsSurfaceMaps>(topNode,"ActsSurfaceMaps");
  std::cout << "got surfmaps" << std::endl;
  if(!tgeometry or !surfmaps)
    {
      std::cout << PHWHERE << "No Acts geometry on node tree. Can't  continue."
		<< std::endl;
    } 
  
  if(Verbosity() > 0) std::cout << std::endl << "original size of cluster map: " << _cluster_map->size() << std::endl;  
  TrkrHitSetContainer::ConstRange hitsetrange = _hitset_map->getHitSets(TrkrDefs::TrkrId::tpcId);
  
  std::vector<TVector3>pos; //position vector in cartesian
  std::vector<int>layer; //cluster layer number
  std::vector<int>i_pair; //vector for pair matching
  std::vector<float>energy;//vector for energy values of clusters
  int nTpcClust = 0;

  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
    {
      auto clusRange = _cluster_map->getClusters(hitsetitr->first);
      for (auto clusiter = clusRange.first; 
	   clusiter != clusRange.second; ++clusiter)
	{
	  TrkrDefs::cluskey cluskey = clusiter->first;
	  auto trkrid = TrkrDefs::getTrkrId(cluskey);
	  if(trkrid != TrkrDefs::tpcId) continue;

	  TrkrCluster *cluster = clusiter->second;
	  ActsTransformations transformer;
	  auto glob = transformer.getGlobalPosition(cluster,surfmaps,tgeometry);

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

	  i_pair.push_back(-1);
	  energy.push_back(cluster->getAdc());
	  nTpcClust++;
	  pos.push_back(TVector3(x,y,z));
	  layer.push_back((int)(TrkrDefs::getLayer(cluskey)));
	  if(Verbosity() > 0)
	      std::cout << ":\t" << x << "\t" << y << "\t" << z <<std::endl;
	}
    }

  TVector3 delta;
  float dphi;
  TVector2 delta2;
  
  for (int i=0; i<nTpcClust;  ++i)
    {
      for (int j=i+1;  j<nTpcClust;  j++)
	{      //redundant to the 'adjacent row' check:  if (i==j) continue; //can't match yourself.
	  if (abs(layer[i]-layer[j])==0)
	    {
	      delta=pos[i]-pos[j];
	      dphi=abs(pos[i].DeltaPhi(pos[j]));
	      if(_histos)
		{
		  hDist->Fill(delta.Mag());
		  hDist2->Fill(dphi);
		  hDistRow->Fill(dphi,layer[i]);
		}
	    }
	  if (abs(layer[i]-layer[j])==1)
	    {  //match those centers to the known/expected stripe positions
	      
	      delta=pos[i]-pos[j];
	      dphi=abs(pos[i].DeltaPhi(pos[j]));
	      if(_histos)
		{
		  hDist2Adj->Fill(dphi);
		  hDistRowAdj->Fill(dphi,layer[i]);
		}
	    }
	}
    }
  
  //now for each cluster, find its nearest partner on an adjacent row:
  const float maxphidist=0.003;//as read off the plots.
  for (int i=0;i<nTpcClust;i++)
    {
      float bestphidist=maxphidist;
      for (int j=0;j<nTpcClust;j++)
	{
	  //redundant to the 'adjacent row' check:  if (i==j) continue; //can't match yourself.
	  if (abs(layer[i]-layer[j])!=1)
	    {
	      continue; //must match to an ADJACENT row.
	    }
	  float newphidist=abs(pos[i].DeltaPhi(pos[j]));
	  if (newphidist<bestphidist)
	    {
	      i_pair[i]=j;
	      bestphidist=newphidist;
	    }
	}  	
    }
  
  // check to see if the cluster pairs each match each other
  vector<bool>goodPair;
  bool allGood=true;
  int  nGood=0;
  
  for (int i=0;i<nTpcClust;i++)
    {
      int myPair=i_pair[i];
      int itsPair=myPair<0?-1:i_pair[myPair];
      if (i!=itsPair){
	goodPair.push_back(false);
	allGood=false;
      } else {
	if (i<myPair) nGood++;
	{
	  goodPair.push_back(true);
	}
      }
    }
  if (allGood)
    {
      printf("All Good!\n");
    } 
  else 
    {
      printf("nGood=%d out of %d\n",nGood,nTpcClust/2);
    }
  
  //build the weighted cluster centers
  //==========================
  vector<float>aveenergy;
  vector<TVector3> avepos;
    
  for (int i=0;i<nTpcClust;++i)
    {
      if(_histos)  hClustE[0]->Fill(energy[i]);

      if (goodPair[i])
	{
	  if (i_pair[i]>i)
	    {
	      if(_histos)  hClustE[1]->Fill(energy[i]+energy[i_pair[i]]);

	      aveenergy.push_back(energy[i]+energy[i_pair[i]]);
	      TVector3 temppos=energy[i]*pos[i];
	      temppos=temppos+(energy[i_pair[i]]*pos[i_pair[i]]);
	      temppos=temppos*(1./(energy[i]+energy[i_pair[i]]));
	      avepos.push_back(temppos);
	    }
	} 
      else 
	{
	  if(_histos)  hClustE[2]->Fill(energy[i]);

	  aveenergy.push_back(energy[i]);
	  avepos.push_back(pos[i]);
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
      
      _corrected_CMcluster_map->addClusterSpecifyKey(iv, cmfc);
      
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
	std::cout << "found cluster " << cmkey << " with adc " << cmclus->getAdc() 
		  << " x " << cmclus->getX() << " y " << cmclus->getY() << " z " << cmclus->getZ() 
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
      fout->cd();
      
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
      
      fout->Close();
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
    cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }


  _cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  if (!_cluster_hit_map)
  {
    cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTERHITASSOC" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }


  _hitset_map = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!_hitset_map)
  {
    cerr << PHWHERE << " ERROR: Can't find node TRKR_HITSET" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
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
      PHIODataNode<PHObject> *TrkrClusterContainerNode =
	new PHIODataNode<PHObject>(_corrected_CMcluster_map, "CORRECTED_CM_CLUSTER", "PHObject");
      DetNode->addNode(TrkrClusterContainerNode);
    }    
  
  return Fun4AllReturnCodes::EVENT_OK;
}


