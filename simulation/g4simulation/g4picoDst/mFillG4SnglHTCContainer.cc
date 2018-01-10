#include "G4SnglHit.h"
#include "G4SnglTower.h"
#include "G4SnglCluster.h"
#include "G4SnglHTCContainer.h"
#include "mFillG4SnglHTCContainer.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4HitContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <phool/PHCompositeNode.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <phool/PHObject.h>

#include <phool/getClass.h>
#include <boost/foreach.hpp>
#include<sstream>
#include <cstdlib>
#include <cmath>


using namespace std;
typedef PHIODataNode<PHObject> PHObjectNode_t;

mFillG4SnglHTCContainer::mFillG4SnglHTCContainer(const std::string &name): SubsysReco( name )
{
 _nodename = "NONE";
 make_hit = false;
 make_twr = true;
 make_clr = true;
}

mFillG4SnglHTCContainer::~mFillG4SnglHTCContainer()
{
}

int
mFillG4SnglHTCContainer::Init( PHCompositeNode* top_node)
{
  cout<< "--------mFillG4SnglHTCContainer::Init( PHCompositeNode* )--------"<<endl;
  //create output node in the node tree
  PHNodeIterator iter(top_node);
  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode","DST"));	
  if(!dstNode) {
      cout << "G4SnglHTCContainer::Init(): failed to find dstNode, exiting now !"<<endl;
      exit(1);
  } else {
      cout << "dstNode is found"<<endl;
  }

  G4SnglHTCContainer *snglhtc = new G4SnglHTCContainer();

  if(snglhtc) {
      if(_nodename != "NONE"){
         PHObjectNode_t *node = new PHIODataNode<PHObject>(snglhtc,_nodename.c_str(),"PHObject");
         dstNode->addNode(node);
         cout << _nodename <<" is added" <<endl;
         } else { 
           cout << "G4SnglHTCContainer::Init() failed to create nodes" <<endl;
         }
   } else {
      cout << ThisName << "G4SnglHTCContainer::Init() failed to create output object, exiting now !"<<endl;
      exit(1);
   }

  nevents = 0;

  return 0;
}

int
mFillG4SnglHTCContainer::process_event( PHCompositeNode* top_node )
{
   if(nevents%100==0){
	cout << "--------Event: " << nevents << "--------" << endl;
   }
   nevents++;

   G4SnglHTCContainer *snglhtcs = findNode::getClass<G4SnglHTCContainer>(top_node,_nodename.c_str());
   if (!snglhtcs) { 
	cout << "mFillG4SnglHTCContainer:: G4SnglHTCContainer not in Node Tree" << endl;
	exit(1);
   }
   snglhtcs->Reset(); 

   // get the primary particle which did this to us....
   PHG4TruthInfoContainer* truthInfoList 
	=  findNode::getClass<PHG4TruthInfoContainer>(top_node , "G4TruthInfo" );

   const PHG4TruthInfoContainer::Range primRange = truthInfoList->GetPrimaryParticleRange();

   snglhtcs->set_PID( primRange.first->second->get_pid() );
   snglhtcs->set_Energy( primRange.first->second->get_e() );
   snglhtcs->set_Theta( atan2(sqrt(pow(primRange.first->second->get_px(),2) + pow(primRange.first->second->get_py(),2)), primRange.first->second->get_pz()) );
   snglhtcs->set_Phi( atan2(primRange.first->second->get_py(), primRange.first->second->get_px()) );
   snglhtcs->set_Px( primRange.first->second->get_px() );
   snglhtcs->set_Py( primRange.first->second->get_py() );
   snglhtcs->set_Pz( primRange.first->second->get_pz() );

   ostringstream hitnode, towernode, towergeomnode, clusternode;
   set<string>::const_iterator iter;
   for (iter = _node_postfix.begin(); iter != _node_postfix.end(); iter++)
   {
	int detid = (_detid.find(*iter))->second;

	/////----- Fill G4 hits into the output files -----/////
	if(make_hit){
	   hitnode.str("");
	   hitnode << "G4HIT_" << *iter;
	   PHG4HitContainer *hits = findNode::getClass<PHG4HitContainer>(top_node, hitnode.str().c_str());
	   G4SnglHit snglhit;
	   process_hit(detid, hits, snglhit, snglhtcs);
	}

	/////----- Fill G4 towers into the output files -----/////
	if(make_twr){
	   string tmpstr(*iter);
	   if( !(tmpstr.compare(0,8,"ABSORBER")==0) ){ 
	   towernode.str("");
	   towernode << "TOWER_CALIB_" << *iter;
	   RawTowerContainer *twrs = findNode::getClass<RawTowerContainer>(top_node, towernode.str().c_str());
	   towergeomnode.str("");
	   towergeomnode << "TOWERGEOM_" << *iter;
	   RawTowerGeomContainer *twrgeom = findNode::getClass<RawTowerGeomContainer>(top_node,towergeomnode.str().c_str());
	   G4SnglTower sngltwr;
	   process_twr(detid, twrs, twrgeom, sngltwr, snglhtcs);
	   }
	}

	/////----- Fill G4 clusters into the output files -----/////
	if(make_clr){
	   string tmpstr(*iter);
	   if( !(tmpstr.compare(0,8,"ABSORBER")==0) ){
	   clusternode.str("");
	   clusternode << "CLUSTER_" << *iter;
	   RawClusterContainer *clrs = findNode::getClass<RawClusterContainer>(top_node, clusternode.str().c_str());
	   G4SnglCluster snglclr;
	   process_clr(detid, clrs, snglclr, snglhtcs);
	   }
	}
   }

   return 0;
}

int
mFillG4SnglHTCContainer::End(PHCompositeNode * top_node)
{
	return 0;
}

int
mFillG4SnglHTCContainer::process_hit(int detid, PHG4HitContainer *hits, G4SnglHit snglhit, G4SnglHTCContainer 
*snglhits)
{

    snglhit.Reset();

 if (hits)   {
        PHG4HitContainer::ConstRange hit_range = hits->getHits();
        for ( PHG4HitContainer::ConstIterator hit_iter = hit_range.first ; hit_iter !=  hit_range.second; hit_iter++ )
         {
           snglhit.set_detid( detid );
           snglhit.set_layer( hit_iter->second->get_layer() );
           snglhit.set_scintid( hit_iter->second->get_scint_id() );
           snglhit.set_trackid( hit_iter->second->get_trkid() );
           snglhit.set_x0( hit_iter->second->get_x(0) );
           snglhit.set_y0( hit_iter->second->get_y(0) );
           snglhit.set_z0( hit_iter->second->get_z(0) );
           snglhit.set_x1( hit_iter->second->get_x(1) );
           snglhit.set_y1( hit_iter->second->get_y(1) );
           snglhit.set_z1( hit_iter->second->get_z(1) );
           snglhit.set_edep( hit_iter->second->get_edep() );

           snglhits->AddG4SnglHit(snglhit);
        }
 }

return 0;
}


int 
mFillG4SnglHTCContainer::process_twr(int detid, RawTowerContainer *twrs, RawTowerGeomContainer *twrgeom, G4SnglTower sngltwr, G4SnglHTCContainer *sngltwrs)
{

    sngltwr.Reset();

 if (twrs)
   {
	RawTowerContainer::ConstRange twr_range = twrs->getTowers();
	for ( RawTowerContainer::ConstIterator twr_iter = twr_range.first ; twr_iter !=  twr_range.second; twr_iter++ )
   	 {
	   sngltwr.set_detid( detid );
	   sngltwr.set_ieta( twr_iter->second->get_bineta() );
	   sngltwr.set_iphi( twr_iter->second->get_binphi() );
	   double eta = twrgeom->get_etacenter(twr_iter->second->get_bineta());
	   double phi = twrgeom->get_phicenter(twr_iter->second->get_binphi());
	   if(phi > M_PI) phi = phi - 2.*M_PI; // convert [0,2pi] to [-pi,pi] for slat geometry
	   sngltwr.set_eta( eta );
	   sngltwr.set_phi( phi );
	   sngltwr.set_edep( twr_iter->second->get_energy() );

	   sngltwrs->AddG4SnglTower(sngltwr);
	}
   }

return 0;
}

int 
mFillG4SnglHTCContainer::process_clr(int detid, RawClusterContainer *clrs, G4SnglCluster snglclr, G4SnglHTCContainer *snglclrs)
{

 snglclr.Reset();

 if (clrs)
   {
	RawClusterContainer::ConstRange clr_range = clrs->getClusters();
	for ( RawClusterContainer::ConstIterator clr_iter = clr_range.first ; clr_iter !=  clr_range.second; clr_iter++ )
   	 {
	   snglclr.set_detid( detid );
	   snglclr.set_ntowers( clr_iter->second->getNTowers() );
//     snglclr.set_eta( clr_iter->second->get_eta() );
     snglclr.set_eta( NAN );
	   snglclr.set_phi( clr_iter->second->get_phi() );
	   snglclr.set_edep( clr_iter->second->get_energy() );

	   snglclrs->AddG4SnglCluster(snglclr);
	}
   }
return 0;
}

void
mFillG4SnglHTCContainer::AddNode(const std::string &name, const int detid)
{
 _node_postfix.insert(name);
 _detid[name] = detid;
 return;
}








