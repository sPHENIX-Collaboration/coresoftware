#include "PHG4TPCClusterizer.h"
#include <vector>
#include "SvtxHitMap.h"
#include "SvtxHit.h"
#include "SvtxClusterMap.h"
#include "SvtxCluster.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <fun4all/getClass.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderCell.h>
#include <g4detectors/PHG4CylinderCellGeom.h>

#include "TMath.h"

#include <iostream>

#include <cstdlib>

using namespace std;


static bool is_local_maximum( std::vector<float> const& amps, unsigned int index )
{
	if( amps[index] == 0. ){return false;}
	if( (index != 0) && (index != (amps.size() - 1) ) )
	{
		if( (amps[index] >= amps[index-1]) && (amps[index] >= amps[index+1]) ){return true;}
	}
	else if( index == 0 )
	{
		if( (amps[index] >= amps[index+1]) && (amps[index] >= amps.back() ) ){return true;}
	}
	else if( index == (amps.size() - 1) )
	{
		if( (amps[index] >= amps[0]) && (amps[index] >= amps[index-1] ) ){return true;}
	}
	return false;
}


static void fit_cluster( std::vector<float>& amps, int& nhits, unsigned int index, PHG4CylinderCellGeom* geo, float& phi, float& e )
{
	nhits -= 1;
	e = amps[index];
	phi = (geo->get_phicenter(index))*amps[index];

	// TODO make this a parameter to be set from the outside
	int span = 3;

	for(int i=-span;i<0;++i)
	{

		int cur_index = (int(index))+i;
		int left_index = cur_index;
		bool left_wrap = false;
		if(cur_index<0){
			left_index += amps.size();
			left_wrap = true;}
		if(amps[left_index]<=0.){continue;}
		nhits -= 1;
		e += amps[left_index];
		if(left_wrap==false){phi += (geo->get_phicenter(left_index))*amps[left_index];}
		else{phi += (geo->get_phicenter(left_index)-2.*TMath::Pi())*amps[left_index];}
		amps[left_index] = 0.;
	}
	for(int i=1;i<=span;++i)
	{
		int cur_index = (int(index))+i;
		int right_index = cur_index;
		bool right_wrap = false;
		if(((unsigned int)cur_index)>=amps.size()){
			right_index -= amps.size();
			right_wrap = true;}
		if(amps[right_index]<=0.){continue;}
		nhits -= 1;
		e += amps[right_index];
		if(right_wrap==false){phi += (geo->get_phicenter(right_index))*amps[right_index];}
		else{phi += (geo->get_phicenter(right_index)+2.*TMath::Pi())*amps[right_index];}
		amps[right_index] = 0.;
	}
	phi /= e;
}


int PHG4TPCClusterizer::process_event(PHCompositeNode *topNode)
{
	PHNodeIterator iter(topNode);

	PHCompositeNode *dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","DST"));
	if (!dstNode) {
		cout << PHWHERE << "DST Node missing, doing nothing." << endl;
		return Fun4AllReturnCodes::ABORTRUN;}

	SvtxHitMap* hits = findNode::getClass<SvtxHitMap>(topNode,"SvtxHitMap");
	if (!hits) {
    	cout << PHWHERE << "ERROR: Can't find node SvtxHitMap" << endl;
    	return Fun4AllReturnCodes::ABORTRUN;}

	PHCompositeNode* svxNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","SVTX"));
	if (!svxNode){
		svxNode = new PHCompositeNode("SVTX");
		dstNode->addNode(svxNode);}

	SvtxClusterMap *svxclusters = findNode::getClass<SvtxClusterMap>(topNode,"SvtxClusterMap");
	if (!svxclusters){
		svxclusters = new SvtxClusterMap();
		PHIODataNode<PHObject> *SvtxClusterMapNode = new PHIODataNode<PHObject>(svxclusters, "SvtxClusterMap", "PHObject");
		svxNode->addNode(SvtxClusterMapNode);}

	PHG4CylinderCellGeomContainer* geom_container = 0;
	PHTypedNodeIterator<PHG4CylinderCellGeomContainer> geomiter(topNode);
	PHIODataNode<PHG4CylinderCellGeomContainer>* PHG4CylinderCellGeomContainerNode = geomiter.find("CYLINDERCELLGEOM_SVTX");
	if(PHG4CylinderCellGeomContainerNode){geom_container = (PHG4CylinderCellGeomContainer*) PHG4CylinderCellGeomContainerNode->getData();}
	if (!geom_container) return Fun4AllReturnCodes::ABORTRUN;;

	PHG4HitContainer* g4hits = 0;
	PHTypedNodeIterator<PHG4HitContainer> g4hititer(topNode);
	PHIODataNode<PHG4HitContainer> *PHG4HitContainerNode = g4hititer.find("G4HIT_SVTX");
	if (PHG4HitContainerNode) {g4hits = (PHG4HitContainer*)PHG4HitContainerNode->getData();}
	if (!g4hits) return Fun4AllReturnCodes::ABORTRUN;;
	
	PHG4CylinderCellContainer* cells = 0;
	PHTypedNodeIterator<PHG4CylinderCellContainer> celliter(topNode);
	PHIODataNode<PHG4CylinderCellContainer>* cell_container_node = celliter.find("G4CELL_SVTX");
	if (cell_container_node) {cells = (PHG4CylinderCellContainer*) cell_container_node->getData();}
	if (!cells) return Fun4AllReturnCodes::ABORTRUN;;

	SvtxClusterMap* clusterlist = 0;
	PHTypedNodeIterator<SvtxClusterMap> clusteriter(topNode);
	PHIODataNode<SvtxClusterMap> *SvtxClusterMapNode = clusteriter.find("SvtxClusterMap");
	if (!SvtxClusterMapNode) {
		cout << PHWHERE << " ERROR: Can't find SvtxClusterMap." << endl;
		return Fun4AllReturnCodes::ABORTRUN;}
		else {clusterlist = (SvtxClusterMap*)SvtxClusterMapNode->getData();}
	clusterlist->Reset();

	std::vector<std::vector<std::vector<float> > > amps;
	std::vector<std::vector<std::vector<int> > > cellids;
	std::vector<std::vector<int> > nhits;
	PHG4CylinderCellGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
	for(PHG4CylinderCellGeomContainer::ConstIterator layeriter = layerrange.first;layeriter != layerrange.second;++layeriter)
	{
    	int nphibins = layeriter->second->get_phibins();
    	int nzbins =  layeriter->second->get_zbins();
    	amps.push_back( std::vector<std::vector<float> >() );
    	amps.back().assign( nzbins, std::vector<float>() );
    	cellids.push_back( std::vector<std::vector<int> >() );
    	cellids.back().assign( nzbins, std::vector<int>() );
    	nhits.push_back( std::vector<int>() );
    	nhits.back().assign( nzbins, 0 );
    	for (int i = 0; i < nzbins; ++i){
    		amps.back()[i].assign( nphibins, 0. );
    		cellids.back()[i].assign( nphibins, 0 );
    	}
	}
	
	for(SvtxHitMap::Iter iter = hits->begin();iter != hits->end();++iter)
	{
		SvtxHit* hit = &iter->second;
		if(hit->get_e() <= 0.){continue;}
		int layer = hit->get_layer();
		PHG4CylinderCell* cell = cells->findCylinderCell(hit->get_cellid());
		int phibin = cell->get_binphi();
		int zbin = cell->get_binz();
		nhits[layer][zbin] += 1;
		amps[layer][zbin][phibin] += hit->get_e();
		cellids[layer][zbin][phibin] = hit->get_cellid();
	}
	for(unsigned int layer=0;layer<amps.size();++layer)
	{
		PHG4CylinderCellGeom* geo = geom_container->GetLayerCellGeom(layer);
		for(unsigned int zbin=0;zbin<amps[layer].size();++zbin)
		{
			int phibin = 0;
			while( nhits[layer][zbin] > 0 )
			{
				if( is_local_maximum( amps[layer][zbin], phibin ) == true )
				{
					float phi,e;
					fit_cluster( amps[layer][zbin], nhits[layer][zbin], phibin, geo, phi, e );
					SvtxCluster clus;
					clus.set_layer( layer );
					clus.set_e(e);
					double radius = geo->get_radius();
					clus.set_position( 0, radius*cos(phi) );
					clus.set_position( 1, radius*sin(phi) );
					clus.set_position( 2, geo->get_zcenter(zbin) );
					clus.insert_hit( cellids[layer][zbin][phibin] );
					clusterlist->insert(clus);
				}
				phibin += 1;
				if(phibin == (int)(amps[layer][zbin].size()))
				{
					phibin = 0;
				}

			}
		}
	}
	return Fun4AllReturnCodes::EVENT_OK;
}




