#include "PHG4TPCClusterizer.h"
#include "SvtxHitMap.h"
#include "SvtxHit.h"
#include "SvtxClusterMap.h"
#include "SvtxClusterMap_v1.h"
#include "SvtxCluster.h"
#include "SvtxCluster_v1.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderCell.h>
#include <g4detectors/PHG4CylinderCellGeom.h>

#include <TMath.h>

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;

static int phi_span=10;
static int z_span=10;


static inline int wrap_bin( int bin, int nbins )
{
	if(bin<0){bin+=nbins;}
	if(bin>=nbins){bin-=nbins;}
	return bin;
}

static bool is_local_maximum( std::vector<std::vector<float> > const& amps, int phi, int z )
{
	if( amps[z][phi] <= 0. ){return false;}
	float cent_val = amps[z][phi];
	bool is_max = true;
	for( int iz=-1;iz<=1;++iz )
	{
		int cz = z+iz;if(cz<0){continue;}if(cz>=(int)(amps.size())){continue;}
		for( int ip=-1;ip<=1;++ip )
		{
			if( (iz==0) && (ip==0) ){continue;}
			int cp = wrap_bin( phi+ip, amps[cz].size() );
			assert (cp >= 0);
			if( amps[cz][cp] > cent_val ){is_max=false;break;}
		}
		if(is_max==false){break;}
	}
	return is_max;
}


static void fit_cluster( std::vector<std::vector<float> >& amps, int& nhits_tot, std::vector<int>& nhits, int phibin, int zbin, PHG4CylinderCellGeom* geo, float& phi, float& z, float& e )
{
	// int phi_span = 3;
	// int z_span = 1;
	e = 0.;
	phi = 0.;
	z = 0.;
	float prop_cut = 0.05;
	float peak = amps[zbin][phibin];

	for( int iz=0;iz<=z_span;++iz )
	{
		int cz = zbin+iz;if(cz<0){continue;}if(cz>=(int)(amps.size())){continue;}
		bool breakout = true;
		for( int ip=1;ip<=phi_span;++ip )
		{
			int cp = wrap_bin( phibin-ip, amps[cz].size() );
			if(amps[cz][cp] <= 0.){break;}
			if( amps[cz][cp] < prop_cut*peak ){break;}
			e += amps[cz][cp];
			phi += amps[cz][cp]*geo->get_phicenter(cp);
			z += amps[cz][cp]*geo->get_zcenter(cz);
			nhits_tot -= 1;
			nhits[cz] -= 1;
			amps[cz][cp] = 0.;
			breakout=false;
		}
		for( int ip=0;ip<=phi_span;++ip )
		{
			int cp = wrap_bin( phibin+ip, amps[cz].size() );
			if(amps[cz][cp] <= 0.){break;}
			if( amps[cz][cp] < prop_cut*peak ){break;}
			e += amps[cz][cp];
			phi += amps[cz][cp]*geo->get_phicenter(cp);
			z += amps[cz][cp]*geo->get_zcenter(cz);
			nhits_tot -= 1;
			nhits[cz] -= 1;
			amps[cz][cp] = 0.;
			breakout=false;
		}
		if(breakout==true){break;}
	}

	for( int iz=1;iz<=z_span;++iz )
	{
		int cz = zbin-iz;if(cz<0){continue;}if(cz>=(int)(amps.size())){continue;}
		bool breakout = true;
		for( int ip=1;ip<=phi_span;++ip )
		{
			int cp = wrap_bin( phibin-ip, amps[cz].size() );
			assert(cp >= 0);
			if(amps[cz][cp] <= 0.){break;}
			if( amps[cz][cp] < prop_cut*peak ){break;}
			e += amps[cz][cp];
			phi += amps[cz][cp]*geo->get_phicenter(cp);
			z += amps[cz][cp]*geo->get_zcenter(cz);
			nhits_tot -= 1;
			nhits[cz] -= 1;
			amps[cz][cp] = 0.;
			breakout=false;
		}
		for( int ip=0;ip<=phi_span;++ip )
		{
			int cp = wrap_bin( phibin+ip, amps[cz].size() );
			if(amps[cz][cp] <= 0.){break;}
			if( amps[cz][cp] < prop_cut*peak ){break;}
			e += amps[cz][cp];
			phi += amps[cz][cp]*geo->get_phicenter(cp);
			z += amps[cz][cp]*geo->get_zcenter(cz);
			nhits_tot -= 1;
			nhits[cz] -= 1;
			amps[cz][cp] = 0.;
			breakout=false;
		}
		if(breakout==true){break;}
	}

	phi /= e;
	z /= e;
}

int PHG4TPCClusterizer::InitRun(PHCompositeNode *topNode)
{
	phi_span = _phi_span;
	z_span = _z_span;

	PHG4CylinderCellGeomContainer* geom_container = 0;
	PHTypedNodeIterator<PHG4CylinderCellGeomContainer> geomiter(topNode);
	PHIODataNode<PHG4CylinderCellGeomContainer>* PHG4CylinderCellGeomContainerNode = geomiter.find("CYLINDERCELLGEOM_SVTX");
	if(PHG4CylinderCellGeomContainerNode){geom_container = (PHG4CylinderCellGeomContainer*) PHG4CylinderCellGeomContainerNode->getData();}
	if (!geom_container)
		{
			cout<<"can't find CYLINDERCELLGEOM_SVTX"<<endl;
			return Fun4AllReturnCodes::ABORTRUN;
		}
	amps.clear();
	cellids.clear();
	nhits.clear();
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
	return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4TPCClusterizer::reset()
{
	for(unsigned int i=0,isize=amps.size();i<isize;i+=1)
	{
		for(unsigned int j=0,jsize=amps[i].size();j<jsize;j+=1)
		{
			for(unsigned int k=0,ksize=amps[i][j].size();k<ksize;k+=1)
			{
				amps[i][j][k] = 0.;
			}
		}
	}
	for(unsigned int i=0,isize=cellids.size();i<isize;i+=1)
	{
		for(unsigned int j=0,jsize=cellids[i].size();j<jsize;j+=1)
		{
			for(unsigned int k=0,ksize=cellids[i][j].size();k<ksize;k+=1)
			{
				cellids[i][j][k] = 0;
			}
		}
	}
	for(unsigned int i=0,isize=nhits.size();i<isize;i+=1)
	{
		for(unsigned int j=0,jsize=nhits[i].size();j<jsize;j+=1)
		{
			nhits[i][j] = 0;
		}
	}
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
		svxclusters = new SvtxClusterMap_v1();
		PHIODataNode<PHObject> *SvtxClusterMapNode = new PHIODataNode<PHObject>(svxclusters, "SvtxClusterMap", "PHObject");
		svxNode->addNode(SvtxClusterMapNode);}

	PHG4CylinderCellGeomContainer* geom_container = 0;
	PHTypedNodeIterator<PHG4CylinderCellGeomContainer> geomiter(topNode);
	PHIODataNode<PHG4CylinderCellGeomContainer>* PHG4CylinderCellGeomContainerNode = geomiter.find("CYLINDERCELLGEOM_SVTX");
	if(PHG4CylinderCellGeomContainerNode){geom_container = (PHG4CylinderCellGeomContainer*) PHG4CylinderCellGeomContainerNode->getData();}
	if (!geom_container) return Fun4AllReturnCodes::ABORTRUN;

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
	
	for(SvtxHitMap::Iter iter = hits->begin();iter != hits->end();++iter)
	{
		SvtxHit* hit = iter->second;
		if(hit->get_e() <= 0.){continue;}
		int layer = hit->get_layer();
		PHG4CylinderCell* cell = cells->findCylinderCell(hit->get_cellid());
		int phibin = cell->get_binphi();
		int zbin = cell->get_binz();
		nhits[layer][zbin] += 1;
		amps[layer][zbin][phibin] += hit->get_e();
		cellids[layer][zbin][phibin] = hit->get_id();
	}
	for(unsigned int layer=0;layer<amps.size();++layer)
	{
		PHG4CylinderCellGeom* geo = geom_container->GetLayerCellGeom(layer);
		int nhits_tot = 0;
		for(int zbin=0;zbin<(int)(nhits[layer].size());++zbin)
		{
			nhits_tot += nhits[layer][zbin];
		}
		while( nhits_tot > 0 )
		{
			for(int zbin=0;zbin<(int)(amps[layer].size());++zbin)
			{
				if(nhits[layer][zbin] <= 0){continue;}
				for( int phibin=0;phibin<(int)(amps[layer][zbin].size());++phibin )
				{
					if( is_local_maximum( amps[layer], phibin, zbin ) == false ){continue;}
					float phi=0.;float z=0.;float e=0.;
					fit_cluster( amps[layer], nhits_tot, nhits[layer], phibin, zbin, geo, phi, z, e );
					SvtxCluster_v1 clus;
					clus.set_layer( layer );
					clus.set_e(e);
					double radius = geo->get_radius();
					clus.set_position( 0, radius*cos(phi) );
					clus.set_position( 1, radius*sin(phi) );
					clus.set_position( 2, z );
					clus.insert_hit( cellids[layer][zbin][phibin] );
					clusterlist->insert(&clus);
				}
			}
		}
	}
	reset();
	return Fun4AllReturnCodes::EVENT_OK;
}




