/*!
 *  \file		PHActsTrkFitter.C
 *  \brief		Refit SvtxTracks with PHActs.
 *  \details	Refit SvtxTracks with PHActs.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#include "PHActsTrkFitter.h"

#include <trackbase/TrkrCluster.h>                  // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/SvtxTrack_v1.h>
#include <trackbase_historic/SvtxVertexMap_v1.h>
#include <trackbase_historic/SvtxVertex_v1.h>
#include <trackbase_historic/SvtxTrackState.h>      // for SvtxTrackState
#include <trackbase_historic/SvtxVertex.h>          // for SvtxVertex
#include <trackbase_historic/SvtxVertexMap.h>       // for SvtxVertexMap

#include <mvtx/MvtxDefs.h>
#include <intt/InttDefs.h>
#include <tpc/TpcDefs.h>

#include <g4detectors/PHG4CylinderGeom.h>           // for PHG4CylinderGeom
#include <g4detectors/PHG4CylinderGeomContainer.h>

//
#include <intt/CylinderGeomIntt.h>

#include <mvtx/CylinderGeom_Mvtx.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Particlev2.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>                    // for PHG4VtxPoint
#include <g4main/PHG4VtxPointv1.h>

#include <phgeom/PHGeomUtility.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <fun4all/SubsysReco.h>                     // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>                           // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>                         // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <phfield/PHFieldUtility.h>
#include <phgeom/PHGeomUtility.h>

/*
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
*/

#include <TRotation.h>
#include <TVector3.h>
#include <TMath.h>                                  // for ATan2
#include <TMatrixT.h>                               // for TMatrixT, operator*
#include <TObject.h>
#include <TGeoManager.h>

#include <cmath>                                   // for sqrt, NAN
#include <iostream>
#include <map>
#include <memory>
#include <utility>
#include <vector>

class PHField;

using namespace std;


/*
 * Constructor
 */
PHActsTrkFitter::PHActsTrkFitter(const string& name)
  : PHTrackFitting(name)
{
  Verbosity(0);

  _event = 0;

}


int PHActsTrkFitter::Setup(PHCompositeNode *topNode)
{
  GetNodes(topNode);
  
  _geomanager = PHGeomUtility::GetTGeoManager(topNode);
  if(!_geomanager )
    {
      cout << PHWHERE << " Did not find TGeoManager, quit! " << endl;
      return false;
    }

  TGeoVolume *topVol = _geomanager->GetTopVolume();
  TObjArray *nodeArray = topVol->GetNodes();

 // recursive search for names containing siactive (INTT) or MVTXSensor (MVTX)
  // This establishes the heirarchy of volumes containing the sensors
  // For the MVTX it is 
  //    av_1_impr_phiindex_MVTXHalfStave_pv_0
  //    MVTXModule_0
  //    MVTXChip_(0 to 8)
  //    MVTXSensor_1
  // For the INTT it is:
  //    ladder_layer_(0 or 1)?_phiindex_posz    (or _negz)
  //    siactive_layer_(0 or 1)?

  TIter iObj(nodeArray); 
  while(TObject *obj = iObj())
    {
      TGeoNode *node = dynamic_cast<TGeoNode*>(obj);
      cout<< " Top Node is " << node->GetName() << " volume name is " << node->GetVolume()->GetName()  << endl;
      cout << " Mother volume name is " << node->GetMotherVolume()->GetName() << endl;
      isActive(node);
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsTrkFitter::isActive(TGeoNode *gnode)
{
  // Recursively searches gnode for silicon sensors, prints out heirarchy

  std::string node_str = gnode->GetName();
  std::string intt_refactive("siactive");
  std::string mvtx_refactive("MVTXSensor");

  if (node_str.compare(0, intt_refactive.length(), intt_refactive) == 0)
    {
      cout << "          ******* Found INTT active volume,  node is " << gnode->GetName() << " volume name is "   << gnode->GetVolume()->GetName() << endl;
      cout << "          Mother volume name is " << gnode->GetMotherVolume()->GetName() << endl;
      return;
    }
  else if (node_str.compare(0, mvtx_refactive.length(), mvtx_refactive) == 0)
    {
      cout << "          ******* Found MVTX active volume,  node is " << gnode->GetName() << " volume name is " << gnode->GetVolume()->GetName() << endl;
      cout << "          Mother volume name is " << gnode->GetMotherVolume()->GetName() << endl;
      return;
    }

  int ndaught = gnode->GetNdaughters();
  for(int i=0; i<ndaught; ++i)
    {
      cout << "     " << gnode->GetVolume()->GetName() << "  daughter " << i << " has name " << gnode->GetDaughter(i)->GetVolume()->GetName() << endl;
      isActive(gnode->GetDaughter(i));      
    }
}
      
int PHActsTrkFitter::Process()
{
  _event++;

  if (Verbosity() > 1)
    std::cout << PHWHERE << "Events processed: " << _event << std::endl;

  cout << "Start PHActsTrkfitter::process_event" << endl;

  // Access the TrkrClusters 
  TrkrClusterContainer::ConstRange clusrange = _clustermap->getClusters();
  for(TrkrClusterContainer::ConstIterator clusiter = clusrange.first; clusiter != clusrange.second; ++clusiter)
    {
      TrkrCluster *cluster = clusiter->second;
      TrkrDefs::cluskey cluskey = clusiter->first;

      // get the cluster parameters in global coordinates
      float x = cluster->getPosition(0);
      float y = cluster->getPosition(1);
      float z = cluster->getPosition(2);

      double radius = sqrt(x*x+y*y);

      // In local coords the covariances are in the  r*phi vs z frame
      // They have been rotated into global coordinates in TrkrCluster

      TMatrixF ERR(3,3);
      for(int i=0; i < 3; ++i)
	for(int j =0; j<3; j++)
	  {
	    ERR[i][j] = cluster->getError(i,j);
	  }

      // extract detector element identifier from cluskey and make Identifier for accessing TGeo element
      unsigned int layer = TrkrDefs::getLayer(cluskey);
      if(Verbosity() > 0) cout << " layer " << layer << endl;

      TVector3 world(x,y,z);
      TVector3 local(0,0,0);
      TMatrixF local_err(3, 3);

      double local_2D[2] = {0};
      double local_err_2D[2][2] = {0};

      unsigned int trkrid = TrkrDefs::getTrkrId(cluskey);  // 0 for MVTX, 1 for INTT, 2 for TPC
      if(trkrid == TrkrDefs::mvtxId)
	{
	  unsigned int staveid = MvtxDefs::getStaveId(cluskey);
	  unsigned int chipid = MvtxDefs::getChipId(cluskey);
	  // make identifier for this chip

	  // transform position back to local coords on chip
	  CylinderGeom_Mvtx *layergeom = dynamic_cast<CylinderGeom_Mvtx *>(_geom_container_mvtx->GetLayerGeom(layer));
	  local = layergeom->get_local_from_world_coords(staveid, 0, 0, chipid, world);

	  // rotate errors back to local coords too
	  double ladder_location[3] = {0.0, 0.0, 0.0};
	  // returns the center of the sensor in world coordinates - used to get the ladder phi location
	  layergeom->find_sensor_center(staveid, 0, 0, chipid, ladder_location);
	  double ladderphi = atan2(ladder_location[1], ladder_location[0]);

	  TMatrixF ROT(3, 3);
	  ROT[0][0] = cos(ladderphi);
	  ROT[0][1] = -1.0 * sin(ladderphi);
	  ROT[0][2] = 0.0;
	  ROT[1][0] = sin(ladderphi);
	  ROT[1][1] = cos(ladderphi);
	  ROT[1][2] = 0.0;
	  ROT[2][0] = 0.0; 
	  ROT[2][1] = 0.0;
	  ROT[2][2] = 1.0;

	  ROT.Invert();

	  TMatrixF ROT_T(3, 3);
	  ROT_T.Transpose(ROT);
	  
	  local_err = ROT * ERR * ROT_T;

	  if(Verbosity() > 0)
	    {
	      for(int i=0;i<3;++i)
		{
		  cout << " i " << i << " local 3D " << local[i] << endl;
		}
	      for(int i=0;i<3;++i)
		for(int j = 0; j<3; ++j)
		  {
		    cout << "  " << i << "    " << j << " local_err 3D " << local_err[i][j] << endl;
		  }
	    }
	  
	  local_2D[0] = local[0];
	  local_2D[1] = local[2];
	  local_err_2D[0][0] = local_err[1][1];
	  local_err_2D[0][1] = local_err[1][2];
	  local_err_2D[1][0] = local_err[2][1];
	  local_err_2D[1][1] = local_err[2][2];

	}
      else if (trkrid == TrkrDefs::inttId)
	{
	  unsigned int ladderzid = InttDefs::getLadderZId(cluskey);
	  unsigned int ladderphiid = InttDefs::getLadderPhiId(cluskey);
	  // make identifier for this sensor

	  // transform position back to local coords on sensor
	  // TBD! convert world to local for INTT
	  CylinderGeomIntt *layergeom = dynamic_cast<CylinderGeomIntt *>(_geom_container_intt->GetLayerGeom(layer));

	  //local = layergeom->get_local_from_world_coords(staveid, 0, 0, chipid, world);

	  // rotate errors back to local coords too	
	  double ladder_location[3] = {0.0, 0.0, 0.0};
	  layergeom->find_segment_center(ladderzid,
				    ladderphiid,
				    ladder_location);
	  double ladderphi = atan2(ladder_location[1], ladder_location[0]);
	  
	  TMatrixF ROT(3, 3);
	  ROT[0][0] = cos(ladderphi);
	  ROT[0][1] = -1.0 * sin(ladderphi);
	  ROT[0][2] = 0.0;
	  ROT[1][0] = sin(ladderphi);
	  ROT[1][1] = cos(ladderphi);
	  ROT[1][2] = 0.0;
	  ROT[2][0] = 0.0; 
	  ROT[2][1] = 0.0;
	  ROT[2][2] = 1.0;

	  ROT.Invert();

	  TMatrixF ROT_T(3, 3);
	  ROT_T.Transpose(ROT);
	  
	  TMatrixF local_err(3, 3);
	  local_err = ROT * ERR * ROT_T;

	  if(Verbosity() > 0)
	    {
	      for(int i=0;i<3;++i)
		{
		  cout << " i " << i << " local 3D " << local[i] << endl;
		}
	      for(int i=0;i<3;++i)
		for(int j = 0; j<3; ++j)
		  {
		    cout << "    " << i << "   " << j << " local_err 3D " << local_err[i][j] << endl;
		  }
	    }

	  local_2D[0] = local[1];
	  local_2D[1] = local[2];
	  local_err_2D[0][0] = local_err[1][1];
	  local_err_2D[0][1] = local_err[1][2];
	  local_err_2D[1][0] = local_err[2][1];
	  local_err_2D[1][1] = local_err[2][2];
	}
      else  // TPC
	{
	  /*
	  unsigned int sectorid = TpcDefs::getSectorId(cluskey);
	  unsigned int side = TpcDefs::getSide(cluskey);
	  unsigned int pad = TpcDefs::getPad(cluskey);
	  unsigned int tbin = TpcDefs::getTBin(cluskey);
	  // make identifier for this TPC layer -- how?
	  */

	  // transform position local coords on cylinder, at center of layer
	  // What do we mean by local coords on a cylinder?
	  // has to be phi and z, right?
	  // so it is just the phi and z part of the global coords

	  double clusphi = atan2(world[1], world[0]);
	  double r_clusphi = radius*clusphi;
	  double ztpc = world[2];

	  // rotate errors back to local coords too	
	  TMatrixF ROT(3, 3);
	  ROT[0][0] = cos(clusphi);
	  ROT[0][1] = -1.0 * sin(clusphi);
	  ROT[0][2] = 0.0;
	  ROT[1][0] = sin(clusphi);
	  ROT[1][1] = cos(clusphi);
	  ROT[1][2] = 0.0;
	  ROT[2][0] = 0.0; 
	  ROT[2][1] = 0.0;
	  ROT[2][2] = 1.0;

	  ROT.Invert();

	  TMatrixF ROT_T(3, 3);
	  ROT_T.Transpose(ROT);
	  
	  TMatrixF local_err(3, 3);
	  local_err = ROT * ERR * ROT_T;

	  if(Verbosity() > 0)
	    {
	      cout << " r " <<  " local 3D " << radius << endl;
	      cout << " r-phi " <<  " local 3D " << r_clusphi << endl;
	      cout << " z " << " local 3D " << ztpc << endl;
	      
	      for(int i=0;i<3;++i)
		for(int j = 0; j<3; ++j)
		  {
		    cout << "   " << i << "   " << j << " local_err 3D " << local_err[i][j] << endl;
		  }
	    }
      
	  local_2D[0] = r_clusphi;
	  local_2D[1] = ztpc;
	  local_err_2D[0][0] = local_err[1][1];
	  local_err_2D[0][1] = local_err[1][2];
	  local_err_2D[1][0] = local_err[2][1];
	  local_err_2D[1][1] = local_err[2][2];
	}

      // local and local_err now contain the position and covariance matrix in local coords
      if(Verbosity() > 0)
	{
	  for(int i=0;i<2;++i)
	    {
	      cout << " i " << i << " local_2D " << local_2D[i]  << endl;
	    }
	  for(int i=0;i<2;++i)
	    for(int j=0;j<2;++j)
	      {
		cout << "   " << i << "   " << j << " cov_2D " << local_err_2D[i][j] << endl;	      	      
	      }
	}

    }
  /*
  // create Acts measurement container
  
  // create the Acts cluster object

  // create Acts measurement
  //==================
  //using SourceLink = MinimalSourceLink;
  using SourceLink = clusiter;  // links from meaurement back to sPHENIX cluster
  
  template <ParID_t... params>
    using MeasurementType = Measurement<SourceLink, params...>;
  
  auto plane = Surface::makeShared<PlaneSurface>(Vector3D(0., 0., 0.),
                                                 Vector3D(1., 0., 0.));

  ActsSymMatrixD<2> covpp;
  covpp << 0.01, 0., 0., 0.02;
  MeasurementType<ParDef::eLOC_0, ParDef::eLOC_1> mpp(
						      plane, {}, std::move(covpp), 0.1, 0.2);
  
  std::vector<FittableMeasurement<SourceLink>> measurements{
    std::move(mc), std::move(mp), std::move(mpp)};
  
  //Acts::ActsSymMatrixD<3> covariance;



  // create Acts track container


  */


  // _trackmap is SvtxTrackMap from the node tree
  // We need to convert to Acts tracks
  for (SvtxTrackMap::Iter iter = _trackmap->begin(); iter != _trackmap->end();
       ++iter)
  {
    SvtxTrack* svtx_track = iter->second;
    if(Verbosity() > 0)
      {
	cout << "   found SVTXTrack " << iter->first << endl;
	svtx_track->identify();
      }
    if (!svtx_track)
      continue;

  }
  return 0;
}

/*
 * End
 */
int PHActsTrkFitter::End(PHCompositeNode* topNode)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * dtor
 */
PHActsTrkFitter::~PHActsTrkFitter()
{

}


int PHActsTrkFitter::CreateNodes(PHCompositeNode* topNode)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * GetNodes():
 *  Get all the all the required nodes off the node tree
 */
int PHActsTrkFitter::GetNodes(PHCompositeNode* topNode)
{
  //DST objects
  /*
  //Truth container
  _truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode,
                                                                "G4TruthInfo");
  if (!_truth_container && _event < 2)
  {
    cout << PHWHERE << " PHG4TruthInfoContainer node not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  */

  _geom_container_mvtx = findNode::getClass<
    PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");
  if (!_geom_container_mvtx)
  {
    cout << PHWHERE << " CYLINDERGEOM_MVTX  node not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _geom_container_intt = findNode::getClass<
    PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  if (!_geom_container_intt)
    {
    cout << PHWHERE << " CYLINDERGEOM_INTT  node not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Input Trkr Clusters
  _clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_clustermap)
  {
    cout << PHWHERE << " TRKR_CLUSTER node not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Input Svtx Tracks
  _trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!_trackmap)
  {
    cout << PHWHERE << " SvtxTrackMap node not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  /*
  // Input Svtx Vertices
  _vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!_vertexmap && _event < 2)
  {
    cout << PHWHERE << " SvtxVertexrMap node not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  */

  /*(
  // Output Svtx Tracks
  if (!(_over_write_svtxtrackmap) || _output_mode == DebugMode)
  {
    _trackmap_refit = findNode::getClass<SvtxTrackMap>(topNode,
                                                       "SvtxTrackMapRefit");
    if (!_trackmap_refit && _event < 2)
    {
      cout << PHWHERE << " SvtxTrackMapRefit node not found on node tree"
           << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  */

  return Fun4AllReturnCodes::EVENT_OK;
}



