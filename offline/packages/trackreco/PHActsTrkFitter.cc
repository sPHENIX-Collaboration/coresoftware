/*!
 *  \file		PHActsTrkFitter.C
 *  \brief		Refit SvtxTracks with PHActs.
 *  \details	Refit SvtxTracks with PHActs.
 *  \author	        Tony Frawley <afrawley@fsu.edu>
 */

#include "PHActsTrkFitter.h"
#include "MakeActsGeometry.h"

#include <trackbase/TrkrCluster.h>                  // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <intt/CylinderGeomIntt.h>
#include <intt/InttDefs.h>

#include <mvtx/CylinderGeom_Mvtx.h>
#include <mvtx/MvtxDefs.h>

#include <tpc/TpcDefs.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <g4detectors/PHG4CylinderGeom.h>           // for PHG4CylinderGeom
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>

#include <phgeom/PHGeomUtility.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/EventData/TrackParameters.hpp>

#include <ACTFW/Detector/IBaseDetector.hpp>
#include <ACTFW/EventData/Track.hpp>
#include <ACTFW/Framework/AlgorithmContext.hpp>
#include <ACTFW/Framework/IContextDecorator.hpp>
#include <ACTFW/Framework/WhiteBoard.hpp>
#include <ACTFW/Plugins/BField/BFieldOptions.hpp>
#include <ACTFW/Geometry/CommonGeometry.hpp>
#include <ACTFW/Options/CommonOptions.hpp>
#include <ACTFW/Plugins/Obj/ObjWriterOptions.hpp>
#include <ACTFW/Fitting/TrkrClusterFittingAlgorithm.hpp>
#include <ACTFW/Utilities/Options.hpp>


#include <TVector3.h>
#include <TMatrixT.h>                               // for TMatrixT, operator*
#include <TObject.h>
#include <TGeoManager.h>
#include <TSystem.h>
#include <TMatrixDSym.h>
#include <cmath>                              // for sqrt, NAN
#include <cstddef>                                              // for size_t
#include <cstdlib>                                              // for atoi
#include <iostream>
#include <map>
#include <memory>
#include <utility>
#include <vector>

using namespace std;

FW::TrkrClusterFittingAlgorithm::Config fitCfg;

/*
 * Constructor
 */
PHActsTrkFitter::PHActsTrkFitter(const string& name)
  : PHTrackFitting(name)
  , _geom_container_mvtx(nullptr)
  , _geom_container_intt(nullptr)
  , _geom_container_tpc(nullptr)
  , _trackmap(nullptr)
  , _clustermap(nullptr)
  , _geomanager(nullptr)
{
  Verbosity(0);

  _event = 0;
}

int PHActsTrkFitter::Setup(PHCompositeNode *topNode)
{
  GetNodes(topNode);

  // run Acts layer builder

  MakeActsGeometry *acts_geo = new MakeActsGeometry();
  acts_geo->BuildAllGeometry(topNode);

  // get the surface maps
  _cluster_surface_map_silicon = acts_geo->getSurfaceMapSilicon();
  _cluster_surface_map_tpc = acts_geo->getSurfaceMapTpc();
  _cluster_node_map = acts_geo->getNodeMap();

  SurfStepZ = acts_geo->getSurfStepZ();
  SurfStepPhi = acts_geo->getSurfStepPhi();
  ModulePhiStart = acts_geo->getModulePhiStart();
  ModuleStepPhi = acts_geo->getModuleStepPhi();
  
  geo_ctxt = acts_geo->getGeoContext();
  contextDecorators = acts_geo->getContextDecorators();

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsTrkFitter::Process()
{
  _event++;

  if (Verbosity() > 1)
    {
      std::cout << PHWHERE << "Events processed: " << _event << std::endl;
            cout << "Start PHActsTrkfitter::process_event" << endl;
    }

  std::map<TrkrDefs::cluskey, unsigned int> cluskey_hitid; 
  unsigned int  hitid = 0;
  
  FW::TrkrClusterSourceLinkContainer sourceLinks;
  TrkrDefs::hitsetkey hsetkey;
  TrkrClusterContainer::ConstRange clusrange = _clustermap->getClusters();
  for(TrkrClusterContainer::ConstIterator clusiter = clusrange.first; clusiter != clusrange.second; ++clusiter)
    {
      TrkrCluster *cluster = clusiter->second;
      TrkrDefs::cluskey cluskey = clusiter->first;

      unsigned int trkrid = TrkrDefs::getTrkrId(cluskey);  

      // map to an arbitrary hitid for later use by Acts 
      cluskey_hitid.insert(std::pair<TrkrDefs::cluskey, unsigned int>(cluskey, hitid));

      // get the cluster parameters in global coordinates
      float x = cluster->getPosition(0);
      float y = cluster->getPosition(1);
      float z = cluster->getPosition(2);
      double radius = sqrt(x*x+y*y);

      // In local coords the covariances are in the  r*phi vs z frame
      // They have been rotated into global coordinates in TrkrCluster
      TMatrixD world_err(3,3);
      for(int i=0; i < 3; ++i)
	for(int j =0; j<3; j++)
	  {
	    world_err[i][j] = cluster->getError(i,j);
	  }

      // extract detector element identifiers from cluskey so we can access the Surface and TGeoNode
      unsigned int layer = TrkrDefs::getLayer(cluskey);
      if(Verbosity() > 0) cout << std::endl << " layer " << layer << endl;

      TVector3 world(x,y,z);
      TVector3 local(0,0,0);
      TMatrixD local_err(3, 3);

      double local_2D[2] = {0};

      TGeoNode *sensor_node;
      std::shared_ptr<const Acts::Surface> surf;

      // Getting the hitsetkey from the cluskey is detector specific
      if(trkrid == TrkrDefs::mvtxId)   	  // MVTX
	{
	  unsigned int staveid = MvtxDefs::getStaveId(cluskey);
	  unsigned int chipid = MvtxDefs::getChipId(cluskey);
	  if(Verbosity() > 0)   cout << "   MVTX cluster with staveid " << staveid << " chipid " << chipid << endl; 

	  // make key for this sensor
	  hsetkey = MvtxDefs::genHitSetKey(layer, staveid, chipid);

	  // get the TGeoNode for it
	  std::map<TrkrDefs::hitsetkey, TGeoNode*>::iterator it;
	  it = _cluster_node_map.find(hsetkey);
	  if(it != _cluster_node_map.end())
	    {
	      sensor_node = it->second;
	      if(Verbosity() > 0) 
		cout << "    Found in _cluster_node_map: layer " << layer << " staveid " << staveid << " chipid " << chipid 
		     <<  " node " << sensor_node->GetName() << endl;
	    }
	  else
	    {
	      cout << PHWHERE << "No entry in TGeo map for cluster: layer " << layer << " staveid " << staveid << " chipid " << chipid  << " - should be impossible!" << endl;
	      return Fun4AllReturnCodes::ABORTEVENT;
	    }
	 
	  // Find Acts surface corresponding to it
	  std::map<TrkrDefs::hitsetkey, std::shared_ptr<const Acts::Surface>>::iterator surf_iter;
	  surf_iter = _cluster_surface_map_silicon.find(hsetkey);  
	  if(surf_iter != _cluster_surface_map_silicon.end())
	    {	      
	      //TrkrDefs::hitsetkey found_key = surf_iter->first;
	      surf = surf_iter->second;
	      if(Verbosity() > 0)
		cout<< "Got surface pair " << surf->name() << " surface type " << surf->type() << std::endl;
	  
	    }
	  else
	    {
	      cout << PHWHERE << "Failed to find associated surface element - should be impossible!" << endl;
	      return Fun4AllReturnCodes::ABORTEVENT;
	    }

	  CylinderGeom_Mvtx *layergeom = dynamic_cast<CylinderGeom_Mvtx *>(_geom_container_mvtx->GetLayerGeom(layer));
	  local = layergeom->get_local_from_world_coords(staveid, 0, 0, chipid, world);
	  local_2D[0] = local[0];
	  local_2D[1] = local[2];

	  if(Verbosity() > 10)
	    {
	      double segcent[3];
	      layergeom->find_sensor_center(staveid, 0, 0, chipid, segcent);
	      cout << "   segment center: " << segcent[0] << " " << segcent[1] << " " << segcent[2] << endl;
	      cout << "   world; " << world[0] << " " << world[1] << " " << world[2] << endl;
	      cout << "   local; " << local[0] << " " << local[1] << " " << local[2] << endl;
	    }

	  // transform covariance matrix back to local coords on chip
	  local_err = GetMvtxCovarLocal(layer, staveid, chipid, world_err);

	}
      else if (trkrid == TrkrDefs::inttId)  	  // INTT
	{
	  unsigned int ladderzid = InttDefs::getLadderZId(cluskey);
	  unsigned int ladderphiid = InttDefs::getLadderPhiId(cluskey);
	  if(Verbosity() > 0) 
	    cout << "    Intt cluster with ladderzid " << ladderzid << " ladderphid " << ladderphiid << endl; 

	  // make identifier for this sensor
	  hsetkey = InttDefs::genHitSetKey(layer, ladderzid, ladderphiid);
	  // get the TGeoNode for it
	  std::map<TrkrDefs::hitsetkey, TGeoNode*>::iterator it;
	  it = _cluster_node_map.find(hsetkey);
	  if(it == _cluster_node_map.end())
	    {
	      cout << PHWHERE << " Did not find entry in TGeo map for this cluster. That should be impossible!" << endl;
	      return Fun4AllReturnCodes::ABORTEVENT; 	      
	    }

	  sensor_node = it->second;
	  if(Verbosity() > 0)
	    cout << "      Found in _cluster_node_map:  layer " << layer << " ladderzid " << ladderzid << " ladderphiid " << ladderphiid 
		 <<  " node " << sensor_node->GetName() << endl;;

	  // Find Acts surface corresponding to this cluster
	  
	  std::map<TrkrDefs::hitsetkey, std::shared_ptr<const Acts::Surface>>::iterator surf_iter;
	  surf_iter = _cluster_surface_map_silicon.find(hsetkey);  
	  if(surf_iter == _cluster_surface_map_silicon.end())
	    {
	      cout << "Failed to find associated surface element - should be impossible " << endl;
	      return Fun4AllReturnCodes::ABORTEVENT;
	    }

	  surf = surf_iter->second;

	  // transform position back to local coords on sensor
	  CylinderGeomIntt *layergeom = dynamic_cast<CylinderGeomIntt *>(_geom_container_intt->GetLayerGeom(layer));
	  local = layergeom->get_local_from_world_coords(ladderzid, ladderphiid, world);
	  local_2D[0] = local[1];     // r*phi
	  local_2D[1] = local[2];    // z

	  if(Verbosity() > 10)
	    {
	      double segcent[3];
	      layergeom->find_segment_center(ladderzid, ladderphiid,segcent);
	      cout << "   segment center: " << segcent[0] << " " << segcent[1] << " " << segcent[2] << endl;
	      cout << "   world; " << world[0] << " " << world[1] << " " << world[2] << endl;
	      cout << "   local; " << local[0] << " " << local[1] << " " << local[2] << endl;
	    }

	  local_err = GetInttCovarLocal(layer, ladderzid, ladderphiid, world_err);

	}
      else        // TPC
	{

	  double clusphi = atan2(world[1], world[0]);
	  double r_clusphi = radius*clusphi;
	  double ztpc = world[2];
	  
	  // figure out which layer, module, side this cluster is in
	  unsigned int tpc_layer = layer;
	  unsigned int sectorid = TpcDefs::getSectorId(cluskey);
	  unsigned int side = TpcDefs::getSide(cluskey);

	  // use the cluster coords to assign the subsurface indices
	  double module_phi_low = ModulePhiStart + (double) sectorid * ModuleStepPhi;
	  unsigned int iphi = (clusphi - module_phi_low) / SurfStepPhi;
	  unsigned int iz = fabs(ztpc) / SurfStepZ;
	  unsigned int i_phi_z = iphi + 100*iz;  // for making map key

	  if(Verbosity() > 0)
	    {
	      double check_surf_rphi_center = radius * (module_phi_low + (double) iphi * SurfStepPhi + SurfStepPhi / 2.0);
	      double check_surf_z_center = (double) iz * SurfStepZ + SurfStepZ / 2.0;
	      if(side == 0)  check_surf_z_center = - check_surf_z_center;
	      
	      std::cout << " xworld " << world[0] << " yworld " << world[1] << " clusphi " << clusphi << std::endl;
	      std::cout << " sectorid " << sectorid << " side " << side << " iphi " << iphi << " iz " << iz << " i_phi_z " << i_phi_z << std::endl;
	      cout <<     "    r_clusphi " << r_clusphi << " check_surf_rphi center " << check_surf_rphi_center 
		   << " ztpc " <<  ztpc  << " check_surf_z_center " << check_surf_z_center 
		   << std::endl;
	    }

	  // get surface
	  TrkrDefs::cluskey surfkey = TpcDefs::genClusKey(tpc_layer, sectorid, side, i_phi_z);
	  std::map<TrkrDefs::cluskey, std::shared_ptr<const Acts::Surface>>::iterator surf_iter;
	  surf_iter = _cluster_surface_map_tpc.find(surfkey);  
	  if(surf_iter == _cluster_surface_map_tpc.end())
	    {
	      std::cout << PHWHERE << "Failed to find surface, should be impossible!" << std::endl;
	      return Fun4AllReturnCodes::ABORTEVENT;
	    }

	  surf = surf_iter->second->getSharedPtr();

	  // transformation of cluster to local surface coords
	  // Coords are r*phi relative to surface r-phi center, and z relative to surface z center

	  Acts::Vector3D center = surf->center(geo_ctxt);
	  double surf_rphi_center = atan2(center[1], center[0]) * radius;
	  double surf_z_center = center[2];
	  if(Verbosity() > 0)
	    {
	      std::cout << "    surface center readback:   x " << center[0] << " y " << center[1] << " phi " << atan2(center[1], center[0]) << " z " << center[2] 
			<< " surf_rphi_center " << surf_rphi_center << " surf_z_center " << surf_z_center 
			<< std::endl;
	    }

	  local_2D[0] = r_clusphi - surf_rphi_center;
	  local_2D[1] = ztpc - surf_z_center;

	  if(Verbosity() > 0)
	    {
	      cout <<     "    r_clusphi " << r_clusphi << " surf_rphi center " << surf_rphi_center 
		   << " ztpc " <<  ztpc  << " surf_z_center " << surf_z_center 
		   << " local rphi " << local_2D[0] << " local z " << local_2D[1]
		   << std::endl;
	    }

	  local_err = TransformCovarToLocal(clusphi, world_err);

	}
      
      //====================================================
      // Finished with detector specific cluster stuff 
      // We have the data needed to construct an Acts  measurement for this cluster
      //====================================================

      // Get the 2D location covariance uncertainty for the cluster
      Acts::BoundMatrix cov = Acts::BoundMatrix::Zero();
      cov(Acts::eLOC_0, Acts::eLOC_0) = local_err[0][0];
      cov(Acts::eLOC_1, Acts::eLOC_0) = local_err[1][0];
      cov(Acts::eLOC_0, Acts::eLOC_1) = local_err[0][1];
      cov(Acts::eLOC_1, Acts::eLOC_1) = local_err[1][1];

      // local and local_err contain the position and covariance matrix in local coords
      if(Verbosity() > 0)
	{
	  std::cout << "    layer " << layer << std::endl;
	  for(int i=0;i<2;++i)
	    {
	      cout << "    i " << i << "   local 2D position " << local_2D[i]  << endl;
	    }
	  
	  std::cout << "    local covariance matrix:" << std::endl;
	  std::cout << cov << std::endl;
	}
 
      // Cluster positions on GeoObject/Surface
      Acts::BoundVector loc = Acts::BoundVector::Zero();     
      loc[Acts::eLOC_0]     = local_2D[0];  
      loc[Acts::eLOC_1]     = local_2D[1];

      if(Verbosity() > 0)
	{      
	  std::cout << "    Layer " << layer << " create measurement for trkrid " << trkrid 
		    << " surface " << surf->name() << " surface type " << surf->type() 
		    << " local x " << loc[Acts::eLOC_0] << " local y " << loc[Acts::eLOC_1]  << endl;
	}
   
      /// TrkrClusterSourceLink creates an Acts::FittableMeasurement
      FW::Data::TrkrClusterSourceLink sourceLink(hitid, surf, loc, cov);
      sourceLinks.emplace_hint(sourceLinks.end(), sourceLink);
      /// Store in map which maps arbitrary hitID to sourceLink. 
      /// hitId can access Clusterkey via cluskey_hitid map
      hitidSourceLink.insert(std::pair<unsigned int, FW::Data::TrkrClusterSourceLink>(hitid, sourceLink));
            
      hitid++;
    }

  /// Construct a perigee surface as the target surface (?)
  auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
      Acts::Vector3D{0.,0.,0.});

  /// Make a vector of source links to fill for each SvtxTrack
  std::vector<FW::Data::TrkrClusterSourceLink> trackSourceLinks;

  /// Setup a context for this event
  FW::WhiteBoard eventStore(Acts::getDefaultLogger("EventStore#" + std::to_string(_event), logLevel));                  
  FW::AlgorithmContext context(0, _event, eventStore);
  
  /// _trackmap is SvtxTrackMap from the node tree
  /// We need to convert to Acts tracks
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

  
      /// Get the necessary parameters and values for the TrackParametersContainer
      Acts::BoundSymMatrix seedCov =  getActsCovMatrix(svtx_track);
      Acts::Vector3D seedPos( svtx_track->get_x() , 
			      svtx_track->get_y() , 
			      svtx_track->get_z() );
      Acts::Vector3D seedMom( svtx_track->get_px() ,
			      svtx_track->get_py() ,
			      svtx_track->get_pz() );
      // Just set to 0?
      double trackTime = 0;
      int trackQ = svtx_track->get_charge();

      FW::TrackParameters trackSeed(seedCov, seedPos, seedMom, trackQ, trackTime);

      /// Loop over clusters for this track and make a list of sourceLinks 
      /// that correspond to this track
      trackSourceLinks.clear();
      for (SvtxTrack::ConstClusterKeyIter iter = svtx_track->begin_cluster_keys();
	   iter != svtx_track->end_cluster_keys();
	   ++iter)
	{
	  TrkrDefs::cluskey cluster_key = *iter;
	  
	  /// Find the corresponding hit index
	  unsigned int hitid = cluskey_hitid.find(cluster_key)->second;

	  if(Verbosity() > 0){
	    cout << "    cluskey " << cluster_key << " has hitid " << hitid << endl;
	  }

	  /// add to the Acts ProtoTrack
	  trackSourceLinks.push_back(hitidSourceLink.find(hitid)->second);
	}

      if(Verbosity() > 0)
	{
	  for(unsigned int i=0;i<trackSourceLinks.size(); ++i)
	    {
	      cout << "   proto_track readback:  hitid " << trackSourceLinks.at(i).hitID()<< endl;
	    }
	}
    

      /// Call KF now. Have a vector of sourceLinks corresponding to clusters
      /// associated to this track and the corresponding track seed which corresponds
      /// to the PHGenFitTrkProp track seeds
      Acts::KalmanFitterOptions kfOptions(context.geoContext, 
					  context.magFieldContext,
					  context.calibContext,
					  &(*pSurface));

      

      /// Run the fitter
      auto result = fitCfg.fit(trackSourceLinks, trackSeed, kfOptions);



      /// Check that the result is okay
      if(result.ok()) {
	const auto& fitOutput = result.value();
	if(fitOutput.fittedParameters){
	  const auto& params = fitOutput.fittedParameters.value();
	  /// Get position, momentum from params
	  if(Verbosity() > 10){
	    std::cout<<"Fitted parameters for track"<<std::endl;
	    std::cout<<" position : " << params.position().transpose()<<std::endl;
	    std::cout<<" momentum : " << params.momentum().transpose()<<std::endl;
	    }
	  }

	}


      /// Add a new track to a container to put on the node tree


    }

  return 0;
}

/**
 * Helper function that puts together the acts covariance matrix from the
 * SvtxTrack covariance matrix
 */

Acts::BoundSymMatrix PHActsTrkFitter::getActsCovMatrix(SvtxTrack *track)
{
  Acts::BoundSymMatrix matrix = Acts::BoundSymMatrix::Zero();
  const double px = track->get_px();
  const double py = track->get_py();
  const double pz = track->get_pz();
  const double p = sqrt(px * px + py * py + pz * pz);

  // Get the track seed covariance matrix
  // These are the variances, so the std devs are sqrt(seed_cov[i][j])
  TMatrixDSym seed_cov(6);
  for(int i = 0; i < 6; i++){
    for(int j= 0; j <6; j++){
      seed_cov[i][j] = track->get_error(i,j);
    }       
  }

  const double sigmap = sqrt(  px * px * seed_cov[3][3]
			     + py * py * seed_cov[4][4] 
			     + pz * pz * seed_cov[5][5] ) / p ;

  // Need to convert seed_cov from x,y,z,px,py,pz basis to Acts basis of
  // x,y,phi/theta of p, qoverp, time
  double phi                  = atan(py / px);
  if(phi < -1 * M_PI)
    phi += 2. * M_PI;
  else if(phi > M_PI)
    phi -= 2. * M_PI;
  const double pxfracerr      = seed_cov[3][3] / (px * px);
  const double pyfracerr      = seed_cov[4][4] / (py * py);
  const double phiPrefactor   = fabs(py)/(fabs(px) * (1 + (py/px)*(py/px) ) );
  const double sigmaPhi       = phi * phiPrefactor * sqrt(pxfracerr + pyfracerr);
  const double theta          = acos(pz / p);
  const double thetaPrefactor = ((fabs(pz)) / ( p * sqrt(1-(pz/p)*(pz/p))));
  const double sigmaTheta     = thetaPrefactor 
    * sqrt(sigmap*sigmap/(p*p) + seed_cov[5][5]/(pz*pz));
  const double sigmaQOverP    = sigmap / (p * p);

  // Just set to 0 for now?
  const double sigmaTime      = 0;

  if(Verbosity() > 10){
    cout << "Track (px,py,pz,p) = (" << px << "," << py 
	 << "," << pz << "," << p << ")" << endl;
    cout << "Track covariance matrix: " << endl;
    for(int i = 0; i < 6; i++){
      for(int j = 0; j < 6; j++){
	cout << seed_cov[i][j] << ", ";
      }
      cout << endl;
    }
    cout << "Corresponding uncertainty calculations: " << endl;
    cout << "perr: " << sigmap << endl;
    cout << "phi: " << phi<< endl;
    cout << "pxfracerr: " << pxfracerr << endl;
    cout << "pyfracerr: " << pyfracerr << endl;
    cout << "phiPrefactor: " << phiPrefactor << endl;
    cout << "sigmaPhi: " << sigmaPhi << endl;
    cout << "theta: " << theta << endl;
    cout << "thetaPrefactor: " << thetaPrefactor << endl;
    cout << "sigmaTheta: " << sigmaTheta << endl;
    cout << "sigmaQOverP: " << sigmaQOverP << endl;

  }

  // seed covariances are already variances, so don't need to square them
  matrix(Acts::eLOC_0, Acts::eLOC_0)  = seed_cov[0][0];
  matrix(Acts::eLOC_1, Acts::eLOC_1)  = seed_cov[1][1];
  matrix(Acts::ePHI, Acts::ePHI )     = sigmaPhi * sigmaPhi;
  matrix(Acts::eTHETA, Acts::eTHETA ) = sigmaTheta * sigmaTheta;
  matrix(Acts::eQOP, Acts::eQOP )     = sigmaQOverP * sigmaQOverP;
  matrix(Acts::eT, Acts::eT )         = sigmaTime;
  
  return matrix;
}
  

// methods for converting TrkrCluster data to what Acts needs

TMatrixD PHActsTrkFitter::GetMvtxCovarLocal(const unsigned int layer, const unsigned int staveid, const unsigned int chipid, TMatrixD world_err)
{
  TMatrixD local_err(3,3);

  // rotate errors back to local coords
  double ladder_location[3] = {0.0, 0.0, 0.0};

  // returns the center of the sensor in world coordinates - used to get the ladder phi location
  CylinderGeom_Mvtx *layergeom = dynamic_cast<CylinderGeom_Mvtx *>(_geom_container_mvtx->GetLayerGeom(layer));
  layergeom->find_sensor_center(staveid, 0, 0, chipid, ladder_location);
  double ladderphi = atan2(ladder_location[1], ladder_location[0]);
  ladderphi += layergeom->get_stave_phi_tilt();

  local_err = TransformCovarToLocal(ladderphi, world_err);
  
  if(Verbosity() > 10)
    {
      for(int i=0;i<3;++i)
	for(int j = 0; j<3; ++j)
	  {
	    cout << "  " << i << "    " << j << " local_err " << local_err[i][j] << endl;
	  }
    }
    
  return local_err;
}



TMatrixD PHActsTrkFitter::GetInttCovarLocal(const unsigned int layer, const unsigned int ladderzid, const unsigned int ladderphiid, TMatrixD world_err)
{
  TMatrixD local_err(3,3);

  // rotate errors back to local coords
  double ladder_location[3] = {0.0, 0.0, 0.0};

  // rotate errors back to local coords 
  CylinderGeomIntt *layergeom = dynamic_cast<CylinderGeomIntt *>(_geom_container_intt->GetLayerGeom(layer));
  layergeom->find_segment_center(ladderzid, ladderphiid, ladder_location);
  double ladderphi = atan2(ladder_location[1], ladder_location[0]);

  local_err = TransformCovarToLocal(ladderphi, world_err);
  
  if(Verbosity() > 10)
    {
      for(int i=0;i<3;++i)
	for(int j = 0; j<3; ++j)
	  {
	    cout << "  INTT: " << i << "    " << j << " local_err " << local_err[i][j] << endl;
	  }
    }
    
  return local_err;
}

TMatrixD PHActsTrkFitter::TransformCovarToLocal(const double ladderphi, TMatrixD world_err)
{
  TMatrixD local_err(3,3);
  
  // this is the matrix that was used to rotate from local to global coords 
  TMatrixD ROT(3, 3);
  ROT[0][0] = cos(ladderphi);
  ROT[0][1] = -1.0 * sin(ladderphi);
  ROT[0][2] = 0.0;
  ROT[1][0] = sin(ladderphi);
  ROT[1][1] = cos(ladderphi);
  ROT[1][2] = 0.0;
  ROT[2][0] = 0.0; 
  ROT[2][1] = 0.0;
  ROT[2][2] = 1.0;
  // we want the inverse rotation
  ROT.Invert();
  
  TMatrixD ROT_T(3, 3);
  ROT_T.Transpose(ROT);
  
  local_err = ROT * world_err * ROT_T;

  return local_err;

}
 
int PHActsTrkFitter::End(PHCompositeNode* topNode)
{
  if(Verbosity() > 10)
    {
      std::cout<<"Finished PHActsTrkFitter"<<std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

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
  _geomanager = PHGeomUtility::GetTGeoManager(topNode);
  if(!_geomanager )
    {
      cout << PHWHERE << " Did not find TGeoManager, quit! " << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  _geom_container_mvtx = findNode::getClass<
    PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");
  if (!_geom_container_mvtx)
  {
    cout << PHWHERE << " CYLINDERGEOM_MVTX  node not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _geom_container_tpc =
    findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!_geom_container_tpc)
    {
      std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
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

  return Fun4AllReturnCodes::EVENT_OK;
}



