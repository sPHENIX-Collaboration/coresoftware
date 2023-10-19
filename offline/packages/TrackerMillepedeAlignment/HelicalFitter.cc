#include "HelicalFitter.h"

#include "Mille.h"

/// Tracking includes
#include <trackbase/TrkrDefs.h>                // for cluskey, getTrkrId, tpcId
#include <trackbase/TpcDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/InttDefs.h>
#include <trackbase/TrkrClusterv3.h>   
#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrClusterContainer.h>   
#include <trackbase/TrkrClusterCrossingAssoc.h>   
#include <trackbase/alignmentTransformationContainer.h>

#include <trackbase_historic/TrackSeed_v1.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/SvtxTrackSeed_v1.h>
#include <trackbase_historic/SvtxVertex.h>     // for SvtxVertex
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxTrack_v4.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>
#include <trackbase_historic/SvtxAlignmentState_v1.h>
#include <trackbase_historic/SvtxAlignmentStateMap_v1.h>
#include <trackbase_historic/SvtxTrackState_v1.h>

#include <Acts/Surfaces/PerigeeSurface.hpp>

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4Particle.h>  // for PHG4Particle
#include <g4main/PHG4HitDefs.h>  // for keytype

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/phool.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TF1.h>
#include <TNtuple.h>
#include <TFile.h>

#include <climits>                            // for UINT_MAX
#include <iostream>                            // for operator<<, basic_ostream
#include <cmath>                              // for fabs, sqrt
#include <set>                                 // for _Rb_tree_const_iterator
#include <utility>                             // for pair
#include <memory>

using namespace std;

//____________________________________________________________________________..
HelicalFitter::HelicalFitter(const std::string &name):
  SubsysReco(name)
  , PHParameterInterface(name)
  , _mille(nullptr)
{
  InitializeParameters();

  vertexPosition(0) = 0;
  vertexPosition(1) = 0;

  vtx_sigma(0) = 0.01;
  vtx_sigma(1) = 0.01;
}

//____________________________________________________________________________..
HelicalFitter::~HelicalFitter()
{

}

//____________________________________________________________________________..
int HelicalFitter::InitRun(PHCompositeNode *topNode)
{
  UpdateParametersWithMacro();

  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = CreateNodes(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  // Instantiate Mille and open output data file
  if(test_output)
    {
      _mille = new Mille(data_outfilename.c_str(), false);   // write text in data files, rather than binary, for debugging only
    }
  else
    {
      _mille = new Mille(data_outfilename.c_str()); 
    }

  // Write the steering file here, and add the data file path to it
  std::ofstream steering_file(steering_outfilename);
  steering_file << data_outfilename << std::endl;
  steering_file.close();

  if(make_ntuple)
    {
      //fout = new TFile("HF_ntuple.root","recreate");
      fout = new TFile(ntuple_outfilename.c_str(),"recreate");
      ntp  = new TNtuple("ntp","HF ntuple","event:trkid:layer:nsilicon:ntpc:nclus:trkrid:sector:side:subsurf:phi:glbl0:glbl1:glbl2:glbl3:glbl4:glbl5:sensx:sensy:sensz:normx:normy:normz:sensxideal:sensyideal:senszideal:normxideal:normyideal:normzideal:xglobideal:yglobideal:zglobideal:R:X0:Y0:Zs:Z0:xglob:yglob:zglob:xfit:yfit:zfit:pcax:pcay:pcaz:tangx:tangy:tangz:X:Y:fitX:fitY:dXdR:dXdX0:dXdY0:dXdZs:dXdZ0:dXdalpha:dXdbeta:dXdgamma:dXdx:dXdy:dXdz:dYdR:dYdX0:dYdY0:dYdZs:dYdZ0:dYdalpha:dYdbeta:dYdgamma:dYdx:dYdy:dYdz");

      track_ntp = new TNtuple("track_ntp","HF track ntuple","track_id:residual_x:residual_y:residualxsigma:residualysigma:dXdR:dXdX0:dXdY0:dXdZs:dXdZ0:dXdx:dXdy:dXdz:dYdR:dYdX0:dYdY0:dYdZs:dYdZ0:dYdx:dYdy:dYdz:xvtx:yvtx:zvtx:event_zvtx:track_phi:perigee_phi");

    }

 
  // print grouping setup to log file:
  std::cout << "HelicalFitter::InitRun: Surface groupings are mvtx " << mvtx_grp << " intt " << intt_grp << " tpc " << tpc_grp << " mms " << mms_grp << std::endl; 
  std::cout << " possible groupings are:" << std::endl
	    << " mvtx " 
	    << AlignmentDefs::mvtxGrp::snsr << "  " 
	    << AlignmentDefs::mvtxGrp::stv << "  " 
	    << AlignmentDefs::mvtxGrp::mvtxlyr << "  "
	    << AlignmentDefs::mvtxGrp::clamshl << "  " << std::endl
	    << " intt " 
	    << AlignmentDefs::inttGrp::chp << "  "
	    << AlignmentDefs::inttGrp::lad << "  "
	    << AlignmentDefs::inttGrp::inttlyr << "  "
	    << AlignmentDefs::inttGrp::inttbrl << "  " << std::endl
	    << " tpc " 
	    << AlignmentDefs::tpcGrp::htst << "  "
	    << AlignmentDefs::tpcGrp::sctr << "  "
	    << AlignmentDefs::tpcGrp::tp << "  " << std::endl
	    << " mms " 
	    << AlignmentDefs::mmsGrp::tl << "  "
	    << AlignmentDefs::mmsGrp::mm << "  " << std::endl;
 
  event=-1;

  return ret;
}

//_____________________________________________________________________
void HelicalFitter::SetDefaultParameters()
{


  return;
}

//____________________________________________________________________________..
int HelicalFitter::process_event(PHCompositeNode*)
{
  // _track_map_tpc contains the TPC seed track stubs
  // _track_map_silicon contains the silicon seed track stubs
  // _svtx_seed_map contains the combined silicon and tpc track seeds

  event++;

  if(Verbosity() > 0)
    cout << PHWHERE 
	 << " TPC seed map size " << _track_map_tpc->size() 
	 << " Silicon seed map size "  << _track_map_silicon->size() 
	 << endl;

  if(_track_map_silicon->size() == 0 && _track_map_tpc->size() == 0)
    return Fun4AllReturnCodes::EVENT_OK;

  // Decide whether we want to make a helical fit for silicon or TPC
  unsigned int maxtracks = 0; 
  unsigned int nsilicon = 0;
  unsigned int ntpc = 0;
  unsigned int nclus = 0;
  std::vector<std::vector<Acts::Vector3>> cumulative_global_vec;
  std::vector<std::vector<TrkrDefs::cluskey>> cumulative_cluskey_vec;
  std::vector<std::vector<float>> cumulative_fitpars_vec;
  std::vector<Acts::Vector3> cumulative_vertex;
  std::vector<TrackSeed> cumulative_someseed;
  std::vector<SvtxTrack_v4> cumulative_newTrack;


  if(fittpc) { maxtracks =  _track_map_tpc->size();  }
  if(fitsilicon)  { maxtracks =  _track_map_silicon->size(); }
  for(unsigned int trackid = 0; trackid < maxtracks; ++trackid)
    {
       TrackSeed* tracklet = nullptr;
      if(fitsilicon) {  tracklet = _track_map_silicon->get(trackid); }
      else if(fittpc) {  tracklet = _track_map_tpc->get(trackid);	 }
      if(!tracklet) { continue; }

      std::vector<Acts::Vector3> global_vec;
      std::vector<TrkrDefs::cluskey> cluskey_vec;

      // Get a vector of cluster keys from the tracklet  
      getTrackletClusterList(tracklet, cluskey_vec);
      // store cluster global positions in a vector global_vec and cluskey_vec
      TrackFitUtils::getTrackletClusters(_tGeometry, _cluster_map, global_vec, cluskey_vec);   
     
      correctTpcGlobalPositions( global_vec, cluskey_vec);

      std::vector<float> fitpars =  TrackFitUtils::fitClusters(global_vec, cluskey_vec);       // do helical fit

      if(fitpars.size() == 0) continue;  // discard this track, not enough clusters to fit

      if(Verbosity() > 1)  
	{ std::cout << " Track " << trackid   << " radius " << fitpars[0] << " X0 " << fitpars[1]<< " Y0 " << fitpars[2]
		    << " zslope " << fitpars[3]  << " Z0 " << fitpars[4] << std::endl; }
      
      //// Create a track map for diagnostics
      SvtxTrack_v4 newTrack;
      newTrack.set_id(trackid);
      if(fitsilicon) { newTrack.set_silicon_seed(tracklet); }
      else if(fittpc) {  newTrack.set_tpc_seed(tracklet); }
      
      // if a full track is requested, get the silicon clusters too and refit
      if(fittpc && fitfulltrack)
	{
	  // this associates silicon clusters and adds them to the vectors
	  ntpc = cluskey_vec.size();
	  nsilicon = TrackFitUtils::addSiliconClusters(fitpars, dca_cut, _tGeometry, _cluster_map, global_vec, cluskey_vec);
	  if(nsilicon < 3) continue;  // discard this TPC seed, did not get a good match to silicon
	  auto trackseed = std::make_unique<TrackSeed_v1>();
	  for(auto& ckey : cluskey_vec)
	    {
	      if(TrkrDefs::getTrkrId(ckey) == TrkrDefs::TrkrId::mvtxId or
		 TrkrDefs::getTrkrId(ckey) == TrkrDefs::TrkrId::inttId)
		{
		  trackseed->insert_cluster_key(ckey);
		}
	    }

	  newTrack.set_silicon_seed(trackseed.get());
	  
	  // fit the full track now
	  fitpars.clear();
	  fitpars =  TrackFitUtils::fitClusters(global_vec, cluskey_vec);       // do helical fit
          if(fitpars.size() == 0) continue;  // discard this track, fit failed

	  if(Verbosity() > 1)  
	    { std::cout << " Full track " << trackid   << " radius " << fitpars[0] << " X0 " << fitpars[1]<< " Y0 " << fitpars[2]
			<< " zslope " << fitpars[3]  << " Z0 " << fitpars[4] << std::endl; }
	}
      else if(fitsilicon)
	{
	  nsilicon = cluskey_vec.size();
	}
      else if(fittpc && !fitfulltrack)
	{
	  ntpc = cluskey_vec.size();
	}
 
      Acts::Vector3 beamline(0,0,0);
      Acts::Vector2 pca2d = TrackFitUtils::get_circle_point_pca(fitpars[0], fitpars[1], fitpars[2], beamline);
      Acts::Vector3 track_vtx (pca2d(0),pca2d(1),fitpars[4]);

      newTrack.set_crossing(tracklet->get_crossing());
      newTrack.set_id(trackid);

      /// use the track seed functions to help get the track trajectory values
      /// in the usual coordinates
  
      TrackSeed_v1 someseed;
      for(auto& ckey : cluskey_vec)
	{ someseed.insert_cluster_key(ckey); }
      someseed.set_qOverR(tracklet->get_charge() / fitpars[0]);
 
      someseed.set_X0(fitpars[1]);
      someseed.set_Y0(fitpars[2]);
      someseed.set_Z0(fitpars[4]);
      someseed.set_slope(fitpars[3]);

      newTrack.set_x(someseed.get_x());
      newTrack.set_y(someseed.get_y());
      newTrack.set_z(someseed.get_z());
      newTrack.set_px(someseed.get_px(_cluster_map,_tGeometry));
      newTrack.set_py(someseed.get_py(_cluster_map,_tGeometry));
      newTrack.set_pz(someseed.get_pz());  
      newTrack.set_charge(tracklet->get_charge());

      nclus = ntpc+nsilicon;

      // some basic track quality requirements
      if(fittpc && ntpc < 35) 
	{
	  if(Verbosity() > 1) { std::cout << " reject this track, ntpc = " << ntpc << std::endl; } 
	  continue;
	}
      if((fitsilicon || fitfulltrack) && nsilicon < 4) 
	{
	  if(Verbosity() > 1) { std::cout << " reject this track, nsilicon = " << nsilicon << std::endl; } 
	  continue; 
	}

      cumulative_global_vec.push_back(global_vec);
      cumulative_cluskey_vec.push_back(cluskey_vec);
      cumulative_vertex.push_back(track_vtx);
      cumulative_fitpars_vec.push_back(fitpars);
      cumulative_someseed.push_back(someseed);
      cumulative_newTrack.push_back(newTrack);
    }

  //terminate loop over tracks
  //Collect fitpars for each track by intializing array of size maxtracks and populaating thorughout the loop
  //Then start new loop over tracks and for each track go over clsutaer 
  // make vector of global_vecs
  float xsum = 0;
  float ysum = 0;
  float zsum = 0;
  unsigned int accepted_tracks = cumulative_fitpars_vec.size();

  //std::cout<<"accepted_tracks: " << accepted_tracks << std::endl;
  for(unsigned int trackid = 0; trackid < accepted_tracks; ++trackid)
    {
      xsum += cumulative_vertex[trackid][0];
      ysum += cumulative_vertex[trackid][1];
      zsum += cumulative_vertex[trackid][2];
    }
  Acts::Vector3 averageVertex (xsum/accepted_tracks,ysum/accepted_tracks,zsum/accepted_tracks); 
      
  for(unsigned int trackid = 0; trackid < accepted_tracks; ++trackid)
    { 
      auto global_vec  = cumulative_global_vec[trackid];
      auto cluskey_vec = cumulative_cluskey_vec[trackid];
      auto fitpars     = cumulative_fitpars_vec[trackid];
      auto someseed    = cumulative_someseed[trackid];
      auto newTrack    = cumulative_newTrack[trackid];
      //std::cout << "trackid " << trackid << " get id " <<newTrack.get_id()<< std::endl;
      SvtxAlignmentStateMap::StateVec statevec;
      
      // get the residuals and derivatives for all clusters
      for(unsigned int ivec=0;ivec<global_vec.size(); ++ivec)
	{
	  auto global  = global_vec[ivec];
	  auto cluskey = cluskey_vec[ivec];
	  auto cluster = _cluster_map->findCluster(cluskey);
	  if(!cluster) { continue;}
	  
	  unsigned int trkrid = TrkrDefs::getTrkrId(cluskey);

	  // What we need now is to find the point on the surface at which the helix would intersect
	  // If we have that point, we can transform the fit back to local coords
	  // we have fitpars for the helix, and the cluster key - from which we get the surface

	  Surface surf = _tGeometry->maps().getSurface(cluskey, cluster);
	  Acts::Vector3 helix_pca(0,0,0);
	  Acts::Vector3 helix_tangent(0,0,0);
	  Acts::Vector3 fitpoint = get_helix_surface_intersection(surf, fitpars, global, helix_pca, helix_tangent);

	  // fitpoint is the point where the helical fit intersects the plane of the surface
	  // Now transform the helix fitpoint to local coordinates to compare with cluster local coordinates
	  Acts::Vector3 fitpoint_local = surf->transform(_tGeometry->geometry().getGeoContext()).inverse() * (fitpoint *  Acts::UnitConstants::cm);
	  
	  fitpoint_local /= Acts::UnitConstants::cm;

	  auto xloc = cluster->getLocalX();  // in cm
	  auto zloc = cluster->getLocalY();	  

	  if(trkrid == TrkrDefs::tpcId) { zloc = convertTimeToZ(cluskey, cluster); }
	  //float yloc = 0.0;   // Because the fitpoint is on the surface, y will always be zero in local coordinates

	  Acts::Vector2 residual(xloc - fitpoint_local(0), zloc - fitpoint_local(1));
 
	  unsigned int layer = TrkrDefs::getLayer(cluskey_vec[ivec]);	  
	  float phi          =  atan2(global(1), global(0));

	  SvtxTrackState_v1 svtxstate(fitpoint.norm());
	  svtxstate.set_x(fitpoint(0));
	  svtxstate.set_y(fitpoint(1));
	  svtxstate.set_z(fitpoint(2));
	  auto tangent = get_helix_tangent(fitpars, global);  
	  svtxstate.set_px(someseed.get_p() * tangent.second.x());
	  svtxstate.set_py(someseed.get_p() * tangent.second.y());
	  svtxstate.set_pz(someseed.get_p() * tangent.second.z());
	  newTrack.insert_state(&svtxstate);
	  
	  if(Verbosity() > 1) {
	    Acts::Vector3 loc_check =  surf->transform(_tGeometry->geometry().getGeoContext()).inverse() * (global *  Acts::UnitConstants::cm);
	    loc_check /= Acts::UnitConstants::cm;
	    std::cout << "    layer " << layer << std::endl
		      << " cluster global " << global(0) << " " << global(1) << " " << global(2) << std::endl
		      << " fitpoint " << fitpoint(0) << " " << fitpoint(1) << " " << fitpoint(2) << std::endl
		      << " fitpoint_local " << fitpoint_local(0) << " " << fitpoint_local(1) << " " << fitpoint_local(2) << std::endl  
		      << " cluster local x " << cluster->getLocalX() << " cluster local y " << cluster->getLocalY() << std::endl
		      << " cluster global to local x " << loc_check(0) << " local y " << loc_check(1) << "  local z " << loc_check(2) << std::endl
		      << " cluster local residual x " << residual(0) << " cluster local residual y " <<residual(1) << std::endl;
	  }

	  if(Verbosity() > 1) 
	    {
	      Acts::Transform3 transform = surf->transform(_tGeometry->geometry().getGeoContext());
	      std::cout << "Transform is:" << std::endl;
	      std::cout <<  transform.matrix() << std::endl;
	      Acts::Vector3 loc_check = surf->transform(_tGeometry->geometry().getGeoContext()).inverse() * (global *  Acts::UnitConstants::cm);
	      loc_check /= Acts::UnitConstants::cm;
	      unsigned int sector = TpcDefs::getSectorId(cluskey_vec[ivec]);	  
	      unsigned int side = TpcDefs::getSide(cluskey_vec[ivec]);	  
	      std::cout << "    layer " << layer << " sector " << sector << " side " << side << " subsurf " <<   cluster->getSubSurfKey() << std::endl
			<< " cluster global " << global(0) << " " << global(1) << " " << global(2) << std::endl
			<< " fitpoint " << fitpoint(0) << " " << fitpoint(1) << " " << fitpoint(2) << std::endl
			<< " fitpoint_local " << fitpoint_local(0) << " " << fitpoint_local(1) << " " << fitpoint_local(2) << std::endl  
			<< " cluster local x " << cluster->getLocalX() << " cluster local y " << cluster->getLocalY() << std::endl
			<< " cluster global to local x " << loc_check(0) << " local y " << loc_check(1) << "  local z " << loc_check(2) << std::endl
			<< " cluster local residual x " << residual(0) << " cluster local residual y " <<residual(1) << std::endl;
	    }
	  
	  // need standard deviation of measurements
	  Acts::Vector2 clus_sigma = getClusterError(cluster, cluskey, global);
	  if(isnan(clus_sigma(0)) || isnan(clus_sigma(1)))  { continue; }

	  int glbl_label[AlignmentDefs::NGL];
	  if(layer < 3)
	    {
	      AlignmentDefs::getMvtxGlobalLabels(surf, glbl_label, mvtx_grp);	      
	    }
	  else if(layer > 2 && layer < 7)
	    {
	      AlignmentDefs::getInttGlobalLabels(surf, glbl_label, intt_grp);
	    }
	  else if (layer < 55)
	    {
	      AlignmentDefs::getTpcGlobalLabels(surf, cluskey, glbl_label, tpc_grp);
	    }
	  else
	    {
	      continue;
	    }
	    
	  // These derivatives are for the local parameters
	  float lcl_derivativeX[AlignmentDefs::NLC];
	  float lcl_derivativeY[AlignmentDefs::NLC];

	  getLocalDerivativesXY(surf, global, fitpars, lcl_derivativeX, lcl_derivativeY, layer);

	  // The global derivs dimensions are [alpha/beta/gamma](x/y/z)
	  float glbl_derivativeX[AlignmentDefs::NGL];
	  float glbl_derivativeY[AlignmentDefs::NGL];
	  getGlobalDerivativesXY(surf, global, fitpoint, fitpars, glbl_derivativeX, glbl_derivativeY, layer);

	  auto alignmentstate = std::make_unique<SvtxAlignmentState_v1>();
	  alignmentstate->set_residual(residual);
	  alignmentstate->set_cluster_key(cluskey);
	  SvtxAlignmentState::GlobalMatrix svtxglob = 
	    SvtxAlignmentState::GlobalMatrix::Zero();
	  SvtxAlignmentState::LocalMatrix svtxloc = 
	    SvtxAlignmentState::LocalMatrix::Zero();
	  for(int i=0; i<AlignmentDefs::NLC; i++)
	    {
	      svtxloc(0,i) = lcl_derivativeX[i];
	      svtxloc(1,i) = lcl_derivativeY[i];
	    }
	  for(int i=0; i<AlignmentDefs::NGL; i++)
	    {
	      svtxglob(0,i) = glbl_derivativeX[i];
	      svtxglob(1,i) = glbl_derivativeY[i];
	    }

	  alignmentstate->set_local_derivative_matrix(svtxloc);
	  alignmentstate->set_global_derivative_matrix(svtxglob);
	  
	  statevec.push_back(alignmentstate.release());
	  
	  for(unsigned int i = 0; i < AlignmentDefs::NGL; ++i) 
	    {
	      if(trkrid == TrkrDefs::mvtxId)
		{
		  // need stave to get clamshell
		  auto stave  = MvtxDefs::getStaveId(cluskey_vec[ivec]);
		  auto clamshell = AlignmentDefs::getMvtxClamshell(layer, stave);
		  if( is_layer_param_fixed(layer, i) || is_mvtx_layer_fixed(layer,clamshell) )
		    {
		      glbl_derivativeX[i] = 0;
		      glbl_derivativeY[i] = 0;
		    }
		}

	      if(trkrid == TrkrDefs::inttId)
		{
		  if( is_layer_param_fixed(layer, i) || is_intt_layer_fixed(layer) )
		    {
		      glbl_derivativeX[i] = 0;
		      glbl_derivativeY[i] = 0;
		    }
		}


	      if(trkrid == TrkrDefs::tpcId)
		{
		  unsigned int sector = TpcDefs::getSectorId(cluskey_vec[ivec]);	  
		  unsigned int side   = TpcDefs::getSide(cluskey_vec[ivec]);	  
		  if(is_layer_param_fixed(layer, i) || is_tpc_sector_fixed(layer, sector, side))
		    {
		      glbl_derivativeX[i] = 0;
		      glbl_derivativeY[i] = 0;
		    }
		}
	    }

	  // Add the measurement separately for each coordinate direction to Mille
	  // set the derivatives non-zero only for parameters we want to be optimized
	  // local parameter numbering is arbitrary:
	  float errinf = 1.0;

	  if(_layerMisalignment.find(layer) != _layerMisalignment.end())
	    {
	      errinf = _layerMisalignment.find(layer)->second;
	    }
	  if(make_ntuple)
	    {
	      // get the local parameters using the ideal transforms
	      alignmentTransformationContainer::use_alignment = false;
	      Acts::Vector3 ideal_center =  surf->center(_tGeometry->geometry().getGeoContext()) * 0.1;  
	      Acts::Vector3 ideal_norm   = -surf->normal(_tGeometry->geometry().getGeoContext());  
	      Acts::Vector3 ideal_local(xloc, zloc, 0.0); // cm
	      Acts::Vector3 ideal_glob = surf->transform(_tGeometry->geometry().getGeoContext())*(ideal_local * Acts::UnitConstants::cm);
	      ideal_glob /= Acts::UnitConstants::cm;
	      alignmentTransformationContainer::use_alignment = true;

	      Acts::Vector3 sensorCenter = surf->center(_tGeometry->geometry().getGeoContext()) * 0.1; // cm
	      Acts::Vector3 sensorNormal = -surf->normal(_tGeometry->geometry().getGeoContext());
	      unsigned int sector  = TpcDefs::getSectorId(cluskey_vec[ivec]);	  
	      unsigned int side    = TpcDefs::getSide(cluskey_vec[ivec]);	  	      
	      unsigned int subsurf = cluster->getSubSurfKey();
	      if(layer < 3)
		{
		  sector  = MvtxDefs::getStaveId(cluskey_vec[ivec]);
		  subsurf = MvtxDefs::getChipId(cluskey_vec[ivec]);
		}
	      else if(layer >2 && layer < 7)
		{
		  sector  = InttDefs::getLadderPhiId(cluskey_vec[ivec]);
		  subsurf = InttDefs::getLadderZId(cluskey_vec[ivec]);
		}
	      float ntp_data[75] = {
		(float) event, (float) trackid,
		(float) layer, (float) nsilicon, (float) ntpc, (float) nclus, (float) trkrid,  (float) sector,  (float) side,
		(float) subsurf, phi,
		(float) glbl_label[0], (float) glbl_label[1], (float) glbl_label[2], (float) glbl_label[3], (float) glbl_label[4], (float) glbl_label[5], 
		(float) sensorCenter(0), (float) sensorCenter(1), (float) sensorCenter(2),
		(float) sensorNormal(0), (float) sensorNormal(1), (float) sensorNormal(2),
		(float) ideal_center(0), (float) ideal_center(1), (float) ideal_center(2),
		(float) ideal_norm(0), (float) ideal_norm(1), (float) ideal_norm(2),
		(float) ideal_glob(0), (float) ideal_glob(1), (float) ideal_glob(2),
		(float) fitpars[0], (float) fitpars[1], (float) fitpars[2], (float) fitpars[3], (float) fitpars[4], 
		(float) global(0), (float) global(1), (float) global(2),
		(float) fitpoint(0), (float) fitpoint(1), (float) fitpoint(2), 
		(float) helix_pca(0), (float) helix_pca(1), (float) helix_pca(2),
		(float) helix_tangent(0), (float) helix_tangent(1), (float) helix_tangent(2),
		xloc,zloc, (float) fitpoint_local(0), (float) fitpoint_local(1), 
		lcl_derivativeX[0],lcl_derivativeX[1],lcl_derivativeX[2],lcl_derivativeX[3],lcl_derivativeX[4],
		glbl_derivativeX[0],glbl_derivativeX[1],glbl_derivativeX[2],glbl_derivativeX[3],glbl_derivativeX[4],glbl_derivativeX[5],
		lcl_derivativeY[0],lcl_derivativeY[1],lcl_derivativeY[2],lcl_derivativeY[3],lcl_derivativeY[4],
		glbl_derivativeY[0],glbl_derivativeY[1],glbl_derivativeY[2],glbl_derivativeY[3],glbl_derivativeY[4],glbl_derivativeY[5] };

	      ntp->Fill(ntp_data);

	      if(Verbosity() > 2)
		{
		  for(int i=0;i<34;++i)
		    {
		      std::cout << ntp_data[i] << "  " ;
		    }
		  std::cout << std::endl;
		}
	    }

	  // add some cluster cuts
	  if(residual(0) > 0.2)  continue;   // 2 mm cut
	  if(residual(1) > 0.2)  continue;   // 2 mm cut

	  if( !isnan(residual(0)) && clus_sigma(0) < 1.0)  // discards crazy clusters
	    {
	      _mille->mille(AlignmentDefs::NLC,lcl_derivativeX,AlignmentDefs::NGL,glbl_derivativeX,glbl_label,residual(0), errinf*clus_sigma(0));
	    }
	  if(!isnan(residual(1)) && clus_sigma(1) < 1.0 && trkrid != TrkrDefs::inttId)
	    {
	      _mille->mille(AlignmentDefs::NLC, lcl_derivativeY,AlignmentDefs::NGL,glbl_derivativeY,glbl_label,residual(1), errinf*clus_sigma(1));
	    }
	}

      m_alignmentmap->insertWithKey(trackid, statevec);
      m_trackmap->insertWithKey(&newTrack, trackid);
	
      //calculate vertex residual with perigee surface 
      Acts::Vector3 event_vtx(0,0,averageVertex(2)); 
      
      // The residual for the vtx case is (event vtx - track vtx) 
      // that is -dca
      float dca3dxy;
      float dca3dz; 
      float dca3dxysigma; 
      float dca3dzsigma;	  
      get_dca(newTrack,dca3dxy,dca3dz,dca3dxysigma,dca3dzsigma,event_vtx);
      Acts::Vector2 vtx_residual(-dca3dxy, -dca3dz);
      
      float lclvtx_derivativeX[AlignmentDefs::NLC];
      float lclvtx_derivativeY[AlignmentDefs::NLC];
      getLocalVtxDerivativesXY(newTrack, event_vtx, fitpars, lclvtx_derivativeX, lclvtx_derivativeY);
      
      // The global derivs dimensions are [alpha/beta/gamma](x/y/z)
      float glblvtx_derivativeX[3];
      float glblvtx_derivativeY[3];
      getGlobalVtxDerivativesXY(newTrack, event_vtx, glblvtx_derivativeX, glblvtx_derivativeY);
  
      if(use_event_vertex)
	{	  
	  if(Verbosity() > 3)
	    {
	      std::cout << "vertex info for track " << trackid << " with charge " << newTrack.get_charge() << std::endl;
	      
	      std::cout << "vertex is " << event_vtx.transpose() << std::endl;
	      std::cout << "vertex residuals " << vtx_residual.transpose() 
			<< std::endl;
	      std::cout << "local derivatives " << std::endl;
	      for(int i=0; i<AlignmentDefs::NLC; i++)
		std::cout << lclvtx_derivativeX[i] << ", ";
	      std::cout << std::endl;
	      for(int i=0; i<AlignmentDefs::NLC; i++)
		std::cout << lclvtx_derivativeY[i] << ", ";
	      std::cout << "global vtx derivaties " << std::endl;
	      for(int i=0; i<3; i++) std::cout << glblvtx_derivativeX[i] << ", ";
	      std::cout << std::endl;
	      for(int i=0; i<3; i++) std::cout << glblvtx_derivativeY[i] << ", ";
	    }

	  // add some track cuts
	  if(fabs(newTrack.get_z() - event_vtx(2)) > 0.2) continue;  // 2 mm cut
	  if(fabs(newTrack.get_x()) > 0.2) continue;  // 2 mm cut
	  if(fabs(newTrack.get_y()) > 0.2) continue;  // 2 mm cut
	  
	  if(!isnan(vtx_residual(0)))
	    {
	      _mille->mille(AlignmentDefs::NLC,lclvtx_derivativeX,AlignmentDefs::NGLVTX,glblvtx_derivativeX,AlignmentDefs::glbl_vtx_label,vtx_residual(0), vtx_sigma(0));
	    }    
	  if(!isnan(vtx_residual(1)))
	    {  
	      _mille->mille(AlignmentDefs::NLC,lclvtx_derivativeY,AlignmentDefs::NGLVTX,glblvtx_derivativeY,AlignmentDefs::glbl_vtx_label,vtx_residual(1), vtx_sigma(1));
	    }
	}

      if(make_ntuple)
	{          
	  Acts::Vector3 mom(newTrack.get_px(),newTrack.get_py(),newTrack.get_pz());
	  Acts::Vector3 r = mom.cross(Acts::Vector3(0.,0.,1.));
	  float perigee_phi       = atan2(r(1), r(0));
	  float track_phi = atan2(newTrack.get_py(), newTrack.get_px());
	  //	  float ntp_data[30] = {(float) trackid,dca3dxy,dca3dz,(float) vtx_sigma(0),(float) vtx_sigma(1),
	  float ntp_data[30] = {(float) trackid,(float) vtx_residual(0),(float) vtx_residual(1),(float) vtx_sigma(0),(float) vtx_sigma(1),
				lclvtx_derivativeX[0],lclvtx_derivativeX[1],lclvtx_derivativeX[2],lclvtx_derivativeX[3],lclvtx_derivativeX[4],
				glblvtx_derivativeX[0],glblvtx_derivativeX[1],glblvtx_derivativeX[2],
				lclvtx_derivativeY[0],lclvtx_derivativeY[1],lclvtx_derivativeY[2],lclvtx_derivativeY[3],lclvtx_derivativeY[4],
				glblvtx_derivativeY[0],glblvtx_derivativeY[1],glblvtx_derivativeY[2],
				newTrack.get_x(), newTrack.get_y(), newTrack.get_z(),(float) event_vtx(2),track_phi, perigee_phi};
	  
	  track_ntp->Fill(ntp_data);
	}
      
      if(Verbosity()>1)
	{
	  std::cout << "vtx_residual xy: " << vtx_residual(0)<< " vtx_residual z: " << vtx_residual(1) << " vtx_sigma xy: " << vtx_sigma(0) << " vtx_sigma z: " << vtx_sigma(1) << std::endl; 
	  std::cout << "track_x" << newTrack.get_x()<<"track_y" << newTrack.get_y()<<"track_z" << newTrack.get_z()<<std::endl;
	}

  
      // close out this track	
      _mille->end();
          
    }  // end loop over tracks
  
  
  return Fun4AllReturnCodes::EVENT_OK;
}

Acts::Vector3 HelicalFitter::get_helix_surface_intersection(Surface surf, std::vector<float>& fitpars, Acts::Vector3 global)
{
  // we want the point where the helix intersects the plane of the surface
  // get the plane of the surface
  Acts::Vector3 sensorCenter      = surf->center(_tGeometry->geometry().getGeoContext()) * 0.1;  // convert to cm
  Acts::Vector3 sensorNormal    = -surf->normal(_tGeometry->geometry().getGeoContext());
  sensorNormal /= sensorNormal.norm();

  // there are analytic solutions for a line-plane intersection.
  // to use this, need to get the vector tangent to the helix near the measurement and a point on it.
  std::pair<Acts::Vector3, Acts::Vector3> line =  get_helix_tangent(fitpars, global);
  Acts::Vector3 pca = line.first;
  Acts::Vector3 tangent = line.second;

  Acts::Vector3 intersection = get_line_plane_intersection(pca, tangent, sensorCenter, sensorNormal);

  return intersection;
}

Acts::Vector3 HelicalFitter::get_helix_surface_intersection(Surface surf, std::vector<float>& fitpars, Acts::Vector3 global, Acts::Vector3& pca, Acts::Vector3& tangent)
{
  // we want the point where the helix intersects the plane of the surface

  // get the plane of the surface
  Acts::Vector3 sensorCenter      = surf->center(_tGeometry->geometry().getGeoContext()) * 0.1;  // convert to cm
  Acts::Vector3 sensorNormal    = -surf->normal(_tGeometry->geometry().getGeoContext());
  sensorNormal /= sensorNormal.norm();

  // there are analytic solutions for a line-plane intersection.
  // to use this, need to get the vector tangent to the helix near the measurement and a point on it.
  std::pair<Acts::Vector3, Acts::Vector3> line =  get_helix_tangent(fitpars, global);
  pca = line.first;
  tangent = line.second;

  Acts::Vector3 intersection = get_line_plane_intersection(pca, tangent, sensorCenter, sensorNormal);

  return intersection;
}


Acts::Vector3 HelicalFitter::get_helix_vtx(Acts::Vector3 event_vtx, const std::vector<float>& fitpars)
{  
  Acts::Vector2 pca2d = TrackFitUtils::get_circle_point_pca(fitpars[0], fitpars[1], fitpars[2], event_vtx);
  Acts::Vector3 helix_vtx (pca2d(0),pca2d(1),fitpars[4]);

  return helix_vtx;
}


Acts::Vector3 HelicalFitter::get_line_plane_intersection(Acts::Vector3 PCA, Acts::Vector3 tangent, Acts::Vector3 sensor_center, Acts::Vector3 sensor_normal)
{
  // get the intersection of the line made by PCA and tangent with the plane of the sensor

  // For a point on the line
  // p = PCA + d * tangent;
  // for a point on the plane
  // (p - sensor_center).sensor_normal = 0

  // The solution is:
  float d = (sensor_center - PCA).dot(sensor_normal) / tangent.dot(sensor_normal);
  Acts::Vector3 intersection = PCA + d * tangent;

  return intersection;
}


std::pair<Acts::Vector3, Acts::Vector3> HelicalFitter::get_helix_tangent(const std::vector<float>& fitpars, Acts::Vector3 global)
{
  auto pair = TrackFitUtils::get_helix_tangent(fitpars, global);
  /*
    save for posterity purposes
  if(Verbosity() > 2)
    {
      // different method for checking:
      // project the circle PCA vector an additional small amount and find the helix PCA to that point 
      
      float projection = 0.25;  // cm
      Acts::Vector3 second_point = pca + projection * pca/pca.norm();
      Acts::Vector2 second_point_pca_circle = TrackFitUtils::get_circle_point_pca(radius, x0, y0, second_point);
      float second_point_pca_z = second_point_pca_circle.norm() * zslope + z0;
      Acts::Vector3 second_point_pca2(second_point_pca_circle(0), second_point_pca_circle(1), second_point_pca_z);
      Acts::Vector3 tangent2 = (second_point_pca2 - pca) /  (second_point_pca2 - pca).norm();
      Acts::Vector3 final_pca2 = getPCALinePoint(global, tangent2, pca);
    
      std::cout << " get_helix_tangent: getting tangent at angle_pca: " << angle_pca * 180.0 / M_PI << std::endl 
		<< " original first point pca                      " << pca(0) << "  " << pca(1) << "  " << pca(2) << std::endl
		<< " original second point pca  " << second_point_pca(0) << "  " << second_point_pca(1) << "  " << second_point_pca(2) << std::endl
		<< " original tangent " << tangent(0) << "  " << tangent(1) << "  " << tangent(2) << std::endl	
		<< " original final pca from line " << final_pca(0) << "  " << final_pca(1) << "  " << final_pca(2) << std::endl;

      if(Verbosity() > 3)
	{
	  std::cout	<< "    Check: 2nd point pca meth 2 "<< second_point_pca2(0)<< "  "<< second_point_pca2(1) << "  "<< second_point_pca2(2) << std::endl
			<< "    check tangent " << tangent2(0) << "  " << tangent2(1) << "  " << tangent2(2) << std::endl	
			<< "    check final pca from line " << final_pca2(0) << "  " << final_pca2(1) << "  " << final_pca2(2) 
			<< std::endl;
	}
      
    }
  */


  return pair;
}
  
int HelicalFitter::End(PHCompositeNode* )
{
  // closes output file in destructor
  delete _mille;

  if(make_ntuple)
    {
      fout->Write();
      fout->Close();
    }

  return Fun4AllReturnCodes::EVENT_OK;
}
int HelicalFitter::CreateNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
    {
      std::cerr << "DST node is missing, quitting" << std::endl;
      throw std::runtime_error("Failed to find DST node in PHActsTrkFitter::createNodes");
    }
  
  PHNodeIterator dstIter(topNode);
  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "SVTX"));
  if (!svtxNode)
    {
      svtxNode = new PHCompositeNode("SVTX");
      dstNode->addNode(svtxNode);
    }
  
  m_trackmap = findNode::getClass<SvtxTrackMap>(topNode, "HelicalFitterTrackMap");
  if(!m_trackmap)
    {
      m_trackmap = new SvtxTrackMap_v2;
      PHIODataNode<PHObject> *node = new PHIODataNode<PHObject>(m_trackmap,"HelicalFitterTrackMap","PHObject");
      svtxNode->addNode(node);
    }

  m_alignmentmap = findNode::getClass<SvtxAlignmentStateMap>(topNode,"HelicalFitterAlignmentStateMap");
  if(!m_alignmentmap)
    {
      m_alignmentmap = new SvtxAlignmentStateMap_v1;
      PHIODataNode<PHObject> *node = new PHIODataNode<PHObject>(m_alignmentmap,"HelicalFitterAlignmentStateMap","PHObject");
      svtxNode->addNode(node);
    }

  return Fun4AllReturnCodes::EVENT_OK;
}
int  HelicalFitter::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get additional objects off the Node Tree
  //---------------------------------

  _track_map_silicon = findNode::getClass<TrackSeedContainer>(topNode, _silicon_track_map_name);
  if (!_track_map_silicon)
    {
      cerr << PHWHERE << " ERROR: Can't find SiliconTrackSeedContainer " << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  _track_map_tpc = findNode::getClass<TrackSeedContainer>(topNode, _track_map_name);
  if (!_track_map_tpc)
    {
      cerr << PHWHERE << " ERROR: Can't find " << _track_map_name.c_str() << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
    {
      std::cout << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  _tGeometry = findNode::getClass<ActsGeometry>(topNode,"ActsGeometry");
  if(!_tGeometry)
    {
      std::cout << PHWHERE << "Error, can't find acts tracking geometry" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
} 

Acts::Vector3 HelicalFitter::get_helix_pca(std::vector<float>& fitpars, Acts::Vector3 global)
{
  return TrackFitUtils::get_helix_pca(fitpars, global);
}


Acts::Vector3 HelicalFitter::getPCALinePoint(Acts::Vector3 global, Acts::Vector3 tangent, Acts::Vector3 posref)
{
  // Approximate track with a straight line consisting of the state position posref and the vector (px,py,pz)   

  // The position of the closest point on the line to global is:
  // posref + projection of difference between the point and posref on the tangent vector
  Acts::Vector3 pca = posref + ( (global - posref).dot(tangent) ) * tangent;

  return pca;
}


float HelicalFitter::convertTimeToZ(TrkrDefs::cluskey cluster_key, TrkrCluster *cluster)
{
  // must convert local Y from cluster average time of arival to local cluster z position
  double drift_velocity = _tGeometry->get_drift_velocity();
  double zdriftlength = cluster->getLocalY() * drift_velocity;
  double surfCenterZ = 52.89; // 52.89 is where G4 thinks the surface center is
  double zloc = surfCenterZ - zdriftlength;   // converts z drift length to local z position in the TPC in north
  unsigned int side = TpcDefs::getSide(cluster_key);
  if(side == 0) zloc = -zloc;
  float z = zloc;  // in cm
 
  return z; 
}

void HelicalFitter::makeTpcGlobalCorrections(TrkrDefs::cluskey cluster_key, short int crossing, Acts::Vector3& global)
{
  // make all corrections to global position of TPC cluster
  unsigned int side = TpcDefs::getSide(cluster_key);
  float z = m_clusterCrossingCorrection.correctZ(global[2], side, crossing);
  global[2] = z;
  
  // apply distortion corrections
  if(_dcc_static) { global = _distortionCorrection.get_corrected_position( global, _dcc_static ); }
  if(_dcc_average) { global = _distortionCorrection.get_corrected_position( global, _dcc_average ); }
  if(_dcc_fluctuation) { global = _distortionCorrection.get_corrected_position( global, _dcc_fluctuation ); }
}

void HelicalFitter::getTrackletClusters(TrackSeed *tracklet, std::vector<Acts::Vector3>& global_vec, std::vector<TrkrDefs::cluskey>& cluskey_vec)
{
  getTrackletClusterList(tracklet, cluskey_vec);
  // store cluster global positions in a vector
  TrackFitUtils::getTrackletClusters(_tGeometry, _cluster_map, global_vec, cluskey_vec);   
}

void HelicalFitter::getTrackletClusterList(TrackSeed *tracklet, std::vector<TrkrDefs::cluskey>& cluskey_vec)
{
  for (auto clusIter = tracklet->begin_cluster_keys();
       clusIter != tracklet->end_cluster_keys();
       ++clusIter)
    {
      auto key = *clusIter;
      auto cluster = _cluster_map->findCluster(key);
      if(!cluster)
	{
	  std::cout << "Failed to get cluster with key " << key << std::endl;
	  continue;
	}	  
      
      /// Make a safety check for clusters that couldn't be attached to a surface
      auto surf = _tGeometry->maps().getSurface(key, cluster);
      if(!surf)  { continue; }

      // drop some bad layers in the TPC completely
      unsigned int layer = TrkrDefs::getLayer(key);
      if(layer == 7 || layer == 22 || layer == 23 || layer == 38 || layer == 39) {continue;}

      cluskey_vec.push_back(key);
      
    } // end loop over clusters for this track 
}

std::vector<float> HelicalFitter::fitClusters(std::vector<Acts::Vector3>& global_vec, std::vector<TrkrDefs::cluskey> cluskey_vec)
{
  return TrackFitUtils::fitClusters(global_vec, cluskey_vec);       // do helical fit
}

Acts::Vector2 HelicalFitter::getClusterError(TrkrCluster *cluster, TrkrDefs::cluskey cluskey, Acts::Vector3& global)
{
  Acts::Vector2 clus_sigma(0,0);

  double clusRadius = sqrt(global[0]*global[0] + global[1]*global[1]);
  auto para_errors = _ClusErrPara.get_clusterv5_modified_error(cluster,clusRadius,cluskey);
  double phierror = sqrt(para_errors.first);
  double zerror = sqrt(para_errors.second);
  clus_sigma(1) = zerror;
  clus_sigma(0) = phierror;
  
  return clus_sigma; 
}

// new one
void HelicalFitter::getLocalDerivativesXY(Surface surf, Acts::Vector3 global, const std::vector<float>& fitpars, float lcl_derivativeX[5], float lcl_derivativeY[5], unsigned int layer)
{
  // Calculate the derivatives of the residual wrt the track parameters numerically
  std::vector<float> temp_fitpars;

  std::vector<float> fitpars_delta;
  fitpars_delta.push_back(0.1);   // radius, cm
  fitpars_delta.push_back(0.1);   // X0, cm
  fitpars_delta.push_back(0.1);   // Y0, cm
  fitpars_delta.push_back(0.1);   // zslope, cm
  fitpars_delta.push_back(0.1);   // Z0, cm

  for(unsigned int ip = 0; ip < fitpars.size(); ++ip)
    {
      temp_fitpars.push_back(fitpars[ip]);
    }

  // calculate projX and projY vectors once for the optimum fit parameters
  if(Verbosity() > 1) std::cout << "Call get_helix_tangent for best fit fitpars" << std::endl;
  std::pair<Acts::Vector3, Acts::Vector3> tangent = get_helix_tangent(fitpars, global);

  Acts::Vector3 projX(0,0,0), projY(0,0,0);
  get_projectionXY(surf, tangent, projX, projY);

  Acts::Vector3 intersection = get_helix_surface_intersection(surf, temp_fitpars, global);

  // loop over the track fit parameters
  for(unsigned int ip = 0; ip < fitpars.size(); ++ip)
    {
      Acts::Vector3 intersection_delta[2];
      for(int ipm = 0; ipm < 2; ++ipm)
	{
	  temp_fitpars[ip]  = fitpars[ip];  // reset to best fit value
	  float deltapm     = pow(-1.0, ipm);
	  temp_fitpars[ip] += deltapm * fitpars_delta[ip];

	  Acts::Vector3 temp_intersection = get_helix_surface_intersection(surf, temp_fitpars, global);
	  intersection_delta[ipm]         = temp_intersection - intersection;
	}
      Acts::Vector3 average_intersection_delta = (intersection_delta[0] - intersection_delta[1]) / (2 * fitpars_delta[ip]);

      if(Verbosity() > 1)
	{std::cout << " average_intersection_delta / delta " << average_intersection_delta(0) << "  " << average_intersection_delta(1) << "  " << average_intersection_delta(2) << std::endl;}

      // calculate the change in fit for X and Y 
      // - note negative sign from ATLAS paper is dropped here because mille wants the derivative of the fit, not the derivative of the residual
      lcl_derivativeX[ip] = average_intersection_delta.dot(projX);
      lcl_derivativeY[ip] = average_intersection_delta.dot(projY);
      if(Verbosity() > 1)
	{std::cout << " layer " << layer << " ip " << ip << "  derivativeX " << lcl_derivativeX[ip] << "  " 
		   << " derivativeY " << lcl_derivativeY[ip] << std::endl;}

      temp_fitpars[ip] = fitpars[ip];
    }
}


void HelicalFitter::getLocalVtxDerivativesXY(SvtxTrack& track, Acts::Vector3 event_vtx, const std::vector<float>& fitpars, float lcl_derivativeX[5], float lcl_derivativeY[5])
{
  // Calculate the derivatives of the residual wrt the track parameters numerically
  std::vector<float> temp_fitpars;
  Acts::Vector3 track_vtx (track.get_x(),track.get_y(),track.get_z());

  std::vector<float> fitpars_delta;
  fitpars_delta.push_back(0.1);   // radius, cm
  fitpars_delta.push_back(0.1);   // X0, cm
  fitpars_delta.push_back(0.1);   // Y0, cm
  fitpars_delta.push_back(0.1);   // zslope, cm
  fitpars_delta.push_back(0.1);   // Z0, cm

  for(unsigned int ip = 0; ip < fitpars.size(); ++ip)
    {
      temp_fitpars.push_back(fitpars[ip]);
    }

  // calculate projX and projY vectors once for the optimum fit parameters
  if(Verbosity() > 1) {std::cout << "Call get_helix_tangent for best fit fitpars" << std::endl;}

  Acts::Vector3 perigeeNormal (track.get_px(),track.get_py(),track.get_pz()) ;  

  // loop over the track fit parameters
  for(unsigned int ip = 0; ip < fitpars.size(); ++ip)
    {
      Acts::Vector3 localPerturb[2];
      Acts::Vector3 paperPerturb[2]; // for local derivative calculation like from the paper

      for(int ipm = 0; ipm < 2; ++ipm)
	{
	  temp_fitpars[ip]  = fitpars[ip];  // reset to best fit value
	  float deltapm     = pow(-1.0, ipm);
	  temp_fitpars[ip] += deltapm * fitpars_delta[ip];

	  Acts::Vector3 temp_track_vtx       = get_helix_vtx(event_vtx, temp_fitpars); // temporary pca
	  paperPerturb[ipm]                  = temp_track_vtx; // for og local derivative calculation 

	  // old method is next two lines
	  Acts::Vector3 localtemp_track_vtx  = globalvtxToLocalvtx(track, event_vtx, temp_track_vtx);
	  localPerturb[ipm]                  =  localtemp_track_vtx  ;
	  
	  if(Verbosity() > 1)
	    {
	      std::cout << "vtx local parameter " << ip << " with ipm " << ipm << " deltapm " << deltapm << " :" << std::endl;
	      std::cout<<" fitpars "<< fitpars[ip]<<" temp_fitpars "<<temp_fitpars[ip]<<std::endl;
	      std::cout << " localtmp_track_vtx: "<< localtemp_track_vtx<<std::endl;
	    }
	}

      Acts::Vector3 projX(0,0,0), projY(0,0,0);
      get_projectionVtxXY(track, event_vtx, projX, projY);

      Acts::Vector3 average_vtxX = (paperPerturb[0] - paperPerturb[1]) / (2 * fitpars_delta[ip]);
      Acts::Vector3 average_vtxY = (paperPerturb[0] - paperPerturb[1]) / (2 * fitpars_delta[ip]);

      // average_vtxX and average_vtxY are the numerical results in global coords for d(fit)/d(par)
      // The ATLAS paper formula is for the derivative of the residual, which is (measurement - fit = event vertex - track vertex)
      // Millepede wants the derivative of the fit, so we drop the minus sign from the paper
      lcl_derivativeX[ip] = average_vtxX.dot(projX); // 
      lcl_derivativeY[ip] = average_vtxY.dot(projY);

      temp_fitpars[ip] = fitpars[ip];
    }
}

void HelicalFitter::getGlobalDerivativesXY(Surface surf, Acts::Vector3 global, Acts::Vector3 fitpoint, const std::vector<float>& fitpars, float glbl_derivativeX[6], float glbl_derivativeY[6], unsigned int layer)
{
  Acts::Vector3 unitx(1, 0, 0);
  Acts::Vector3 unity(0, 1, 0);
  Acts::Vector3 unitz(0, 0, 1);

  // calculate projX and projY vectors once for the optimum fit parameters
  std::pair<Acts::Vector3, Acts::Vector3> tangent = get_helix_tangent(fitpars, global);  // should this be global, not fitpoint?

  Acts::Vector3 projX(0,0,0), projY(0,0,0);
  get_projectionXY(surf, tangent, projX, projY);
  
  // translations
  glbl_derivativeX[3] = unitx.dot(projX);
  glbl_derivativeX[4] = unity.dot(projX);
  glbl_derivativeX[5] = unitz.dot(projX);

  glbl_derivativeY[3] = unitx.dot(projY);
  glbl_derivativeY[4] = unity.dot(projY);
  glbl_derivativeY[5] = unitz.dot(projY);

  /*
  // note: the global derivative sign should be reversed from the ATLAS paper 
  // because mille wants the derivative of the fit, while the ATLAS paper gives the derivative of the residual.
  // But this sign reversal does NOT work. 
  // Verified that not reversing the sign here produces the correct sign of the prediction of the residual..
  for(unsigned int i = 3; i < 6; ++i)
  {
  glbl_derivativeX[i] *= -1.0; 
  glbl_derivativeY[i] *= -1.0; 
  }
  */
  // rotations
  // need center of sensor to intersection point
  Acts::Vector3 sensorCenter = surf->center(_tGeometry->geometry().getGeoContext()) / Acts::UnitConstants::cm;  // convert to cm
  Acts::Vector3 OM           = fitpoint - sensorCenter;   // this effectively reverses the sign from the ATLAS paper

  glbl_derivativeX[0] = (unitx.cross(OM)).dot(projX);
  glbl_derivativeX[1] = (unity.cross(OM)).dot(projX);
  glbl_derivativeX[2] = (unitz.cross(OM)).dot(projX);

  glbl_derivativeY[0] = (unitx.cross(OM)).dot(projY);
  glbl_derivativeY[1] = (unity.cross(OM)).dot(projY);
  glbl_derivativeY[2] = (unitz.cross(OM)).dot(projY);
  

  if(Verbosity() > 1)
    {
      for(int ip=0;ip<6;++ip)
	{
	  std::cout << " layer " << layer << " ip " << ip << "  glbl_derivativeX " << glbl_derivativeX[ip] << "  " 
		    << " glbl_derivativeY " << glbl_derivativeY[ip] << std::endl;
	}
    }
  
}


void HelicalFitter::getGlobalVtxDerivativesXY(SvtxTrack& track, Acts::Vector3 event_vtx, float glbl_derivativeX[3], float glbl_derivativeY[3])
{
  Acts::Vector3 unitx(1, 0, 0);
  Acts::Vector3 unity(0, 1, 0);
  Acts::Vector3 unitz(0, 0, 1);

  Acts::Vector3 track_vtx (track.get_x(),track.get_y(),track.get_z());
  Acts::Vector3 mom(track.get_px(),track.get_py(),track.get_pz());

  // calculate projX and projY vectors once for the optimum fit parameters
  Acts::Vector3 projX(0,0,0), projY(0,0,0);
  get_projectionVtxXY(track, event_vtx, projX, projY);

  // translations
  glbl_derivativeX[0] = unitx.dot(projX);
  glbl_derivativeX[1] = unity.dot(projX);
  glbl_derivativeX[2] = unitz.dot(projX);
  glbl_derivativeY[0] = unitx.dot(projY);
  glbl_derivativeY[1] = unity.dot(projY);
  glbl_derivativeY[2] = unitz.dot(projY);

  // The derivation in the ATLAS paper used above gives the derivative of the residual (= measurement - fit)
  // pede wants the derivative of the fit, so we reverse that - valid if our residual is (event vertex - track vertex)

  // Verified that reversing these signs produces the correct sign and magnitude for the prediction of the residual.
  // tested this by offsetting the simulated event vertex with zero misalignments. Pede fit reproduced simulated (xvtx, yvtx) within 7%.
  //   -- test gave zero for zvtx param, since this is determined relative to the measured event z vertex. 
  for(int i = 0; i<3;++i)
    {
      glbl_derivativeX[i] *= -1.0;
      glbl_derivativeY[i] *= -1.0;
    }

}

void HelicalFitter::get_projectionXY(Surface surf, std::pair<Acts::Vector3, Acts::Vector3> tangent, Acts::Vector3& projX, Acts::Vector3& projY)
{
  // we only need the direction part of the tangent
  Acts::Vector3 tanvec = tangent.second;
  // get the plane of the surface
  Acts::Vector3 sensorCenter = surf->center(_tGeometry->geometry().getGeoContext()) / Acts::UnitConstants::cm;  // convert to cm
  // sensorNormal is the Z vector
  Acts::Vector3 Z = -surf->normal(_tGeometry->geometry().getGeoContext()) / Acts::UnitConstants::cm; 
  // get surface X and Y unit vectors in global frame
  // transform Xlocal = 1.0 to global, subtract the surface center, normalize to 1
  Acts::Vector3 xloc(1.0,0.0,0.0); //local coord unit vector in x
  Acts::Vector3 xglob =  surf->transform(_tGeometry->geometry().getGeoContext()) * (xloc *  Acts::UnitConstants::cm);
  xglob /=  Acts::UnitConstants::cm;
  Acts::Vector3 yloc(0.0,1.0,0.0);
  Acts::Vector3 yglob =  surf->transform(_tGeometry->geometry().getGeoContext()) * (yloc *  Acts::UnitConstants::cm);
  yglob /=  Acts::UnitConstants::cm;
  Acts::Vector3 X = (xglob-sensorCenter) / (xglob-sensorCenter).norm();
  Acts::Vector3 Y = (yglob-sensorCenter) / (yglob-sensorCenter).norm();
  // see equation 31 of the ATLAS paper (and discussion) for this
  projX = X - (tanvec.dot(X) / tanvec.dot(Z)) * Z;
  projY = Y - (tanvec.dot(Y) / tanvec.dot(Z)) * Z;
  return;
}

void HelicalFitter::get_projectionVtxXY(SvtxTrack& track, Acts::Vector3 event_vtx, Acts::Vector3& projX, Acts::Vector3& projY)
{
  Acts::Vector3 tanvec(track.get_px(),track.get_py(),track.get_pz());
  Acts::Vector3 normal(track.get_px(),track.get_py(),0); 

  tanvec /= tanvec.norm();
  normal /= normal.norm();

  // get surface X and Y unit vectors in global frame
  Acts::Vector3 xloc(1.0,0.0,0.0);
  Acts::Vector3 yloc(0.0,0.0,1.0); // local y 
  Acts::Vector3 xglob = localvtxToGlobalvtx(track, event_vtx, xloc);
  Acts::Vector3 yglob = yloc + event_vtx;
  Acts::Vector3 X     = (xglob-event_vtx) / (xglob-event_vtx).norm(); // local unit vector transformed to global coordinates
  Acts::Vector3 Y     = (yglob-event_vtx) / (yglob-event_vtx).norm();

  // see equation 31 of the ATLAS paper (and discussion) for this
  projX = X - (tanvec.dot(X) / tanvec.dot(normal)) * normal;
  projY = Y - (tanvec.dot(Y) / tanvec.dot(normal)) * normal;

  return;
}

unsigned int HelicalFitter::addSiliconClusters(std::vector<float>& fitpars, std::vector<Acts::Vector3>& global_vec,  std::vector<TrkrDefs::cluskey>& cluskey_vec)
{

  return TrackFitUtils::addSiliconClusters(fitpars, dca_cut, _tGeometry, _cluster_map, global_vec, cluskey_vec);
}

bool HelicalFitter::is_intt_layer_fixed(unsigned int layer)
{
  bool ret = false;
  auto it = fixed_intt_layers.find(layer);
  if(it != fixed_intt_layers.end()) 
    ret = true;

  return ret;
}

bool HelicalFitter::is_mvtx_layer_fixed(unsigned int layer, unsigned int clamshell)
{
  bool ret = false;

  std::pair pair = std::make_pair(layer,clamshell);
  auto it = fixed_mvtx_layers.find(pair);
  if(it != fixed_mvtx_layers.end()) 
    ret = true;

  return ret;
}

void HelicalFitter::set_intt_layer_fixed(unsigned int layer)
{
  fixed_intt_layers.insert(layer);
}

void HelicalFitter::set_mvtx_layer_fixed(unsigned int layer, unsigned int clamshell)
{
  fixed_mvtx_layers.insert(std::make_pair(layer,clamshell));
}

 
bool HelicalFitter::is_layer_param_fixed(unsigned int layer, unsigned int param)
{
  bool ret = false;
  std::pair<unsigned int, unsigned int> pair = std::make_pair(layer, param);
  auto it = fixed_layer_params.find(pair);
  if(it != fixed_layer_params.end()) 
    ret = true;

  return ret;
}

void HelicalFitter::set_layer_param_fixed(unsigned int layer, unsigned int param)
{
  std::pair<unsigned int, unsigned int> pair = std::make_pair(layer, param);
  fixed_layer_params.insert(pair);
}

void HelicalFitter::set_tpc_sector_fixed(unsigned int region, unsigned int sector, unsigned int side)
{
  // make a combined subsector index
  unsigned int subsector = region * 24 + side * 12 + sector;
  fixed_sectors.insert(subsector);
}

bool HelicalFitter::is_tpc_sector_fixed(unsigned int layer, unsigned int sector, unsigned int side)
{
  bool ret = false;
  unsigned int region = AlignmentDefs::getTpcRegion(layer);
  unsigned int subsector = region * 24 + side * 12 + sector;
  auto it = fixed_sectors.find(subsector);
  if(it != fixed_sectors.end()) 
    ret = true;
  
  return ret;
}

void HelicalFitter::correctTpcGlobalPositions(std::vector<Acts::Vector3> global_vec,  std::vector<TrkrDefs::cluskey> cluskey_vec)
{
  for(unsigned int iclus=0;iclus<cluskey_vec.size();++iclus)
    {
      auto cluskey = cluskey_vec[iclus];
      auto global = global_vec[iclus];
      const unsigned int trkrId = TrkrDefs::getTrkrId(cluskey);	  
      if(trkrId == TrkrDefs::tpcId) 
	{  
	  // have to add corrections for TPC clusters after transformation to global
	  int crossing = 0;  // for now
	  makeTpcGlobalCorrections(cluskey, crossing, global); 
	  global_vec[iclus] = global;
	}
    }
}

float HelicalFitter::getVertexResidual(Acts::Vector3 vtx)
{
  float phi = atan2(vtx(1),vtx(0));
  float r   = vtx(0)/cos(phi);
  float test_r = sqrt(vtx(0)*vtx(0)+vtx(1)*vtx(1));

  if(Verbosity() > 1)
    {
      std::cout << "my method position: " << vtx << std::endl;
      std::cout << "r " << r <<" phi: " << phi*180/M_PI<<" test_r"<< test_r << std::endl;
    }
  return r;
}

void HelicalFitter::get_dca(SvtxTrack& track,float& dca3dxy, float& dca3dz, float& dca3dxysigma, float& dca3dzsigma, Acts::Vector3 event_vertex)
{
  //give trackseed 
  dca3dxy = NAN;
  Acts::Vector3 track_vtx(track.get_x(),track.get_y(),track.get_z());
  Acts::Vector3 mom(track.get_px(),track.get_py(),track.get_pz());

  track_vtx -= event_vertex; // difference between track_vertex and event_vtx
  
  Acts::ActsSquareMatrix<3> posCov;
  for(int i = 0; i < 3; ++i)
    {
      for(int j = 0; j < 3; ++j)
	{
	  posCov(i, j) = track.get_error(i, j);
	} 
    }
  
  Acts::Vector3 r = mom.cross(Acts::Vector3(0.,0.,1.));

  float phi       = atan2(r(1), r(0));
  Acts::RotationMatrix3 rot;
  Acts::RotationMatrix3 rot_T;
  phi *= -1;
  rot(0,0) = cos(phi);
  rot(0,1) = -sin(phi);
  rot(0,2) = 0;
  rot(1,0) = sin(phi);
  rot(1,1) = cos(phi);
  rot(1,2) = 0;
  rot(2,0) = 0;
  rot(2,1) = 0;
  rot(2,2) = 1;
  rot_T    = rot.transpose();

  Acts::Vector3 pos_R           = rot * track_vtx;
  Acts::ActsSquareMatrix<3> rotCov = rot * posCov * rot_T;
  dca3dxy      = pos_R(0);
  dca3dz       = pos_R(2);
  dca3dxysigma = sqrt(rotCov(0,0));
  dca3dzsigma  = sqrt(rotCov(2,2));
	
  if(Verbosity()>1)
    {
      std::cout << " momentum X z: "<<r<< " phi: " << phi*180/M_PI << std::endl;
      std::cout << "dca3dxy " << dca3dxy << " dca3dz: " << dca3dz << " pos_R(1): "<<pos_R(1)<<" dca3dxysigma " << dca3dxysigma << " dca3dzsigma " << dca3dzsigma << std::endl;
    }
  
}

Acts::Vector3 HelicalFitter::globalvtxToLocalvtx(SvtxTrack& track, Acts::Vector3 event_vertex)
{

  Acts::Vector3 track_vtx(track.get_x(),track.get_y(),track.get_z());
  Acts::Vector3 mom(track.get_px(),track.get_py(),track.get_pz());
  track_vtx -= event_vertex; // difference between track_vertex and event_vtx

  Acts::Vector3 r = mom.cross(Acts::Vector3(0.,0.,1.));
  float phi       = atan2(r(1), r(0));
  Acts::RotationMatrix3 rot;
  Acts::RotationMatrix3 rot_T;
  phi     *= -1;
  rot(0,0) = cos(phi);
  rot(0,1) = -sin(phi);
  rot(0,2) = 0;
  rot(1,0) = sin(phi);
  rot(1,1) = cos(phi);
  rot(1,2) = 0;
  rot(2,0) = 0;
  rot(2,1) = 0;
  rot(2,2) = 1;
  rot_T    = rot.transpose();

  Acts::Vector3 pos_R = rot * track_vtx;

  if(Verbosity()>1)
    {
      std::cout << " momentum X z: "<<r<< " phi: " << phi*180/M_PI << std::endl;
      std::cout << " pos_R(0): "<<pos_R(0)<<" pos_R(1): "<<pos_R(1) << std::endl;
    }
  return pos_R; 
}

Acts::Vector3 HelicalFitter::globalvtxToLocalvtx(SvtxTrack& track, Acts::Vector3 event_vertex, Acts::Vector3 PCA)
{
  Acts::Vector3 mom(track.get_px(),track.get_py(),track.get_pz());
  PCA -= event_vertex; // difference between track_vertex and event_vtx

  Acts::Vector3 r = mom.cross(Acts::Vector3(0.,0.,1.));
  float phi       = atan2(r(1), r(0));
  Acts::RotationMatrix3 rot;
  Acts::RotationMatrix3 rot_T;
  phi     *= -1;
  rot(0,0) = cos(phi);
  rot(0,1) = -sin(phi);
  rot(0,2) = 0;
  rot(1,0) = sin(phi);
  rot(1,1) = cos(phi);
  rot(1,2) = 0;
  rot(2,0) = 0;
  rot(2,1) = 0;
  rot(2,2) = 1;
  rot_T    = rot.transpose();

  Acts::Vector3 pos_R = rot * PCA;

  if(Verbosity()>1)
    {
      std::cout << " momentum X z: "<<r<< " phi: " << phi*180/M_PI << std::endl;
      std::cout << " pos_R(0): "<<pos_R(0)<<" pos_R(1): "<<pos_R(1) << std::endl;
    }
  return pos_R; 
}


Acts::Vector3 HelicalFitter::localvtxToGlobalvtx(SvtxTrack& track, Acts::Vector3 event_vtx,  Acts::Vector3 local)
{

  //Acts::Vector3 track_vtx = local;
  Acts::Vector3 mom(track.get_px(),track.get_py(),track.get_pz());
  //std::cout << "first pos: " << pos << " mom: " << mom << std::endl; 
  //local -= event_vertex; // difference between track_vertex and event_vtx

  Acts::Vector3 r = mom.cross(Acts::Vector3(0.,0.,1.));
  float phi       = atan2(r(1), r(0));
  Acts::RotationMatrix3 rot;
  Acts::RotationMatrix3 rot_T;
  // phi *= -1;
  rot(0,0) = cos(phi);
  rot(0,1) = -sin(phi);
  rot(0,2) = 0;
  rot(1,0) = sin(phi);
  rot(1,1) = cos(phi);
  rot(1,2) = 0;
  rot(2,0) = 0;
  rot(2,1) = 0;
  rot(2,2) = 1;
  
  rot_T = rot.transpose();

  Acts::Vector3 pos_R = rot * local;
  pos_R += event_vtx;
  if(Verbosity()>1)
    {
      std::cout << " momentum X z: "<<r<< " phi: " << phi*180/M_PI << std::endl;
      std::cout << " pos_R(0): "<<pos_R(0)<<" pos_R(1): "<<pos_R(1)<<"pos_R(2): "<<pos_R(2) << std::endl;
    }
  return pos_R; 
}

