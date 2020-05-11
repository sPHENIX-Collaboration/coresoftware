/*!
 *  \file		PHG4TrackFastSim.cc
 *  \brief		Kalman Filter based on smeared truth PHG4Hit
 *  \details	Kalman Filter based on smeared truth PHG4Hit
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#include "PHG4TrackFastSim.h"

#include <phgenfit/Fitter.h>
#include <phgenfit/Measurement.h>  // for Measurement
#include <phgenfit/PlanarMeasurement.h>
#include <phgenfit/SpacepointMeasurement.h>
#include <phgenfit/Track.h>
#include <phgeom/PHGeomUtility.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxTrackState.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/SvtxTrack_FastSim.h>
#include <trackbase_historic/SvtxVertex.h>     // for SvtxVertex
#include <trackbase_historic/SvtxVertexMap.h>  // for SvtxVertexMap
#include <trackbase_historic/SvtxVertexMap_v1.h>
#include <trackbase_historic/SvtxVertex_v1.h>

#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <phfield/PHFieldUtility.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>  // for PHG4HitContainer
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <GenFit/AbsMeasurement.h>
#include <GenFit/EventDisplay.h>
#include <GenFit/MeasuredStateOnPlane.h>
#include <GenFit/RKTrackRep.h>

#include <GenFit/FitStatus.h>              // for FitStatus
#include <GenFit/GFRaveTrackParameters.h>  // for GFRaveTrackParameters
#include <GenFit/GFRaveVertex.h>
#include <GenFit/GFRaveVertexFactory.h>
#include <GenFit/Track.h>

#include <TMath.h>
#include <TMatrixDSymfwd.h>  // for TMatrixDSym
#include <TMatrixTSym.h>     // for TMatrixTSym
#include <TMatrixTUtils.h>   // for TMatrixTRow
#include <TVector3.h>        // for TVector3, operator*
#include <TVectorDfwd.h>     // for TVectorD
#include <TVectorT.h>        // for TVectorT

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <cassert>  // for assert
#include <cmath>
#include <iostream>  // for operator<<, basic_...
#include <map>
#include <memory>  // for unique_ptr, alloca...
#include <utility>

class PHField;
class TGeoManager;
namespace genfit
{
class AbsTrackRep;
}  // namespace genfit

#define LogDebug(exp) \
  if (Verbosity()) std::cout << "PHG4TrackFastSim (DEBUG): " << __FILE__ << ": " << __LINE__ << ": " << exp << "\n"
#define LogError(exp) \
  std::cout << "PHG4TrackFastSim (ERROR): " << __FILE__ << ": " << __LINE__ << ": " << exp << "\n"
#define LogWarning(exp) \
  std::cout << "PHG4TrackFastSim (WARNING): " << __FILE__ << ": " << __LINE__ << ": " << exp << "\n"

using namespace std;

PHG4TrackFastSim::PHG4TrackFastSim(const std::string& name)
  : SubsysReco(name)
  , _event(-1)
  , m_SmearingFlag(true)
  , _truth_container(nullptr)
  , _sub_top_node_name("SVTX")
  , _trackmap_out_name("SvtxTrackMap")
  , _trackmap_out(nullptr)
  , _vertexmap(nullptr)
  , _fitter(nullptr)
  , _vertex_finder(nullptr)
  , _vertexing_method("kalman-smoothing:1")
  , _vertex_min_ndf(10)
  , _do_vertexing(false)
  , _fit_alg_name("DafRef")  // was ("KalmanFitterRefTrack")
  , _primary_assumption_pid(211)
  , m_DoEvtDisplayFlag(false)
  , _use_vertex_in_fitting(true)
  , _vertex_xy_resolution(50E-4)
  , _vertex_z_resolution(50E-4)
  , _primary_tracking(1)
{
  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  m_RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(m_RandomGenerator, seed);
}

PHG4TrackFastSim::~PHG4TrackFastSim()
{
  delete _fitter;
  delete _vertex_finder;
  gsl_rng_free(m_RandomGenerator);
}

int PHG4TrackFastSim::InitRun(PHCompositeNode* topNode)
{
  _event = -1;

  int ret = CreateNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }
  ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  TGeoManager* tgeo_manager = PHGeomUtility::GetTGeoManager(topNode);
  PHField* field = PHFieldUtility::GetFieldMapNode(nullptr, topNode);

  _fitter = PHGenFit::Fitter::getInstance(tgeo_manager,
                                          field, _fit_alg_name, "RKTrackRep",
                                          m_DoEvtDisplayFlag);

  if (!_fitter)
  {
    cerr << PHWHERE << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  _fitter->set_verbosity(Verbosity());

  // tower geometry for track states

  for (map<string,double>::iterator iter = m_ZStateMap.begin(); iter != m_ZStateMap.end(); ++iter)
  {
// if the plane is given do not look for the geometry object
    if (isfinite(iter->second))
    {
      continue;
    }
// this is for the legacy finding of the Babar forward calorimeters which have their geometry on the node tree
      string towergeonodename = "TOWERGEOM_" + iter->first;
      RawTowerGeomContainer* towergeo = findNode::getClass<RawTowerGeomContainer>(topNode, towergeonodename);
      if (!towergeo)
      {
        cerr << PHWHERE << " ERROR: Can't find node " << towergeonodename << endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }

      // Grab the first tower, us it to get the
      // location along the beamline
      RawTowerGeomContainer::ConstRange twr_range = towergeo->get_tower_geometries();
      RawTowerGeomContainer::ConstIterator twr_iter = twr_range.first;
      RawTowerGeom* temp_geo = twr_iter->second;

      //Changed by Barak on 12/10/19
      _state_location.push_back(temp_geo->get_center_z());
      iter->second=temp_geo->get_center_z();
    }

  for (map<string,double>::iterator iter = m_CylinderStateMap.begin(); iter != m_CylinderStateMap.end(); ++iter)
  {
    if (isfinite(iter->second))
    {
      continue;
    }
// this is for the legacy finding of the Babar sPHENIX calorimeters which have their geometry on the node tree
      string nodename = "TOWERGEOM_" + iter->first;
      RawTowerGeomContainer* geo = findNode::getClass<RawTowerGeomContainer>(topNode, nodename);
      if (geo)
      {
        _state_location.push_back(geo->get_radius());
	iter->second=geo->get_radius();
      }
  }

/*
  for (unsigned int i = 0; i < _state_names.size(); i++)
  {

    if ((_state_names[i] == "FHCAL") || (_state_names[i] == "FEMC") || (_state_names[i] == "EEMC"))
    {
      // Get the z-location of the detector plane

      string towergeonodename = "TOWERGEOM_" + _state_names[i];
      RawTowerGeomContainer* towergeo = findNode::getClass<RawTowerGeomContainer>(topNode, towergeonodename.c_str());
      if (!towergeo)
      {
        cerr << PHWHERE << " ERROR: Can't find node " << towergeonodename << endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }

      // Grab the first tower, us it to get the
      // location along the beamline
      RawTowerGeomContainer::ConstRange twr_range = towergeo->get_tower_geometries();
      RawTowerGeomContainer::ConstIterator twr_iter = twr_range.first;
      RawTowerGeom* temp_geo = twr_iter->second;

      //Changed by Barak on 12/10/19
      _state_location.push_back(temp_geo->get_center_z());
    }
    else if ((_state_names[i] == "CEMC") || (_state_names[i] == "HCALIN") || (_state_names[i] == "OHCAL"))
    {
      // Get the calorimeter radius

      string nodename = "TOWERGEOM_" + _state_names[i];
      RawTowerGeomContainer* geo = findNode::getClass<RawTowerGeomContainer>(topNode, nodename.c_str());
      if (geo)
      {
        _state_location.push_back(geo->get_radius());
      }
      else
      {
        cerr << PHWHERE << " ERROR: Can't find node " << nodename << endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }
    }
    else
    {
      cerr << PHWHERE << " ERROR: Unrecognized detector name for state projection:  " << _state_names[i] << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
*/
  if (_do_vertexing)
  {
    _vertex_finder = new genfit::GFRaveVertexFactory(Verbosity(), true);
    //_vertex_finder->setMethod("kalman-smoothing:1"); //! kalman-smoothing:1 is the defaul method
    _vertex_finder->setMethod(_vertexing_method.data());
    //_vertex_finder->setBeamspot();

    //_vertex_finder = new PHRaveVertexFactory(Verbosity());

    if (!_vertex_finder)
    {
      cout << PHWHERE << " no Vertex Finder" << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TrackFastSim::End(PHCompositeNode* topNode)
{
  if (m_DoEvtDisplayFlag && _fitter)
  {
    _fitter->displayEvent();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TrackFastSim::process_event(PHCompositeNode* topNode)
{
  _event++;

  if (Verbosity() >= 2)
    std::cout << "PHG4TrackFastSim::process_event: " << _event << ".\n";

  //	if(_clustermap_out)
  //		_clustermap_out->empty();
  //	else {
  //		LogError("_clustermap_out not found!");
  //		return Fun4AllReturnCodes::ABORTRUN;
  //	}

  if (_trackmap_out)
  {
    _trackmap_out->empty();
  }
  else
  {
    LogError("_trackmap_out not found!");
    return Fun4AllReturnCodes::ABORTRUN;
  }

  vector<PHGenFit::Track*> rf_tracks;

  PHG4VtxPoint* truthVtx = _truth_container->GetPrimaryVtx(_truth_container->GetPrimaryVertexIndex());
  TVector3 vtxPoint(truthVtx->get_x(), truthVtx->get_y(), truthVtx->get_z());
  // Smear the vertex ONCE for all particles in the event
  vtxPoint.SetX(vtxPoint.x() + gsl_ran_gaussian(m_RandomGenerator, _vertex_xy_resolution));
  vtxPoint.SetY(vtxPoint.y() + gsl_ran_gaussian(m_RandomGenerator, _vertex_xy_resolution));
  vtxPoint.SetZ(vtxPoint.z() + gsl_ran_gaussian(m_RandomGenerator, _vertex_z_resolution));

  PHG4TruthInfoContainer::ConstRange itr_range;
  if (_primary_tracking)
  {
    // Tracking for primaries only
    itr_range = _truth_container->GetPrimaryParticleRange();
  }
  else
  {
    // Check ALL particles
    itr_range = _truth_container->GetParticleRange();
  }

  GenFitTrackMap gf_track_map;
  // Now we can loop over the particles

  for (PHG4TruthInfoContainer::ConstIterator itr = itr_range.first;
       itr != itr_range.second; ++itr)
  {
    PHG4Particle* particle = itr->second;

    TVector3 seed_pos(vtxPoint.x(), vtxPoint.y(), vtxPoint.z());
    TVector3 seed_mom(0, 0, 0);
    TMatrixDSym seed_cov(6);

    //! Create measurements
    std::vector<PHGenFit::Measurement*> measurements;

    PHGenFit::Measurement* vtx_meas = nullptr;

    if (_use_vertex_in_fitting)
    {
      vtx_meas = VertexMeasurement(TVector3(vtxPoint.x(),
                                            vtxPoint.y(),
                                            vtxPoint.z()),
                                   _vertex_xy_resolution,
                                   _vertex_z_resolution);
      measurements.push_back(vtx_meas);
    }

    PseudoPatternRecognition(particle, measurements, seed_pos, seed_mom,
                             seed_cov);

    if (measurements.size() < 3)
    {
      if (Verbosity() >= 2)
      {
        //LogWarning("measurements.size() < 3");
        std::cout << "event: " << _event << " : measurements.size() < 3"
                  << "\n";
      }
      // Delete the measurements
      // We need to also delete the underlying genfit::AbsMeasurement object
      for (unsigned int im = 0; im < measurements.size(); im++)
      {
        delete measurements[im]->getMeasurement();
        delete measurements[im];
      }
      continue;
    }

    //! Build TrackRep from particle assumption
    /*!
	   * mu+:	-13
	   * mu-:	13
	   * pi+:	211
	   * pi-:	-211
	   * e-:	11
	   * e+:	-11
	   */
    //int pid = 13; //
    //SMART(genfit::AbsTrackRep) rep = NEW(genfit::RKTrackRep)(pid);
    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(_primary_assumption_pid);

    //rep->setDebugLvl(1); //DEBUG

    //! Initialize track with seed from pattern recognition

    PHGenFit::Track* track = new PHGenFit::Track(rep, seed_pos, seed_mom, seed_cov);

    // NOTE: We need to keep a list of tracks so they can
    // all be cleanly deleted at the end
    rf_tracks.push_back(track);

    //LogDEBUG;
    //! Add measurements to track
    track->addMeasurements(measurements);

    //LogDEBUG;
    //! Fit the track
    int fitting_err = _fitter->processTrack(track, false);

    if (fitting_err != 0)
    {
      if (Verbosity() >= 2)
      {
        //LogWarning("measurements.size() < 3");
        std::cout << "event: " << _event
                  << " : fitting_err != 0, next track."
                  << "\n";
      }
      continue;
    }

    TVector3 vtx(vtxPoint.x(), vtxPoint.y(), vtxPoint.z());
    SvtxTrack* svtx_track_out = MakeSvtxTrack(track,
                                              particle->get_track_id(),
                                              measurements.size(), vtx);
    if (Verbosity() > 1)
    {
      svtx_track_out->identify();
    }

    if (svtx_track_out)
    {
      //      track -> output container

      const unsigned int track_id = _trackmap_out->insert(svtx_track_out)->get_id();
      gf_track_map.insert({track->getGenFitTrack(), track_id});

      delete svtx_track_out;  // insert makes a clone
    }

  }  // Loop all primary particles

  //vertex finding
  if (_do_vertexing)
  {
    if (!_vertex_finder)
    {
      cout << __PRETTY_FUNCTION__ << "Failed to init vertex finder" << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    if (!_vertexmap)
    {
      cout << __PRETTY_FUNCTION__ << "Failed to init vertex map" << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    //    genfit::GFRaveVertexFactory* _vertex_finder = new genfit::GFRaveVertexFactory(10, true);
    //    _vertex_finder->setMethod("kalman-smoothing:1");
    //    _vertex_finder->setBeamspot();

    vector<genfit::GFRaveVertex*> rave_vertices;
    if (rf_tracks.size() >= 2)
    {
      try
      {
        vector<genfit::Track*> rf_gf_tracks;
        for (std::vector<PHGenFit::Track*>::iterator it = rf_tracks.begin(); it != rf_tracks.end(); ++it)
        {
          genfit::Track* track = (*it)->getGenFitTrack();

          if (Verbosity())
          {
            TVector3 pos, mom;
            TMatrixDSym cov;

            track->getFittedState().getPosMomCov(pos, mom, cov);

            cout << "Track getCharge = " << track->getFitStatus()->getCharge() << " getChi2 = " << track->getFitStatus()->getChi2() << " getNdf = " << track->getFitStatus()->getNdf() << endl;
            pos.Print();
            mom.Print();
            cov.Print();
          }
          if (track->getFitStatus()->getNdf() > _vertex_min_ndf)
            rf_gf_tracks.push_back(track);
        }
        _vertex_finder->findVertices(&rave_vertices, rf_gf_tracks);
      }
      catch (...)
      {
        if (Verbosity() > 1)
          std::cout << PHWHERE << "GFRaveVertexFactory::findVertices failed!";
      }
    }

    if (Verbosity())
    {
      cout << __PRETTY_FUNCTION__ << __LINE__ << " rf_tracks = " << rf_tracks.size() << " rave_vertices = " << rave_vertices.size() << endl;
    }
    FillSvtxVertexMap(rave_vertices, gf_track_map);
  }

  //! add tracks to event display
  if (m_DoEvtDisplayFlag)
  {
    vector<genfit::Track*> rf_gf_tracks;
    for (std::vector<PHGenFit::Track*>::iterator it = rf_tracks.begin(); it != rf_tracks.end(); ++it)
    {
      rf_gf_tracks.push_back((*it)->getGenFitTrack());
    }
    _fitter->getEventDisplay()->addEvent(rf_gf_tracks);
  }
  else
  {
    for (std::vector<PHGenFit::Track*>::iterator it = rf_tracks.begin(); it != rf_tracks.end(); ++it)
    {
      delete (*it);
    }
    rf_tracks.clear();
  }

  //	if(_trackmap_out->get(0)) {
  //		_trackmap_out->get(0)->identify();
  //		std::cout<<"DEBUG : "<< _trackmap_out->get(0)->get_px() <<"\n";
  //		std::cout<<"DEBUG : "<< _trackmap_out->get(0)->get_truth_track_id() <<"\n";
  //	}

  return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * Fill SvtxVertexMap from GFRaveVertexes and Tracks
 */
bool PHG4TrackFastSim::FillSvtxVertexMap(
    const std::vector<genfit::GFRaveVertex*>& rave_vertices,
    const GenFitTrackMap& gf_track_map)
{
  for (genfit::GFRaveVertex* rave_vtx : rave_vertices)
  {
    if (!rave_vtx)
    {
      cerr << PHWHERE << endl;
      return false;
    }

    std::shared_ptr<SvtxVertex> svtx_vtx(new SvtxVertex_v1());

    svtx_vtx->set_chisq(rave_vtx->getChi2());
    svtx_vtx->set_ndof(rave_vtx->getNdf());
    svtx_vtx->set_position(0, rave_vtx->getPos().X());
    svtx_vtx->set_position(1, rave_vtx->getPos().Y());
    svtx_vtx->set_position(2, rave_vtx->getPos().Z());

    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        svtx_vtx->set_error(i, j, rave_vtx->getCov()[i][j]);

    for (unsigned int i = 0; i < rave_vtx->getNTracks(); i++)
    {
      //TODO improve speed
      const genfit::Track* rave_track =
          rave_vtx->getParameters(i)->getTrack();
      //      for(auto iter : gf_track_map) {
      //        if (iter.second == rave_track)
      //          svtx_vtx->insert_track(iter.first);
      //      }
      auto iter = gf_track_map.find(rave_track);
      if (iter != gf_track_map.end())
      {
        svtx_vtx->insert_track(iter->second);
      }
    }

    if (_vertexmap)
    {
      _vertexmap->insert_clone(svtx_vtx.get());
    }
    else
    {
      LogError("!_vertexmap");
    }

  }  //loop over RAVE vertices

  return true;
}

int PHG4TrackFastSim::CreateNodes(PHCompositeNode* topNode)
{
  // create nodes...
  PHNodeIterator iter(topNode);

  PHCompositeNode* dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cout << PHWHERE << " DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  PHNodeIterator iter_dst(dstNode);

  // Create the FGEM node
  PHCompositeNode* tb_node = dynamic_cast<PHCompositeNode*>(iter_dst.findFirst(
      "PHCompositeNode", _sub_top_node_name.c_str()));
  if (!tb_node)
  {
    tb_node = new PHCompositeNode(_sub_top_node_name.c_str());
    dstNode->addNode(tb_node);
    if (Verbosity() > 0)
    {
      cout << _sub_top_node_name.c_str() << " node added" << endl;
    }
  }

  //	_clustermap_out = new SvtxClusterMap_v1;
  //
  //	PHIODataNode<PHObject>* clusters_node = new PHIODataNode<PHObject>(
  //			_clustermap_out, _clustermap_out_name.c_str(), "PHObject");
  //	tb_node->addNode(clusters_node);
  //	if (Verbosity() > 0)
  //		cout << _clustermap_out_name.c_str() <<" node added" << endl;

  _trackmap_out = new SvtxTrackMap_v1;

  PHIODataNode<PHObject>* tracks_node = new PHIODataNode<PHObject>(_trackmap_out, _trackmap_out_name.c_str(), "PHObject");
  tb_node->addNode(tracks_node);
  if (Verbosity() > 0)
  {
    cout << _trackmap_out_name.c_str() << " node added" << endl;
  }

  if (_do_vertexing)
  {
    _vertexmap = new SvtxVertexMap_v1;
    PHIODataNode<PHObject>* vertexes_node = new PHIODataNode<PHObject>(_vertexmap, "SvtxVertexMap", "PHObject");
    tb_node->addNode(vertexes_node);
    if (Verbosity() > 0)
    {
      cout << "Svtx/SvtxVertexMap node added" << endl;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TrackFastSim::GetNodes(PHCompositeNode* topNode)
{
  //DST objects
  //Truth container
  _truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_truth_container)
  {
    cout << PHWHERE << " PHG4TruthInfoContainer node not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  for (unsigned int i = 0; i < _phg4hits_names.size(); i++)
  {
    PHG4HitContainer* phg4hit = findNode::getClass<PHG4HitContainer>(topNode, _phg4hits_names[i]);
    if (!phg4hit)
    {
      cout << PHWHERE << _phg4hits_names[i].c_str()
           << " node not found on node tree" << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    if (Verbosity() > 0)
    {
      cout << "PHG4TrackFastSim::GetNodes - node added: " << _phg4hits_names[i].c_str() << endl;
    }
    _phg4hits.push_back(phg4hit);
  }

  //checks
  assert(_phg4hits_names.size() == _phg4hits.size());
  assert(_phg4_detector_type.size() == _phg4hits.size());
  assert(_phg4_detector_radres.size() == _phg4hits.size());
  assert(_phg4_detector_phires.size() == _phg4hits.size());
  assert(_phg4_detector_lonres.size() == _phg4hits.size());
  assert(_phg4_detector_hitfindeff.size() == _phg4hits.size());
  assert(_phg4_detector_noise.size() == _phg4hits.size());

  //	_clustermap_out = findNode::getClass<SvtxClusterMap>(topNode,
  //			_clustermap_out_name.c_str());
  //	if (!_clustermap_out && _event < 2) {
  //		cout << PHWHERE << _clustermap_out_name.c_str() << " node not found on node tree"
  //				<< endl;
  //		return Fun4AllReturnCodes::ABORTEVENT;
  //	}

  _trackmap_out = findNode::getClass<SvtxTrackMap>(topNode, _trackmap_out_name.c_str());
  if (!_trackmap_out && _event < 2)
  {
    cout << PHWHERE << _trackmap_out_name.c_str()
         << " node not found on node tree" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TrackFastSim::PseudoPatternRecognition(const PHG4Particle* particle,
                                               std::vector<PHGenFit::Measurement*>& meas_out, TVector3& seed_pos,
                                               TVector3& seed_mom, TMatrixDSym& seed_cov, const bool do_smearing)
{
  seed_cov.ResizeTo(6, 6);

  seed_pos.SetXYZ(0, 0, 0);
  // reset the seed resolution to the approximate position resolution of the last detector
  seed_cov[0][0] = .1 * .1;
  seed_cov[1][1] = .1 * .1;
  seed_cov[2][2] = 30 * 30;
  //  for (int i = 0; i < 3; i++)
  //  {
  //    seed_cov[i][i] = _phi_resolution * _phi_resolution;
  //  }

  seed_mom.SetXYZ(0, 0, 10);
  for (int i = 3; i < 6; i++)
  {
    seed_cov[i][i] = 10;
  }

  if (particle)
  {
    TVector3 True_mom(particle->get_px(), particle->get_py(),
                      particle->get_pz());

    seed_mom.SetXYZ(particle->get_px(), particle->get_py(),
                    particle->get_pz());
    if (do_smearing)
    {
      const double momSmear = 3. / 180. * TMath::Pi();  // rad
      const double momMagSmear = 0.1;                   // relative

      seed_mom.SetMag(
          True_mom.Mag() + gsl_ran_gaussian(m_RandomGenerator,
                                            momMagSmear * True_mom.Mag()));
      seed_mom.SetTheta(True_mom.Theta() + gsl_ran_gaussian(m_RandomGenerator, momSmear));
      seed_mom.SetPhi(True_mom.Phi() + gsl_ran_gaussian(m_RandomGenerator, momSmear));
    }
  }

  if (Verbosity())
  {
    std::cout << "PHG4TrackFastSim::PseudoPatternRecognition - DEBUG: "
              << "searching for hits from  " << _phg4hits.size() << " PHG4Hit nodes" << endl;
  }

  // order measurement with g4hit time via stl multimap
  multimap<double, PHGenFit::Measurement*> ordered_measurements;

  for (unsigned int ilayer = 0; ilayer < _phg4hits.size(); ilayer++)
  {
    if (!_phg4hits[ilayer])
    {
      LogError("No _phg4hits[i] found!");
      continue;
    }

    int dettype = _phg4_detector_type[ilayer];
    float detradres = _phg4_detector_radres[ilayer];
    float detphires = _phg4_detector_phires[ilayer];
    float detlonres = _phg4_detector_lonres[ilayer];
    float dethiteff = _phg4_detector_hitfindeff[ilayer];
    float detnoise = _phg4_detector_noise[ilayer];
    if (Verbosity())
    {
      std::cout << "PHG4TrackFastSim::PseudoPatternRecognition - DEBUG: "
                << "ilayer: "
                << ilayer << ",  " << _phg4hits_names[ilayer]
                << " with nsublayers: " << _phg4hits[ilayer]->num_layers()
                << ", detradres = " << detradres
                << ", detphires = " << detphires
                << ", detlonres = " << detlonres
                << ", dethiteff = " << dethiteff
                << ", detnoise = " << detnoise
                << " \n";
    }
    for (PHG4HitContainer::LayerIter layerit =
             _phg4hits[ilayer]->getLayers().first;
         layerit != _phg4hits[ilayer]->getLayers().second; layerit++)
    {
      for (PHG4HitContainer::ConstIterator itr =
               _phg4hits[ilayer]->getHits(*layerit).first;
           itr != _phg4hits[ilayer]->getHits(*layerit).second; ++itr)
      {
        PHG4Hit* hit = itr->second;
        if (!hit)
        {
          LogDebug("No PHG4Hit Found!");
          continue;
        }

        if (hit->get_trkid() == particle->get_track_id() || gsl_ran_binomial(m_RandomGenerator, detnoise, 1) > 0)
        {
          if (gsl_ran_binomial(m_RandomGenerator, dethiteff, 1) > 0)
          {
            PHGenFit::Measurement* meas = nullptr;
            if (dettype == Vertical_Plane)
            {
              if (Verbosity())
              {
                std::cout << "PHG4TrackFastSim::PseudoPatternRecognition -adding vertical plane hit ilayer: "
                          << ilayer << "; detphires: " << detphires << "; detradres: " << detradres << " \n";
                hit->identify();
              }
              meas = PHG4HitToMeasurementVerticalPlane(hit,
                                                       detphires, detradres);
            }
            else if (dettype == Cylinder)
            {
              if (Verbosity())
              {
                std::cout << "PHG4TrackFastSim::PseudoPatternRecognition -adding cylinder hit ilayer: "
                          << ilayer << "; detphires: " << detphires << "; detlonres : " << detlonres << " \n";
                hit->identify();
              }
              meas = PHG4HitToMeasurementCylinder(hit,
                                                  detphires, detlonres);
            }
            else
            {
              LogError("Type not implemented!");
              return Fun4AllReturnCodes::ABORTEVENT;
            }
            //            meas_out.push_back(meas);
            ordered_measurements.insert(make_pair(hit->get_avg_t(), meas));

            //meas->getMeasurement()->Print(); //DEBUG
          }
        }
      }
    } /*Loop layers within one detector layer*/
  }   /*Loop detector layers*/

  for (auto& pair : ordered_measurements)
  {
    meas_out.push_back(pair.second);

    if (Verbosity())
    {
      std::cout << "PHG4TrackFastSim::PseudoPatternRecognition - measruement at t =  " << pair.first << " ns: ";
      pair.second->getMeasurement()->Print();
    }
  }

  if (Verbosity())
  {
    std::cout << "PHG4TrackFastSim::PseudoPatternRecognition - meas_out.size = " << meas_out.size() << " for "
              << "particle: "
              << endl;
    particle->identify();

    std::cout << "PHG4TrackFastSim::PseudoPatternRecognition - seed_pos = "
              << seed_pos.x() << ", " << seed_pos.y() << ". " << seed_pos.z() << endl;
    std::cout << "PHG4TrackFastSim::PseudoPatternRecognition - seed_pos = "
              << seed_mom.x() << ", " << seed_mom.y() << ". " << seed_mom.z() << endl;
    std::cout << "PHG4TrackFastSim::PseudoPatternRecognition - seed_cov = "
              << sqrt(seed_cov[0][0]) << ", " << sqrt(seed_cov[1][1]) << ". " << sqrt(seed_cov[2][2])
              << ","
              << sqrt(seed_cov[3][3]) << ", " << sqrt(seed_cov[4][4]) << ". " << sqrt(seed_cov[5][5]) << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

SvtxTrack* PHG4TrackFastSim::MakeSvtxTrack(const PHGenFit::Track* phgf_track,
                                           const unsigned int truth_track_id,
                                           const unsigned int nmeas,
                                           const TVector3& vtx)
{
  double chi2 = phgf_track->get_chi2();
  double ndf = phgf_track->get_ndf();

  double pathlenth_from_first_meas = -999999;
  unique_ptr<genfit::MeasuredStateOnPlane> gf_state(new genfit::MeasuredStateOnPlane());

  //  if (_detector_type == Vertical_Plane)
  //  {
  //    pathlenth_orig_from_first_meas = phgf_track->extrapolateToPlane(*gf_state, vtx,
  //                                                                    TVector3(0., 0., 1.), 0);
  //  }
  //  else if (_detector_type == Cylinder)
  //    pathlenth_orig_from_first_meas = phgf_track->extrapolateToLine(*gf_state, vtx,
  //                                                                   TVector3(0., 0., 1.));
  //  else
  //  {
  //    LogError("Detector Type NOT implemented!");
  //    return nullptr;
  //  }

  // always extrapolate to a z-line through the vertex
  double pathlenth_orig_from_first_meas = phgf_track->extrapolateToLine(*gf_state, vtx,
                                                                        TVector3(0., 0., 1.));

  if (Verbosity() > 1)
  {
    cout << __PRETTY_FUNCTION__ << __LINE__ << " pathlenth_orig_from_first_meas = " << pathlenth_orig_from_first_meas << endl;
  }
  if (pathlenth_orig_from_first_meas < -999990)
  {
    LogError("Extraction faild!");
    return nullptr;
  }

  TVector3 mom = gf_state->getMom();
  TVector3 pos = gf_state->getPos();
  TMatrixDSym cov = gf_state->get6DCov();

  SvtxTrack_FastSim* out_track = new SvtxTrack_FastSim();
  out_track->set_truth_track_id(truth_track_id);
  /*!
	 * TODO: check the definition
	 *  1/p, u'/z', v'/z', u, v
	 *  u is defined as mom X beam line at POCA
	 *  so u is the dca2d direction
	 */
  double dca2d = gf_state->getState()[3];
  out_track->set_dca2d(dca2d);
  out_track->set_dca2d_error(gf_state->getCov()[3][3]);
  double dca3d = sqrt(dca2d * dca2d + gf_state->getState()[4] * gf_state->getState()[4]);
  out_track->set_dca(dca3d);

  out_track->set_chisq(chi2);
  out_track->set_ndf(ndf);
  out_track->set_charge(phgf_track->get_charge());

  out_track->set_num_measurements(nmeas);

  out_track->set_px(mom.Px());
  out_track->set_py(mom.Py());
  out_track->set_pz(mom.Pz());

  out_track->set_x(pos.X());
  out_track->set_y(pos.Y());
  out_track->set_z(pos.Z());

  for (int i = 0; i < 6; i++)
  {
    for (int j = i; j < 6; j++)
    {
      out_track->set_error(i, j, cov[i][j]);
    }
  }
  // State Projections
  map<string,double> state_name_path_map;
  for (map<string,double>::iterator iter = m_ZStateMap.begin(); iter != m_ZStateMap.end(); ++iter)
  {
    // Project to a plane at fixed z
    cout << "extrapolate to " << iter->second << endl;
    pathlenth_from_first_meas = phgf_track->extrapolateToPlane(*gf_state, TVector3(0., 0.,iter->second),
							       TVector3(1., 0., iter->second), 0);
    if (pathlenth_from_first_meas < -999990) // don't add for failure
    {
      continue;
    }
//    state_name_path_map.insert(make_pair(iter->first,pathlenth_from_first_meas));
    SvtxTrackState* state = new SvtxTrackState_v1(pathlenth_from_first_meas - pathlenth_orig_from_first_meas);
    state->set_x(gf_state->getPos().x());
    state->set_y(gf_state->getPos().y());
    state->set_z(gf_state->getPos().z());

    state->set_px(gf_state->getMom().x());
    state->set_py(gf_state->getMom().y());
    state->set_pz(gf_state->getMom().z());

    state->set_name(iter->first);
    for (int i = 0; i < 6; i++)
    {
      for (int j = i; j < 6; j++)
      {
        state->set_error(i, j, gf_state->get6DCov()[i][j]);
      }
    }
    out_track->insert_state(state);
    // the state is cloned on insert_state, so delete this copy here!
    delete state;

  }
  for (map<string,double>::iterator iter = m_CylinderStateMap.begin(); iter != m_CylinderStateMap.end(); ++iter)
  {
    pathlenth_from_first_meas = phgf_track->extrapolateToCylinder(*gf_state, iter->second, TVector3(0., 0., 0.),
								  TVector3(0., 0., 1.), 0);
    if (pathlenth_from_first_meas < -999990) // don't add for failure
    {
      continue;
    }
//    state_name_path_map.insert(make_pair(iter->first,pathlenth_from_first_meas));
    SvtxTrackState* state = new SvtxTrackState_v1(pathlenth_from_first_meas - pathlenth_orig_from_first_meas);
    state->set_x(gf_state->getPos().x());
    state->set_y(gf_state->getPos().y());
    state->set_z(gf_state->getPos().z());

    state->set_px(gf_state->getMom().x());
    state->set_py(gf_state->getMom().y());
    state->set_pz(gf_state->getMom().z());

    state->set_name(iter->first);
    for (int i = 0; i < 6; i++)
    {
      for (int j = i; j < 6; j++)
      {
        state->set_error(i, j, gf_state->get6DCov()[i][j]);
      }
    }
    out_track->insert_state(state);
    // the state is cloned on insert_state, so delete this copy here!
    delete state;
  }
/*  
  for (unsigned int i = 0; i < _state_names.size(); i++)
  {
    if ((_state_names[i] == "FHCAL") || (_state_names[i] == "FEMC") || (_state_names[i] == "EEMC"))
    {
      // Project to a plane at fixed z
      pathlenth_from_first_meas = phgf_track->extrapolateToPlane(*gf_state, TVector3(0., 0., _state_location[i]),
                                                                 TVector3(1., 0., _state_location[i]), 0);
      cout << "extracting path for " << _state_names[i] << " at " << _state_location[i] 
	   << " : " << pathlenth_from_first_meas << endl;
    }
    else if ((_state_names[i] == "CEMC") || (_state_names[i] == "HCALIN") || (_state_names[i] == "OHCAL"))
    {
      // Project to a cylinder at fixed r
      pathlenth_from_first_meas = phgf_track->extrapolateToCylinder(*gf_state, _state_location[i], TVector3(0., 0., 0.),
                                                                    TVector3(0., 0., 1.), 0);
      cout << "extracting path for " << _state_names[i] << " at " << _state_location[i] 
	   << " : " << pathlenth_from_first_meas << endl;
    }
    else
    {
      LogError("Unrecognized detector name for state projection");
      continue;
    }
*/
/*
    // if projection fails, bail out
  for (map<string,double>::iterator iter = state_name_path_map.begin(); iter != state_name_path_map.end(); ++iter)
  {
    SvtxTrackState* state = new SvtxTrackState_v1(iter->second - pathlenth_orig_from_first_meas);
    state->set_x(gf_state->getPos().x());
    state->set_y(gf_state->getPos().y());
    state->set_z(gf_state->getPos().z());

    state->set_px(gf_state->getMom().x());
    state->set_py(gf_state->getMom().y());
    state->set_pz(gf_state->getMom().z());

    state->set_name(iter->first);
    for (int i = 0; i < 6; i++)
    {
      for (int j = i; j < 6; j++)
      {
        state->set_error(i, j, gf_state->get6DCov()[i][j]);
      }
    }
    out_track->insert_state(state);
    // the state is cloned on insert_state, so delete this copy here!
    delete state;
  }
*/
  //  // State Projections
  //  {
  //    // Project to a cylinder at fixed r
  //    pathlenth_from_first_meas = phgf_track->extrapolateToCylinder(*gf_state, 2., TVector3(0., 0., 0.),
  //                                                                  TVector3(0., 0., 1.), 0, -1);
  //
  //    SvtxTrackState* state = new SvtxTrackState_v1(pathlenth_from_first_meas - pathlenth_orig_from_first_meas);
  //    state->set_x(gf_state->getPos().x());
  //    state->set_y(gf_state->getPos().y());
  //    state->set_z(gf_state->getPos().z());
  //
  //    state->set_px(gf_state->getMom().x());
  //    state->set_py(gf_state->getMom().y());
  //    state->set_pz(gf_state->getMom().z());
  //
  //    //    state->set_name(string("BeamPipe"));
  //
  //    for (int i = 0; i < 6; i++)
  //    {
  //      for (int j = i; j < 6; j++)
  //      {
  //        state->set_error(i, j, gf_state->get6DCov()[i][j]);
  //      }
  //    }
  //    out_track->insert_state(state);
  //
  //    cout << __PRETTY_FUNCTION__ << __LINE__;
  //    state->identify();
  //
  //    // the state is cloned on insert_state, so delete this copy here!
  //    delete state;
  //  }
  //
  //  cout << __PRETTY_FUNCTION__ << __LINE__;
  //  out_track->identify();

  return static_cast<SvtxTrack*>(out_track);
}

PHGenFit::PlanarMeasurement* PHG4TrackFastSim::PHG4HitToMeasurementVerticalPlane(
    const PHG4Hit* g4hit, const double phi_resolution,
    const double r_resolution)
{
  TVector3 pos(g4hit->get_avg_x(), g4hit->get_avg_y(), g4hit->get_avg_z());

  TVector3 v(pos.X(), pos.Y(), 0);
  v = 1 / v.Mag() * v;

  TVector3 u = v.Cross(TVector3(0, 0, 1));
  u = 1 / u.Mag() * u;

  double u_smear = 0.;
  double v_smear = 0.;
  if (m_SmearingFlag)
  {
   u_smear = gsl_ran_gaussian(m_RandomGenerator, phi_resolution);
   v_smear = gsl_ran_gaussian(m_RandomGenerator, r_resolution);
  }
  pos.SetX(g4hit->get_avg_x() + u_smear * u.X() + v_smear * v.X());
  pos.SetY(g4hit->get_avg_y() + u_smear * u.Y() + v_smear * v.Y());

  PHGenFit::PlanarMeasurement* meas = new PHGenFit::PlanarMeasurement(pos, u, v, phi_resolution,
                                                                      r_resolution);

  //	std::cout<<"------------\n";
  //	pos.Print();
  //	std::cout<<"at "<<istation<<" station, "<<ioctant << " octant \n";
  //	u.Print();
  //	v.Print();

  //dynamic_cast<PHGenFit::PlanarMeasurement*> (meas)->getMeasurement()->Print();

  return meas;
}

PHGenFit::PlanarMeasurement* PHG4TrackFastSim::PHG4HitToMeasurementCylinder(
    const PHG4Hit* g4hit, const double phi_resolution,
    const double z_resolution)
{
  TVector3 pos(g4hit->get_avg_x(), g4hit->get_avg_y(), g4hit->get_avg_z());

  TVector3 v(0, 0, 1);

  TVector3 u = v.Cross(TVector3(pos.X(), pos.Y(), 0));
  u = 1 / u.Mag() * u;

  double u_smear = 0.;
  double v_smear = 0.;
  if (m_SmearingFlag)
  {
  u_smear = gsl_ran_gaussian(m_RandomGenerator, phi_resolution);
  v_smear = gsl_ran_gaussian(m_RandomGenerator, z_resolution);
  }
  pos.SetX(g4hit->get_avg_x() + u_smear * u.X());
  pos.SetY(g4hit->get_avg_y() + u_smear * u.Y());
  pos.SetZ(g4hit->get_avg_z() + v_smear);

  PHGenFit::PlanarMeasurement* meas = new PHGenFit::PlanarMeasurement(pos, u, v, phi_resolution,
                                                                      z_resolution);

  //	std::cout<<"------------\n";
  //	pos.Print();
  //	std::cout<<"at "<<istation<<" station, "<<ioctant << " octant \n";
  //	u.Print();
  //	v.Print();

  //dynamic_cast<PHGenFit::PlanarMeasurement*> (meas)->getMeasurement()->Print();

  return meas;
}

PHGenFit::Measurement* PHG4TrackFastSim::VertexMeasurement(const TVector3& vtx, double dxy, double dz)
{
  TMatrixDSym cov(3);
  cov.Zero();
  cov(0, 0) = dxy * dxy;
  cov(1, 1) = dxy * dxy;
  cov(2, 2) = dz * dz;

  TVector3 pos = vtx;
  pos.SetX(vtx.X());
  pos.SetY(vtx.Y());
  pos.SetZ(vtx.Z());

  PHGenFit::Measurement* meas = new PHGenFit::SpacepointMeasurement(pos, cov);

  return meas;
}

void PHG4TrackFastSim::DisplayEvent() const
{
  if (m_DoEvtDisplayFlag && _fitter)
  {
    _fitter->displayEvent();
  }
  return;
}

void PHG4TrackFastSim::add_state_name(const std::string& stateName)
{
    _state_names.push_back(stateName);
    if (stateName == "FEMC" || stateName == "FHCAL" || stateName == "EEMC")
    {
      m_ZStateMap.insert(make_pair(stateName,NAN));
    }
    else if (stateName == "CEMC" || stateName == "HCALIN" || stateName == "HCALOUT")
    {
      m_CylinderStateMap.insert(make_pair(stateName,NAN));
    }
    return;
}
