/*!
 *  \file		PHGenFitTrkFitter.C
 *  \brief		Refit SvtxTracks with PHGenFit.
 *  \details	Refit SvtxTracks with PHGenFit.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#include "PHGenFitTrkFitter.h"

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

#include <phgenfit/Fitter.h>
#include <phgenfit/Measurement.h>                   // for Measurement
#include <phgenfit/PlanarMeasurement.h>
#include <phgenfit/SpacepointMeasurement.h>
#include <phgenfit/Track.h>

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

#include <GenFit/AbsMeasurement.h>                  // for AbsMeasurement
#include <GenFit/EventDisplay.h>                    // for EventDisplay
#include <GenFit/Exception.h>                       // for Exception
#include <GenFit/GFRaveConverters.h>
#include <GenFit/GFRaveTrackParameters.h>           // for GFRaveTrackParame...
#include <GenFit/GFRaveVertex.h>
#include <GenFit/GFRaveVertexFactory.h>
#include <GenFit/KalmanFitterInfo.h>
#include <GenFit/MeasuredStateOnPlane.h>
#include <GenFit/RKTrackRep.h>
#include <GenFit/Track.h>
#include <GenFit/TrackPoint.h>                      // for TrackPoint

//Rave
#include <rave/ConstantMagneticField.h>
#include <rave/VacuumPropagator.h>                  // for VacuumPropagator
#include <rave/VertexFactory.h>

#include <TClonesArray.h>
#include <TRotation.h>
#include <TTree.h>
#include <TVector3.h>
#include <TMath.h>                                  // for ATan2
#include <TMatrixDSymfwd.h>                         // for TMatrixDSym
#include <TMatrixFfwd.h>                            // for TMatrixF
#include <TMatrixT.h>                               // for TMatrixT, operator*
#include <TMatrixTSym.h>                            // for TMatrixTSym
#include <TMatrixTUtils.h>                          // for TMatrixTRow
#include <TVectorDfwd.h>                            // for TVectorD
#include <TVectorT.h>                               // for TVectorT

#include <cmath>                                   // for sqrt, NAN
#include <iostream>
#include <map>
#include <memory>
#include <utility>
#include <vector>

class PHField;
class TGeoManager;
namespace genfit { class AbsTrackRep; }

#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl
#define LogError(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl
#define LogWarning(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl

static constexpr float WILD_FLOAT = -9999;

#define _DEBUG_MODE_ 0

//#define _DEBUG_

using namespace std;

//______________________________________________________
namespace {

  // square
  template< class T > T square( T x ) { return x*x; }

  // convert gf state to SvtxTrackState_v1
  SvtxTrackState_v1 create_track_state( float pathlength, const genfit::MeasuredStateOnPlane* gf_state )
  {

    SvtxTrackState_v1 out( pathlength );
    out.set_x(gf_state->getPos().x());
    out.set_y(gf_state->getPos().y());
    out.set_z(gf_state->getPos().z());

    out.set_px(gf_state->getMom().x());
    out.set_py(gf_state->getMom().y());
    out.set_pz(gf_state->getMom().z());

    for (int i = 0; i < 6; i++)
    {
      for (int j = i; j < 6; j++)
      { out.set_error(i, j, gf_state->get6DCov()[i][j]); }
    }

    return out;

  }

}

/*
 * Constructor
 */
PHGenFitTrkFitter::PHGenFitTrkFitter(const string& name)
  : SubsysReco(name)
{
  Verbosity(0);
  _cluster_eval_tree_x = WILD_FLOAT;
  _cluster_eval_tree_y = WILD_FLOAT;
  _cluster_eval_tree_z = WILD_FLOAT;
  _cluster_eval_tree_gx = WILD_FLOAT;
  _cluster_eval_tree_gy = WILD_FLOAT;
  _cluster_eval_tree_gz = WILD_FLOAT;
}

/*
 * Init
 */
int PHGenFitTrkFitter::Init(PHCompositeNode* topNode)
{
  //	CreateNodes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * Init run
 */
int PHGenFitTrkFitter::InitRun(PHCompositeNode* topNode)
{
  CreateNodes(topNode);

  TGeoManager* tgeo_manager = PHGeomUtility::GetTGeoManager(topNode);
  PHField* field = PHFieldUtility::GetFieldMapNode(nullptr, topNode);

  //_fitter = new PHGenFit::Fitter("sPHENIX_Geo.root","sPHENIX.2d.root", 1.4 / 1.5);
  _fitter.reset( PHGenFit::Fitter::getInstance(tgeo_manager,
                                          field, _track_fitting_alg_name,
                                          "RKTrackRep", _do_evt_display) );
  _fitter->set_verbosity(Verbosity());

  if (!_fitter)
  {
    cerr << PHWHERE << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  //LogDebug(genfit::FieldManager::getInstance()->getFieldVal(TVector3(0, 0, 0)).Z());

  _vertex_finder.reset( new genfit::GFRaveVertexFactory(Verbosity()) );
  _vertex_finder->setMethod(_vertexing_method.data());


  if (!_vertex_finder)
  {
    cerr << PHWHERE << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (_do_eval)
  {
    if (Verbosity() >= 1)
      cout << PHWHERE << " Openning file: " << _eval_outname << endl;
    PHTFileServer::get().open(_eval_outname, "RECREATE");
    init_eval_tree();
  }

  // print disabled layers
  for( const auto& layer:_disabled_layers )
  { std::cout << PHWHERE << " Layer " << layer << " is disabled." << std::endl; }

  return Fun4AllReturnCodes::EVENT_OK;
}
/*
 * process_event():
 *  Call user instructions for every event.
 *  This function contains the analysis structure.
 *
 */
int PHGenFitTrkFitter::process_event(PHCompositeNode* topNode)
{
  _event++;

  if (Verbosity() > 1)
    std::cout << PHWHERE << "Events processed: " << _event << std::endl;
  //	if (_event % 1000 == 0)
  //		cout << PHWHERE << "Events processed: " << _event << endl;

  //cout << "Start PHGenFitTrkfitter::process_event" << endl;

  GetNodes(topNode);

  //! stands for Refit_GenFit_Tracks
  vector<genfit::Track*> rf_gf_tracks;
  vector<std::shared_ptr<PHGenFit::Track> > rf_phgf_tracks;

  map<unsigned int, unsigned int> svtxtrack_genfittrack_map;

  if (_trackmap_refit)
    _trackmap_refit->empty();

  // _trackmap is SvtxTrackMap from the node tree
  for (SvtxTrackMap::Iter iter = _trackmap->begin(); iter != _trackmap->end();
       ++iter)
  {
    SvtxTrack* svtx_track = iter->second;
    if(Verbosity() > 10)
      {
  cout << "   process SVTXTrack " << iter->first << endl;
  svtx_track->identify();
      }
    if (!svtx_track)
      continue;
    if (!(svtx_track->get_pt() > _fit_min_pT))
      continue;

    // This is the final track (re)fit. It does not include the collision vertex. If fit_primary_track is set, a refit including the vertex is done below.
    //! rf_phgf_track stands for Refit_PHGenFit_Track
    std::shared_ptr<PHGenFit::Track> rf_phgf_track = ReFitTrack(topNode, svtx_track);
    if (rf_phgf_track)
      {
  svtxtrack_genfittrack_map[svtx_track->get_id()] =  rf_phgf_tracks.size();
  rf_phgf_tracks.push_back(rf_phgf_track);
  if (rf_phgf_track->get_ndf() > _vertex_min_ndf)
    rf_gf_tracks.push_back(rf_phgf_track->getGenFitTrack());
  if(Verbosity() > 10) cout << "Done refitting input track" << svtx_track->get_id() << " or rf_phgf_track " <<   rf_phgf_tracks.size() << endl;
      }

  }

  /*
   * add tracks to event display
   * needs to make copied for smart ptrs will be destroied even
   * there are still references in TGeo::EventView
   */
  if (_do_evt_display)
  {
    vector<genfit::Track*> copy;
    for (genfit::Track* t : rf_gf_tracks)
    {
      copy.push_back(new genfit::Track(*t));
    }
    _fitter->getEventDisplay()->addEvent(copy);
  }

  //! find vertices using final tracks
  std::vector<genfit::GFRaveVertex*> rave_vertices;

  if (rf_gf_tracks.size() >= 2)
  {
    if(Verbosity() > 10)
      {cout << "Call Rave vertex finder" << endl;
  cout << " rf_gf_tracks size " << rf_gf_tracks.size() << endl;
  for(long unsigned int i=0;i<rf_gf_tracks.size();++i)
    {
      cout << "   track " << i << " num points " << rf_gf_tracks[i]->getNumPointsWithMeasurement() << endl;
    }
      }
    try
    {
      _vertex_finder->findVertices(&rave_vertices, rf_gf_tracks);
    }
    catch (...)
    {
      if (Verbosity() > 1)
        std::cout << PHWHERE << "GFRaveVertexFactory::findVertices failed!";
    }
  }

  if(Verbosity() > 10 && rave_vertices.size() == 0)
    {
      cout << PHWHERE << " Rave found no vertices - SvtxVertexMapRefit will be empty" << endl;
    }

  // Creates new SvtxVertex objects and copies in the rave vertex info and associated rave tracks. Then adds the vertices to _svtxmap
  // also makes a map of (trackid/vertex), _rave_vertex_gf_track_map for later use
  FillSvtxVertexMap(rave_vertices, rf_gf_tracks);

  // Finds the refitted rf_phgf_track corresponding to each SvtxTrackMap entry
  // Converts it to an SvtxTrack in MakeSvtxTrack
  // MakeSvtxTrack takes a vertex that it gets from the map made in FillSvtxVertex
  // If the refit was succesful, the track on the node tree is replaced with the new one
  // If not, the track is erased from the node tree
  for (SvtxTrackMap::Iter iter = _trackmap->begin(); iter != _trackmap->end();)
  {
    std::shared_ptr<PHGenFit::Track> rf_phgf_track;

    // find the genfit track that corresponds to this one on the node tree
    unsigned int itrack =0;
    if (svtxtrack_genfittrack_map.find(iter->second->get_id()) != svtxtrack_genfittrack_map.end())
    {
      itrack =
  svtxtrack_genfittrack_map[iter->second->get_id()];
      rf_phgf_track = rf_phgf_tracks[itrack];
    }

    if (rf_phgf_track)
    {
      SvtxVertex* vertex = nullptr;

      unsigned int ivert = 0;
      ivert = _rave_vertex_gf_track_map[itrack];

      if (_vertexmap_refit->size() > 0)
  {
    vertex = _vertexmap_refit->get(ivert);

    if(Verbosity() > 20) cout << PHWHERE << "     gf track " << itrack << " will add to track: _vertexmap_refit vertex " << ivert
            << " with position x,y,z = " << vertex->get_x() << "  " << vertex->get_y() << "  " << vertex->get_z() << endl;
  }
      std::shared_ptr<SvtxTrack> rf_track = MakeSvtxTrack(iter->second, rf_phgf_track,
                vertex);

#ifdef _DEBUG_
      cout << __LINE__ << endl;
#endif
      if (!rf_track)
  {
    //if (_output_mode == OverwriteOriginalNode)
#ifdef _DEBUG_
    LogDebug("!rf_track, continue.");
#endif
    if (_over_write_svtxtrackmap)
      {
        auto key = iter->first;
        ++iter;
        _trackmap->erase(key);
        continue;
      } else {
        ++iter;
        continue;
      }
  }

      //			delete vertex;//DEBUG

      //			rf_phgf_tracks.push_back(rf_phgf_track);
      //			rf_gf_tracks.push_back(rf_phgf_track->getGenFitTrack());

      if (!(_over_write_svtxtrackmap) || _output_mode == DebugMode)
        if (_trackmap_refit)
    {
          _trackmap_refit->insert(rf_track.get());
          //					delete rf_track;
        }

      if (_over_write_svtxtrackmap || _output_mode == DebugMode)
      {
        *(dynamic_cast<SvtxTrack_v1*>(iter->second)) =
            *(dynamic_cast<SvtxTrack_v1*>(rf_track.get()));
//				delete rf_track;
      }
    }
    else
    {
      if (_over_write_svtxtrackmap)
      {
        auto key = iter->first;
        ++iter;
        _trackmap->erase(key);
        continue;
      }
    }

    ++iter;
  }

#ifdef _DEBUG_
  cout << __LINE__ << endl;
#endif

  // Need to keep tracks if _do_evt_display
  if (!_do_evt_display)
  {
    rf_phgf_tracks.clear();
  }

#ifdef _DEBUG_
  cout << __LINE__ << endl;
#endif
  /*!
   * Optionally fit track as primary track by including collision vertex,
   This part need to be called after FillSvtxVertexMap
   */
  if (_fit_primary_tracks && rave_vertices.size() > 0)
  {
    _primary_trackmap->empty();

    //FIXME figure out which vertex to use.
    SvtxVertex* vertex = nullptr;
    if (_vertexmap_refit->size() > 0)
      vertex = _vertexmap_refit->get(0);

    // fix this, have to get vertex ID from track

    if (vertex)
    {
      for (SvtxTrackMap::ConstIter iter = _trackmap->begin();
           iter != _trackmap->end(); ++iter)
      {
        SvtxTrack* svtx_track = iter->second;
        if (!svtx_track)
          continue;
        if (!(svtx_track->get_pt() > _fit_min_pT))
          continue;
        /*!
         * rf_phgf_track stands for Refit_PHGenFit_Track
         */
        std::shared_ptr<PHGenFit::Track> rf_phgf_track = ReFitTrack(topNode, svtx_track,
                                                                    vertex);
        if (rf_phgf_track)
        {
          //					//FIXME figure out which vertex to use.
          //					SvtxVertex* vertex = nullptr;
          //					if (_vertexmap_refit->size() > 0)
          //						vertex = _vertexmap_refit->get(0);

          std::shared_ptr<SvtxTrack> rf_track = MakeSvtxTrack(svtx_track,
                                                              rf_phgf_track, vertex);
          //delete rf_phgf_track;
          if (!rf_track)
          {
#ifdef _DEBUG_
            LogDebug("!rf_track, continue.");
#endif
            continue;
          }
          _primary_trackmap->insert(rf_track.get());
        }
      }
    }
    else
    {
      LogError("No vertex in SvtxVertexMapRefit!");
    }
  }
#ifdef _DEBUG_
  cout << __LINE__ << endl;
#endif
  for (genfit::GFRaveVertex* vertex : rave_vertices)
  {
    delete vertex;
  }
  rave_vertices.clear();

  if (_do_eval)
  {
    fill_eval_tree(topNode);
  }
#ifdef _DEBUG_
  cout << __LINE__ << endl;
#endif
  return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * End
 */
int PHGenFitTrkFitter::End(PHCompositeNode* topNode)
{
  if (_do_eval)
  {
    if (Verbosity() >= 1)
      cout << PHWHERE << " Writing to file: " << _eval_outname << endl;
    PHTFileServer::get().cd(_eval_outname);
    _eval_tree->Write();
    _cluster_eval_tree->Write();
  }

  if (_do_evt_display)
    _fitter->displayEvent();

  return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * fill_eval_tree():
 */
void PHGenFitTrkFitter::fill_eval_tree(PHCompositeNode* topNode)
{
  //! Make sure to reset all the TTree variables before trying to set them.
  reset_eval_variables();

  int i = 0;
  for (PHG4TruthInfoContainer::ConstIterator itr =
           _truth_container->GetPrimaryParticleRange().first;
       itr != _truth_container->GetPrimaryParticleRange().second; ++itr)
    new ((*_tca_particlemap)[i++])(PHG4Particlev2)(
        *dynamic_cast<PHG4Particlev2*>(itr->second));

  i = 0;
  for (PHG4TruthInfoContainer::ConstVtxIterator itr =
           _truth_container->GetPrimaryVtxRange().first;
       itr != _truth_container->GetPrimaryVtxRange().second; ++itr)
    new ((*_tca_vtxmap)[i++])(PHG4VtxPointv1)(
        *dynamic_cast<PHG4VtxPointv1*>(itr->second));

  i = 0;
  for (SvtxTrackMap::ConstIter itr = _trackmap->begin();
       itr != _trackmap->end(); ++itr)
    new ((*_tca_trackmap)[i++])(SvtxTrack_v1)(
        *dynamic_cast<SvtxTrack_v1*>(itr->second));

  i = 0;
  if (_vertexmap)
    for (SvtxVertexMap::ConstIter itr = _vertexmap->begin();
         itr != _vertexmap->end(); ++itr)
      new ((*_tca_vertexmap)[i++])(SvtxVertex_v1)(
          *dynamic_cast<SvtxVertex_v1*>(itr->second));

  if (_trackmap_refit)
  {
    i = 0;
    for (SvtxTrackMap::ConstIter itr = _trackmap_refit->begin();
         itr != _trackmap_refit->end(); ++itr)
      new ((*_tca_trackmap_refit)[i++])(SvtxTrack_v1)(
          *dynamic_cast<SvtxTrack_v1*>(itr->second));
  }

  if (_fit_primary_tracks)
  {
    i = 0;
    for (SvtxTrackMap::ConstIter itr = _primary_trackmap->begin();
         itr != _primary_trackmap->end(); ++itr)
      new ((*_tca_primtrackmap)[i++])(SvtxTrack_v1)(
          *dynamic_cast<SvtxTrack_v1*>(itr->second));
  }

  if (_vertexmap_refit)
  {
    i = 0;
    for (SvtxVertexMap::ConstIter itr = _vertexmap_refit->begin();
         itr != _vertexmap_refit->end(); ++itr)
      new ((*_tca_vertexmap_refit)[i++])(SvtxVertex_v1)(
          *dynamic_cast<SvtxVertex_v1*>(itr->second));
  }

  _eval_tree->Fill();

  return;
}

/*
 * init_eval_tree
 */
void PHGenFitTrkFitter::init_eval_tree()
{
  if (!_tca_particlemap)
    _tca_particlemap = new TClonesArray("PHG4Particlev2");
  if (!_tca_vtxmap)
    _tca_vtxmap = new TClonesArray("PHG4VtxPointv1");

  if (!_tca_trackmap)
    _tca_trackmap = new TClonesArray("SvtxTrack_v1");
  if (!_tca_vertexmap)
    _tca_vertexmap = new TClonesArray("SvtxVertex_v1");
  if (!_tca_trackmap_refit)
    _tca_trackmap_refit = new TClonesArray("SvtxTrack_v1");
  if (_fit_primary_tracks)
    if (!_tca_primtrackmap)
      _tca_primtrackmap = new TClonesArray("SvtxTrack_v1");
  if (!_tca_vertexmap_refit)
    _tca_vertexmap_refit = new TClonesArray("SvtxVertex_v1");

  //! create TTree
  _eval_tree = new TTree("T", "PHGenFitTrkFitter Evaluation");

  _eval_tree->Branch("PrimaryParticle", _tca_particlemap);
  _eval_tree->Branch("TruthVtx", _tca_vtxmap);

  _eval_tree->Branch("SvtxTrack", _tca_trackmap);
  _eval_tree->Branch("SvtxVertex", _tca_vertexmap);
  _eval_tree->Branch("SvtxTrackRefit", _tca_trackmap_refit);
  if (_fit_primary_tracks)
    _eval_tree->Branch("PrimSvtxTrack", _tca_primtrackmap);
  _eval_tree->Branch("SvtxVertexRefit", _tca_vertexmap_refit);

  _cluster_eval_tree = new TTree("cluster_eval", "cluster eval tree");
  _cluster_eval_tree->Branch("x", &_cluster_eval_tree_x, "x/F");
  _cluster_eval_tree->Branch("y", &_cluster_eval_tree_y, "y/F");
  _cluster_eval_tree->Branch("z", &_cluster_eval_tree_z, "z/F");
  _cluster_eval_tree->Branch("gx", &_cluster_eval_tree_gx, "gx/F");
  _cluster_eval_tree->Branch("gy", &_cluster_eval_tree_gy, "gy/F");
  _cluster_eval_tree->Branch("gz", &_cluster_eval_tree_gz, "gz/F");
}

/*
 * reset_eval_variables():
 *  Reset all the tree variables to their default values.
 *  Needs to be called at the start of every event
 */
void PHGenFitTrkFitter::reset_eval_variables()
{
  _tca_particlemap->Clear();
  _tca_vtxmap->Clear();

  _tca_trackmap->Clear();
  _tca_vertexmap->Clear();
  _tca_trackmap_refit->Clear();
  if (_fit_primary_tracks)
    _tca_primtrackmap->Clear();
  _tca_vertexmap_refit->Clear();

  _cluster_eval_tree_x = WILD_FLOAT;
  _cluster_eval_tree_y = WILD_FLOAT;
  _cluster_eval_tree_z = WILD_FLOAT;
  _cluster_eval_tree_gx = WILD_FLOAT;
  _cluster_eval_tree_gy = WILD_FLOAT;
  _cluster_eval_tree_gz = WILD_FLOAT;
}

int PHGenFitTrkFitter::CreateNodes(PHCompositeNode* topNode)
{
  // create nodes...
  PHNodeIterator iter(topNode);

  PHCompositeNode* dstNode = static_cast<PHCompositeNode*>(iter.findFirst(
      "PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cerr << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  PHNodeIterator iter_dst(dstNode);

  // Create the SVTX node
  PHCompositeNode* tb_node = dynamic_cast<PHCompositeNode*>(iter_dst.findFirst(
      "PHCompositeNode", "SVTX"));
  if (!tb_node)
  {
    tb_node = new PHCompositeNode("SVTX");
    dstNode->addNode(tb_node);
    if (Verbosity() > 0)
      cout << "SVTX node added" << endl;
  }

  if (!(_over_write_svtxtrackmap) || _output_mode == DebugMode)
  {
    _trackmap_refit = new SvtxTrackMap_v1;
    PHIODataNode<PHObject>* tracks_node = new PHIODataNode<PHObject>(
        _trackmap_refit, "SvtxTrackMapRefit", "PHObject");
    tb_node->addNode(tracks_node);
    if (Verbosity() > 0)
      cout << "Svtx/SvtxTrackMapRefit node added" << endl;
  }

  if (_fit_primary_tracks)
  {
    _primary_trackmap = new SvtxTrackMap_v1;
    PHIODataNode<PHObject>* primary_tracks_node =
        new PHIODataNode<PHObject>(_primary_trackmap, "PrimaryTrackMap",
                                   "PHObject");
    tb_node->addNode(primary_tracks_node);
    if (Verbosity() > 0)
      cout << "Svtx/PrimaryTrackMap node added" << endl;
  }

  // always write final vertex results to SvtxVertexMapRefit
  _vertexmap_refit = new SvtxVertexMap_v1;
  PHIODataNode<PHObject>* vertexes_node = new PHIODataNode<PHObject>(
                     _vertexmap_refit, "SvtxVertexMapRefit", "PHObject");
  tb_node->addNode(vertexes_node);
  if (Verbosity() > 0)
    cout << "Svtx/SvtxVertexMapRefit node added" << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//______________________________________________________
void PHGenFitTrkFitter::disable_layer( int layer, bool disabled )
{
  if( disabled ) _disabled_layers.insert( layer );
  else _disabled_layers.erase( layer );
}

//______________________________________________________
void PHGenFitTrkFitter::set_disabled_layers( const std::set<int>& layers )
{ _disabled_layers = layers; }

//______________________________________________________
void PHGenFitTrkFitter::clear_disabled_layers()
{ _disabled_layers.clear(); }

//______________________________________________________
const std::set<int>& PHGenFitTrkFitter::get_disabled_layers() const
{ return _disabled_layers; }

/*
 * GetNodes():
 *  Get all the all the required nodes off the node tree
 */
int PHGenFitTrkFitter::GetNodes(PHCompositeNode* topNode)
{
  //DST objects
  //Truth container
  _truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode,
                                                                "G4TruthInfo");
  if (!_truth_container && _event < 2)
  {
    cout << PHWHERE << " PHG4TruthInfoContainer node not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Input Svtx Clusters
  //_clustermap = findNode::getClass<SvtxClusterMap>(topNode, "SvtxClusterMap");
  _clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_clustermap && _event < 2)
  {
    cout << PHWHERE << " TRKR_CLUSTER node not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Input Svtx Tracks
  _trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!_trackmap && _event < 2)
  {
    cout << PHWHERE << " SvtxTrackMap node not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Input Svtx Vertices
  _vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!_vertexmap && _event < 2)
  {
    cout << PHWHERE << " SvtxVertexrMap node not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

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

  // Output Primary Svtx Tracks
  if (_fit_primary_tracks)
  {
    _primary_trackmap = findNode::getClass<SvtxTrackMap>(topNode,
                                                         "PrimaryTrackMap");
    if (!_primary_trackmap && _event < 2)
    {
      cout << PHWHERE << " PrimaryTrackMap node not found on node tree"
           << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  // Output Svtx Vertices
  _vertexmap_refit = findNode::getClass<SvtxVertexMap>(topNode,
                   "SvtxVertexMapRefit");
  if (!_vertexmap_refit && _event < 2)
    {
      cout << PHWHERE << " SvtxVertexMapRefit node not found on node tree"
           << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

//struct CompMeasurementByR {
//  bool operator() (PHGenFit::Measurement *m1,PHGenFit::Measurement *m2) {
//	  float x1 = m1->getMeasurement()
//
//	  return (i<j);}
//} myobject;

/*
 * fit track with SvtxTrack as input seed.
 * \param intrack Input SvtxTrack
 * \param invertex Input Vertex, if fit track as a primary vertex
 */
//PHGenFit::Track* PHGenFitTrkFitter::ReFitTrack(PHCompositeNode *topNode, const SvtxTrack* intrack,
std::shared_ptr<PHGenFit::Track> PHGenFitTrkFitter::ReFitTrack(PHCompositeNode* topNode, const SvtxTrack* intrack,
                                                               const SvtxVertex* invertex)
{
  //std::shared_ptr<PHGenFit::Track> empty_track(nullptr);
  if (!intrack)
  {
    cerr << PHWHERE << " Input SvtxTrack is nullptr!" << endl;
    return nullptr;
  }

  PHG4CylinderGeomContainer* geom_container_intt = findNode::getClass<
      PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");

  PHG4CylinderGeomContainer* geom_container_mvtx = findNode::getClass<
      PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");

  // prepare seed
  TVector3 seed_mom(100, 0, 0);
  TVector3 seed_pos(0, 0, 0);
  TMatrixDSym seed_cov(6);
  for (int i = 0; i < 6; i++)
  {
    for (int j = 0; j < 6; j++)
    {
      seed_cov[i][j] = 100.;
    }
  }

  // Create measurements
  std::vector<PHGenFit::Measurement*> measurements;

  /*!
   * if fit track as a primary track
   */

  /*
  if(invertex and Verbosity() >= 2)
    {
      LogDebug(invertex->size_tracks());
      LogDebug(invertex->get_chisq());
      LogDebug(invertex->get_ndof());
      for (unsigned int i = 0; i < 3; i++)
  for (unsigned int j = 0; j < 3; j++)
    {
      LogDebug(invertex->get_error(i,j));
    }
    }
  */

  /*!
   *
   */
#if _DEBUG_MODE_ == 1
  if (invertex
      //			and invertex->size_tracks() == 1
  )
  {
    TRandom3 rand(0);
    double dxy = 0.0007;  //7 um
    double dz = 0.003;    //30 um

    TVector3 pos(invertex->get_x(), invertex->get_y(), invertex->get_z());
    TMatrixDSym cov(3);

    // Use smeared position instead of reco'd one.
    double x = rand.Gaus(0, dxy);
    double y = rand.Gaus(0, dxy);
    double z = rand.Gaus(0, dz);
    pos.SetXYZ(x, y, z);

    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        cov[i][j] = 0;

    cov[0][0] = dxy * dxy;
    cov[1][1] = dxy * dxy;
    cov[2][2] = dz * dz;

    PHGenFit::Measurement* meas = new PHGenFit::SpacepointMeasurement(
        pos, cov);
    measurements.push_back(meas);
  }
#else
  //! 1000 is a arbitrary number for now
  const double vertex_chi2_over_dnf_cut = 1000;
  const double vertex_cov_element_cut = 10000;  //arbitrary cut cm*cm

  if (invertex and invertex->size_tracks() > 1 and invertex->get_chisq() / invertex->get_ndof() < vertex_chi2_over_dnf_cut)
  {
    TVector3 pos(invertex->get_x(), invertex->get_y(), invertex->get_z());
    TMatrixDSym cov(3);
    cov.Zero();
    bool is_vertex_cov_sane = true;
    for (unsigned int i = 0; i < 3; i++)
      for (unsigned int j = 0; j < 3; j++)
      {
        cov(i, j) = invertex->get_error(i, j);

        if (i == j)
        {
          if (!(invertex->get_error(i, j) > 0 and invertex->get_error(i, j) < vertex_cov_element_cut))
            is_vertex_cov_sane = false;
        }
      }

    if (is_vertex_cov_sane)
    {
      PHGenFit::Measurement* meas = new PHGenFit::SpacepointMeasurement(
          pos, cov);
      measurements.push_back(meas);
      if(Verbosity() >= 2)
  {
    meas->getMeasurement()->Print();
  }
    }
  }
#endif

  // sort clusters with radius before fitting
  if(Verbosity() > 10)   intrack->identify();
  std::map<float, TrkrDefs::cluskey> m_r_cluster_id;
  for (auto iter = intrack->begin_cluster_keys();
       iter != intrack->end_cluster_keys(); ++iter)
  {
    TrkrDefs::cluskey cluster_key = *iter;
    TrkrCluster* cluster = _clustermap->findCluster(cluster_key);
    float x = cluster->getPosition(0);
    float y = cluster->getPosition(1);
    float r = sqrt(square(x) + square(y));
    m_r_cluster_id.insert(std::pair<float, TrkrDefs::cluskey>(r, cluster_key));
    int layer_out = TrkrDefs::getLayer(cluster_key);
    if(Verbosity() > 10) cout << "    Layer " << layer_out << " cluster " << cluster_key << " radius " << r << endl;
  }

  for (auto iter = m_r_cluster_id.begin();
       iter != m_r_cluster_id.end();
       ++iter)
  {
    TrkrDefs::cluskey cluster_key = iter->second;
    const int layer = TrkrDefs::getLayer(cluster_key);

    // skip disabled layers
    if( _disabled_layers.find( layer ) != _disabled_layers.end() )
    { continue; }

    TrkrCluster* cluster = _clustermap->findCluster(cluster_key);
    if (!cluster)
    {
      LogError("No cluster Found!");
      continue;
    }

#ifdef _DEBUG_
    int output_layer = TrkrDefs::getLayer(cluster_key);
    cout
        << __LINE__
        << ": ID: " << cluster_key
        << ": layer: " << output_layer
        << endl;
#endif

    TVector3 pos(cluster->getPosition(0), cluster->getPosition(1), cluster->getPosition(2));

    seed_mom.SetPhi(pos.Phi());
    seed_mom.SetTheta(pos.Theta());

    //TODO use u, v explicitly?
    TVector3 n(cluster->getPosition(0), cluster->getPosition(1), 0);

    //------------------------------
    // new

    // Replace n for the silicon subsystems

    // get the trkrid
    unsigned int trkrid = TrkrDefs::getTrkrId(cluster_key);

    if(trkrid == TrkrDefs::mvtxId)
      {
  int stave_index = MvtxDefs::getStaveId(cluster_key);
  int chip_index = MvtxDefs::getChipId(cluster_key);

  double ladder_location[3] = {0.0, 0.0, 0.0};
  CylinderGeom_Mvtx* geom =
          dynamic_cast<CylinderGeom_Mvtx*>(geom_container_mvtx->GetLayerGeom(layer));
  // returns the center of the sensor in world coordinates - used to get the ladder phi location
  geom->find_sensor_center(stave_index, 0,
         0, chip_index, ladder_location);

  //cout << " MVTX stave phi tilt = " <<  geom->get_stave_phi_tilt()
  //   << " seg.X " << ladder_location[0] << " seg.Y " << ladder_location[1] << " seg.Z " << ladder_location[2] << endl;
  n.SetXYZ(ladder_location[0], ladder_location[1], 0);
  n.RotateZ(geom->get_stave_phi_tilt());
      }
    else if(trkrid == TrkrDefs::inttId)
      {
  CylinderGeomIntt* geom =
          dynamic_cast<CylinderGeomIntt*>(geom_container_intt->GetLayerGeom(layer));
  double hit_location[3] = {0.0, 0.0, 0.0};
  geom->find_segment_center(InttDefs::getLadderZId(cluster_key),
          InttDefs::getLadderPhiId(cluster_key), hit_location);

  //cout << " Intt strip phi tilt = " <<  geom->get_strip_phi_tilt()
  //   << " seg.X " << hit_location[0] << " seg.Y " << hit_location[1] << " seg.Z " << hit_location[2] << endl;
  n.SetXYZ(hit_location[0], hit_location[1], 0);
  n.RotateZ(geom->get_strip_phi_tilt());
      }
    // end new
    //-----------------

    PHGenFit::Measurement* meas = new PHGenFit::PlanarMeasurement(pos, n,
                  cluster->getRPhiError(), cluster->getZError());

    if(Verbosity() > 10)
      {
  cout << "Add meas layer " << layer << " cluskey " << cluster_key
       << endl
       << " pos.X " << pos.X() << " pos.Y " << pos.Y() << " pos.Z " << pos.Z()
       << "  n.X " <<  n.X() << " n.Y " << n.Y()
       << " RPhiErr " << cluster->getRPhiError()
       << " ZErr " << cluster->getZError()
       << endl;
      }
      measurements.push_back(meas);
  }

  /*!
   * mu+:	-13
   * mu-:	13
   * pi+:	211
   * pi-:	-211
   * e-:	11
   * e+:	-11
   */
  //TODO Add multiple TrackRep choices.
  //int pid = 211;
  genfit::AbsTrackRep* rep = new genfit::RKTrackRep(_primary_pid_guess);
  std::shared_ptr<PHGenFit::Track> track(new PHGenFit::Track(rep, seed_pos, seed_mom,
                                                             seed_cov));

  //TODO unsorted measurements, should use sorted ones?
  track->addMeasurements(measurements);

  /*!
   *  Fit the track
   *  ret code 0 means 0 error or good status
   */
  if (_fitter->processTrack(track.get(), false) != 0)
  {
    if (Verbosity() >= 1)
      {
  LogWarning("Track fitting failed");
  cout << " track->getChisq() " << track->get_chi2() << " get_ndf " << track->get_ndf()
       << " mom.X " << track->get_mom().X()
       << " mom.Y " << track->get_mom().Y()
       << " mom.Z " << track->get_mom().Z()
       << endl;
      }
    //delete track;
    return nullptr;
  }

  if(Verbosity() > 10)
    cout << " track->getChisq() " << track->get_chi2() << " get_ndf " << track->get_ndf()
   << " mom.X " << track->get_mom().X()
   << " mom.Y " << track->get_mom().Y()
   << " mom.Z " << track->get_mom().Z()
   << endl;

  return track;
}

/*
 * Make SvtxTrack from PHGenFit::Track and SvtxTrack
 */
//SvtxTrack* PHGenFitTrkFitter::MakeSvtxTrack(const SvtxTrack* svtx_track,
std::shared_ptr<SvtxTrack> PHGenFitTrkFitter::MakeSvtxTrack(const SvtxTrack* svtx_track,
                                                            const std::shared_ptr<PHGenFit::Track>& phgf_track, const SvtxVertex* vertex)
{
  double chi2 = phgf_track->get_chi2();
  double ndf = phgf_track->get_ndf();

  TVector3 vertex_position(0, 0, 0);
  TMatrixF vertex_cov(3, 3);
  double dvr2 = 0;
  double dvz2 = 0;

  if (_use_truth_vertex)
  {
    PHG4VtxPoint* first_point = _truth_container->GetPrimaryVtx(_truth_container->GetPrimaryVertexIndex());
    vertex_position.SetXYZ(first_point->get_x(), first_point->get_y(), first_point->get_z());
    if (Verbosity() > 1)
    {
      cout << PHWHERE << "Using: truth vertex: {" << vertex_position.X() << ", " << vertex_position.Y() << ", " << vertex_position.Z() << "} " << endl;
    }
  }
  else if (vertex)
  {
    vertex_position.SetXYZ(vertex->get_x(), vertex->get_y(),
                           vertex->get_z());
    dvr2 = vertex->get_error(0, 0) + vertex->get_error(1, 1);
    dvz2 = vertex->get_error(2, 2);

    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        vertex_cov[i][j] = vertex->get_error(i, j);
  }

  std::unique_ptr<genfit::MeasuredStateOnPlane> gf_state_beam_line_ca;
  try
  {
    gf_state_beam_line_ca.reset(phgf_track->extrapolateToLine(vertex_position,
                                                              TVector3(0., 0., 1.)));
  }
  catch (...)
  {
    if (Verbosity() >= 2)
      LogWarning("extrapolateToLine failed!");
  }
  if (!gf_state_beam_line_ca) return nullptr;

  /*!
   *  1/p, u'/z', v'/z', u, v
   *  u is defined as momentum X beam line at POCA of the beam line
   *  v is alone the beam line
   *  so u is the dca2d direction
   */

  double u = gf_state_beam_line_ca->getState()[3];
  double v = gf_state_beam_line_ca->getState()[4];

  double du2 = gf_state_beam_line_ca->getCov()[3][3];
  double dv2 = gf_state_beam_line_ca->getCov()[4][4];
  //cout << PHWHERE << "        u " << u << " v " << v << " du2 " << du2 << " dv2 " << dv2 << " dvr2 " << dvr2 << endl;
  //delete gf_state_beam_line_ca;

  // create new track
  std::shared_ptr<SvtxTrack_v1> out_track( new SvtxTrack_v1() );

  // assign all members from old
  /* we need asignment operator instead of copy constructor, because the former makes a deep copy of the track state, as opposed to the latter */
  *out_track = *static_cast<const SvtxTrack_v1*>(svtx_track);

  // clear states and insert empty one for vertex position
  out_track->clear_states();
  {
    /*
    insert first, dummy state, as done in SvtxTrack_v1 constructor,
    so that the track state list is never empty. Note that insert_state, despite taking a pointer as argument,
    does not take ownership of the state
    */
    SvtxTrackState_v1 first(0.0);
    out_track->insert_state( &first );
  }

  out_track->set_dca2d(u);
  out_track->set_dca2d_error(sqrt(du2 + dvr2));

  std::unique_ptr<genfit::MeasuredStateOnPlane> gf_state_vertex_ca;
  try
  {
    gf_state_vertex_ca.reset( phgf_track->extrapolateToPoint(vertex_position) );
  }
  catch (...)
  {
    if (Verbosity() >= 2)
      LogWarning("extrapolateToPoint failed!");
  }
  if (!gf_state_vertex_ca)
  {
    //delete out_track;
    return nullptr;
  }

  TVector3 mom = gf_state_vertex_ca->getMom();
  TVector3 pos = gf_state_vertex_ca->getPos();
  TMatrixDSym cov = gf_state_vertex_ca->get6DCov();

  //	genfit::MeasuredStateOnPlane* gf_state_vertex_ca =
  //			phgf_track->extrapolateToLine(vertex_position,
  //					TVector3(0., 0., 1.));

  u = gf_state_vertex_ca->getState()[3];
  v = gf_state_vertex_ca->getState()[4];

  du2 = gf_state_vertex_ca->getCov()[3][3];
  dv2 = gf_state_vertex_ca->getCov()[4][4];

  double dca3d = sqrt(u * u + v * v);
  double dca3d_error = sqrt(du2 + dv2 + dvr2 + dvz2);

  out_track->set_dca(dca3d);
  out_track->set_dca_error(dca3d_error);

  /*!
   * dca3d_xy, dca3d_z
   */

  /*
  // Rotate from u,v,n to r: n X Z, Z': n X r, n using 5D state/cov
  // commented on 2017-10-09

  TMatrixF pos_in(3,1);
  TMatrixF cov_in(3,3);
  pos_in[0][0] = gf_state_vertex_ca->getState()[3];
  pos_in[1][0] = gf_state_vertex_ca->getState()[4];
  pos_in[2][0] = 0.;

  cov_in[0][0] = gf_state_vertex_ca->getCov()[3][3];
  cov_in[0][1] = gf_state_vertex_ca->getCov()[3][4];
  cov_in[0][2] = 0.;
  cov_in[1][0] = gf_state_vertex_ca->getCov()[4][3];
  cov_in[1][1] = gf_state_vertex_ca->getCov()[4][4];
  cov_in[1][2] = 0.;
  cov_in[2][0] = 0.;
  cov_in[2][1] = 0.;
  cov_in[2][2] = 0.;

  TMatrixF pos_out(3,1);
  TMatrixF cov_out(3,3);

  TVector3 vu = gf_state_vertex_ca->getPlane().get()->getU();
  TVector3 vv = gf_state_vertex_ca->getPlane().get()->getV();
  TVector3 vn = vu.Cross(vv);

  pos_cov_uvn_to_rz(vu, vv, vn, pos_in, cov_in, pos_out, cov_out);

  //! vertex cov in (u',v',n')
  TMatrixF vertex_cov_out(3,3);

  get_vertex_error_uvn(vu,vv,vn, vertex_cov, vertex_cov_out);

  float dca3d_xy = pos_out[0][0];
  float dca3d_z  = pos_out[1][0];

  float dca3d_xy_error = sqrt(cov_out[0][0] + vertex_cov_out[0][0]);
  float dca3d_z_error  = sqrt(cov_out[1][1] + vertex_cov_out[1][1]);

    //Begin DEBUG
//	LogDebug("rotation debug---------- ");
//	gf_state_vertex_ca->Print();
//	LogDebug("dca rotation---------- ");
//	pos_out = pos_in;
//	cov_out = cov_in;
//	pos_in.Print();
//	cov_in.Print();
//	pos_out.Print();
//	cov_out.Print();
//	cout
//		<<"dca3d_xy: "<<dca3d_xy <<" +- "<<dca3d_xy_error*dca3d_xy_error
//		<<"; dca3d_z: "<<dca3d_z<<" +- "<< dca3d_z_error*dca3d_z_error
//		<<"\n";
//	gf_state_vertex_ca->get6DCov().Print();
//	LogDebug("vertex rotation---------- ");
//	vertex_position.Print();
//	vertex_cov.Print();
//	vertex_cov_out.Print();
  //End DEBUG
  */

  //
  // in: X, Y, Z; out; r: n X Z, Z X r, Z

  float dca3d_xy = NAN;
  float dca3d_z = NAN;
  float dca3d_xy_error = NAN;
  float dca3d_z_error = NAN;

  try
  {
    TMatrixF pos_in(3, 1);
    TMatrixF cov_in(3, 3);
    TMatrixF pos_out(3, 1);
    TMatrixF cov_out(3, 3);

    TVectorD state6(6);      // pos(3), mom(3)
    TMatrixDSym cov6(6, 6);  //

    gf_state_vertex_ca->get6DStateCov(state6, cov6);

    TVector3 vn(state6[3], state6[4], state6[5]);

    // mean of two multivariate gaussians Pos - Vertex
    pos_in[0][0] = state6[0] - vertex_position.X();
    pos_in[1][0] = state6[1] - vertex_position.Y();
    pos_in[2][0] = state6[2] - vertex_position.Z();

    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        cov_in[i][j] = cov6[i][j] + vertex_cov[i][j];
      }
    }

    // vn is momentum vector, pos_in is position vector (of what?)
    pos_cov_XYZ_to_RZ(vn, pos_in, cov_in, pos_out, cov_out);

    if(Verbosity() > 30)
      {
  cout << " vn.X " << vn.X() << " vn.Y " << vn.Y() << " vn.Z " << vn.Z() << endl;
  cout << " pos_in.X " << pos_in[0][0] << " pos_in.Y " << pos_in[1][0] << " pos_in.Z " << pos_in[2][0] << endl;
  cout << " pos_out.X " << pos_out[0][0] << " pos_out.Y " << pos_out[1][0] << " pos_out.Z " << pos_out[2][0] << endl;
      }


    dca3d_xy = pos_out[0][0];
    dca3d_z = pos_out[2][0];
    dca3d_xy_error = sqrt(cov_out[0][0]);
    dca3d_z_error = sqrt(cov_out[2][2]);

#ifdef _DEBUG_
    cout << __LINE__ << ": Vertex: ----------------" << endl;
    vertex_position.Print();
    vertex_cov.Print();

    cout << __LINE__ << ": State: ----------------" << endl;
    state6.Print();
    cov6.Print();

    cout << __LINE__ << ": Mean: ----------------" << endl;
    pos_in.Print();
    cout << "===>" << endl;
    pos_out.Print();

    cout << __LINE__ << ": Cov: ----------------" << endl;
    cov_in.Print();
    cout << "===>" << endl;
    cov_out.Print();

    cout << endl;
#endif
  }
  catch (...)
  {
    if (Verbosity() > 0)
      LogWarning("DCA calculationfailed!");
  }

  out_track->set_dca3d_xy(dca3d_xy);
  out_track->set_dca3d_z(dca3d_z);
  out_track->set_dca3d_xy_error(dca3d_xy_error);
  out_track->set_dca3d_z_error(dca3d_z_error);

  //if(gf_state_vertex_ca) delete gf_state_vertex_ca;

  out_track->set_chisq(chi2);
  out_track->set_ndf(ndf);
  out_track->set_charge(phgf_track->get_charge());

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

  //	for (SvtxTrack::ConstClusterIter iter = svtx_track->begin_clusters();
  //			iter != svtx_track->end_clusters(); ++iter) {
  //		unsigned int cluster_id = *iter;
  //		SvtxCluster* cluster = _clustermap->get(cluster_id);
  //		if (!cluster) {
  //			LogError("No cluster Found!");
  //			continue;
  //		}
  //		//cluster->identify; //DEBUG
  //
  //		//unsigned int l = cluster->get_layer();
  //
  //		TVector3 pos(cluster->get_x(), cluster->get_y(), cluster->get_z());
  //
  //		double radius = pos.Pt();
  //
  //		std::shared_ptr<genfit::MeasuredStateOnPlane> gf_state = nullptr;
  //		try {
  //			gf_state = std::shared_ptr < genfit::MeasuredStateOnPlane
  //					> (phgf_track->extrapolateToCylinder(radius,
  //							TVector3(0, 0, 0), TVector3(0, 0, 1), 0));
  //		} catch (...) {
  //			if (Verbosity() >= 2)
  //				LogWarning("Exrapolation failed!");
  //		}
  //		if (!gf_state) {
  //			if (Verbosity() > 1)
  //				LogWarning("Exrapolation failed!");
  //			continue;
  //		}
  //
  //		//SvtxTrackState* state = new SvtxTrackState_v1(radius);
  //		std::shared_ptr<SvtxTrackState> state = std::shared_ptr<SvtxTrackState> (new SvtxTrackState_v1(radius));
  //		state->set_x(gf_state->getPos().x());
  //		state->set_y(gf_state->getPos().y());
  //		state->set_z(gf_state->getPos().z());
  //
  //		state->set_px(gf_state->getMom().x());
  //		state->set_py(gf_state->getMom().y());
  //		state->set_pz(gf_state->getMom().z());
  //
  //		//gf_state->getCov().Print();
  //
  //		for (int i = 0; i < 6; i++) {
  //			for (int j = i; j < 6; j++) {
  //				state->set_error(i, j, gf_state->get6DCov()[i][j]);
  //			}
  //		}
  //
  //		out_track->insert_state(state.get());
  //
  //#ifdef _DEBUG_
  //		cout
  //		<<__LINE__
  //		<<": " << radius <<" => "
  //		<<sqrt(state->get_x()*state->get_x() + state->get_y()*state->get_y())
  //		<<endl;
  //#endif
  //	}

#ifdef _DEBUG_
  cout << __LINE__ << endl;
#endif

  const genfit::Track* gftrack = phgf_track->getGenFitTrack();
  const genfit::AbsTrackRep* rep = gftrack->getCardinalRep();
  for (unsigned int id = 0; id < gftrack->getNumPointsWithMeasurement(); ++id)
  {
    genfit::TrackPoint* trpoint = gftrack->getPointWithMeasurementAndFitterInfo(id, gftrack->getCardinalRep());

    if (!trpoint)
    {
      if (Verbosity() > 1)
        LogWarning("!trpoint");
      continue;
    }

    genfit::KalmanFitterInfo* kfi = static_cast<genfit::KalmanFitterInfo*>(trpoint->getFitterInfo(rep));
    if (!kfi)
    {
      if (Verbosity() > 1)
        LogWarning("!kfi");
      continue;
    }

    const genfit::MeasuredStateOnPlane* gf_state = nullptr;
    try
    {
      // this works because KalmanFitterInfo returns a const reference to internal object and not a temporary object
      gf_state = &kfi->getFittedState(true);
    }
    catch (...)
    {
      if (Verbosity() > 1)
        LogWarning("Exrapolation failed!");
    }
    if (!gf_state)
    {
      if (Verbosity() > 1)
        LogWarning("Exrapolation failed!");
      continue;
    }
    genfit::MeasuredStateOnPlane temp;
    float pathlength = -phgf_track->extrapolateToPoint(temp, vertex_position, id);

    // create new svtx state and add to track
    auto state = create_track_state( pathlength, gf_state );
    out_track->insert_state( &state );

#ifdef _DEBUG_
    cout
        << __LINE__
        << ": " << id
        << ": " << pathlength << " => "
        << sqrt(state->get_x() * state->get_x() + state->get_y() * state->get_y())
        << endl;
#endif
  }

  // loop over clusters, check if layer is disabled, include extrapolated SvtxTrackState
  if( !_disabled_layers.empty() )
  {

    unsigned int id_min = 0;
    for (auto iter = svtx_track->begin_cluster_keys(); iter != svtx_track->end_cluster_keys(); ++iter)
    {

      auto cluster_key = *iter;
      auto cluster = _clustermap->findCluster(cluster_key);
      const int layer = TrkrDefs::getLayer(cluster_key);

      // skip enabled layers
      if( _disabled_layers.find( layer ) == _disabled_layers.end() )
      { continue; }

      // get position
      TVector3 pos(cluster->getPosition(0), cluster->getPosition(1), cluster->getPosition(2));
      float r_cluster = std::sqrt( square(pos[0]) + square(pos[1]) );

      // loop over states
      /* find first state whose radius is larger than that of cluster if any */
      unsigned int id = id_min;
      for( ; id < gftrack->getNumPointsWithMeasurement(); ++id)
      {

        auto trpoint = gftrack->getPointWithMeasurementAndFitterInfo(id, rep);
        if (!trpoint) continue;

        auto kfi = static_cast<genfit::KalmanFitterInfo*>(trpoint->getFitterInfo(rep));
        if (!kfi) continue;

        const genfit::MeasuredStateOnPlane* gf_state = nullptr;
        try
        {

          gf_state = &kfi->getFittedState(true);

        } catch (...) {

          if (Verbosity() > 1)
          { LogWarning("Failed to get kf fitted state"); }

        }

        if( !gf_state ) continue;

        float r_track = std::sqrt( square( gf_state->getPos().x() ) + square( gf_state->getPos().y() ) );
        if( r_track > r_cluster ) break;

      }

      // forward extrapolation
      genfit::MeasuredStateOnPlane gf_state;
      float pathlength = 0;

      // first point is previous, if valid
      if( id > 0 )  id_min = id-1;

      // extrapolate forward
      {
        auto trpoint = gftrack->getPointWithMeasurementAndFitterInfo(id_min, rep);
        auto kfi = static_cast<genfit::KalmanFitterInfo*>(trpoint->getFitterInfo(rep));
        gf_state = *kfi->getForwardUpdate();
        pathlength = gf_state.extrapolateToPoint( pos );
        auto tmp = *kfi->getBackwardUpdate();
        pathlength -= tmp.extrapolateToPoint( vertex_position );
      }

      // also extrapolate backward from next state if any
      // and take the weighted average between both points
      if( id > 0 && id < gftrack->getNumPointsWithMeasurement() )
      {
        auto trpoint = gftrack->getPointWithMeasurementAndFitterInfo(id, rep);
        auto kfi = static_cast<genfit::KalmanFitterInfo*>(trpoint->getFitterInfo(rep));
        genfit::KalmanFittedStateOnPlane gf_state_backward = *kfi->getBackwardUpdate();
        gf_state_backward.extrapolateToPlane( gf_state.getPlane() );
        gf_state = genfit::calcAverageState( gf_state, gf_state_backward );
      }

      // create new svtx state and add to track
      auto state = create_track_state(  pathlength, &gf_state );
      out_track->insert_state( &state );

    }

  }

  return out_track;
}

/*
 * Fill SvtxVertexMap from GFRaveVertexes and Tracks
 */
bool PHGenFitTrkFitter::FillSvtxVertexMap(
    const std::vector<genfit::GFRaveVertex*>& rave_vertices,
    const std::vector<genfit::Track*>& gf_tracks)
{

  if(Verbosity() > 0) cout << "Rave vertices size " << rave_vertices.size() << endl;
  if(rave_vertices.size() > 0)
    {
      for (unsigned int ivtx = 0; ivtx < rave_vertices.size(); ++ivtx)
  {
    genfit::GFRaveVertex* rave_vtx = rave_vertices[ivtx];

    if (!rave_vtx)
      {
        cerr << PHWHERE << endl;
        return false;
      }

    if(Verbosity() > 0) cout << "   ivtx " << ivtx << " has  Z = " << rave_vtx->getPos().Z() << endl;

    SvtxVertex_v1 svtx_vtx;
    svtx_vtx.set_chisq(rave_vtx->getChi2());
    svtx_vtx.set_ndof(rave_vtx->getNdf());
    svtx_vtx.set_position(0, rave_vtx->getPos().X());
    svtx_vtx.set_position(1, rave_vtx->getPos().Y());
    svtx_vtx.set_position(2, rave_vtx->getPos().Z());

    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        svtx_vtx.set_error(i, j, rave_vtx->getCov()[i][j]);

    for (unsigned int i = 0; i < rave_vtx->getNTracks(); i++)
      {
        //TODO Assume id's are sync'ed between _trackmap_refit and gf_tracks, need to change?
        const genfit::Track* rave_track =
    rave_vtx->getParameters(i)->getTrack();
        for (unsigned int j = 0; j < gf_tracks.size(); j++)
    {
      if (rave_track == gf_tracks[j])
        {
          svtx_vtx.insert_track(j);
          _rave_vertex_gf_track_map.insert(std::pair<unsigned int, unsigned int>(j, ivtx));
          if(Verbosity() > 0) cout << " rave vertex " << ivtx << " at Z " << svtx_vtx.get_position(2) << " rave track " << i  << " genfit track ID " << j << endl;
        }
    }
      }

    if (_vertexmap_refit)
      {
        if(Verbosity() > 0) cout << "insert svtx_vtx into _vertexmap_refit " << endl;
        _vertexmap_refit->insert_clone( &svtx_vtx );
        if(Verbosity() > 10) _vertexmap_refit->identify();
      }
    else
      {
        LogError("!_vertexmap_refit");
      }
  }  //loop over RAVE vertices
    }

  return true;
}

//bool PHGenFitTrkFitter::pos_cov_uvn_to_rz(const TVector3 u, const TVector3 v,
//		const TVector3 n, const TMatrixF pos_in, const TMatrixF cov_in,
//		TMatrixF& pos_out, TMatrixF& cov_out) const {
//
//	if(pos_in.GetNcols() != 1 || pos_in.GetNrows() != 3) {
//		if(Verbosity() > 0) LogWarning("pos_in.GetNcols() != 1 || pos_in.GetNrows() != 3");
//		return false;
//	}
//
//	if(cov_in.GetNcols() != 3 || cov_in.GetNrows() != 3) {
//		if(Verbosity() > 0) LogWarning("cov_in.GetNcols() != 3 || cov_in.GetNrows() != 3");
//		return false;
//	}
//
//	TVector3 up = TVector3(0., 0., 1.).Cross(n);
//	if(up.Mag() < 0.00001){
//		if(Verbosity() > 0) LogWarning("n is parallel to z");
//		return false;
//	}
//
//	TMatrixF R(3, 3);
//	TMatrixF R_inv(3,3);
//	TMatrixF R_inv_T(3,3);
//
//	try {
//		TMatrixF ROT1(3, 3);
//		TMatrixF ROT2(3, 3);
//		TMatrixF ROT3(3, 3);
//
//		// rotate n along z to xz plane
//		float phi = -TMath::ATan2(n.Y(), n.X());
//		ROT1[0][0] = cos(phi);
//		ROT1[0][1] = -sin(phi);
//		ROT1[0][2] = 0;
//		ROT1[1][0] = sin(phi);
//		ROT1[1][1] = cos(phi);
//		ROT1[1][2] = 0;
//		ROT1[2][0] = 0;
//		ROT1[2][1] = 0;
//		ROT1[2][2] = 1;
//
//		// rotate n along y to z
//		TVector3 n1(n);
//		n1.RotateZ(phi);
//		float theta = -TMath::ATan2(n1.X(), n1.Z());
//		ROT2[0][0] = cos(theta);
//		ROT2[0][1] = 0;
//		ROT2[0][2] = sin(theta);
//		ROT2[1][0] = 0;
//		ROT2[1][1] = 1;
//		ROT2[1][2] = 0;
//		ROT2[2][0] = -sin(theta);
//		ROT2[2][1] = 0;
//		ROT2[2][2] = cos(theta);
//
//		// rotate u along z to x
//		TVector3 u2(u);
//		u2.RotateZ(phi);
//		u2.RotateY(theta);
//		float phip = -TMath::ATan2(u2.Y(), u2.X());
//		phip -= -TMath::ATan2(up.Y(), up.X());
//		ROT3[0][0] = cos(phip);
//		ROT3[0][1] = -sin(phip);
//		ROT3[0][2] = 0;
//		ROT3[1][0] = sin(phip);
//		ROT3[1][1] = cos(phip);
//		ROT3[1][2] = 0;
//		ROT3[2][0] = 0;
//		ROT3[2][1] = 0;
//		ROT3[2][2] = 1;
//
//		// R: rotation from u,v,n to (z X n), v', z
//		R = ROT3 * ROT2 * ROT1;
//		R_inv = R.Invert();
//		R_inv_T.Transpose(R_inv);
//
//	} catch (...) {
//		if (Verbosity() > 0)
//			LogWarning("Can't get rotation matrix");
//
//		return false;
//	}
//
//	pos_out.ResizeTo(3, 1);
//	cov_out.ResizeTo(3, 3);
//
//	pos_out = R_inv * pos_in;
//	cov_out = R_inv * cov_in * R_inv_T;
//
//	return true;
//}

bool PHGenFitTrkFitter::pos_cov_uvn_to_rz(const TVector3& u, const TVector3& v,
                                          const TVector3& n, const TMatrixF& pos_in, const TMatrixF& cov_in,
                                          TMatrixF& pos_out, TMatrixF& cov_out) const
{
  if (pos_in.GetNcols() != 1 || pos_in.GetNrows() != 3)
  {
    if (Verbosity() > 0) LogWarning("pos_in.GetNcols() != 1 || pos_in.GetNrows() != 3");
    return false;
  }

  if (cov_in.GetNcols() != 3 || cov_in.GetNrows() != 3)
  {
    if (Verbosity() > 0) LogWarning("cov_in.GetNcols() != 3 || cov_in.GetNrows() != 3");
    return false;
  }

  TVector3 Z_uvn(u.Z(), v.Z(), n.Z());
  TVector3 up_uvn = TVector3(0., 0., 1.).Cross(Z_uvn);  // n_uvn X Z_uvn

  if (up_uvn.Mag() < 0.00001)
  {
    if (Verbosity() > 0) LogWarning("n is parallel to z");
    return false;
  }

  // R: rotation from u,v,n to n X Z, nX(nXZ), n
  TMatrixF R(3, 3);
  TMatrixF R_T(3, 3);

  try
  {
    // rotate u along z to up
    float phi = -TMath::ATan2(up_uvn.Y(), up_uvn.X());
    R[0][0] = cos(phi);
    R[0][1] = -sin(phi);
    R[0][2] = 0;
    R[1][0] = sin(phi);
    R[1][1] = cos(phi);
    R[1][2] = 0;
    R[2][0] = 0;
    R[2][1] = 0;
    R[2][2] = 1;

    R_T.Transpose(R);
  }
  catch (...)
  {
    if (Verbosity() > 0)
      LogWarning("Can't get rotation matrix");

    return false;
  }

  pos_out.ResizeTo(3, 1);
  cov_out.ResizeTo(3, 3);

  pos_out = R * pos_in;
  cov_out = R * cov_in * R_T;

  return true;
}

bool PHGenFitTrkFitter::get_vertex_error_uvn(const TVector3& u,
                                             const TVector3& v, const TVector3& n, const TMatrixF& cov_in,
                                             TMatrixF& cov_out) const
{
  /*!
   * Get matrix that rotates frame (u,v,n) to (x,y,z)
   * or the matrix that rotates vector defined in (x,y,z) to defined (u,v,n)
   */

  TMatrixF R = get_rotation_matrix(u, v, n);
  //
  //	LogDebug("PHGenFitTrkFitter::get_vertex_error_uvn::R = ");
  //	R.Print();
  //	cout<<"R.Determinant() = "<<R.Determinant()<<"\n";

  if (!(abs(R.Determinant() - 1) < 0.01))
  {
    if (Verbosity() > 0)
      LogWarning("!(abs(R.Determinant()-1)<0.0001)");
    return false;
  }

  if (R.GetNcols() != 3 || R.GetNrows() != 3)
  {
    if (Verbosity() > 0)
      LogWarning("R.GetNcols() != 3 || R.GetNrows() != 3");
    return false;
  }

  if (cov_in.GetNcols() != 3 || cov_in.GetNrows() != 3)
  {
    if (Verbosity() > 0)
      LogWarning("cov_in.GetNcols() != 3 || cov_in.GetNrows() != 3");
    return false;
  }

  TMatrixF R_T(3, 3);

  R_T.Transpose(R);

  cov_out.ResizeTo(3, 3);

  cov_out = R * cov_in * R_T;

  return true;
}

bool PHGenFitTrkFitter::pos_cov_XYZ_to_RZ(
    const TVector3& n, const TMatrixF& pos_in, const TMatrixF& cov_in,
    TMatrixF& pos_out, TMatrixF& cov_out) const
{
  if (pos_in.GetNcols() != 1 || pos_in.GetNrows() != 3)
  {
    if (Verbosity() > 0) LogWarning("pos_in.GetNcols() != 1 || pos_in.GetNrows() != 3");
    return false;
  }

  if (cov_in.GetNcols() != 3 || cov_in.GetNrows() != 3)
  {
    if (Verbosity() > 0) LogWarning("cov_in.GetNcols() != 3 || cov_in.GetNrows() != 3");
    return false;
  }

  // produces a vector perpendicular to both the momentum vector and beam line - i.e. in the direction of the dca_xy
  // only the angle of r will be used, not the magnitude
  TVector3 r = n.Cross(TVector3(0., 0., 1.));
  if (r.Mag() < 0.00001)
  {
    if (Verbosity() > 0) LogWarning("n is parallel to z");
    return false;
  }

  // R: rotation from u,v,n to n X Z, nX(nXZ), n
  TMatrixF R(3, 3);
  TMatrixF R_T(3, 3);

  try
  {
    // rotate u along z to up
    float phi = -TMath::ATan2(r.Y(), r.X());
    R[0][0] = cos(phi);
    R[0][1] = -sin(phi);
    R[0][2] = 0;
    R[1][0] = sin(phi);
    R[1][1] = cos(phi);
    R[1][2] = 0;
    R[2][0] = 0;
    R[2][1] = 0;
    R[2][2] = 1;

    R_T.Transpose(R);
  }
  catch (...)
  {
    if (Verbosity() > 0)
      LogWarning("Can't get rotation matrix");

    return false;
  }

  pos_out.ResizeTo(3, 1);
  cov_out.ResizeTo(3, 3);

  pos_out = R * pos_in;
  cov_out = R * cov_in * R_T;

  return true;
}

/*!
 * Get 3D Rotation Matrix that rotates frame (x,y,z) to (x',y',z')
 * Default rotate local to global, or rotate vector in global to local representation
 */
TMatrixF PHGenFitTrkFitter::get_rotation_matrix(const TVector3 x,
                                                const TVector3 y, const TVector3 z, const TVector3 xp, const TVector3 yp,
                                                const TVector3 zp) const
{
  TMatrixF R(3, 3);

  TVector3 xu = x.Unit();
  TVector3 yu = y.Unit();
  TVector3 zu = z.Unit();

  const float max_diff = 0.01;

  if (!(
          abs(xu * yu) < max_diff and
          abs(xu * zu) < max_diff and
          abs(yu * zu) < max_diff))
  {
    if (Verbosity() > 0)
      LogWarning("input frame error!");
    return R;
  }

  TVector3 xpu = xp.Unit();
  TVector3 ypu = yp.Unit();
  TVector3 zpu = zp.Unit();

  if (!(
          abs(xpu * ypu) < max_diff and
          abs(xpu * zpu) < max_diff and
          abs(ypu * zpu) < max_diff))
  {
    if (Verbosity() > 0)
      LogWarning("output frame error!");
    return R;
  }

  /*!
   * Decompose x',y',z' in x,y,z and call them u,v,n
   * Then the question will be rotate the standard X,Y,Z to u,v,n
   */

  TVector3 u(xpu.Dot(xu), xpu.Dot(yu), xpu.Dot(zu));
  TVector3 v(ypu.Dot(xu), ypu.Dot(yu), ypu.Dot(zu));
  TVector3 n(zpu.Dot(xu), zpu.Dot(yu), zpu.Dot(zu));

  try
  {
    std::shared_ptr<TRotation> rotation(new TRotation());
    //TRotation *rotation = new TRotation();

    //! Rotation that rotate standard (X, Y, Z) to (u, v, n)
    rotation->RotateAxes(u, v, n);

    R[0][0] = rotation->XX();
    R[0][1] = rotation->XY();
    R[0][2] = rotation->XZ();
    R[1][0] = rotation->YX();
    R[1][1] = rotation->YY();
    R[1][2] = rotation->YZ();
    R[2][0] = rotation->ZX();
    R[2][1] = rotation->ZY();
    R[2][2] = rotation->ZZ();
    //
    //		LogDebug("PHGenFitTrkFitter::get_rotation_matrix: TRotation:");
    //		R.Print();
    //		cout<<"R.Determinant() = "<<R.Determinant()<<"\n";

    //delete rotation;

    //		TMatrixF ROT1(3, 3);
    //		TMatrixF ROT2(3, 3);
    //		TMatrixF ROT3(3, 3);
    //
    //		// rotate n along z to xz plane
    //		float phi = -TMath::ATan2(n.Y(), n.X());
    //		ROT1[0][0] = cos(phi);
    //		ROT1[0][1] = -sin(phi);
    //		ROT1[0][2] = 0;
    //		ROT1[1][0] = sin(phi);
    //		ROT1[1][1] = cos(phi);
    //		ROT1[1][2] = 0;
    //		ROT1[2][0] = 0;
    //		ROT1[2][1] = 0;
    //		ROT1[2][2] = 1;
    //
    //		// rotate n along y to z
    //		TVector3 n1(n);
    //		n1.RotateZ(phi);
    //		float theta = -TMath::ATan2(n1.X(), n1.Z());
    //		ROT2[0][0] = cos(theta);
    //		ROT2[0][1] = 0;
    //		ROT2[0][2] = sin(theta);
    //		ROT2[1][0] = 0;
    //		ROT2[1][1] = 1;
    //		ROT2[1][2] = 0;
    //		ROT2[2][0] = -sin(theta);
    //		ROT2[2][1] = 0;
    //		ROT2[2][2] = cos(theta);
    //
    //		// rotate u along z to x
    //		TVector3 u2(u);
    //		u2.RotateZ(phi);
    //		u2.RotateY(theta);
    //		float phip = -TMath::ATan2(u2.Y(), u2.X());
    //		ROT3[0][0] = cos(phip);
    //		ROT3[0][1] = -sin(phip);
    //		ROT3[0][2] = 0;
    //		ROT3[1][0] = sin(phip);
    //		ROT3[1][1] = cos(phip);
    //		ROT3[1][2] = 0;
    //		ROT3[2][0] = 0;
    //		ROT3[2][1] = 0;
    //		ROT3[2][2] = 1;
    //
    //		// R: rotation from u,v,n to (z X n), v', z
    //		R = ROT3 * ROT2 * ROT1;
    //
    //		R.Invert();
    //		LogDebug("PHGenFitTrkFitter::get_rotation_matrix: Home Brew:");
    //		R.Print();
    //		cout<<"R.Determinant() = "<<R.Determinant()<<"\n";
  }
  catch (...)
  {
    if (Verbosity() > 0)
      LogWarning("Can't get rotation matrix");

    return R;
  }

  return R;
}
