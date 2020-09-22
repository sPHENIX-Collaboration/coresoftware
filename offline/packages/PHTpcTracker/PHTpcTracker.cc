/*!
 *  \file       PHTpcTracker.cc
 *  \brief      
 *  \author     Dmitry Arkhipkin <arkhipkin@gmail.com>
 */

#include "PHTpcTracker.h"

#include "PHTpcEventExporter.h"
#include "PHTpcLookup.h"
#include "PHTpcSeedFinder.h"
#include "PHTpcTrackFollower.h"
#include "PHTpcVertexFinder.h"

#include "Fitter.h"
#include "Track.h"  // for Track
#include "externals/kdfinder.hpp"

#include <phfield/PHField.h>
#include <phfield/PHFieldConfigv2.h>
#include <phfield/PHFieldUtility.h>

#include <phgeom/PHGeomUtility.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>  // for cluskey

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack_v1.h>

#include <trackreco/PHTrackSeeding.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHLog.h>
#include <phool/getClass.h>

#include <log4cpp/CategoryStream.hh>  // for CategoryStream

#include <CLHEP/Units/SystemOfUnits.h>

#include <TMatrixDSymfwd.h>  // for TMatrixDSym
#include <TMatrixTSym.h>     // for TMatrixTSym
#include <TMatrixTUtils.h>   // for TMatrixTRow
#include <TVector3.h>        // for TVector3
#include <TVectorDfwd.h>     // for TVectorD
#include <TVectorT.h>        // for TVectorT

#include <iostream>  // for operator<<, endl, basic...
#include <limits>    // for numeric_limits
#include <memory>    // for shared_ptr, __shared_ptr
#include <vector>    // for vector

class PHCompositeNode;

using namespace std;

PHTpcTracker::PHTpcTracker(const std::string& name)
  : PHTrackSeeding(name)
  ,
  //		mSeedFinder(nullptr), mTrackFollower(nullptr),
  //		mVertexFinder(nullptr), mEventExporter(nullptr), mLookup(nullptr),
  mFitter(nullptr)
  , mTGeoManager(nullptr)
  , mField(nullptr)
  , mB(1.4)
  , mEnableVertexing(false)
  , mEnableJsonExport(false)
{
  //	if ( !mSeedFinder ) {
  mSeedFinder = new PHTpcSeedFinder();
  //	}
  //	if ( !mTrackFollower ) {
  mTrackFollower = new PHTpcTrackFollower();
  //	}
  //	if ( !mVertexFinder ) {
  mVertexFinder = new PHTpcVertexFinder();
  //	}
  //	if ( !mEventExporter ) {
  mEventExporter = new PHTpcEventExporter();
  //	}
  //	if ( !mLookup ) {
  mLookup = new PHTpcLookup();
  //	}
}

PHTpcTracker::~PHTpcTracker()
{
  delete mSeedFinder;
  delete mTrackFollower;
  delete mVertexFinder;
  delete mEventExporter;
  delete mLookup;
  delete mFitter;
  delete mField;
}

int PHTpcTracker::Setup(PHCompositeNode* topNode)
{
  int ret = PHTrackSeeding::Setup(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTpcTracker::Process(PHCompositeNode* topNode)
{
  LOG_INFO("tracking.PHTpcTracker.process_event") << "---- process event started -----";

  // ----- magnetic field -----
  if (!mField)
  {
    mField = getMagField(topNode, mB);
  }

  // ----- Setup Geometry and Fitter -----
  // FIXME: get TGeoManager only once per file?
  if (!mTGeoManager)
  {
    mTGeoManager = PHGeomUtility::GetTGeoManager(topNode);
    LOG_ERROR_IF("tracking.PHTpcTracker.process_event", !mTGeoManager) << "Cannot find TGeoManager, track propagation will fail";
  }

  if (!mFitter && mField && mTGeoManager)
  {
    mFitter = new PHGenFit2::Fitter(mTGeoManager, mField);
  }

  // ----- Seed finding -----
  TrkrClusterContainer* cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  std::vector<kdfinder::TrackCandidate<double>*> candidates;
  candidates = mSeedFinder->findSeeds(cluster_map, mB);
  LOG_INFO("tracking.PHTpcTracker.process_event") << "SeedFinder produced " << candidates.size() << " track seeds";

  // ----- Track Following -----
  mLookup->init(cluster_map);
  std::vector<PHGenFit2::Track*> gtracks;
  gtracks = mTrackFollower->followTracks(cluster_map, candidates, mField, mLookup, mFitter);
  LOG_INFO("tracking.PHTpcTracker.process_event") << "TrackFollower reconstructed " << gtracks.size() << " tracks";
  // write tracks to fun4all server
  //  cout<<  gtracks[1]->get_vertex_id() << endl;
  for (int i = 0, ilen = gtracks.size(); i < ilen; i++)
  {
    //  for (auto it = gtracks.begin(); it != gtracks.end(); ++it)
    std::shared_ptr<SvtxTrack_v1> svtx_track(new SvtxTrack_v1());
    ////// from here:

    svtx_track->Reset();
    svtx_track->set_id(1);
    // cout << gtracks[i]->get_vertex_id() << endl;
    TVectorD state = gtracks[i]->getGenFitTrack()->getStateSeed();
    TVector3 pos(state(0), state(1), state(2));
    TVector3 mom(state(3), state(4), state(5));
    TMatrixDSym cov = gtracks[i]->getGenFitTrack()->getCovSeed();
    //cout<< "pt: " << pos.Perp() << endl;

    //double charge =  gtracks[i]->get_charge();
    double charge =  -gtracks[i]->get_charge();   // kludge to get Acts tracking to work - disagreement about magfield sign?

    svtx_track->set_charge(charge);

    for (int k = 0; k < 6; k++)
    {
      for (int j = 0; j < 6; j++)
      {
        svtx_track->set_error(k, j, cov[k][j]);
      }
    }

    svtx_track->set_px(mom.Px());
    svtx_track->set_py(mom.Py());
    svtx_track->set_pz(mom.Pz());

    svtx_track->set_x(pos.X());
    svtx_track->set_y(pos.Y());
    svtx_track->set_z(pos.Z());

    for (TrkrDefs::cluskey cluster_key : gtracks[i]->get_cluster_keys())
    {
      svtx_track->insert_cluster_key(cluster_key);
    }

    _track_map->insert(svtx_track.get());
  }

  // ----- cleanup -----

  // candidates cleanup
  for (auto it = candidates.begin(); it != candidates.end(); ++it)
  {
    delete (*it);
  }

  // genfit tracks cleanup
  for (auto it = gtracks.begin(); it != gtracks.end(); ++it)
  {
    delete (*it);
  }

  // hit lookup cleanup
  mLookup->clear();

  LOG_INFO("tracking.PHTpcTracker.process_event") << "---- process event finished -----";

  return Fun4AllReturnCodes::EVENT_OK;
}

PHField* PHTpcTracker::getMagField(PHCompositeNode* topNode, double& B)
{
  PHField* field = nullptr;
  PHFieldConfigv2 bconfig(0, 0, B);
  // ----- check file -----
  PHField* field_file = PHFieldUtility::GetFieldMapNode(nullptr, topNode);
  if (field_file)
  {
    const double point[] = {0, 0, 0, 0};  // x,y,z,t
    double Bx = 0, By = 0, Bz = 0;
    double Bfield[] = {std::numeric_limits<double>::signaling_NaN(),
                       std::numeric_limits<double>::signaling_NaN(),
                       std::numeric_limits<double>::signaling_NaN(),
                       std::numeric_limits<double>::signaling_NaN(),
                       std::numeric_limits<double>::signaling_NaN(),
                       std::numeric_limits<double>::signaling_NaN()};
    field_file->GetFieldValue(point, Bfield);
    Bx = Bfield[0];
    By = Bfield[1];
    Bz = Bfield[2];
    B = Bz / CLHEP::tesla;
    LOG_DEBUG("tracking.PHTpcTracker.process_event") << "Importing B field from file, Bx,By,Bz Tesla = " << Bx / CLHEP::tesla << "," << By / CLHEP::tesla << "," << Bz / CLHEP::tesla;
    bconfig.set_field_mag_z(B);
  }
  else
  {
    LOG_WARN("tracking.PHTpcTracker.process_event") << "No field found in file, using default Bz value = " << B << " Tesla";
  }
  field = PHFieldUtility::BuildFieldMap(&bconfig, 1);
  return field;
}

void PHTpcTracker::set_seed_finder_options(double maxdistance1, double tripletangle1, size_t minhits1,
                                           double maxdistance2, double tripletangle2, size_t minhits2, size_t nthreads)
{
  mSeedFinder->set_options(maxdistance1, tripletangle1, minhits1,
                           maxdistance2, tripletangle2, minhits2, nthreads);
}

void PHTpcTracker::set_seed_finder_optimization_remove_loopers(bool opt, double minr, double maxr)
{
  mSeedFinder->set_optimization_remove_loopers(opt, minr, maxr);
}

void PHTpcTracker::set_track_follower_optimization_helix(bool opt)
{
  mTrackFollower->set_optimization_helix(opt);
}

void PHTpcTracker::set_track_follower_optimization_precise_fit(bool opt)
{
  mTrackFollower->set_optimization_precise_fit(opt);
}
