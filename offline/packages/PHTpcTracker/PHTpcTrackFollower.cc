
/*!
 *  \file       PHTpcTrackFollower.cc
 *  \brief      
 *  \author     Dmitry Arkhipkin <arkhipkin@gmail.com>
 */

#include "PHTpcTrackFollower.h"
#include "PHTpcConst.h"

#include "Fitter.h"       // for Fitter
#include "PHTpcLookup.h"  // for PHTpcLookup
#include "SpacepointMeasurement2.h"
#include "Track.h"
#include "externals/kdfinder.hpp"  // for TrackCandidate, Helix, TVector

#include <trackbase/TrkrDefs.h>  // for cluskey

#include <phool/PHLog.h>

#include <GenFit/FieldManager.h>
#include <GenFit/KalmanFitStatus.h>
#include <GenFit/KalmanFitterInfo.h>
#include <GenFit/MeasuredStateOnPlane.h>
#include <GenFit/RKTrackRep.h>
#include <GenFit/Track.h>
#include <GenFit/TrackPoint.h>

#include <TMatrixDSymfwd.h>  // for TMatrixDSym
#include <TMatrixTSym.h>     // for TMatrixTSym
#include <TVector3.h>        // for TVector3

#include <log4cpp/CategoryStream.hh>  // for CategoryStream

#include <algorithm>
#include <cmath>    // for isnan, sqrt
#include <cstddef>  // for size_t
#include <cstdint>  // for uint64_t, int64_t
#include <limits>   // for numeric_limits
#include <unordered_set>

class PHField;
class TrkrClusterContainer;

namespace PHGenFit
{
  class Measurement;
}
namespace genfit
{
  class AbsTrackRep;
}

PHTpcTrackFollower::PHTpcTrackFollower()
  : mOptHelix(true)
  , mOptPreciseFit(true)
{
}

std::vector<PHGenFit2::Track*> PHTpcTrackFollower::followTracks(std::vector<kdfinder::TrackCandidate<double>*>& candidates, PHField* B, PHTpcLookup* lookup, PHGenFit2::Fitter* fitter)
{
  LOG_DEBUG("tracking.PHTpcTrackFollower.followTracks") << "start";

  // ----- sort track seeds by n hits desc -----
  std::sort(candidates.begin(), candidates.end(), [](kdfinder::TrackCandidate<double>* a, kdfinder::TrackCandidate<double>* b) {
    return a->nhits() > b->nhits();
  });

  // ----- process seeds in a loop -----
  std::vector<PHGenFit2::Track*> gtracks;

  std::unordered_set<uint64_t> usedClusterKeys;
  int reused_hits = 0;

  const int NCANDIDATES = candidates.size();
  for (int i = 0; i < NCANDIDATES; i++)
  {
    LOG_DEBUG("tracking.PHTpcTrackFollower.followTracks") << "+++++ following track " << i << " +++++";

    // ---- check if candidate hits belong to some already reco-ed track -----
    const int NHITS = candidates[i]->nhits();
    reused_hits = 0;
    for (int j = 0; j < NHITS; j++)
    {
      std::vector<double>& hit = candidates[i]->getHit(j);
      auto it = usedClusterKeys.find(*((int64_t*) &hit[3]));
      if (it != usedClusterKeys.end())
      {
        ++reused_hits;
      }
      if (reused_hits > 1)
      {
        break;
      }
    }
    // skip seed if it has at least 2 previously used hits
    if (reused_hits > 1)
    {
      LOG_DEBUG("tracking.PHTpcTrackFollower.followTracks") << "more than one hit already used, skipping seed";
      continue;
    }

    // ----- seed looks ok, proceed -----
    PHGenFit2::Track* gtrack = nullptr;
    try
    {
      gtrack = propagateTrack(candidates[i], lookup, fitter);
    }
    catch (...)
    {
      LOG_DEBUG("tracking.PHTpcTrackFollower.followTracks") << "exception caught during propagation, skipping track";
      continue;
    }
    LOG_DEBUG("tracking.PHTpcTrackFollower.followTracks") << "track propagated";

    if (!gtrack || gtrack->getGenFitTrack()->getNumPoints() < 10)
    {
      LOG_DEBUG("tracking.PHTpcTrackFollower.followTracks") << "no track or track has less than 10 hit points, skipping";
      delete gtrack;
      continue;
    }

    genfit::KalmanFitStatus* fstatus = gtrack->getGenFitTrack()->getKalmanFitStatus();
    double chi2 = fstatus->getChi2(),
           ndf = fstatus->getNdf();
    int nfailedpts = fstatus->getNFailedPoints();

    LOG_DEBUG("tracking.PHTpcTrackFollower.propagateTrack") << "chi2: " << chi2 << ", ndf: " << ndf << ", chi2/ndf: " << chi2 / ndf
                                                            << ", failed: " << nfailedpts << ", pts: " << gtrack->getGenFitTrack()->getNumPoints();

    if (nfailedpts > 0)
    {
      LOG_DEBUG("tracking.PHTpcTrackFollower.followTracks") << "track has " << nfailedpts << " failed points, skipping";
      delete gtrack;
      continue;
    }

    gtracks.push_back(gtrack);
    LOG_DEBUG("tracking.PHTpcTrackFollower.followTracks") << "track added";
    // ----- mark hits as used in track reco -----
    const std::vector<TrkrDefs::cluskey>& cluskeys = gtrack->get_cluster_keys();
    for (TrkrDefs::cluskey cluskey : cluskeys)
    {
      usedClusterKeys.insert(cluskey);
    }
  }

  return gtracks;
}

PHGenFit2::Track* PHTpcTrackFollower::propagateTrack(kdfinder::TrackCandidate<double>* candidate,
                                                     PHTpcLookup* lookup, PHGenFit2::Fitter* fitter)
{
  //----- candidate hits are pre-sorted, inner R -> outer R -----
  std::vector<std::vector<double>>& hits = candidate->getHits();

  LOG_DEBUG("tracking.PHTpcTrackFollower.propagateTrack") << "ctrack nhits: " << hits.size();

  // TEST: check pos, mom of candidate vs pos, mom of gtrack
  LOG_DEBUG("tracking.PHTpcTrackFollower.propagateTrack") << "ctrack mom: " << candidate->momentum();

  std::vector<double> cmom = candidate->getMomForHit(0);
  std::vector<double> cpos = candidate->getPosForHit(0);
  LOG_DEBUG("tracking.PHTpcTrackFollower.propagateTrack")
      << "ctrack pos: " << cpos[0] << ", " << cpos[1] << ", " << cpos[2] << " | "
      << "mom: " << cmom[0] << ", " << cmom[1] << ", " << cmom[2];

  // ----- make genfit track out of a track candidate -----
  PHGenFit2::Track* gtrack = candidate_to_genfit(candidate);
  try
  {
    gtrack->getGenFitTrack()->checkConsistency();
  }
  catch (...)
  {
    LOG_DEBUG("tracking.PHTpcTrackFollower.propagateTrack") << "Found inconsistent track, skipping";
    return nullptr;
  }
  // ----- make initial fit of the track -----
  try
  {
    if (mOptPreciseFit)
    {
      fitter->processTrack5(gtrack);
    }
    else
    {
      fitter->processTrack1(gtrack);
    }
  }
  catch (...)
  {
    // genfit exception thrown
    LOG_DEBUG("tracking.PHTpcTrackFollower.propagateTrack") << "GenFit exception on track fit, skipping";
    return nullptr;
  }

  if (!gtrack->getGenFitTrack()->hasFitStatus())
  {
    LOG_DEBUG("tracking.PHTpcTrackFollower.propagateTrack") << "Track has not been fit, skipping";
    return nullptr;
  }

  genfit::KalmanFitStatus* fstatus = gtrack->getGenFitTrack()->getKalmanFitStatus();
  double chi2 = fstatus->getChi2(), ndf = fstatus->getNdf();
  int nfailedpts = fstatus->getNFailedPoints();

  if (nfailedpts > 0)
  {
    LOG_DEBUG("tracking.PHTpcTrackFollower.propagateTrack") << "Track has failed points, skipping";
    return nullptr;
  }

  LOG_DEBUG("tracking.PHTpcTrackFollower.propagateTrack") << "SEED chi2: " << chi2 << ", ndf: " << ndf << ", chi2/ndf: " << chi2 / ndf
                                                          << ", failed: " << nfailedpts << ", pts: " << gtrack->getGenFitTrack()->getNumPoints();

  //----- check if track has max possible hits already -----
  if (hits.size() >= PHTpcConst::TPC_LAYERS_MAX)
  {
    LOG_DEBUG("tracking.PHTpcTrackFollower.propagateTrack") << "track has max possible points from seed, no propagation needed";
    return gtrack;
  }

  //----- propagate track outwards using KDTree lookups -----
  LOG_DEBUG("tracking.PHTpcTrackFollower.propagateTrack") << "follow track outwards";
  int rc1 = 0;
  try
  {
    rc1 = followTrack(gtrack, lookup, fitter, 1);
  }
  catch (...)
  {
    // bad track?
    LOG_DEBUG("tracking.PHTpcTrackFollower.propagateTrack") << "Errors when following track outwards, skipping";
    return nullptr;
  }

  // ---- refit track again ----
  if (rc1 >= 2)
  {
    LOG_DEBUG("tracking.PHTpcTrackFollower.propagateTrack") << "new hits added during outwards following, refitting track";
    try
    {
      // FIXME: 1-pass imprecise fit or partial fit to keep inwards prop precise?
      fitter->processTrack1(gtrack);
    }
    catch (...)
    {
    }
  }

  //----- propagate track inwards using KDTree lookups ------
  LOG_DEBUG("tracking.PHTpcTrackFollower.propagateTrack") << "follow track inwards";
  int rc2 = 0;
  try
  {
    rc2 = followTrack(gtrack, lookup, fitter, -1);
  }
  catch (...)
  {
    LOG_DEBUG("tracking.PHTpcTrackFollower.propagateTrack") << "Errors when following track inwards";
  }

  //----- refit completed track -----
  if (rc1 || rc2)
  {
    LOG_DEBUG("tracking.PHTpcTrackFollower.propagateTrack") << "new hits added during inwards following, refitting track";
    try
    {
      if (mOptPreciseFit)
      {
        fitter->processTrack5(gtrack);
      }
      else
      {
        fitter->processTrack1(gtrack);
      }
    }
    catch (...)
    {
    }
  }

  try
  {
    gtrack->getGenFitTrack()->checkConsistency();
  }
  catch (...)
  {
    LOG_DEBUG("tracking.PHTpcTrackFollower.propagateTrack") << "got inconsistent track, skipping";
    return nullptr;
  }

  return gtrack;
}

PHGenFit2::Track* PHTpcTrackFollower::candidate_to_genfit(kdfinder::TrackCandidate<double>* candidate)
{
  // ----- set RungeKutta track rep with the most common particle hypothesis pi- -----
  genfit::AbsTrackRep* rep = new genfit::RKTrackRep(candidate->sign() > 0 ? 211 : -211 /* Pion PDG ID, as the most common */);

  // ----- pos, mom from the innermost hit of the track candidate -----
  std::vector<double> cmom = candidate->getMomForHit(0);
  std::vector<double> cpos = candidate->getPosForHit(0);
  TVector3 seed_pos(cpos[0], cpos[1], cpos[2]);  // x, y, z
  TVector3 seed_mom(cmom[0], cmom[1], cmom[2]);  // px, py, pz

  // ----- set up cov matrix -----
  TMatrixDSym seed_cov(6);
  seed_cov(0, 0) = 0.1 * 0.1;                      // dx * dx
  seed_cov(1, 1) = 0.1 * 0.1;                      // dy * dy
  seed_cov(2, 2) = 0.2 * 0.2;                      // dz * dz
  seed_cov(3, 3) = cmom[0] * 0.3 * cmom[0] * 0.3;  // dpx * dpx
  seed_cov(4, 4) = cmom[1] * 0.3 * cmom[1] * 0.3;  // dpy * dpy
  seed_cov(5, 5) = cmom[2] * 0.5 * cmom[2] * 0.5;  // dpz * dpz

  PHGenFit2::Track* track = new PHGenFit2::Track(rep, seed_pos, seed_mom, seed_cov);

  // ----- add track candidate hits as PHGenFit measurements -----
  const int NHITS = candidate->nhits();
  std::vector<PHGenFit::Measurement*> measurements;
  measurements.reserve(NHITS);
  for (int i = 0; i < NHITS; i++)
  {
    std::vector<double>& hit = candidate->getHit(i);
    PHGenFit::SpacepointMeasurement2* meas = hit_to_measurement(hit);
    measurements.push_back(meas);
  }
  track->addMeasurements(measurements);

  return track;
}

PHGenFit::SpacepointMeasurement2* PHTpcTrackFollower::hit_to_measurement(std::vector<double>& hit)
{
  TVector3 pos(hit[0], hit[1], hit[2]);
  TVector3 res(0.1, 0.1, 0.2);
  PHGenFit::SpacepointMeasurement2* meas = new PHGenFit::SpacepointMeasurement2(pos, res);
  meas->set_cluster_key(*((int64_t*) &hit[3]));
  return meas;
}

int PHTpcTrackFollower::get_track_layer(PHGenFit2::Track* track, int dir)
{
  return get_track_layer(track->getGenFitTrack(), dir);
}

int PHTpcTrackFollower::get_track_layer(genfit::Track* gftrack, int dir)
{
  int layer = -1;  // negative = no layer assigned yet
  unsigned int GFNUMPOINTS = gftrack->getNumPoints();
  if (GFNUMPOINTS <= 0)
  {
    LOG_ERROR("tracking.PHTpcTrackFollower.get_track_layer") << "track has no points";
    return layer;
  }
  if (!gftrack->hasKalmanFitStatus())
  {
    LOG_ERROR("tracking.PHTpcTrackFollower.get_track_layer") << "track has not been fitted";
    return layer;
  }
  try
  {
    LOG_DEBUG("tracking.PHTpcTrackFollower.get_track_layer") << "starting layer calc";
    genfit::TrackPoint* pt = gftrack->getPointWithMeasurementAndFitterInfo(dir == 1 ? GFNUMPOINTS - 1 : 0);
    if (!pt)
    {
      LOG_DEBUG("tracking.PHTpcTrackFollower.get_track_layer") << "track does not have a point with measurement and fitter info";
      return -1;
    }
    LOG_DEBUG("tracking.PHTpcTrackFollower.get_track_layer") << "got point";
    genfit::KalmanFitterInfo* kfinfo = pt->getKalmanFitterInfo();
    if (!kfinfo)
    {
      LOG_DEBUG("tracking.PHTpcTrackFollower.get_track_layer") << "track point does not have kalman fitter info";
      return -1;
    }
    LOG_DEBUG("tracking.PHTpcTrackFollower.get_track_layer") << "got kfinfo";
    const genfit::MeasuredStateOnPlane& kfstate = kfinfo->getFittedState();
    LOG_DEBUG("tracking.PHTpcTrackFollower.get_track_layer") << "got kfstate";
    TVector3 pos = kfstate.getPos();
    TVector3 mom = kfstate.getMom();
    double ptradius = pos.Perp();
    for (int i = 0; i < PHTpcConst::TPC_LAYERS_MAX; i++)
    {
      if (ptradius > (PHTpcConst::TPC_LAYERS_RADIUS[i] - PHTpcConst::TPC_LAYERS_DELTAR[i] / 2.0) &&
          ptradius <= (PHTpcConst::TPC_LAYERS_RADIUS[i] + PHTpcConst::TPC_LAYERS_DELTAR[i] / 2.0))
      {
        layer = i + dir;
        break;
      }
    }
  }
  catch (...)
  {
    LOG_DEBUG("tracking.PHTpcTrackFollower.get_track_layer") << "error during layer calculation";
    return -1;
  }
  return layer;
}

std::pair<genfit::MeasuredStateOnPlane*, double> PHTpcTrackFollower::get_projected_coordinate(genfit::Track* gftrack, int dir, double radius)
{
  // project from last fitted point to radius
  // TODO: for constant field is faster to use helix crossing cylinder - analytical
  // kdfinder::Helix helix( const TVector<T>& p, const TVector<T>& o, T B, T q )
  double pathlength = 0;
  genfit::MeasuredStateOnPlane* state = 0;
  try
  {
    genfit::TrackPoint* pt = gftrack->getPointWithMeasurementAndFitterInfo(dir == 1 ? -1 : 0);
    if (!pt)
    {
      LOG_DEBUG("tracking.PHTpcTrackFollower.get_projected_coordinate") << "no trackpoint with measurement and fitter info";
      return std::make_pair(nullptr, 0);
    }
    genfit::KalmanFitterInfo* kfinfo = pt->getKalmanFitterInfo();
    if (!kfinfo)
    {
      LOG_DEBUG("tracking.PHTpcTrackFollower.get_projected_coordinate") << "trackpoint does not have Kalman Fitter info";
      return std::make_pair(nullptr, 0);
    }
    const genfit::MeasuredStateOnPlane& kfstate = kfinfo->getFittedState();

    if (mOptHelix)
    {
      // make const-field helix projection first to save time if loopers present
      // FIXME: check if tracks are still correct ( i.e. only loopers are split, but not primaries! )
      TVector3 pos, mom;
      kfstate.getPosMom(pos, mom);
      double charge = kfstate.getCharge();
      double posX = pos.X(), posY = pos.Y(), posZ = pos.Z(), Bx, By, Bz;
      genfit::FieldManager::getInstance()->getFieldVal(posX, posY, posZ, Bx, By, Bz);
      kdfinder::Helix<double> h(kdfinder::TVector<double>(mom.X(), mom.Y(), mom.Z()),
                                kdfinder::TVector<double>(posX, posY, posZ), Bz, charge < 0 ? -1 : 1);
      std::pair<double, double> s = h.pathLength(radius);
      if ((!std::isnan(s.first) && s.first < 999.9) || (!std::isnan(s.second) && s.second < 999.9))
      {
        LOG_DEBUG("tracking.PHTpcTrackFollower.get_projected_coordinate") << "helix projects to cylinder: " << s.first << ", " << s.second;
        state = kfstate.clone();
        pathlength = state->extrapolateToCylinder(radius,
                                                  TVector3(0., 0., 0.),
                                                  TVector3(0., 0., 1.),
                                                  false,
                                                  false);
      }
      else
      {
        LOG_DEBUG("tracking.PHTpcTrackFollower.get_projected_coordinate") << "helix does not project to cylinder";
        // does not project to cylinder :(
        if (state)
        {
          delete state;
          state = 0;
        }
        return std::make_pair(nullptr, 0);
      }
    }
    else
    {  // mOptHelix = false

      // ALT: no-helix projection, painfully slow for loopers that change direction
      // GenFit does up to x1000 iterations to see if track projects on helix
      state = kfstate.clone();
      pathlength = state->extrapolateToCylinder(radius,
                                                TVector3(0., 0., 0.),
                                                TVector3(0., 0., 1.),
                                                false,
                                                false);
    }
  }
  catch (...)
  {
    // can't project to helix
    if (state)
    {
      delete state;
      state = 0;
    }
    return std::make_pair(nullptr, 0);
  }
  LOG_DEBUG("tracking.PHTpcTrackFollower.get_projected_coordinate") << "pathlength is: " << pathlength;
  return std::make_pair(state, pathlength);
}

std::pair<genfit::MeasuredStateOnPlane*, double> PHTpcTrackFollower::get_projected_coordinate(genfit::Track* gftrack, int dir, const TVector3& point)
{
  // project from last fitted point to radius
  double pathlength = 0;
  genfit::MeasuredStateOnPlane* state = nullptr;

  try
  {
    genfit::TrackPoint* pt = gftrack->getPointWithMeasurementAndFitterInfo(dir == 1 ? -1 : 0);
    if (!pt)
    {
      LOG_DEBUG("tracking.PHTpcTrackFollower.get_projected_coordinate") << "no trackpoint with measurement and fitter info";
      return std::make_pair(nullptr, 0);
    }
    genfit::KalmanFitterInfo* kfinfo = pt->getKalmanFitterInfo();
    if (!kfinfo)
    {
      LOG_DEBUG("tracking.PHTpcTrackFollower.get_projected_coordinate") << "trackpoint does not have Kalman Fitter info";
      return std::make_pair(nullptr, 0);
    }
    const genfit::MeasuredStateOnPlane& kfstate = kfinfo->getFittedState();

    state = kfstate.clone();
    pathlength = state->extrapolateToPoint(
        point,
        false,
        false);
  }
  catch (...)
  {
    // material effects exception suppression
  }

  LOG_DEBUG("tracking.PHTpcTrackFollower.get_projected_coordinate") << "pathlength is: " << pathlength;
  return std::make_pair(state, pathlength);
}

int PHTpcTrackFollower::followTrack(PHGenFit2::Track* track, PHTpcLookup* lookup, PHGenFit2::Fitter* fitter, int dir)
{
  if (!track)
  {
    LOG_ERROR("tracking.PHTpcTrackFollower.follow_track") << "cannot follow track, empty pointer!";
    return 0;
  }

  genfit::Track* gftrack = track->getGenFitTrack();

  int nHitsAdded = 0;

  LOG_DEBUG("tracking.PHTpcTrackFollower.follow_track") << "about to get layer info";
  int layer = -1;
  try
  {
    layer = get_track_layer(gftrack, dir);
  }
  catch (...)
  {
    LOG_DEBUG("tracking.PHTpcTrackFollower.follow_track") << "cannot get track layer, skipping propagation";
    return 0;
  }
  LOG_DEBUG("tracking.PHTpcTrackFollower.follow_track") << "got layer info: " << layer;

  if (!(layer >= 0 && layer < PHTpcConst::TPC_LAYERS_MAX))
  {
    LOG_DEBUG("tracking.PHTpcTrackFollower.followTrack") << "track goes beyond min/max layer, stop propagation. Layer: " << layer;
    return 0;
  }
  LOG_DEBUG("tracking.PHTpcTrackFollower.follow_track") << "layer info is good";

  while (1)
  {
    LOG_DEBUG("tracking.PHTpcTrackFollower.followTrack") << "dir: " << dir << ", propagating track to layer: " << layer << ", " << PHTpcConst::TPC_LAYERS_RADIUS[layer] << " cm";
    std::pair<genfit::MeasuredStateOnPlane*, double> p = get_projected_coordinate(gftrack, dir, PHTpcConst::TPC_LAYERS_RADIUS[layer]);
    if (!p.first)
    {
      break;
    }  // can't project to cylinder, likely curler track
    TVector3 pos = p.first->getPos();
    delete p.first;

    LOG_DEBUG("tracking.PHTpcTrackFollower.followTrack") << "projected position: " << pos.X() << ", " << pos.Y() << ", " << pos.Z() << ", radius: " << pos.Perp();

    size_t nMatches = 0;
    std::vector<std::vector<double>*> matches = lookup->find(pos.X(), pos.Y(), pos.Z(), PHTpcConst::TPC_LAYERS_DELTAR[layer] / 2.0, nMatches);
    LOG_DEBUG("tracking.PHTpcTrackFollower.followTrack") << "lookup returned " << nMatches << " hits nearby, size: " << matches.size();

    // options: no hits found, one hit found, several hits found
    if (nMatches > 0)
    {
      if (nMatches == 1)
      {
        LOG_DEBUG("tracking.PHTpcTrackFollower.followTrack") << "processing one hit";
        std::vector<double>* hit = matches[0];
        LOG_DEBUG("tracking.PHTpcTrackFollower.followTrack") << "hit x: " << (*hit)[0] << ", y: " << (*hit)[1] << ", z: " << (*hit)[2];
        // add hit as a measurement to the track, check chi2, remove hit if chi2 = bad
        std::pair<genfit::MeasuredStateOnPlane*, double> p2 = get_projected_coordinate(gftrack, dir, TVector3((*hit)[0], (*hit)[1], (*hit)[2]));
        if (!p2.first)
        {
          LOG_DEBUG("tracking.PHTpcTrackFollower.followTrack") << "cannot project to point, skipping";
          break;
        }
        TVector3 pos2 = p2.first->getPos();
        delete p2.first;
        LOG_DEBUG("tracking.PHTpcTrackFollower.followTrack") << "projected point: " << pos2.X() << ", " << pos2.Y() << ", " << pos2.Z() << ", radius: " << pos.Perp();
        LOG_DEBUG("tracking.PHTpcTrackFollower.followTrack") << "distance to hit: " << std::sqrt(std::pow(pos2.X() - (*hit)[0], 2) + std::pow(pos2.Y() - (*hit)[1], 2) + std::pow(pos2.Z() - (*hit)[2], 2));

        //FIXME: convert hit to measurement and add measurement to track
        PHGenFit::SpacepointMeasurement2* meas = hit_to_measurement(*hit);
        LOG_DEBUG_IF("tracking.PHTpcTrackFollower.followTrack", !meas) << "bad, no measurement";
        if (dir == 1)
        {
          track->addMeasurement(meas, -1);
          try
          {
            fitter->processTrackPartially(gftrack, -2, -1);
          }
          catch (...)
          {
          }
        }
        else
        {
          track->addMeasurement(meas, 0);
          try
          {
            fitter->processTrackPartially(gftrack, 1, 0);
          }
          catch (...)
          {
          }
        }
        ++nHitsAdded;
      }
      else
      {
        LOG_DEBUG("tracking.PHTpcTrackFollower.followTrack") << "processing multiple hits";
        // project track to point, sort matches by pathlength ASC, try hits one by one, choose best hit
        int ind = -1;
        double dist = std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < nMatches; i++)
        {
          std::vector<double>* hit = matches[i];
          LOG_DEBUG("tracking.PHTpcTrackFollower.followTrack") << "hit x: " << (*hit)[0] << ", y: " << (*hit)[1] << ", z: " << (*hit)[2];
          std::pair<genfit::MeasuredStateOnPlane*, double> p2 = get_projected_coordinate(gftrack, dir, TVector3((*hit)[0], (*hit)[1], (*hit)[2]));
          if (!p2.first)
          {
            continue;
          }
          TVector3 pos2 = p2.first->getPos();
          delete p2.first;
          double dist2 = std::sqrt(std::pow(pos2.X() - (*hit)[0], 2) + std::pow(pos2.Y() - (*hit)[1], 2) + std::pow(pos2.Z() - (*hit)[2], 2));
          LOG_DEBUG("tracking.PHTpcTrackFollower.followTrack") << "projected point: " << pos2.X() << ", " << pos2.Y() << ", " << pos2.Z() << ", radius: " << pos.Perp();
          LOG_DEBUG("tracking.PHTpcTrackFollower.followTrack") << "distance to hit: " << dist2;
          if (dist2 < dist)
          {
            ind = i;
            dist = dist2;
          }
        }
        if (ind != -1)
        {
          PHGenFit::SpacepointMeasurement2* meas = hit_to_measurement(*matches[ind]);
          if (dir == 1)
          {
            track->addMeasurement(meas, -1);
            try
            {
              fitter->processTrackPartially(gftrack, -2, -1);
            }
            catch (...)
            {
            }
          }
          else
          {
            track->addMeasurement(meas, 0);
            try
            {
              fitter->processTrackPartially(gftrack, 1, 0);
            }
            catch (...)
            {
            }
          }
          ++nHitsAdded;  // FIXME: only if hit is really added
        }
      }
    }  // else => no hits, scan next layer

    layer += dir;
    if (layer < 0 || layer >= PHTpcConst::TPC_LAYERS_MAX)
    {
      LOG_DEBUG("tracking.PHTpcTrackFollower.followTrack") << "done with " << dir << " following, layer: " << layer;
      break;
    }
  }

  return nHitsAdded;  // number of hits added to the track
}
