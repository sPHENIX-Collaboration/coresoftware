#include "ALICEKF.h"

#include "GPUTPCTrackLinearisation.h"
#include "GPUTPCTrackParam.h"

#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase_historic/ActsTransformations.h>

#include <Geant4/G4SystemOfUnits.hh>

#include <TMatrixFfwd.h>
#include <TMatrixT.h>
#include <TMatrixTUtils.h>

#include <omp.h>

//#define _DEBUG_

#if defined(_DEBUG_)
#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp
#else
#define LogDebug(exp) (void) 0
#endif

#define LogError(exp) \
  if (Verbosity() > 0) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogWarning(exp) \
  if (Verbosity() > 0) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp

using keylist = std::vector<TrkrDefs::cluskey>;

// anonymous namespace for local functions
namespace
{
  // square
  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }
}  // namespace

bool ALICEKF::checknan(double val, const std::string& name, int num) const
{
  if (std::isnan(val))
  {
    if (Verbosity() > 0)
    {
      std::cout << "WARNING: " << name << " is NaN for seed " << num << ". Aborting this seed.\n";
    }
  }
  return std::isnan(val);
}

double ALICEKF::get_Bz(double x, double y, double z) const
{
  //  if(_use_const_field) return _const_field;
  if (_use_const_field || fabs(z) > 105.5)
  {
    return _const_field;
  }
  double p[4] = {x * cm, y * cm, z * cm, 0. * cm};
  double bfield[3];

  // check thread number. Use uncached field accessor for all but thread 0.
  if( omp_get_thread_num() == 0 )
  {
    _B->GetFieldValue(p, bfield);
  } else {
    _B->GetFieldValue_nocache(p, bfield);
  }

  return bfield[2] / tesla;
}

double ALICEKF::getClusterError(TrkrCluster* c, TrkrDefs::cluskey key, Acts::Vector3 global, int i, int j) const
{
  if (_use_fixed_clus_error)
  {
    if (i == j)
    {
      return _fixed_clus_error.at(i) * _fixed_clus_error.at(i);
    }
    else
    {
      return 0.;
    }
  }
  else
  {
    TMatrixF localErr(3, 3);

    double clusRadius = sqrt(global[0] * global[0] + global[1] * global[1]);
    auto para_errors = _ClusErrPara->get_clusterv5_modified_error(c, clusRadius, key);

    localErr[0][0] = 0.;
    localErr[0][1] = 0.;
    localErr[0][2] = 0.;
    localErr[1][0] = 0.;
    localErr[1][1] = para_errors.first;
    localErr[1][2] = 0.;
    localErr[2][0] = 0.;
    localErr[2][1] = 0.;
    localErr[2][2] = para_errors.second;

    float clusphi = atan2(global(1), global(0));
    TMatrixF ROT(3, 3);
    ROT[0][0] = cos(clusphi);
    ROT[0][1] = -sin(clusphi);
    ROT[0][2] = 0.0;
    ROT[1][0] = sin(clusphi);
    ROT[1][1] = cos(clusphi);
    ROT[1][2] = 0.0;
    ROT[2][0] = 0.0;
    ROT[2][1] = 0.0;
    ROT[2][2] = 1.0;
    TMatrixF ROT_T(3, 3);
    ROT_T.Transpose(ROT);

    TMatrixF err(3, 3);
    err = ROT * localErr * ROT_T;

    return err[i][j];
  }
}

bool ALICEKF::TransportAndRotate(double old_radius, double new_radius, double& phi, GPUTPCTrackParam& kftrack, GPUTPCTrackParam::GPUTPCTrackFitParam& fp) const
{
  const double transport_spacing = .05;
  const int Ndivisions = floor(fabs(new_radius - old_radius) / transport_spacing);
  if (Verbosity() > 2)
  {
    std::cout << "old_radius: " << old_radius << ", new_radius: " << new_radius << std::endl;
  }
  if (Verbosity() > 2)
  {
    std::cout << "Ndivisions: " << Ndivisions << std::endl;
  }
  for (int i = 1; i <= Ndivisions + 1; i++)
  {
    if (std::isnan(kftrack.GetX()) ||
        std::isnan(kftrack.GetY()) ||
        std::isnan(kftrack.GetZ()))
    {
      if (Verbosity() > 0)
      {
        std::cout << "position is NaN, exiting" << std::endl;
      }
      return false;
    }

    if (new_radius > 78.)
    {
      if (Verbosity() > 1)
      {
        std::cout << "outside TPC, exiting" << std::endl;
      }
      return false;
    }

    double r_div;
    if (i == Ndivisions + 1)
    {
      r_div = new_radius;
    }
    else if (old_radius < new_radius)
    {
      r_div = old_radius + transport_spacing * i;
    }
    else
    {
      r_div = old_radius - transport_spacing * i;
    }

    // track state position relative to radial direction
    const double tX = kftrack.GetX();
    const double tY = kftrack.GetY();
    const double tz = kftrack.GetZ();

    // track state global position (including tz above)
    const double tx = tX * cos(phi) - tY * sin(phi);
    const double ty = tX * sin(phi) + tY * cos(phi);

    const double Bz = _Bzconst * get_Bz(tx, ty, tz);

    kftrack.CalculateFitParameters(fp);

    // transport to radius
    if (!kftrack.TransportToXWithMaterial(r_div, fp, Bz, 1.))
    {
      if (Verbosity() > 1)
      {
        std::cout << "transport failed" << std::endl;
      }
      return false;
    }

    // rotate track state reference frame
    const double new_tX = kftrack.GetX();
    const double new_tY = kftrack.GetY();
    const double new_tx = new_tX * cos(phi) - new_tY * sin(phi);
    const double new_ty = new_tX * sin(phi) + new_tY * cos(phi);
    const double new_phi = atan2(new_ty, new_tx);
    const double alpha = new_phi - phi;
    if (!kftrack.Rotate(alpha, 1.))
    {
      if (Verbosity() > 1)
      {
        std::cout << "rotate failed" << std::endl;
      }
      return false;
    }
    phi = new_phi;
  }

  const double final_tX = kftrack.GetX();
  const double final_tY = kftrack.GetY();
  const double final_tx = final_tX * cos(phi) - final_tY * sin(phi);
  const double final_ty = final_tX * sin(phi) + final_tY * cos(phi);
  const double final_tz = kftrack.GetZ();

  if (Verbosity() > 1)
  {
    std::cout << "track position after transport: (" << final_tx << ", " << final_ty << ", " << final_tz << ")" << std::endl;
  }
  return true;
}

bool ALICEKF::FilterStep(TrkrDefs::cluskey ckey, keylist& keys, double& current_phi, GPUTPCTrackParam& kftrack, GPUTPCTrackParam::GPUTPCTrackFitParam& fp, const PositionMap& globalPositions) const
{
  // give up if position vector has NaN for any component
  if (std::isnan(kftrack.GetX()) ||
      std::isnan(kftrack.GetY()) ||
      std::isnan(kftrack.GetZ()))
  {
    if (Verbosity() > 0)
    {
      std::cout << "position is NaN, exiting" << std::endl;
    }
    return false;
  }

  // track state position relative to radial direction
  const double tX = kftrack.GetX();
  const double tY = kftrack.GetY();
  const double tz = kftrack.GetZ();

  // track state global position (including tz above)
  const double tx = tX * cos(current_phi) - tY * sin(current_phi);
  const double ty = tX * sin(current_phi) + tY * cos(current_phi);

  if (Verbosity() > 1)
  {
    std::cout << std::endl;
    std::cout << "current phi: " << current_phi << std::endl;
    std::cout << "track position: (" << tx << ", " << ty << ", " << tz << ")" << std::endl;
  }

  const auto& cluster_pos = globalPositions.at(ckey);
  const double cx = cluster_pos(0);
  const double cy = cluster_pos(1);
  const double cz = cluster_pos(2);

  if (Verbosity() > 1)
  {
    std::cout << "cluster position: (" << cx << ", " << cy << ", " << cz << ")" << std::endl;
  }

  const double cX = sqrt(cx * cx + cy * cy);

  if (Verbosity() > 1)
  {
    std::cout << "transporting to radius: " << cX << std::endl;
  }

  if (fabs(tX - cX) < 0.1 || !TransportAndRotate(tX, cX, current_phi, kftrack, fp))
  {
    if (std::isnan(kftrack.GetX()) ||
        std::isnan(kftrack.GetY()) ||
        std::isnan(kftrack.GetZ()))
    {
      return false;
    }

    if (Verbosity() > 1)
    {
      std::cout << "track turned around near X=" << kftrack.GetX() << std::endl;
    }

    // TransportAndRotate failure indicates that track has turned around before it reaches next layer
    // The track parameters don't update in that last step, so we're likely somewhere between two layers
    // So, first, we turn the track around ourselves, setting it to its next intersection at its current radius

    // basically circle project in xy, linear project in z
    double pt = 1. / fabs(kftrack.GetQPt());
    double end_tx = kftrack.GetX() * cos(current_phi) - kftrack.GetY() * sin(current_phi);
    double end_ty = kftrack.GetX() * sin(current_phi) + kftrack.GetY() * cos(current_phi);
    double end_tz = kftrack.GetZ();
    if (Verbosity() > 1)
    {
      std::cout << "current parameters: pt=" << pt << ", (x, y, z) = (" << end_tx << ", " << end_ty << ", " << end_tz << ")" << std::endl;
    }
    // pt[GeV] = 0.3 B[T] R[m]
    double R = 100. * pt / (0.3 * get_Bz(end_tx, end_ty, end_tz));
    if (Verbosity() > 2)
    {
      std::cout << "R=" << R << std::endl;
    }
    double pX = pt * kftrack.GetCosPhi();
    double pY = pt * kftrack.GetSinPhi();
    if (Verbosity() > 2)
    {
      std::cout << "(pX, pY) = (" << pX << ", " << pY << ")" << std::endl;
    }
    double px = pX * cos(current_phi) - pY * sin(current_phi);
    double py = pX * sin(current_phi) + pY * cos(current_phi);
    double tangent_phi = atan2(py, px);
    double center_phi;
    if (kftrack.GetQPt() > 0)
    {
      center_phi = tangent_phi + M_PI / 2.;
    }
    else
    {
      center_phi = tangent_phi - M_PI / 2.;
    }
    if (center_phi > M_PI)
    {
      center_phi -= 2. * M_PI;
    }
    if (center_phi < -M_PI)
    {
      center_phi += 2. * M_PI;
    }
    double xc = end_tx - R * cos(center_phi);
    double yc = end_ty - R * sin(center_phi);
    if (Verbosity() > 2)
    {
      std::cout << "(xc, yc) = (" << xc << ", " << yc << ")" << std::endl;
    }
    auto circle_output = TrackFitUtils::circle_circle_intersection(sqrt(pow(kftrack.GetX(), 2.) + pow(kftrack.GetY(), 2.)), R, xc, yc);
    // pick the furthest point from current track position
    double new_tx;
    double new_ty;
    double xplus = std::get<0>(circle_output);
    double yplus = std::get<1>(circle_output);
    double xminus = std::get<2>(circle_output);
    double yminus = std::get<3>(circle_output);
    if (Verbosity() > 1)
    {
      std::cout << "circle-circle intersection: (" << xplus << ", " << yplus << "), (" << xminus << ", " << yminus << ")" << std::endl;
    }

    if (sqrt(pow(end_tx - xplus, 2.) + pow(end_ty - yplus, 2.)) > sqrt(pow(end_tx - xminus, 2.) + pow(end_ty - yminus, 2.)))
    {
      new_tx = xplus;
      new_ty = yplus;
    }
    else
    {
      new_tx = xminus;
      new_ty = yminus;
    }
    double rot_phi = atan2(new_ty, new_tx);
    //    double rot_alpha = rot_phi - current_phi;

    // new track point is existing track point rotated by alpha
    double new_tX = new_tx * cos(rot_phi) + new_ty * sin(rot_phi);
    double new_tY = -new_tx * sin(rot_phi) + new_ty * cos(rot_phi);
    double new_centerphi = atan2(new_ty - yc, new_tx - xc);
    double dcenterphi = new_centerphi - center_phi;
    if (dcenterphi > M_PI)
    {
      dcenterphi = 2. * M_PI - dcenterphi;
    }
    if (dcenterphi < -M_PI)
    {
      dcenterphi = 2. * M_PI + dcenterphi;
    }
    double ds = R * fabs(dcenterphi);
    double dz = kftrack.GetDzDs() * ds;

    current_phi = rot_phi;
    kftrack.SetX(new_tX);
    kftrack.SetY(new_tY);
    kftrack.SetZ(end_tz + dz);
    // no change to sinPhi
    // no change to DzDs
    // no change to Q/pt

    kftrack.SetSignCosPhi(-kftrack.GetSignCosPhi());

    // now finish transport
    if (!TransportAndRotate(kftrack.GetX(), cX, current_phi, kftrack, fp))
    {
      return false;
    }
  }

  if (Verbosity() > 1)
  {
    const double new_tX = kftrack.GetX();
    const double new_tY = kftrack.GetY();
    const double new_tx = new_tX * cos(current_phi) - new_tY * sin(current_phi);
    const double new_ty = new_tX * sin(current_phi) + new_tY * cos(current_phi);
    const double new_tz = kftrack.GetZ();
    std::cout << "cluster position: (" << cx << ", " << cy << ", " << cz << ")" << std::endl;
    std::cout << "current phi: " << current_phi << std::endl;
    std::cout << "new track position: (" << new_tx << ", " << new_ty << ", " << new_tz << ")" << std::endl;

    const double tYerr = sqrt(kftrack.GetCov(0));
    const double txerr = fabs(tYerr * sin(current_phi));
    const double tyerr = fabs(tYerr * cos(current_phi));
    const double tzerr = sqrt(kftrack.GetCov(5));
    std::cout << "track position error: (" << txerr << ", " << tyerr << ", " << tzerr << ")" << std::endl;
  }

  TrkrCluster* cluster = _cluster_map->findCluster(ckey);
  const double cxerr = sqrt(getClusterError(cluster, ckey, cluster_pos, 0, 0));
  const double cyerr = sqrt(getClusterError(cluster, ckey, cluster_pos, 1, 1));
  const double czerr = sqrt(getClusterError(cluster, ckey, cluster_pos, 2, 2));

  if (Verbosity() > 1)
  {
    std::cout << "cluster position error: (" << cxerr << ", " << cyerr << ", " << czerr << ")" << std::endl;
  }

  const double cY = -cx * sin(current_phi) + cy * cos(current_phi);
  const double cxycov2 = getClusterError(cluster, ckey, cluster_pos, 0, 1);
  const double cYerr2 = cxerr * cxerr * sin(current_phi) * sin(current_phi) + cxycov2 * sin(current_phi) * cos(current_phi) + cyerr * cyerr * cos(current_phi) * cos(current_phi);
  const double czerr2 = czerr * czerr;

  if (Verbosity() > 1)
  {
    std::cout << "Filtering cluster with Y=" << cY << ", z=" << cz << ", Yerr2=" << cYerr2 << ", zerr2=" << czerr2 << std::endl;
  }

  bool filter_success = kftrack.Filter(cY, cz, cYerr2, czerr2, _max_sin_phi);
  if (Verbosity() > 1)
  {
    std::cout << "position after filter: (" << kftrack.GetX() * cos(current_phi) - kftrack.GetY() * sin(current_phi) << ", " << kftrack.GetX() * sin(current_phi) + kftrack.GetY() * cos(current_phi) << ", " << kftrack.GetZ() << std::endl;
    std::cout << "track current parameters:" << std::endl;
    std::cout << "(X, Y, Z) = (" << kftrack.GetX() << ", " << kftrack.GetY() << ", " << kftrack.GetZ() << ")" << std::endl;
    std::cout << "pt = " << 1. / fabs(kftrack.GetQPt()) << std::endl;
    std::cout << "QPt = " << kftrack.GetQPt() << std::endl;
    std::cout << "sinPhi = " << kftrack.GetSinPhi() << std::endl;
    std::cout << "cosPhi = " << kftrack.GetCosPhi() << std::endl;
    std::cout << "dzds = " << kftrack.GetDzDs() << std::endl;
  }
  if (!filter_success)
  {
    keys.erase(std::find(keys.begin(), keys.end(), ckey));
  }
  return true;
}

TrackSeedAliceSeedMap ALICEKF::ALICEKalmanFilter(const std::vector<keylist>& trackSeedKeyLists, bool use_nhits_limit, const PositionMap& globalPositions, std::vector<float>& trackChi2) const
{
  //  TFile* f = new TFile("/sphenix/u/mjpeters/macros_hybrid/detectors/sPHENIX/pull.root", "RECREATE");
  //  TNtuple* ntp = new TNtuple("pull","pull","cx:cy:cz:xerr:yerr:zerr:tx:ty:tz:layer:xsize:ysize:phisize:phierr:zsize");
  std::vector<TrackSeed_v2> seeds_vector;
  std::vector<GPUTPCTrackParam> alice_seeds_vector;
  int nseeds = 0;
  int ncandidates = -1;
  if (Verbosity() > 0)
  {
    std::cout << "min clusters per track: " << _min_clusters_per_track << "\n";
  }
  for (auto trackKeyChain : trackSeedKeyLists)
  {
    ++ncandidates;

    if (trackKeyChain.size() < 2)
    {
      continue;
    }
    if (use_nhits_limit && trackKeyChain.size() < _min_clusters_per_track)
    {
      continue;
    }
    //    if(TrkrDefs::getLayer(trackKeyChain.front())<TrkrDefs::getLayer(trackKeyChain.back())) { std::reverse(trackKeyChain.begin(),trackKeyChain.end()); }
    // get starting cluster from key
    // Transform sPHENIX coordinates into ALICE-compatible coordinates
    const auto& globalpos = globalPositions.at(trackKeyChain.at(0));
    double x0 = globalpos(0);
    double y0 = globalpos(1);
    double z0 = globalpos(2);
    ;
    LogDebug("Initial (x,y,z): (" << x0 << "," << y0 << "," << z0 << ")" << std::endl);
    // ALICE x coordinate = distance from beampipe
    double alice_x0 = sqrt(x0 * x0 + y0 * y0);
    double alice_y0 = 0;
    double alice_z0 = z0;
    // Initialize track and linearisation
    GPUTPCTrackParam trackSeed;
    trackSeed.setNeonFraction(Ne_frac);
    trackSeed.setArgonFraction(Ar_frac);
    trackSeed.setCF4Fraction(CF4_frac);
    trackSeed.setNitrogenFraction(N2_frac);
    trackSeed.setIsobutaneFraction(isobutane_frac);
    trackSeed.InitParam();
    trackSeed.SetX(alice_x0);
    trackSeed.SetY(alice_y0);
    trackSeed.SetZ(alice_z0);
    double x = x0;
    double y = y0;
#if defined(_DEBUG_)
    double z = z0;
    double alice_x = sqrt(x0 * x0 + y0 * y0);
#endif
    // Pre-set momentum-based parameters to improve numerical stability
    const auto& secondpos = globalPositions.at(trackKeyChain.at(1));

    const double second_x = secondpos(0);
    const double second_y = secondpos(1);
    const double second_z = secondpos(2);
    const double first_phi = atan2(y0, x0);
    const double second_alice_x = second_x * std::cos(first_phi) + second_y * std::sin(first_phi);
    const double delta_alice_x = second_alice_x - alice_x0;
    // double second_alice_y = (second_x/cos(first_phi)-second_y/sin(first_phi))/(sin(first_phi)/cos(first_phi)+cos(first_phi)/sin(first_phi));
    const double second_alice_y = -second_x * std::sin(first_phi) + second_y * std::cos(first_phi);
    double init_SinPhi = second_alice_y / std::sqrt(square(delta_alice_x) + square(second_alice_y));
    const double delta_z = second_z - z0;
    double init_DzDs = delta_z / std::sqrt(square(delta_alice_x) + square(second_alice_y));
    if (delta_alice_x < 0.)
    {
      init_SinPhi *= -1.;
      init_DzDs *= -1.;
    }
    trackSeed.SetSinPhi(init_SinPhi);
    // trackSeed.SetSignCosPhi(delta_alice_x / std::sqrt(square(delta_alice_x) + square(second_alice_y)));
    LogDebug("Set initial SinPhi to " << init_SinPhi << std::endl);
    trackSeed.SetDzDs(init_DzDs);
    LogDebug("Set initial DzDs to " << init_DzDs << std::endl);

    // get initial pt estimate
    std::vector<std::pair<double, double>> pts;
    std::transform(trackKeyChain.begin(), trackKeyChain.end(), std::back_inserter(pts), [&globalPositions](const TrkrDefs::cluskey& key)
                   {
      const auto& clpos = globalPositions.at(key);
      return std::make_pair(clpos(0),clpos(1)); });

    const auto [R, x_center, y_center] = TrackFitUtils::circle_fit_by_taubin(pts);
    if (Verbosity() > 1)
    {
      std::cout << std::endl
                << "candidate " << ncandidates << " seed " << nseeds << " circle fit parameters: R=" << R << ", X0=" << x_center << ", Y0=" << y_center << std::endl;
    }

    // check circle fit success
    /* failed fit will result in infinite momentum for the track, which in turn will break the kalman filter */
    if (std::isnan(R))
    {
      continue;
    }

    double init_QPt = 1. / (0.3 * R / 100. * get_Bz(x0, y0, z0));
    // determine charge
    double phi_first = atan2(y0, x0);
    if (Verbosity() > 1)
    {
      std::cout << "phi_first: " << phi_first << std::endl;
    }
    double phi_second = atan2(second_y, second_x);
    if (Verbosity() > 1)
    {
      std::cout << "phi_second: " << phi_second << std::endl;
    }
    double dphi = phi_second - phi_first;
    if (Verbosity() > 1)
    {
      std::cout << "dphi: " << dphi << std::endl;
    }
    if (dphi > M_PI)
    {
      dphi = 2 * M_PI - dphi;
    }
    if (dphi < -M_PI)
    {
      dphi = 2 * M_PI + dphi;
    }
    if (Verbosity() > 1)
    {
      std::cout << "corrected dphi: " << dphi << std::endl;
    }
    if ((dphi > 0. && x0 * x0 + y0 * y0 < second_x * second_x + second_y * second_y) ||
        (dphi < 0. && x0 * x0 + y0 * y0 > second_x * second_x + second_y * second_y))
    {
      init_QPt = -1 * init_QPt;
    }
    LogDebug("initial QPt: " << init_QPt << std::endl);
    trackSeed.SetQPt(init_QPt);
    /*
        if (trackSeed.GetSignCosPhi() < 0) {
          trackSeed.SetSignCosPhi(-trackSeed.GetSignCosPhi());
          trackSeed.SetSinPhi(-trackSeed.GetSinPhi());
          trackSeed.SetDzDs(-trackSeed.GetDzDs());
          trackSeed.SetQPt(-trackSeed.GetQPt());
          trackSeed.SetCov(3,-trackSeed.GetCov(3));
          trackSeed.SetCov(4,-trackSeed.GetCov(4));
          trackSeed.SetCov(6,-trackSeed.GetCov(6));
          trackSeed.SetCov(7,-trackSeed.GetCov(7));
          trackSeed.SetCov(10,-trackSeed.GetCov(10));
          trackSeed.SetCov(11,-trackSeed.GetCov(11));
        }
    */
    GPUTPCTrackLinearisation trackLine(trackSeed);
    GPUTPCTrackParam::GPUTPCTrackFitParam fp{};
    trackSeed.CalculateFitParameters(fp);

    LogDebug(std::endl
             << std::endl
             << "------------------------" << std::endl
             << "seed size: " << trackKeyChain.size() << std::endl
             << std::endl
             << std::endl);
    double current_phi = phi_first;
    bool filter_failed = false;
    keylist outputKeyChain = trackKeyChain;
    for (auto clusterkey = std::next(trackKeyChain.begin()); clusterkey != trackKeyChain.end(); ++clusterkey)
    {
      if (!FilterStep(*clusterkey, outputKeyChain, current_phi, trackSeed, fp, globalPositions))
      {
        if (Verbosity() > 0)
        {
          std::cout << "Kalman filter failed, exiting" << std::endl;
        }
        filter_failed = true;
        break;
      }
    }
    if (filter_failed)
    {
      continue;
    }
    if (Verbosity() > 0)
    {
      std::cout << "finished track\n";
    }

    double track_phi = atan2(y, x);

    if (Verbosity() > 1)
    {
      std::cout << "final QPt = " << trackSeed.GetQPt() << std::endl;
    }

    double track_pt = fabs(1. / trackSeed.GetQPt());
#if defined(_DEBUG_)
    double track_pY = track_pt * trackSeed.GetSinPhi();
    double track_pX = sqrt(track_pt * track_pt - track_pY * track_pY);
    double track_px = track_pX * cos(track_phi) - track_pY * sin(track_phi);
    double track_py = track_pX * sin(track_phi) + track_pY * cos(track_phi);
    double track_pz = track_pt * trackSeed.GetDzDs();
#endif
    double track_pterr = sqrt(trackSeed.GetErr2QPt()) / (trackSeed.GetQPt() * trackSeed.GetQPt());
    // If Kalman filter doesn't do its job (happens often with short seeds), use the circle-fit estimate as the central value
    // if(trackKeyChain.size()<10) {track_pt = fabs(1./init_QPt);}
    LogDebug("track pt = " << track_pt << " +- " << track_pterr << std::endl);
    LogDebug("track ALICE p = (" << track_pX << ", " << track_pY << ", " << track_pz << ")" << std::endl);
    LogDebug("track p = (" << track_px << ", " << track_py << ", " << track_pz << ")" << std::endl);

    /*
        if(cluster_ctr!=1 && !trackSeed.CheckNumericalQuality())
        {
          std::cout << "ERROR: Track seed failed numerical quality check before conversion to sPHENIX coordinates! Skipping this one.\n";
          aborted = true;
          continue;
        }
    */
    //    pt:z:dz:phi:dphi:c:dc
    // Fill NT with track parameters
    // double StartEta = -log(tan(atan(z0/sqrt(x0*x0+y0*y0))));
    //    if(aborted) continue;
    //    double track_pt = fabs( 1./(trackSeed.GetQPt()));
    if (checknan(track_pt, "pT", nseeds))
    {
      continue;
    }
    //    double track_pterr = sqrt(trackSeed.GetErr2QPt())/(trackSeed.GetQPt()*trackSeed.GetQPt());
    if (checknan(track_pterr, "pT err", nseeds))
    {
      continue;
    }
    LogDebug("Track pterr = " << track_pterr << std::endl);
    double track_x = trackSeed.GetX() * cos(track_phi) - trackSeed.GetY() * sin(track_phi);
    double track_y = trackSeed.GetX() * sin(track_phi) + trackSeed.GetY() * cos(track_phi);
    double track_z = trackSeed.GetZ();
    if (checknan(track_z, "z", nseeds))
    {
      continue;
    }
    double track_zerr = sqrt(trackSeed.GetErr2Z());
    if (checknan(track_zerr, "zerr", nseeds))
    {
      continue;
    }
    auto lcluster = _cluster_map->findCluster(trackKeyChain.back());
    const auto& lclusterglob = globalPositions.at(trackKeyChain.back());
    const float lclusterrad = sqrt(lclusterglob(0) * lclusterglob(0) + lclusterglob(1) * lclusterglob(1));
    double last_cluster_phierr = lcluster->getRPhiError() / lclusterrad;
    ;

    // phi error assuming error in track radial coordinate is zero
    double track_phierr = sqrt(pow(last_cluster_phierr, 2) + (pow(trackSeed.GetX(), 2) * trackSeed.GetErr2Y()) /
                                                                 pow(pow(trackSeed.GetX(), 2) + pow(trackSeed.GetY(), 2), 2));
    if (checknan(track_phierr, "phierr", nseeds))
    {
      continue;
    }
    LogDebug("Track phi = " << atan2(track_py, track_px) << std::endl);
    LogDebug("Track phierr = " << track_phierr << std::endl);
    double track_curvature = trackSeed.GetKappa(_Bzconst * get_Bz(track_x, track_y, track_z));
    if (checknan(track_curvature, "curvature", nseeds))
    {
      continue;
    }
    double track_curverr = sqrt(trackSeed.GetErr2QPt()) * _Bzconst * get_Bz(track_x, track_y, track_z);
    if (checknan(track_curverr, "curvature error", nseeds))
    {
      continue;
    }
    TrackSeed_v2 track;
    //    track.set_vertex_id(_vertex_ids[best_vtx]);
    for (unsigned long j : outputKeyChain)
    {
      track.insert_cluster_key(j);
    }

    int track_charge = 0;
    if (trackSeed.GetQPt() < 0)
    {
      track_charge = -1 * trackSeed.GetSignCosPhi();  // * _fieldDir;
    }
    else
    {
      track_charge = 1 * trackSeed.GetSignCosPhi();  // * _fieldDir;
    }

    double sinphi = sin(track_phi);  // who had the idea to use s here?????
    double c = cos(track_phi);
    double p = trackSeed.GetSinPhi();

    /// Shows the transformation between ALICE and sPHENIX coordinates
    // track.set_x(trackSeed.GetX()*c-trackSeed.GetY()*s);//_vertex_x[best_vtx]);  //track.set_x(cl->getX());
    // track.set_y(trackSeed.GetX()*s+trackSeed.GetY()*c);//_vertex_y[best_vtx]);  //track.set_y(cl->getY());
    // track.set_z(trackSeed.GetZ());//_vertex_z[best_vtx]);  //track.set_z(cl->getZ());
    // if(Verbosity()>0) std::cout << "x " << track.get_x() << "\n";
    // if(Verbosity()>0) std::cout << "y " << track.get_y() << "\n";
    // if(Verbosity()>0) std::cout << "z " << track.get_z() << "\n";
    // if(checknan(p,"ALICE sinPhi",nseeds)) continue;
    double d = trackSeed.GetDzDs();
    if (checknan(d, "ALICE dz/ds", nseeds))
    {
      continue;
    }

    /// Shows the transformation between ALICE and sPHENIX coordinates
    // double pY = track_pt*p;
    // double pX = sqrt(track_pt*track_pt-pY*pY);
    /// We set the qoverR to get the good charge estimate from the KF
    /// which helps the Acts fit
    track.set_qOverR(trackSeed.GetQPt() * (0.3 * _const_field) / 100.);
    // track.set_px(pX*c-pY*s);
    // track.set_py(pX*s+pY*c);
    // track.set_pz(track_pt * trackSeed.GetDzDs());
    const double* cov = trackSeed.GetCov();
    bool cov_nan = false;
    for (int i = 0; i < 15; i++)
    {
      if (checknan(cov[i], "covariance element " + std::to_string(i), nseeds))
      {
        cov_nan = true;
      }
    }
    if (cov_nan)
    {
      continue;
    }
    // make this into an actual Eigen matrix
    Eigen::Matrix<double, 5, 5> ecov;
    ecov(0, 0) = cov[0];
    ecov(0, 1) = cov[1];
    ecov(0, 2) = cov[2];
    ecov(0, 3) = cov[3];
    ecov(0, 4) = cov[4];
    ecov(1, 1) = cov[5];
    ecov(1, 2) = cov[6];
    ecov(1, 3) = cov[7];
    ecov(1, 4) = cov[8];
    ecov(2, 2) = cov[9];
    ecov(2, 3) = cov[10];
    ecov(2, 4) = cov[11];
    ecov(3, 3) = cov[12];
    ecov(3, 4) = cov[13];
    ecov(4, 4) = cov[14];
    // symmetrize
    ecov(1, 0) = ecov(0, 1);
    ecov(2, 0) = ecov(0, 2);
    ecov(3, 0) = ecov(0, 3);
    ecov(4, 0) = ecov(0, 4);
    ecov(2, 1) = ecov(1, 2);
    ecov(3, 1) = ecov(1, 3);
    ecov(4, 1) = ecov(1, 4);
    ecov(3, 2) = ecov(2, 3);
    ecov(4, 2) = ecov(2, 4);
    ecov(4, 3) = ecov(3, 4);
    // make rotation matrix based on the following:
    // x = X*cos(track_phi) - Y*sin(track_phi)
    // y = X*sin(track_phi) + Y*cos(track_phi)
    // z = Z
    // pY = pt*sinphi
    // pX = sqrt(pt**2 - pY**2)
    // px = pX*cos(track_phi) - pY*sin(track_phi)
    // py = pX*sin(track_phi) + pY*cos(track_phi)
    // pz = pt*(dz/ds)
    Eigen::Matrix<double, 6, 5> J;
    J(0, 0) = -sinphi;  // dx/dY
    J(0, 1) = 0.;       // dx/dZ
    J(0, 2) = 0.;       // dx/d(sinphi)
    J(0, 3) = 0.;       // dx/d(dz/ds)
    J(0, 4) = 0.;       // dx/d(Q/pt)

    J(1, 0) = c;   // dy/dY
    J(1, 1) = 0.;  // dy/dZ
    J(1, 2) = 0.;  // dy/d(sinphi)
    J(1, 3) = 0.;  // dy/d(dz/ds)
    J(1, 4) = 0.;  // dy/d(Q/pt)

    J(2, 0) = 0.;  // dz/dY
    J(2, 1) = 1.;  // dz/dZ
    J(2, 2) = 0.;  // dz/d(sinphi)
    J(2, 3) = 0.;  // dz/d(dz/ds)
    J(2, 4) = 0.;  // dz/d(Q/pt)

    J(3, 0) = 0.;                                                                       // dpx/dY
    J(3, 1) = 0.;                                                                       // dpx/dZ
    J(3, 2) = -track_pt * (p * c / sqrt(1 - p * p) + sinphi);                           // dpx/d(sinphi)
    J(3, 3) = 0.;                                                                       // dpx/d(dz/ds)
    J(3, 4) = track_pt * track_pt * track_charge * (p * sinphi - c * sqrt(1 - p * p));  // dpx/d(Q/pt)

    J(4, 0) = 0.;                                                                        // dpy/dY
    J(4, 1) = 0.;                                                                        // dpy/dZ
    J(4, 2) = track_pt * (c - p * sinphi / sqrt(1 - p * p));                             // dpy/d(sinphi)
    J(4, 3) = 0.;                                                                        // dpy/d(dz/ds)
    J(4, 4) = -track_pt * track_pt * track_charge * (p * c + sinphi * sqrt(1 - p * p));  // dpy/d(Q/pt)

    J(5, 0) = 0.;                                       // dpz/dY
    J(5, 1) = 0.;                                       // dpz/dZ
    J(5, 2) = 0.;                                       // dpz/d(sinphi)
    J(5, 3) = track_pt;                                 // dpz/d(dz/ds)
    J(5, 4) = -track_pt * track_pt * track_charge * d;  // dpz/d(Q/pt)
                                                        /*    bool cov_rot_nan = false;
                                                            for(int i=0;i<6;i++)
                                                            {
                                                              for(int j=0;j<5;j++)
                                                              {
                                                                if(checknan(J(i,j),"covariance rotator element ("+std::to_string(i)+","+std::to_string(j)+")",nseeds))
                                                                {
                                                                  cov_rot_nan = true;
                                                                  continue;
                                                                }
                                                              }
                                                            }
                                                            if(cov_rot_nan) continue;
                                                        */
    // the heavy lifting happens here
    Eigen::Matrix<double, 6, 6> scov = J * ecov * J.transpose();
    if (!covIsPosDef(scov))
    {
      repairCovariance(scov);
    }
    /*
    // Derived from:
    // 1) Taking the Jacobian of the conversion from (Y,Z,SinPhi,DzDs,Q/Pt) to (x,y,z,px,py,pz)
    // 2) Computing (Jacobian)*(ALICE covariance matrix)*(transpose of Jacobian)
    track.set_error(0, 0, cov[0]*s*s);
    track.set_error(0, 1, -cov[0]*c*s);
    track.set_error(0, 2, -cov[1]*s);
    track.set_error(0, 3, cov[2]*s*s/q-cov[4]*s*(-c/(q*q)+p*s/(q*q)));
    track.set_error(0, 4, -cov[2]*c*s/q-cov[4]*s*(-c*p/(q*q)-s/(q*q)));
    track.set_error(0, 5, cov[4]*d*s/(q*q)-cov[3]*s/q);
    track.set_error(1, 1, cov[0]*c*c);
    track.set_error(1, 2, cov[1]*c);
    track.set_error(1, 3, -cov[2]*c*s/q+cov[4]*c*(-c/(q*q)+p*s/(q*q)));
    track.set_error(1, 4, cov[2]*c*c/q+cov[4]*c*(-c*p/(q*q)-s/(q*q)));
    track.set_error(1, 5, cov[4]*d*c/(q*q)+cov[3]*c/q);
    track.set_error(2, 2, cov[5]);
    track.set_error(2, 3, -cov[6]*s/q+cov[8]*(-c/(q*q)+p*s/(q*q)));
    track.set_error(2, 4, cov[6]*c/q+cov[8]*(-c*p/(q*q)-s/(q*q)));
    track.set_error(2, 5, -cov[8]*d/(q*q)+cov[7]/q);
    track.set_error(3, 3, cov[9]*s*s/(q*q)-cov[11]*(-c/(q*q*q)+p*s/(q*q*q)) + (-c/(q*q)+p*s/(q*q))*(-cov[11]*s/q+cov[14]*(-c/(q*q)+p*s/(q*q))));
    track.set_error(3, 4, -cov[9]*c*s/(q*q)+cov[11]*(-c/(q*q*q)+p*s/(q*q*q)) + (-c*p/(q*q)-s/(q*q))*(-cov[11]*s/q+cov[14]*(-c/(q*q)+p*s/(q*q))));
    track.set_error(3, 5, -cov[10]*s/(q*q)+cov[13]/q*(-c/(q*q)+p*s/(q*q))-d/(q*q)*(-cov[11]*s/q+cov[14]*(-c/(q*q)+p*s/(q*q))));
    track.set_error(4, 4, c/q*(c/q*cov[9]+cov[11]*(-c*p/(q*q)-s/(q*q)))+(-c*p/(q*q)-s/(q*q))*(c/q*cov[11]+cov[14]*(-c*p/(q*q)-s/(q*q))));
    track.set_error(4, 5, cov[10]*c/(q*q)+cov[13]/q*(-c*p/(q*q)-s/(q*q))-d/(q*q)*(c/q*cov[11]+cov[14]*(-c*p/(q*q)-s/(q*q))));
    track.set_error(5, 5, -d/(q*q)*(-d*cov[14]/(q*q)+cov[13]/q)-d*cov[13]/(q*q*q)+cov[12]/(q*q));
    // symmetrize covariance
    track.set_error(1, 0, track.get_error(0, 1));
    track.set_error(2, 0, track.get_error(0, 2));
    track.set_error(3, 0, track.get_error(0, 3));
    track.set_error(4, 0, track.get_error(0, 4));
    track.set_error(5, 0, track.get_error(0, 5));
    track.set_error(2, 1, track.get_error(1, 2));
    track.set_error(3, 1, track.get_error(1, 3));
    track.set_error(4, 1, track.get_error(1, 4));
    track.set_error(5, 1, track.get_error(1, 5));
    track.set_error(3, 2, track.get_error(2, 3));
    track.set_error(4, 2, track.get_error(2, 4));
    track.set_error(5, 2, track.get_error(2, 5));
    track.set_error(4, 3, track.get_error(3, 4));
    track.set_error(5, 3, track.get_error(3, 5));
    track.set_error(5, 4, track.get_error(4, 5));
*/

    /*
        for(int w=0;w<cx.size();w++)
        {
          ntp->Fill(cx[w],cy[w],cz[w],xerr[w],yerr[w],zerr[w],tx[w],ty[w],tz[w],layer[w],xsize[w],ysize[w],phisize[w],phierr[w],zsize[w]);
        }
        cx.clear();
        cy.clear();
        cz.clear();
        tx.clear();
        ty.clear();
        tz.clear();
        xerr.clear();
        yerr.clear();
        zerr.clear();
        layer.clear();
        xsize.clear();
        ysize.clear();
        phisize.clear();
        phierr.clear();
        zsize.clear();
    */
    seeds_vector.push_back(track);
    alice_seeds_vector.push_back(trackSeed);
    trackChi2.push_back(trackSeed.GetChi2() / trackSeed.GetNDF());

    ++nseeds;
  }
  //  f->cd();
  //  ntp->Write();
  //  f->Close();
  if (Verbosity() > 0)
  {
    std::cout << "number of seeds: " << nseeds << "\n";
  }

  return std::make_pair(seeds_vector, alice_seeds_vector);
}

bool ALICEKF::covIsPosDef(Eigen::Matrix<double, 6, 6>& cov) const
{
  // attempt Cholesky decomposition
  Eigen::LLT<Eigen::Matrix<double, 6, 6>> chDec(cov);
  // if Cholesky decomposition does not exist, matrix is not positive definite
  return (chDec.info() != Eigen::NumericalIssue);
}

void ALICEKF::repairCovariance(Eigen::Matrix<double, 6, 6>& cov) const
{
  Eigen::Matrix<double, 6, 6> repaircov = cov;
  // find closest positive definite matrix
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 6, 6>> solver(repaircov);
  const Eigen::Matrix<double, 6, 1>& D = solver.eigenvalues();
  Eigen::Matrix<double, 6, 6> Q = solver.eigenvectors();
  Eigen::Matrix<double, 6, 1> Dp = D.cwiseMax(1e-15);
  Eigen::Matrix<double, 6, 6> Z = Q * Dp.asDiagonal() * Q.transpose();
  // updates covariance matrix
  for (int i = 0; i < 6; i++)
  {
    for (int j = 0; j < 6; j++)
    {
      cov(i, j) = Z(i, j);
    }
  }
}

std::vector<double> ALICEKF::GetCircleClusterResiduals(const std::vector<std::pair<double, double>>& points, double R, double X0, double Y0) const
{
  std::vector<double> residues;
  std::transform(points.begin(), points.end(), std::back_inserter(residues), [R, X0, Y0](const std::pair<double, double>& point)
                 {
    double x = point.first;
    double y = point.second;

    // The shortest distance of a point from a circle is along the radial; line from the circle center to the point
    return std::sqrt( square(x-X0) + square(y-Y0) )  -  R; });
  return residues;
}

std::vector<double> ALICEKF::GetLineClusterResiduals(const std::vector<std::pair<double, double>>& points, double A, double B) const
{
  std::vector<double> residues;
  // calculate cluster residuals from the fitted circle
  std::transform(points.begin(), points.end(), std::back_inserter(residues), [A, B](const std::pair<double, double>& point)
                 {
    double r = point.first;
    double z = point.second;

    // The shortest distance of a point from a circle is along the radial; line from the circle center to the point

    double a = -A;
    double b = 1.0;
    double c = -B;
    return std::abs(a*r+b*z+c)/sqrt(square(a)+square(b)); });
  return residues;
}
