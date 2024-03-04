#include "TrackFitUtils.h"

#include <trackbase/MvtxDefs.h>
#include "ActsGeometry.h"
#include "TpcDefs.h"
#include "TrkrClusterContainerv4.h"
#include "TrkrClusterv5.h"
#include "TrkrDefs.h"  // for cluskey, getTrkrId, tpcId

#include <cmath>

namespace
{

  //! convenience square method
  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }
}  // namespace

std::pair<Acts::Vector3, Acts::Vector3> TrackFitUtils::get_helix_tangent(const std::vector<float>& fitpars, Acts::Vector3& global)
{
  // no analytic solution for the coordinates of the closest approach of a helix to a point
  // Instead, we get the PCA in x and y to the circle, and the PCA in z to the z vs R line at the R of the PCA

  float radius = fitpars[0];
  float x0 = fitpars[1];
  float y0 = fitpars[2];
  float zslope = fitpars[3];
  float z0 = fitpars[4];

  Acts::Vector2 pca_circle = TrackFitUtils::get_circle_point_pca(radius, x0, y0, global);

  // The radius of the PCA determines the z position:
  float pca_circle_radius = pca_circle.norm();  // radius of the PCA of the circle to the point
  float pca_z = pca_circle_radius * zslope + z0;
  Acts::Vector3 pca(pca_circle(0), pca_circle(1), pca_z);

  // now we want a second point on the helix so we can get a local straight line approximation to the track
  // Get the angle of the PCA relative to the fitted circle center
  float angle_pca = atan2(pca_circle(1) - y0, pca_circle(0) - x0);
  // calculate coords of a point at a slightly larger angle
  float d_angle = 0.005;
  float newx = radius * std::cos(angle_pca + d_angle) + x0;
  float newy = radius * std::sin(angle_pca + d_angle) + y0;
  float newz = std::sqrt(newx * newx + newy * newy) * zslope + z0;
  Acts::Vector3 second_point_pca(newx, newy, newz);

  // pca and second_point_pca define a straight line approximation to the track
  Acts::Vector3 tangent = (second_point_pca - pca) / (second_point_pca - pca).norm();

  // get the PCA of the cluster to that line
  // Approximate track with a straight line consisting of the state position posref and the vector (px,py,pz)

  // The position of the closest point on the line to global is:
  // posref + projection of difference between the point and posref on the tangent vector
  Acts::Vector3 final_pca = pca + ((global - pca).dot(tangent)) * tangent;

  auto line = std::make_pair(final_pca, tangent);

  return line;
}

//_________________________________________________________________________________
TrackFitUtils::circle_fit_output_t TrackFitUtils::circle_fit_by_taubin(const TrackFitUtils::position_vector_t& positions)
{
  // Compute x- and y- sample means
  double meanX = 0;
  double meanY = 0;
  double weight = 0;

  for (const auto& [x, y] : positions)
  {
    meanX += x;
    meanY += y;
    ++weight;
  }
  meanX /= weight;
  meanY /= weight;

  //     computing moments

  double Mxx = 0;
  double Myy = 0;
  double Mxy = 0;
  double Mxz = 0;
  double Myz = 0;
  double Mzz = 0;

  for (auto& [x, y] : positions)
  {
    double Xi = x - meanX;  //  centered x-coordinates
    double Yi = y - meanY;  //  centered y-coordinates
    double Zi = square(Xi) + square(Yi);

    Mxy += Xi * Yi;
    Mxx += Xi * Xi;
    Myy += Yi * Yi;
    Mxz += Xi * Zi;
    Myz += Yi * Zi;
    Mzz += Zi * Zi;
  }
  Mxx /= weight;
  Myy /= weight;
  Mxy /= weight;
  Mxz /= weight;
  Myz /= weight;
  Mzz /= weight;

  //  computing coefficients of the characteristic polynomial

  const double Mz = Mxx + Myy;
  const double Cov_xy = Mxx * Myy - Mxy * Mxy;
  const double Var_z = Mzz - Mz * Mz;
  const double A3 = 4 * Mz;
  const double A2 = -3 * Mz * Mz - Mzz;
  const double A1 = Var_z * Mz + 4 * Cov_xy * Mz - Mxz * Mxz - Myz * Myz;
  const double A0 = Mxz * (Mxz * Myy - Myz * Mxy) + Myz * (Myz * Mxx - Mxz * Mxy) - Var_z * Cov_xy;
  const double A22 = A2 + A2;
  const double A33 = A3 + A3 + A3;

  //    finding the root of the characteristic polynomial
  //    using Newton's method starting at x=0
  //    (it is guaranteed to converge to the right root)
  static constexpr int iter_max = 99;
  double x = 0;
  double y = A0;

  // usually, 4-6 iterations are enough
  for (int iter = 0; iter < iter_max; ++iter)
  {
    const double Dy = A1 + x * (A22 + A33 * x);
    const double xnew = x - y / Dy;
    if ((xnew == x) || (!std::isfinite(xnew)))
    {
      break;
    }

    const double ynew = A0 + xnew * (A1 + xnew * (A2 + xnew * A3));
    if (std::abs(ynew) >= std::abs(y))
    {
      break;
    }

    x = xnew;
    y = ynew;
  }

  //  computing parameters of the fitting circle
  const double DET = square(x) - x * Mz + Cov_xy;
  const double Xcenter = (Mxz * (Myy - x) - Myz * Mxy) / DET / 2;
  const double Ycenter = (Myz * (Mxx - x) - Mxz * Mxy) / DET / 2;

  //  assembling the output
  double X0 = Xcenter + meanX;
  double Y0 = Ycenter + meanY;
  double R = std::sqrt(square(Xcenter) + square(Ycenter) + Mz);
  return std::make_tuple(R, X0, Y0);
}

//_________________________________________________________________________________
TrackFitUtils::circle_fit_output_t TrackFitUtils::circle_fit_by_taubin(const std::vector<Acts::Vector3>& positions)
{
  position_vector_t positions_2d;
  for (const auto& position : positions)
  {
    positions_2d.emplace_back(position.x(), position.y());
  }

  return circle_fit_by_taubin(positions_2d);
}

//_________________________________________________________________________________
TrackFitUtils::line_fit_output_t TrackFitUtils::line_fit(const TrackFitUtils::position_vector_t& positions)
{
  double xsum = 0;
  double x2sum = 0;
  double ysum = 0;
  double xysum = 0;
  double y2sum = 0;
  for (const auto& [r, z] : positions)
  {
    xsum = xsum + r;            // calculate sigma(xi)
    ysum = ysum + z;            // calculate sigma(yi)
    x2sum = x2sum + square(r);  // calculate sigma(x^2i)
    y2sum = y2sum + square(z);  // calculate sigma(y^2i)
    xysum = xysum + r * z;      // calculate sigma(xi*yi)
  }

  const auto npts = positions.size();
  const double denominator = (x2sum * npts - square(xsum));
  const double a = (xysum * npts - xsum * ysum) / denominator;   // calculate slope
  const double b = (x2sum * ysum - xsum * xysum) / denominator;  // calculate intercept

  //! The fit fails for the case where the line is close to vertical, because
  //! there is little dependence on x and thus the denominator is very small
  //! we can swap the x-y points to find a y-x horizontal line in this case
  //! then check which line minimizes the sums of squares of residuals
  //! that is the best fit, and should be returned
  const double denominatoryx = (y2sum * npts - square(ysum));
  const double inva = (xysum * npts - xsum * ysum) / denominatoryx;
  const double invb = (y2sum * xsum - ysum * xysum) / denominatoryx;

  //! now determine which one minimizes the sum of residuals to return
  float sumresid = 0;
  float suminvresid = 0;
  for (const auto& [r, z] : positions)
  {
    sumresid += square(z - (a * r + b));
    suminvresid += square(r - (inva * z + invb));
  }

  if (sumresid < suminvresid)
  {
    return std::make_tuple(a, b);
  }

  //! the best fit is swapped in y-x, so find the perpendicular line
  return std::make_tuple(1. / inva, -invb / inva);
}

//_________________________________________________________________________________
TrackFitUtils::line_fit_output_t TrackFitUtils::line_fit(const std::vector<Acts::Vector3>& positions)
{
  position_vector_t positions_2d;
  for (const auto& position : positions)
  {
    positions_2d.emplace_back(std::sqrt(square(position.x()) + square(position.y())), position.z());
  }

  return line_fit(positions_2d);
}

//_________________________________________________________________________________
TrackFitUtils::circle_circle_intersection_output_t TrackFitUtils::circle_circle_intersection(double r1, double r2, double x2, double y2)
{
  const double D = square(r1) - square(r2) + square(x2) + square(y2);
  const double a = 1.0 + square(x2 / y2);
  const double b = -D * x2 / square(y2);
  const double c = square(D / (2. * y2)) - square(r1);
  const double delta = square(b) - 4. * a * c;

  const double sqdelta = std::sqrt(delta);

  const double xplus = (-b + sqdelta) / (2. * a);
  const double xminus = (-b - sqdelta) / (2. * a);

  const double yplus = -(2 * x2 * xplus - D) / (2. * y2);
  const double yminus = -(2 * x2 * xminus - D) / (2. * y2);

  return std::make_tuple(xplus, yplus, xminus, yminus);
}

unsigned int TrackFitUtils::addClustersOnLine(TrackFitUtils::line_fit_output_t& fitpars,
                                              const bool& isXY,
                                              const double& dca_cut,
                                              ActsGeometry* tGeometry,
                                              TrkrClusterContainer* clusterContainer,
                                              std::vector<Acts::Vector3>& global_vec,
                                              std::vector<TrkrDefs::cluskey>& cluskey_vec,
                                              const unsigned int& startLayer,
                                              const unsigned int& endLayer)
{
  float slope = std::get<0>(fitpars);
  float intercept = std::get<1>(fitpars);

  int nclusters = 0;
  std::set<TrkrDefs::cluskey> keys_to_add;
  std::set<TrkrDefs::TrkrId> detectors = {TrkrDefs::TrkrId::mvtxId,
                                          TrkrDefs::TrkrId::inttId,
                                          TrkrDefs::TrkrId::tpcId,
                                          TrkrDefs::TrkrId::micromegasId};

  if (startLayer > 2)
  {
    detectors.erase(TrkrDefs::TrkrId::mvtxId);
    if (startLayer > 6)
    {
      detectors.erase(TrkrDefs::TrkrId::inttId);
      if (startLayer > 54)
      {
        detectors.erase(TrkrDefs::TrkrId::tpcId);
      }
    }
  }
  if (endLayer < 56)
  {
    detectors.erase(TrkrDefs::TrkrId::micromegasId);
    if (endLayer < 7)
    {
      detectors.erase(TrkrDefs::TrkrId::tpcId);
      if (endLayer < 3)
      {
        detectors.erase(TrkrDefs::TrkrId::inttId);
      }
    }
  }
  for (const auto& det : detectors)
  {
    for (const auto& hitsetkey : clusterContainer->getHitSetKeys(det))
    {
      if (TrkrDefs::getLayer(hitsetkey) < startLayer ||
          TrkrDefs::getLayer(hitsetkey) > endLayer)
      {
        continue;
      }

      auto range = clusterContainer->getClusters(hitsetkey);
      for (auto clusIter = range.first; clusIter != range.second; ++clusIter)
      {
        TrkrDefs::cluskey cluskey = clusIter->first;

        TrkrCluster* cluster = clusIter->second;

        auto global = tGeometry->getGlobalPosition(cluskey, cluster);
        float x, y;
        if (isXY)
        {
          x = global.x();
          y = global.y();
        }
        else
        {
          //! use r-z
          x = global.z();
          y = std::sqrt(square(global.x()) + square(global.y()));
          if (global.y() < 0)
          {
            y *= -1;
          }
        }

        //! Need to find the point on the line closest to the cluster position
        //! The line connecting the cluster and the line is necessarily
        //! perpendicular, so it has the opposite and reciprocal slope
        //! The cluster is on the line so we use it to find the intercept
        float perpSlope = -1 / slope;
        float perpIntercept = y - perpSlope * x;

        //! Now we have two lines, so solve the system of eqns to get the intersection
        float pcax = (perpIntercept - intercept) / (slope - perpSlope);
        float pcay = slope * pcax + intercept;
        //! Take the difference to find the distance

        float dcax = pcax - x;
        float dcay = pcay - y;
        float dca = std::sqrt(square(dcax) + square(dcay));

        if (dca < dca_cut)
        {
          keys_to_add.insert(cluskey);
        }
      }
    }
  }

  for (auto& key : keys_to_add)
  {
    cluskey_vec.push_back(key);
    auto clus = clusterContainer->findCluster(key);
    auto global = tGeometry->getGlobalPosition(key, clus);
    global_vec.push_back(global);
    nclusters++;
  }

  return nclusters;
}

//_________________________________________________________________________________
unsigned int TrackFitUtils::addClusters(std::vector<float>& fitpars,
                                        double dca_cut,
                                        ActsGeometry* _tGeometry,
                                        TrkrClusterContainer* _cluster_map,
                                        std::vector<Acts::Vector3>& global_vec,
                                        std::vector<TrkrDefs::cluskey>& cluskey_vec,
                                        unsigned int startLayer,
                                        unsigned int endLayer)
{
  // project the fit of the TPC clusters to each silicon layer, and find the nearest silicon cluster
  // iterate over the cluster map and find silicon clusters that match this track fit

  unsigned int nclusters = 0;

  // We want the best match in each layer
  std::set<TrkrDefs::cluskey> keysToAdd;
  std::set<TrkrDefs::TrkrId> detectors = {TrkrDefs::TrkrId::mvtxId,
                                          TrkrDefs::TrkrId::inttId,
                                          TrkrDefs::TrkrId::tpcId,
                                          TrkrDefs::TrkrId::micromegasId};

  if (startLayer > 2)
  {
    detectors.erase(TrkrDefs::TrkrId::mvtxId);
    if (startLayer > 6)
    {
      detectors.erase(TrkrDefs::TrkrId::inttId);
      if (startLayer > 54)
      {
        detectors.erase(TrkrDefs::TrkrId::tpcId);
      }
    }
  }
  if (endLayer < 56)
  {
    detectors.erase(TrkrDefs::TrkrId::micromegasId);
    if (endLayer < 7)
    {
      detectors.erase(TrkrDefs::TrkrId::tpcId);
      if (endLayer < 3)
      {
        detectors.erase(TrkrDefs::TrkrId::inttId);
      }
    }
  }

  for (const auto& det : detectors)
  {
    for (const auto& hitsetkey : _cluster_map->getHitSetKeys(det))
    {
      if (TrkrDefs::getLayer(hitsetkey) < startLayer ||
          TrkrDefs::getLayer(hitsetkey) > endLayer)
      {
        continue;
      }

      auto range = _cluster_map->getClusters(hitsetkey);
      for (auto clusIter = range.first; clusIter != range.second; ++clusIter)
      {
        TrkrDefs::cluskey cluskey = clusIter->first;

        TrkrCluster* cluster = clusIter->second;

        auto global = _tGeometry->getGlobalPosition(cluskey, cluster);

        Acts::Vector3 pca = get_helix_pca(fitpars, global);
        float dca = (pca - global).norm();

        Acts::Vector2 global_xy(global(0), global(1));
        Acts::Vector2 pca_xy(pca(0), pca(1));
        Acts::Vector2 pca_xy_residual = pca_xy - global_xy;
        dca = pca_xy_residual.norm();
        if (dca < dca_cut)
        {
          keysToAdd.insert(cluskey);
        }
      }  // end cluster iteration
    }    // end hitsetkey iteration
  }

  for (auto& key : keysToAdd)
  {
    cluskey_vec.push_back(key);
    auto clus = _cluster_map->findCluster(key);
    auto global = _tGeometry->getGlobalPosition(key, clus);
    global_vec.push_back(global);
    nclusters++;
  }

  return nclusters;
}

//_________________________________________________________________________________
Acts::Vector3 TrackFitUtils::get_helix_pca(std::vector<float>& fitpars,
                                           const Acts::Vector3& global)
{
  // no analytic solution for the coordinates of the closest approach of a helix to a point
  // Instead, we get the PCA in x and y to the circle, and the PCA in z to the z vs R line at the R of the PCA

  float radius = fitpars[0];
  float x0 = fitpars[1];
  float y0 = fitpars[2];
  float zslope = fitpars[3];
  float z0 = fitpars[4];

  Acts::Vector2 pca_circle = get_circle_point_pca(radius, x0, y0, global);

  // The radius of the PCA determines the z position:
  float pca_circle_radius = pca_circle.norm();
  float pca_z = pca_circle_radius * zslope + z0;
  Acts::Vector3 pca(pca_circle(0), pca_circle(1), pca_z);

  // now we want a second point on the helix so we can get a local straight line approximation to the track
  // project the circle PCA vector an additional small amount and find the helix PCA to that point
  float projection = 0.25;  // cm
  Acts::Vector3 second_point = pca + projection * pca / pca.norm();
  Acts::Vector2 second_point_pca_circle = get_circle_point_pca(radius, x0, y0, second_point);
  float second_point_pca_z = pca_circle_radius * zslope + z0;
  Acts::Vector3 second_point_pca(second_point_pca_circle(0), second_point_pca_circle(1), second_point_pca_z);

  // pca and second_point_pca define a straight line approximation to the track
  Acts::Vector3 tangent = (second_point_pca - pca) / (second_point_pca - pca).norm();

  // get the PCA of the cluster to that line
  Acts::Vector3 final_pca = getPCALinePoint(global, tangent, pca);

  return final_pca;
}

//_________________________________________________________________________________
Acts::Vector2 TrackFitUtils::get_circle_point_pca(float radius, float x0, float y0, Acts::Vector3 global)
{
  // get the PCA of a cluster (x,y) position to a circle
  // draw a line from the origin of the circle to the point
  // the intersection of the line with the circle is at the distance radius from the origin along that line

  Acts::Vector2 origin(x0, y0);
  Acts::Vector2 point(global(0), global(1));

  Acts::Vector2 pca = origin + radius * (point - origin) / (point - origin).norm();

  return pca;
}

//_________________________________________________________________________________
std::vector<float> TrackFitUtils::fitClusters(std::vector<Acts::Vector3>& global_vec,
                                              std::vector<TrkrDefs::cluskey> cluskey_vec)
{
  std::vector<float> fitpars;

  // make the helical fit using TrackFitUtils
  if (global_vec.size() < 3)
  {
    return fitpars;
  }
  std::tuple<double, double, double> circle_fit_pars = TrackFitUtils::circle_fit_by_taubin(global_vec);

  // It is problematic that the large errors on the INTT strip z values are not allowed for - drop the INTT from the z line fit
  std::vector<Acts::Vector3> global_vec_noINTT;
  for (unsigned int ivec = 0; ivec < global_vec.size(); ++ivec)
  {
    unsigned int trkrid = TrkrDefs::getTrkrId(cluskey_vec[ivec]);

    if (trkrid != TrkrDefs::inttId and cluskey_vec[ivec] != 0)
    {
      global_vec_noINTT.push_back(global_vec[ivec]);
    }
  }

  if (global_vec_noINTT.size() < 3)
  {
    return fitpars;
  }
  std::tuple<double, double> line_fit_pars = TrackFitUtils::line_fit(global_vec_noINTT);

  fitpars.push_back(std::get<0>(circle_fit_pars));
  fitpars.push_back(std::get<1>(circle_fit_pars));
  fitpars.push_back(std::get<2>(circle_fit_pars));
  fitpars.push_back(std::get<0>(line_fit_pars));
  fitpars.push_back(std::get<1>(line_fit_pars));

  return fitpars;
}

//_________________________________________________________________________________
void TrackFitUtils::getTrackletClusters(ActsGeometry* _tGeometry,
                                        TrkrClusterContainer* _cluster_map,
                                        std::vector<Acts::Vector3>& global_vec,
                                        std::vector<TrkrDefs::cluskey>& cluskey_vec)
{
  for (unsigned long key : cluskey_vec)
  {
    auto cluster = _cluster_map->findCluster(key);
    if (!cluster)
    {
      std::cout << "Failed to get cluster with key " << key << std::endl;
      continue;
    }

    Acts::Vector3 global = _tGeometry->getGlobalPosition(key, cluster);

    /*
    const unsigned int trkrId = TrkrDefs::getTrkrId(key);
    // have to add corrections for TPC clusters after transformation to global
    if(trkrId == TrkrDefs::tpcId)
      {
        int crossing = 0;  // for now
        makeTpcGlobalCorrections(key, crossing, global);
      }
    */

    // add the global positions to a vector to return
    global_vec.push_back(global);

  }  // end loop over clusters for this track
}

//_________________________________________________________________________________
Acts::Vector3 TrackFitUtils::getPCALinePoint(const Acts::Vector3& global, const Acts::Vector3& tangent, const Acts::Vector3& posref)
{
  // Approximate track with a straight line consisting of the state position posref and the vector (px,py,pz)

  // The position of the closest point on the line to global is:
  // posref + projection of difference between the point and posref on the tangent vector
  Acts::Vector3 pca = posref + ((global - posref).dot(tangent)) * tangent;
  /*
  if( (pca-posref).norm() > 0.001)
    {
      std::cout << " getPCALinePoint: old pca " << posref(0) << "  " << posref(1) << "  " << posref(2) << std::endl;
      std::cout << " getPCALinePoint: new pca " << pca(0) << "  " << pca(1) << "  " << pca(2) << std::endl;
      std::cout << " getPCALinePoint: delta pca " << pca(0) - posref(0) << "  " << pca(1)-posref(1) << "  " << pca(2) -posref(2) << std::endl;
    }
  */

  return pca;
}

Acts::Vector3 TrackFitUtils::get_helix_surface_intersection(const Surface& surf,
                                                            std::vector<float>& fitpars, Acts::Vector3& global, ActsGeometry* _tGeometry)
{
  // we want the point where the helix intersects the plane of the surface
  // get the plane of the surface
  Acts::Vector3 sensorCenter = surf->center(_tGeometry->geometry().getGeoContext()) * 0.1;  // convert to cm
  Acts::Vector3 sensorNormal = -surf->normal(_tGeometry->geometry().getGeoContext());
  sensorNormal /= sensorNormal.norm();

  // there are analytic solutions for a line-plane intersection.
  // to use this, need to get the vector tangent to the helix near the measurement and a point on it.
  std::pair<Acts::Vector3, Acts::Vector3> line = get_helix_tangent(fitpars, global);
  Acts::Vector3 pca = line.first;
  Acts::Vector3 tangent = line.second;

  Acts::Vector3 intersection = get_line_plane_intersection(pca, tangent, sensorCenter, sensorNormal);

  return intersection;
}
Acts::Vector3 TrackFitUtils::get_line_plane_intersection(const Acts::Vector3& pca, const Acts::Vector3& tangent,
                                                         const Acts::Vector3& sensorCenter, const Acts::Vector3& sensorNormal)
{
  // get the intersection of the line made by PCA and tangent with the plane of the sensor

  // For a point on the line
  // p = PCA + d * tangent;
  // for a point on the plane
  // (p - sensor_center).sensor_normal = 0

  // The solution is:
  float d = (sensorCenter - pca).dot(sensorNormal) / tangent.dot(sensorNormal);
  return pca + d * tangent;
}
std::vector<double> TrackFitUtils::getLineClusterResiduals(TrackFitUtils::position_vector_t& rz_pts, float slope, float intercept)
{
  std::vector<double> residuals;
  // calculate cluster residuals from the fitted circle
  std::transform(rz_pts.begin(), rz_pts.end(),
                 std::back_inserter(residuals), [slope, intercept](const std::pair<double, double>& point)
                 {
    double r = point.first;
    double z = point.second;
    
    // The shortest distance of a point from a circle is along the radial; line from the circle center to the point
    
    double a = -slope;
    double b = 1.0;
    double c = -intercept;
    return std::abs(a*r+b*z+c)/sqrt(square(a)+square(b)); });
  return residuals;
}

std::vector<double> TrackFitUtils::getCircleClusterResiduals(TrackFitUtils::position_vector_t& xy_pts, float R, float X0, float Y0)
{
  std::vector<double> residuals;
  std::transform(xy_pts.begin(), xy_pts.end(),
                 std::back_inserter(residuals), [R, X0, Y0](const std::pair<double, double>& point)
                 {
    double x = point.first;
    double y = point.second;

    // The shortest distance of a point from a circle is along the radial; line from the circle center to the point
    return std::sqrt( square(x-X0) + square(y-Y0) )  -  R; });

  return residuals;
}
