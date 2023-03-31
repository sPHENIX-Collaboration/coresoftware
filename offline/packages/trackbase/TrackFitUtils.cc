#include "TrackFitUtils.h"

#include <trackbase_historic/TrackSeed_v1.h>
#include "ActsGeometry.h"
#include "TrkrDefs.h"                // for cluskey, getTrkrId, tpcId
#include "TpcDefs.h"
#include <trackbase/MvtxDefs.h>
#include "TrkrClusterv5.h"
#include "TrkrClusterContainerv4.h"

#include <cmath>

namespace
{

  //! convenience square method
  template<class T>
    inline constexpr T square( const T& x ) { return x*x; }
}

//_________________________________________________________________________________
TrackFitUtils::circle_fit_output_t TrackFitUtils::circle_fit_by_taubin( const TrackFitUtils::position_vector_t& positions )
{

  // Compute x- and y- sample means
  double meanX = 0;
  double meanY = 0;
  double weight = 0;

  for( const auto& [x,y] : positions)
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

  for (auto& [x,y] : positions)
  {
    double Xi = x - meanX;   //  centered x-coordinates
    double Yi = y - meanY;   //  centered y-coordinates
    double Zi = square(Xi) + square(Yi);

    Mxy += Xi*Yi;
    Mxx += Xi*Xi;
    Myy += Yi*Yi;
    Mxz += Xi*Zi;
    Myz += Yi*Zi;
    Mzz += Zi*Zi;
  }
  Mxx /= weight;
  Myy /= weight;
  Mxy /= weight;
  Mxz /= weight;
  Myz /= weight;
  Mzz /= weight;

  //  computing coefficients of the characteristic polynomial

  const double Mz = Mxx + Myy;
  const double Cov_xy = Mxx*Myy - Mxy*Mxy;
  const double Var_z = Mzz - Mz*Mz;
  const double A3 = 4*Mz;
  const double A2 = -3*Mz*Mz - Mzz;
  const double A1 = Var_z*Mz + 4*Cov_xy*Mz - Mxz*Mxz - Myz*Myz;
  const double A0 = Mxz*(Mxz*Myy - Myz*Mxy) + Myz*(Myz*Mxx - Mxz*Mxy) - Var_z*Cov_xy;
  const double A22 = A2 + A2;
  const double A33 = A3 + A3 + A3;

  //    finding the root of the characteristic polynomial
  //    using Newton's method starting at x=0
  //    (it is guaranteed to converge to the right root)
  static constexpr int iter_max = 99;
  double x = 0;
  double y = A0;

  // usually, 4-6 iterations are enough
  for( int iter=0; iter<iter_max; ++iter)
  {
    const double Dy = A1 + x*(A22 + A33*x);
    const double xnew = x - y/Dy;
    if ((xnew == x)||(!std::isfinite(xnew))) break;

    const double ynew = A0 + xnew*(A1 + xnew*(A2 + xnew*A3));
    if (std::abs(ynew)>=std::abs(y))  break;

    x = xnew;  y = ynew;

  }

  //  computing parameters of the fitting circle
  const double DET = square(x) - x*Mz + Cov_xy;
  const double Xcenter = (Mxz*(Myy - x) - Myz*Mxy)/DET/2;
  const double Ycenter = (Myz*(Mxx - x) - Mxz*Mxy)/DET/2;

  //  assembling the output
  double X0 = Xcenter + meanX;
  double Y0 = Ycenter + meanY;
  double R = std::sqrt( square(Xcenter) + square(Ycenter) + Mz);
  return std::make_tuple( R, X0, Y0 );

}

//_________________________________________________________________________________
TrackFitUtils::circle_fit_output_t TrackFitUtils::circle_fit_by_taubin( const std::vector<Acts::Vector3>& positions )
{
  position_vector_t positions_2d;
  for( const auto& position:positions )
  { positions_2d.emplace_back( position.x(), position.y() ); }

  return circle_fit_by_taubin( positions_2d );
}

//_________________________________________________________________________________
TrackFitUtils::line_fit_output_t TrackFitUtils::line_fit(const TrackFitUtils::position_vector_t& positions)
{
  double xsum=0;
  double x2sum=0;
  double ysum=0;
  double xysum=0;
  for( const auto& [r,z]:positions )
  {
    xsum=xsum+r;                        //calculate sigma(xi)
    ysum=ysum+z;                        //calculate sigma(yi)
    x2sum=x2sum+square(r);              //calculate sigma(x^2i)
    xysum=xysum+r*z;                    //calculate sigma(xi*yi)
  }

  const auto npts = positions.size();
  const double denominator = (x2sum*npts-square(xsum));
  const double a= (xysum*npts-xsum*ysum)/denominator;            //calculate slope
  const double b= (x2sum*ysum-xsum*xysum)/denominator;           //calculate intercept
  return std::make_tuple( a, b );
}

//_________________________________________________________________________________
TrackFitUtils::line_fit_output_t TrackFitUtils::line_fit( const std::vector<Acts::Vector3>& positions )
{
  position_vector_t positions_2d;
  for( const auto& position:positions )
  { positions_2d.emplace_back( std::sqrt( square(position.x())+square(position.y())), position.z() ); }

  return line_fit( positions_2d );
}

//_________________________________________________________________________________
TrackFitUtils::circle_circle_intersection_output_t TrackFitUtils::circle_circle_intersection(double r1, double r2, double x2, double y2 )
{
  
  const double D = square(r1) - square(r2) + square(x2) + square(y2);
  const double a = 1.0 + square(x2/y2);
  const double b = - D * x2/square(y2);
  const double c = square(D/(2.*y2)) - square(r1);
  const double delta = square(b)-4.*a*c;
  
  const double sqdelta = std::sqrt( delta );
  
  const double xplus = (-b + sqdelta ) / (2. * a);
  const double xminus = (-b - sqdelta ) / (2. * a);
  
  const double yplus  = -(2*x2*xplus - D) / (2.*y2);
  const double yminus = -(2*x2*xminus - D) / (2.*y2);

  
  return std::make_tuple( xplus, yplus, xminus, yminus );

}

//_________________________________________________________________________________
unsigned int TrackFitUtils::addSiliconClusters(std::vector<float>& fitpars, 
					       double dca_cut,
					       ActsGeometry* _tGeometry, 
					       TrkrClusterContainer * _cluster_map,
					       std::vector<Acts::Vector3>& global_vec,  
					       std::vector<TrkrDefs::cluskey>& cluskey_vec)
{
  // project the fit of the TPC clusters to each silicon layer, and find the nearest silicon cluster
  // iterate over the cluster map and find silicon clusters that match this track fit

  unsigned int nsilicon = 0;

  // We want the best match in each layer
  std::vector<float> best_layer_dca;
  best_layer_dca.assign(7, 999.0);
  std::vector<TrkrDefs::cluskey> best_layer_cluskey;
  best_layer_cluskey.assign(7, 0);

  for(const auto& hitsetkey:_cluster_map->getHitSetKeys())
    {
      auto range = _cluster_map->getClusters(hitsetkey);
      for( auto clusIter = range.first; clusIter != range.second; ++clusIter )
	{
	  TrkrDefs::cluskey cluskey = clusIter->first;
	  unsigned int layer = TrkrDefs::getLayer(cluskey);
	  unsigned int trkrid = TrkrDefs::getTrkrId(cluskey);
	  
	  if(trkrid != TrkrDefs::mvtxId && trkrid != TrkrDefs::inttId)  continue;
	  
	  TrkrCluster* cluster = clusIter->second;
	  auto global = _tGeometry->getGlobalPosition(cluskey, cluster);

	  Acts::Vector3 pca = get_helix_pca(fitpars, global);
	  float dca = (pca - global).norm();
	  if(trkrid == TrkrDefs::inttId || trkrid == TrkrDefs::mvtxId)
	    {
	      Acts::Vector2 global_xy(global(0), global(1));
	      Acts::Vector2 pca_xy(pca(0), pca(1));
	      Acts::Vector2 pca_xy_residual = pca_xy - global_xy;
	      dca = pca_xy_residual.norm();
	    }

	  if(dca < best_layer_dca[layer])
	    {
	      best_layer_dca[layer] = dca;
	      best_layer_cluskey[layer] = cluskey;
	    }
	}  // end cluster iteration
    } // end hitsetkey iteration

  for(unsigned int layer = 0; layer < 7; ++layer)
    {
      if(best_layer_dca[layer] < dca_cut)
	{
	  if(best_layer_cluskey[layer] != 0)
          {
            cluskey_vec.push_back(best_layer_cluskey[layer]);
	    auto clus =  _cluster_map->findCluster(best_layer_cluskey[layer]);
	    auto global = _tGeometry->getGlobalPosition(best_layer_cluskey[layer], clus);
	    global_vec.push_back(global);
	    nsilicon++;
          }
	}
    }

  return nsilicon;
}

//_________________________________________________________________________________
Acts::Vector3 TrackFitUtils::get_helix_pca(std::vector<float>& fitpars, 
					   Acts::Vector3 global)
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
  Acts::Vector3 second_point = pca + projection * pca/pca.norm();
  Acts::Vector2 second_point_pca_circle = get_circle_point_pca(radius, x0, y0, second_point);
  float second_point_pca_z = pca_circle_radius * zslope + z0;
  Acts::Vector3 second_point_pca(second_point_pca_circle(0), second_point_pca_circle(1), second_point_pca_z);

  // pca and second_point_pca define a straight line approximation to the track
  Acts::Vector3 tangent = (second_point_pca - pca) /  (second_point_pca - pca).norm();

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
      if(global_vec.size() < 3)  
	{ return fitpars; }
      std::tuple<double, double, double> circle_fit_pars = TrackFitUtils::circle_fit_by_taubin(global_vec);

      // It is problematic that the large errors on the INTT strip z values are not allowed for - drop the INTT from the z line fit
      std::vector<Acts::Vector3> global_vec_noINTT;
      for(unsigned int ivec=0;ivec<global_vec.size(); ++ivec)
	{
	  unsigned int trkrid = TrkrDefs::getTrkrId(cluskey_vec[ivec]);
	  if(trkrid != TrkrDefs::inttId) { global_vec_noINTT.push_back(global_vec[ivec]); }
	}      
      if(global_vec_noINTT.size() < 3) 
	{ return fitpars; }
     std::tuple<double,double> line_fit_pars = TrackFitUtils::line_fit(global_vec_noINTT);

     fitpars.push_back( std::get<0>(circle_fit_pars));
     fitpars.push_back( std::get<1>(circle_fit_pars));
     fitpars.push_back( std::get<2>(circle_fit_pars));
     fitpars.push_back( std::get<0>(line_fit_pars));
     fitpars.push_back( std::get<1>(line_fit_pars));

     return fitpars; 
}

//_________________________________________________________________________________
void TrackFitUtils::getTrackletClusters( ActsGeometry * _tGeometry, 
					TrkrClusterContainer * _cluster_map,
					std::vector<Acts::Vector3>& global_vec, 
					std::vector<TrkrDefs::cluskey>& cluskey_vec)
{
  for (unsigned int iclus=0;  iclus < cluskey_vec.size(); ++iclus)
    {
      auto key = cluskey_vec[iclus];
      auto cluster = _cluster_map->findCluster(key);
      if(!cluster)
	{
	  std::cout << "Failed to get cluster with key " << key << std::endl;
	  continue;
	}	  
            
      Acts::Vector3 global  = _tGeometry->getGlobalPosition(key, cluster);	  
      
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
      
    } // end loop over clusters for this track 
}

//_________________________________________________________________________________
Acts::Vector3 TrackFitUtils::getPCALinePoint(Acts::Vector3 global, Acts::Vector3 tangent, Acts::Vector3 posref)
{
  // Approximate track with a straight line consisting of the state position posref and the vector (px,py,pz)   

  // The position of the closest point on the line to global is:
  // posref + projection of difference between the point and posref on the tangent vector
  Acts::Vector3 pca = posref + ( (global - posref).dot(tangent) ) * tangent;

  return pca;
}

