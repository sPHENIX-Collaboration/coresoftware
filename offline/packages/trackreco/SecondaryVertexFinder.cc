#include "SecondaryVertexFinder.h"

/// Tracking includes

#include <trackbase/TrkrCluster.h>            // for TrkrCluster
#include <trackbase/TrkrDefs.h>               // for cluskey, getLayer, TrkrId
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TpcDefs.h>

#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack, SvtxTrack::C...
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex_v1.h>
#include <trackbase_historic/SvtxVertexMap_v1.h>
#include <trackbase_historic/TrackSeed_v1.h>
#include <trackbase_historic/ActsTransformations.h>

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>

#include <phool/PHCompositeNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <cmath>                              // for sqrt, fabs, atan2, cos
#include <iostream>                           // for operator<<, basic_ostream
#include <map>                                // for map
#include <set>                                // for _Rb_tree_const_iterator
#include <utility>                            // for pair, make_pair
#include <iomanip>

#include <vector>
#include <cassert>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <functional>

#include <Eigen/Dense>

//____________________________________________________________________________..
SecondaryVertexFinder::SecondaryVertexFinder(const std::string &name)
  : SubsysReco(name)
{

}

//____________________________________________________________________________..
SecondaryVertexFinder::~SecondaryVertexFinder()
{

}

//____________________________________________________________________________..
int SecondaryVertexFinder::InitRun(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return ret;
}

//____________________________________________________________________________..
int SecondaryVertexFinder::process_event(PHCompositeNode */*topNode*/)
{

  if(Verbosity() > 0)
    std::cout << PHWHERE << " track map size " << _track_map->size()  << std::endl;

  // reset maps
  _vertex_track_map.clear();
  _track_pair_map.clear();
  _track_pair_pca_map.clear();
  _vertex_position_map.clear();
  _vertex_covariance_map.clear();
  _vertex_set.clear();
  
  // Loop over tracks and check for close DCA match with all other tracks
  for(auto tr1_it = _track_map->begin(); tr1_it != _track_map->end(); ++tr1_it)
    {
      auto id1 = tr1_it->first;
      auto tr1 = tr1_it->second;
      if(tr1->get_quality() > _qual_cut) continue;

      // look for close DCA matches with all other such tracks
      for(auto tr2_it = std::next(tr1_it); tr2_it != _track_map->end(); ++tr2_it)
	{
	  auto id2 = tr2_it->first;
	  auto tr2 = tr2_it->second;
	  if(tr2->get_quality() > _qual_cut) continue;

	  // find DCA and PCA of these two tracks
	  if(Verbosity() > 3) std::cout << "Check DCA for tracks " << id1 << " and  " << id2 << std::endl;
	  Eigen::Vector3d PCA1(0,0,0), PCA2(0,0,0);
	  double dca = findTwoTrackPCA(tr1, tr2, PCA1, PCA2);

	  // check for close pair
	  if(dca > _two_track_dcacut) { continue; }

	  std::cout << "Found a decay pair: id1" << tr1_it->first << " id2 " << tr1_it->first 
		    << " PCA1 " << PCA1(0) << "  " << PCA1(1) << "  " << PCA1(2)
		    << " PCA2 " << PCA2(0) << "  " << PCA2(1) << "  " << PCA2(2)
		    << " pair dca " << dca << std::endl;  
	}

    }

  return Fun4AllReturnCodes::EVENT_OK;
}

void SecondaryVertexFinder::updateSvtxTrack(SvtxTrack* track, 
					     const Acts::BoundTrackParameters& params)
{
  auto position = params.position(_tGeometry->geometry().getGeoContext());
  
  if(Verbosity() > 2)
    {
      std::cout << "Updating position track parameters from " << track->get_x()
		<< ", " << track->get_y() << ", " << track->get_z() << " to " 
		<< position.transpose() / 10.
		<< std::endl;
    }

  track->set_x(position(0) / Acts::UnitConstants::cm);
  track->set_y(position(1) / Acts::UnitConstants::cm);
  track->set_z(position(2) / Acts::UnitConstants::cm);

  ActsTransformations rotater;
  rotater.setVerbosity(Verbosity());
  if(params.covariance())
    {
      auto rotatedCov = rotater.rotateActsCovToSvtxTrack(params);
    
      /// Update covariance
      for(int i = 0; i < 3; i++) {
	for(int j = 0; j < 3; j++) {
	  track->set_error(i,j, rotatedCov(i,j));
	}
      }
    }
}

BoundTrackParamPtrResult SecondaryVertexFinder::propagateTrack(
		         const Acts::BoundTrackParameters& params,
			 const Eigen::Vector3d PCA)
{
  
  /// create perigee surface
  auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(PCA);

  using Stepper = Acts::EigenStepper<>;
  using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
  
  Stepper stepper(_tGeometry->geometry().magField);
  Acts::Navigator::Config cfg{_tGeometry->geometry().tGeometry};
  Acts::Navigator navigator(cfg);
  Propagator propagator(stepper, navigator);
  
  Acts::Logging::Level logLevel = Acts::Logging::INFO;
  if(Verbosity() > 3)
    { logLevel = Acts::Logging::VERBOSE; }
  
  auto logger = Acts::getDefaultLogger("PHActsVertexPropagator", 
				       logLevel);

  Acts::PropagatorOptions<> options(_tGeometry->geometry().getGeoContext(),
				    _tGeometry->geometry().magFieldContext,
				    Acts::LoggerWrapper{*logger});
  
  auto result = propagator.propagate(params, *perigee, 
				     options);
  if(result.ok())
    { return std::move((*result).endParameters); }
  
  return result.error();

}

double SecondaryVertexFinder::findTwoTrackPCA(SvtxTrack *track1, SvtxTrack *track2, Eigen::Vector3d &PCA1, Eigen::Vector3d &PCA2)
{
  // For secondary vertex finding we cannot assume that the vertex is close to the beam line
  // we start by fitting circles to the TPC clusters and determining the circle-circle intersection - roughly the decay vertex
  // then we project the tracks to that point and get the momentum vector there
  // then we call the line-line DCA/PCA method to get the precise result

  double dca = 999.0;

  TrackSeed *tr1 = track1->get_tpc_seed();
  TrackSeed *tr2 = track2->get_tpc_seed();

  std::vector<float> circle_fitpars1 = fitClusters(tr1);
  if(circle_fitpars1.size() == 0) return dca;  // discard this track, not enough clusters to fit
  std::vector<float> circle_fitpars2 = fitClusters(tr2);
  if(circle_fitpars2.size() == 0) return dca;  // discard this track, not enough clusters to fit

  // get intersection point for these two tracks in x,y plane
  std::vector<double> intersections =  circle_circle_intersection(
								  circle_fitpars1[0], circle_fitpars1[1], circle_fitpars1[2],
								  circle_fitpars2[0], circle_fitpars2[1], circle_fitpars2[2] );

  // which intersection solution is the one we want?
  // The correct intersection is the one in the direction of the track vector
  double x0 = intersections[0];
  double y0 = intersections[1];
  double x1 = intersections[2];
  double y1 = intersections[3];

  Eigen::Vector2d intersect(0,0);
  // check the direction relative to the beam line
  if( (x0 - track1->get_x()) / track1->get_px() > 0 && (y0-track1->get_y()) / track1->get_py() > 0)  
    {
      intersect(0) = x0;
      intersect(1) = y0;
    }
  else  
    {
      intersect(0) = x1;
      intersect(1) = y1;
    }

  // Project both tracks to this point and get the momentum vectors

  // track 1  
  auto boundParams1 = makeTrackParams(track1);
  auto propresult1 = propagateTrack(boundParams1, PCA1);
  if(propresult1.ok())
    {
      auto paramsAtVertex = std::move(**propresult1);
      updateSvtxTrack(track1,paramsAtVertex);
    }

  // track 2  
  auto boundParams2 = makeTrackParams(track2);
  auto propresult2 = propagateTrack(boundParams2, PCA2);
  if(propresult2.ok())
    {
      auto paramsAtVertex = std::move(**propresult2);
      updateSvtxTrack(track2,paramsAtVertex);
    }
  
  // Now make the final two-track PCA determination using the updated tracks
  Eigen::Vector3d a1(track1->get_x(), track1->get_y(), track1->get_z());
  Eigen::Vector3d b1(track1->get_px() / track1->get_p(), track1->get_py() / track1->get_p(), track1->get_pz() / track1->get_p());
  Eigen::Vector3d a2(track2->get_x(), track2->get_y(), track2->get_z());
  Eigen::Vector3d b2(track2->get_px() / track2->get_p(), track2->get_py() / track2->get_p(), track2->get_pz() / track2->get_p());

  dca =  dcaTwoLines(a1,b1,a2, b2, PCA1, PCA2);

  return dca;
}

 Acts::Vector3 SecondaryVertexFinder::getVertex(SvtxTrack* track)
{
  auto vertexId = track->get_vertex_id();
  const SvtxVertex* svtxVertex = _svtx_vertex_map->get(vertexId);
  Acts::Vector3 vertex = Acts::Vector3::Zero();
  if (svtxVertex)
  {
    vertex(0) = svtxVertex->get_x() * Acts::UnitConstants::cm;
    vertex(1) = svtxVertex->get_y() * Acts::UnitConstants::cm;
    vertex(2) = svtxVertex->get_z() * Acts::UnitConstants::cm;
  }

  return vertex;
}

Acts::BoundTrackParameters SecondaryVertexFinder::makeTrackParams(SvtxTrack* track)
{
  Acts::Vector3 momentum(track->get_px(),
                         track->get_py(),
                         track->get_pz());

  auto actsVertex = getVertex(track);
  auto perigee =
      Acts::Surface::makeShared<Acts::PerigeeSurface>(actsVertex);
  auto actsFourPos =
      Acts::Vector4(track->get_x() * Acts::UnitConstants::cm,
                    track->get_y() * Acts::UnitConstants::cm,
                    track->get_z() * Acts::UnitConstants::cm,
                    10 * Acts::UnitConstants::ns);

  ActsTransformations transformer;

  Acts::BoundSymMatrix cov = transformer.rotateSvtxTrackCovToActs(track);

  return ActsExamples::TrackParameters::create(perigee, _tGeometry->geometry().getGeoContext(),
                                               actsFourPos, momentum,
                                               track->get_charge() / track->get_p(),
                                               cov).value();
}

std::vector<double> SecondaryVertexFinder::circle_circle_intersection(double r0, double x0, double y0, double r1, double x1, double y1 )
{
  std::vector<double> ret;
 
  Eigen::Vector2d p0(x0,y0);
  Eigen::Vector2d p1(x1,y1);

  double d = (p0 - p1).norm();

  if(d > r0+r1) return ret;
  if(d < fabs(r1-r0)) return ret;
     
  double a=(r0*r0-r1*r1+d*d)/(2*d);
  double h = sqrt(r0*r0 - a*a);

  double x2= x0 + a*(x1-x0)/d;   
  double y2=y0+a*(y1-y0)/d;

  double x3a=x2+h*(y1-y0)/d;       // also x3=x2-h*(y1-y0)/d
  double y3a=y2-h*(x1-x0)/d;       // also y3=y2+h*(x1-x0)/d

  double x3b=x2-h*(y1-y0)/d;       // also x3=x2-h*(y1-y0)/d
  double y3b=y2+h*(x1-x0)/d;       // also y3=y2+h*(x1-x0)/d

    ret.push_back(x3a);
    ret.push_back(y3a);
    ret.push_back(x3b);
    ret.push_back(y3b);

  
    return ret;

}

double SecondaryVertexFinder::dcaTwoLines(const Eigen::Vector3d &a1,const Eigen::Vector3d &b1,
					 const Eigen::Vector3d &a2,const Eigen::Vector3d &b2,
					 Eigen::Vector3d &PCA1, Eigen::Vector3d &PCA2)
{
  auto bcrossb = b1.cross(b2);
  auto mag_bcrossb = bcrossb.norm();
  // a2-a1 is the vector joining any arbitrary points on the two lines
  auto aminusa = a2-a1;

  // The DCA of these two lines is the projection of a2-a1 along the direction of the perpendicular to both 
  // remember that a2-a1 is longer than (or equal to) the dca by definition
  double dca = 999;
  if( mag_bcrossb != 0)
    dca = bcrossb.dot(aminusa) / mag_bcrossb;
  else
    return dca;   // same track, skip combination
  
  double X =  b1.dot(b2) - b1.dot(b1) * b2.dot(b2) / b2.dot(b1);
  double Y =  (a2.dot(b2) - a1.dot(b2)) - (a2.dot(b1) - a1.dot(b1)) * b2.dot(b2) / b2.dot(b1) ;
  double c = Y/X;

  double F = b1.dot(b1) / b2.dot(b1); 
  double G = - (a2.dot(b1) - a1.dot(b1)) / b2.dot(b1);
  double d = c * F + G;

  // then the points of closest approach are:
  PCA1 = a1+c*b1;
  PCA2 = a2+d*b2;

  return dca;
}

std::vector<float> SecondaryVertexFinder::fitClusters(TrackSeed *tracklet)
{
     std::vector<float> fitpars;

     std::vector<Eigen::Vector3d> global_vec;
     std::vector<TrkrDefs::cluskey> cluskey_vec;
     getTrackletClusters(tracklet, global_vec, cluskey_vec);   // store cluster corrected global positions in a vector
     
      // make the helical fit using TrackFitUtils
      if(global_vec.size() < 3)  
	if(Verbosity() > 0) {  std::cout << " track has too few clusters for circle fit, skip it" << std::endl; return fitpars; }
      std::tuple<double, double, double> circle_fit_pars = TrackFitUtils::circle_fit_by_taubin(global_vec);

      // It is problematic that the large errors on the INTT strip z values are not allowed for - drop the INTT from the z line fit
      std::vector<Eigen::Vector3d> global_vec_noINTT;
      for(unsigned int ivec=0;ivec<global_vec.size(); ++ivec)
	{
	  unsigned int trkrid = TrkrDefs::getTrkrId(cluskey_vec[ivec]);
	  if(trkrid != TrkrDefs::inttId) { global_vec_noINTT.push_back(global_vec[ivec]); }
	}      
      if(global_vec_noINTT.size() < 3) 
	if(Verbosity() > 0) { std::cout << " track has too few non-INTT clusters for z fit, skip it" << std::endl; return fitpars; }
     std::tuple<double,double> line_fit_pars = TrackFitUtils::line_fit(global_vec_noINTT);

     fitpars.push_back( std::get<0>(circle_fit_pars));
     fitpars.push_back( std::get<1>(circle_fit_pars));
     fitpars.push_back( std::get<2>(circle_fit_pars));
     fitpars.push_back( std::get<0>(line_fit_pars));
     fitpars.push_back( std::get<1>(line_fit_pars));

     return fitpars; 
}

void SecondaryVertexFinder::getTrackletClusters(TrackSeed *tracklet, std::vector<Eigen::Vector3d>& global_vec, std::vector<TrkrDefs::cluskey>& cluskey_vec)
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
      if(!surf)
	{ continue; }
      
      Eigen::Vector3d global  = _tGeometry->getGlobalPosition(key, cluster);	  
      
      const unsigned int trkrId = TrkrDefs::getTrkrId(key);	  
      
      // have to add corrections for TPC clusters after transformation to global
      if(trkrId == TrkrDefs::tpcId) 
	{  
	  int crossing = 0;  // for now
	  makeTpcGlobalCorrections(key, crossing, global); 
	}
      
      // add the global positions to a vector to give to the helical fitter
      global_vec.push_back(global);
      cluskey_vec.push_back(key);
      
    } // end loop over clusters for this track 
}

int SecondaryVertexFinder::End(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}


void SecondaryVertexFinder::makeTpcGlobalCorrections(TrkrDefs::cluskey cluster_key, short int crossing, Eigen::Vector3d& global)
{
  // make all corrections to global position of TPC cluster
  unsigned int side = TpcDefs::getSide(cluster_key);
  float z = _clusterCrossingCorrection.correctZ(global[2], side, crossing);
  global[2] = z;
  
  // apply distortion corrections
  if(_dcc_static) { global = _distortionCorrection.get_corrected_position( global, _dcc_static ); }
  if(_dcc_average) { global = _distortionCorrection.get_corrected_position( global, _dcc_average ); }
  if(_dcc_fluctuation) { global = _distortionCorrection.get_corrected_position( global, _dcc_fluctuation ); }
}

int SecondaryVertexFinder::GetNodes(PHCompositeNode* topNode)
{

  _track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!_track_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find SvtxTrackMap: " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
    {
      std::cout << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  _svtx_vertex_map = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMap");  
  if(!_svtx_vertex_map)
    {
      std::cout << PHWHERE << "No vertex node, quit!" << std::endl;
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

/*
void SecondaryVertexFinder::findDcaTwoTracks(SvtxTrack *tr1, SvtxTrack *tr2)
{
  // assumes vertex is very close to beam line - not good for secondary vertex finding 

  if(tr1->get_pt() < _track_pt_cut) return;
  if(tr2->get_pt() < _track_pt_cut) return;

  unsigned int id1 = tr1->get_id();
  unsigned int id2 = tr2->get_id();

  // get the line equations for the tracks
  
  Eigen::Vector3d a1(tr1->get_x(), tr1->get_y(), tr1->get_z());
  Eigen::Vector3d b1(tr1->get_px() / tr1->get_p(), tr1->get_py() / tr1->get_p(), tr1->get_pz() / tr1->get_p());
  Eigen::Vector3d a2(tr2->get_x(), tr2->get_y(), tr2->get_z());
  Eigen::Vector3d b2(tr2->get_px() / tr2->get_p(), tr2->get_py() / tr2->get_p(), tr2->get_pz() / tr2->get_p());
  
  Eigen::Vector3d PCA1(0,0,0);
  Eigen::Vector3d PCA2(0,0,0);  
  double dca = dcaTwoLines(a1, b1, a2,  b2, PCA1, PCA2);

  // check dca cut is satisfied, and that PCA is close to beam line
  if( fabs(dca) < _dcacut && (fabs(PCA1.x()) < _beamline_xy_cut && fabs(PCA1.y()) < _beamline_xy_cut) )
    {
      if(Verbosity() > 3)
	{
	  std::cout << " good match for tracks " << tr1->get_id() << " and " << tr2->get_id() << " with pT " << tr1->get_pt()  << " and " << tr2->get_pt() << std::endl;
	  std::cout << "    a1.x " << a1.x() << " a1.y " << a1.y() << " a1.z " << a1.z() << std::endl;
	  std::cout << "    a2.x  " << a2.x()  << " a2.y " << a2.y() << " a2.z " << a2.z() << std::endl;
	  std::cout << "    PCA1.x() " << PCA1.x() << " PCA1.y " << PCA1.y() << " PCA1.z " << PCA1.z() << std::endl;
	  std::cout << "    PCA2.x() " << PCA2.x() << " PCA2.y " << PCA2.y() << " PCA2.z " << PCA2.z() << std::endl;      
	  std::cout << "    dca " << dca << std::endl;
	}  

      // capture the results for successful matches
      _track_pair_map.insert(std::make_pair(id1,std::make_pair(id2, dca)));
      _track_pair_pca_map.insert( std::make_pair(id1, std::make_pair(id2, std::make_pair(PCA1, PCA2))) );
    }

  return;
}
*/

/*
double SecondaryVertexFinder::dcaTwoLines(const Eigen::Vector3d &a1,const Eigen::Vector3d &b1,
					 const Eigen::Vector3d &a2,const Eigen::Vector3d &b2,
					 Eigen::Vector3d &PCA1, Eigen::Vector3d &PCA2)
{
  // The shortest distance between two skew lines described by
  //  a1 + c * b1
  //  a2 + d * b2
  // where a1, a2, are vectors representing points on the lines, b1, b2 are direction vectors, and c and d are scalars
  // is:
  // dca = (b1 x b2) .(a2-a1) / |b1 x b2|

  // bcrossb/mag_bcrossb is a unit vector perpendicular to both direction vectors b1 and b2
  auto bcrossb = b1.cross(b2);
  auto mag_bcrossb = bcrossb.norm();
  // a2-a1 is the vector joining any arbitrary points on the two lines
  auto aminusa = a2-a1;

  // The DCA of these two lines is the projection of a2-a1 along the direction of the perpendicular to both 
  // remember that a2-a1 is longer than (or equal to) the dca by definition
  double dca = 999;
  if( mag_bcrossb != 0)
    dca = bcrossb.dot(aminusa) / mag_bcrossb;
  else
    return dca;   // same track, skip combination
  
  // get the points at which the normal to the lines intersect the lines, where the lines are perpendicular
  // Assume the shortest distance occurs at points A on line 1, and B on line 2, call the line AB
  //    AB = a2+d*b2 - (a1+c*b1)
  // we need to find c and d where AB is perpendicular to both lines. so AB.b1 = 0 and AB.b2 = 0
  //    ( (a2 -a1) + d*b2 -c*b1 ).b1 = 0
  //    ( (a2 -a1) + d*b2 -c*b1 ).b2 = 0
  // so we have two simultaneous equations in 2 unknowns
  //    (a2-a1).b1 + d*b2.b1 -c*b1.b1 = 0 => d*b2.b1 = c*b1.b1 - (a2-a1).b1 => d = (1/b2.b1) * (c*b1.b1 - (a2-a1).b1) 
  //    (a2-a1).b2 + d*b2.b2 - c*b1.b2 = 0 => c*b1.b2 =  (a2-a1).b2 + [(1/b2.b1) * (c*b1*b1 -(a2-a1).b1)}*b2.b2
  //    c*b1.b2 - (1/b2.b1) * c*b1.b1*b2.b2 = (a2-a1).b2 - (1/b2.b1)*(a2-a1).b1*b2.b2
  //    c *(b1.b2 - (1/b2.b1)*b1.b1*b2.b2)  = (a2-a1).b2 - (1/b2.b1)*(a2-a1).b1*b2.b2
  // call this: c*X = Y
  // plug back into the d equation
  //     d = c*b1.b1 / b2.b1 - (a2-a1).b1 / b2.b1
  // and call the d equation: d = c * F - G 

  double X =  b1.dot(b2) - b1.dot(b1) * b2.dot(b2) / b2.dot(b1);
  double Y =  (a2.dot(b2) - a1.dot(b2)) - (a2.dot(b1) - a1.dot(b1)) * b2.dot(b2) / b2.dot(b1) ;
  double c = Y/X;

  double F = b1.dot(b1) / b2.dot(b1); 
  double G = - (a2.dot(b1) - a1.dot(b1)) / b2.dot(b1);
  double d = c * F + G;

  // then the points of closest approach are:
  PCA1 = a1+c*b1;
  PCA2 = a2+d*b2;

  return dca;


}
*/

/*
std::vector<std::set<unsigned int>> SecondaryVertexFinder::findConnectedTracks()
{
 std::vector<std::set<unsigned int>> connected_tracks;
  std::set<unsigned int> connected;
  std::set<unsigned int> used;
  for(auto it : _track_pair_map)
    {
      unsigned int id1 = it.first;
      unsigned int id2 = it.second.first;

      if( (used.find(id1) != used.end()) && (used.find(id2) != used.end()) )
	{
	  if(Verbosity() > 3) std::cout << " tracks " << id1 << " and " << id2 << " are both in used , skip them" << std::endl;
	  continue;
	}
      else if( (used.find(id1) == used.end()) && (used.find(id2) == used.end()))
	{
	  if(Verbosity() > 3) std::cout << " tracks " << id1 << " and " << id2 << " are both not in used , start a new connected set" << std::endl;
	  // close out and start a new connections set
	  if(connected.size() > 0)
	    {
	      connected_tracks.push_back(connected);
	      connected.clear();
	      if(Verbosity() > 3) std::cout << "           closing out set " << std::endl;
	    }
	}

      // get everything connected to id1 and id2
      connected.insert(id1);
      used.insert(id1);
      connected.insert(id2);
      used.insert(id2);
      for(auto cit :  _track_pair_map)
	{
	  unsigned int id3 = cit.first;
	  unsigned int id4 = cit.second.first;
	  if( (connected.find(id3) != connected.end()) || (connected.find(id4) != connected.end()) )
	    {
	      if(Verbosity() > 3) std::cout << " found connection to " << id3 << " and " << id4 << std::endl;
	      connected.insert(id3);
	      used.insert(id3);
	      connected.insert(id4);
	      used.insert(id4);
	    }
	}
    }
  
  // close out the last set
  if(connected.size() > 0)
    {
      connected_tracks.push_back(connected);
      connected.clear();
      if(Verbosity() > 3) std::cout << "           closing out last set " << std::endl;
    }
    
  if(Verbosity() > 3)   std::cout << "connected_tracks size " << connected_tracks.size() << std::endl;

  return connected_tracks;
}
*/

/*
void SecondaryVertexFinder::removeOutlierTrackPairs()
{
  //  Note: std::multimap<unsigned int, std::pair<unsigned int, std::pair<Eigen::Vector3d,  Eigen::Vector3d>>>  _track_pair_pca_map
 
  for(auto it : _vertex_set) 
    {
      unsigned int vtxid = it;
      if(Verbosity() > 1) std::cout << "calculate average position for vertex " << vtxid << std::endl; 

      // we need the median values of the x and y positions 
      std::vector<double> vx;
      std::vector<double> vy;
      std::vector<double> vz;

      double pca_median_x = 0.;
      double pca_median_y = 0.;
      double pca_median_z = 0.;

      Eigen::Vector3d new_pca_avge(0.,0.,0.);
      double new_wt = 0.0;
      
      auto ret = _vertex_track_map.equal_range(vtxid);
      
      // Start by getting the positions for this vertex into vectors for the median calculation
      for (auto cit=ret.first; cit!=ret.second; ++cit)
	{	  
	  unsigned int tr1id = cit->second;
	  if(Verbosity() > 2) std::cout << "   vectors: get entries for track " << tr1id << " for vertex " << vtxid << std::endl; 
	  
	  // find all pairs for this vertex with tr1id
	  auto pca_range = _track_pair_pca_map.equal_range(tr1id);
	  for (auto pit=pca_range.first; pit!=pca_range.second; ++pit)
	    {
	      unsigned int tr2id = pit->second.first;
	      
	      Eigen::Vector3d PCA1 = pit->second.second.first;
	      Eigen::Vector3d PCA2 = pit->second.second.second;
	      
	      if(Verbosity() > 2)
		std::cout << " vectors: tr1id " << tr1id << " tr2id " << tr2id
			  << " PCA1 " << PCA1.x() << "  " << PCA1.y() << "  " << PCA1.z()
			  << " PCA2 " << PCA2.x() << "  " << PCA2.y() << "  " << PCA2.z()
			  << std::endl;
	      
	      vx.push_back(PCA1.x());
	      vx.push_back(PCA2.x());
	      vy.push_back(PCA1.y());
	      vy.push_back(PCA2.y());
	      vz.push_back(PCA1.z());
	      vz.push_back(PCA2.z());
	    }
	}
      
      // Get the medians for this vertex
      // Using the median as a reference for rejecting outliers only makes sense for more than 2 tracks
      if(vx.size() < 3)
	{
	  new_pca_avge.x() = getAverage(vx);
	  new_pca_avge.y() = getAverage(vy);
	  new_pca_avge.z() = getAverage(vz);
	  _vertex_position_map.insert(std::make_pair(vtxid, new_pca_avge));
	  if(Verbosity() > 1) 
	    std::cout << " Vertex has only 2 tracks, use average for PCA: " << new_pca_avge.x() << "  " << new_pca_avge.y() << "  " << new_pca_avge.z() << std::endl; 

	  // done with this vertex
	  continue;
	}
      
      pca_median_x = getMedian(vx);
      pca_median_y = getMedian(vy);
      pca_median_z = getMedian(vz);
      if(Verbosity() > 1) std::cout << "Median values: x " << pca_median_x << " y " << pca_median_y << " z : " << pca_median_z << std::endl;
      
      // Make the average vertex position with outlier rejection wrt the median
      for (auto cit=ret.first; cit!=ret.second; ++cit)
	{	  
	  unsigned int tr1id = cit->second;
	  if(Verbosity() > 2) std::cout << "   average: get entries for track " << tr1id << " for vertex " << vtxid << std::endl; 
	  
	  // find all pairs for this vertex with tr1id
	  auto pca_range = _track_pair_pca_map.equal_range(tr1id);
	  for (auto pit=pca_range.first; pit!=pca_range.second; ++pit)
	    {
	      unsigned int tr2id = pit->second.first;
	      
	      Eigen::Vector3d PCA1 = pit->second.second.first;
	      Eigen::Vector3d PCA2 = pit->second.second.second;
	      
	      if(
		 fabs(PCA1.x() - pca_median_x) < _outlier_cut &&
		 fabs(PCA1.y() - pca_median_y) < _outlier_cut &&
		 fabs(PCA2.x() - pca_median_x) < _outlier_cut &&
		 fabs(PCA2.y() - pca_median_y) < _outlier_cut 
		 )
		{
		  // good track pair, add to new average
		  
		  new_pca_avge += PCA1;
		  new_wt++;
		  new_pca_avge += PCA2;
		  new_wt++;	  
		}
	      else
		{
		  if(Verbosity() > 1) std::cout << "Reject pair with tr1id " << tr1id << " tr2id " << tr2id << std::endl;
		}
	    }
	}
      if(new_wt > 0.0)      
	new_pca_avge = new_pca_avge / new_wt;
      else
	{
	  // There were no pairs that survived the track cuts, use the median values
	  new_pca_avge.x() = pca_median_x;
	  new_pca_avge.y() = pca_median_y;
	  new_pca_avge.z() = pca_median_z;
	}

      _vertex_position_map.insert(std::make_pair(vtxid, new_pca_avge));

    }
  
  return;
 }
*/

/*  
double SecondaryVertexFinder::getMedian(std::vector<double> &v)
{
  double median = 0.0;

  if( (v.size() % 2) == 0)
    {
      // even number of entries
      // we want the average of the middle two numbers, v.size()/2 and v.size()/2-1
      auto m1 = v.begin() + v.size()/2;
      std::nth_element(v.begin(), m1, v.end());
      double median1 =  v[v.size()/2]; 

      auto m2 = v.begin() + v.size()/2 - 1;
      std::nth_element(v.begin(), m2, v.end());
      double median2 =  v[v.size()/2 - 1]; 

      median = (median1 + median2) / 2.0; 
      if(Verbosity() > 2) std::cout << "The vector size is " << v.size() 
				    << " element m is " << v.size() / 2  << " = " << v[v.size()/2] 
				    << " element m-1 is " << v.size() / 2 -1 << " = " << v[v.size()/2-1] 
				    <<  std::endl;
    } 
  else
    {
      // odd number of entries
      auto m = v.begin() + v.size()/2;
      std::nth_element(v.begin(), m, v.end());
      median =  v[v.size()/2];
      if(Verbosity() > 2) std::cout << "The vector size is " << v.size() << " element m is " << v.size() / 2 << " = " << v[v.size()/2] <<  std::endl;
    }

    return median ;
}
double SecondaryVertexFinder::getAverage(std::vector<double> &v)
{
  double avge = 0.0;
  double wt = 0.0;
  for(auto it : v)
    {
      avge += it;
      wt++;
    }

  avge /= wt;
  if(Verbosity() > 2)
    {
      std::cout << " average = " << avge << std::endl;
    }
  
  return avge;
}
*/

