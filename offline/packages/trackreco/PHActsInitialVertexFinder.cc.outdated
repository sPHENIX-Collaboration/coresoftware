#include "PHActsInitialVertexFinder.h"

#include <trackbase_historic/ActsTransformations.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHRandomSeed.h>

#if __cplusplus < 201402L
#include <boost/make_unique.hpp>
#endif

#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertex_v1.h>
#include <trackbase_historic/SvtxVertexMap_v1.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <ActsExamples/Plugins/BField/BFieldOptions.hpp>
#include <ActsExamples/Plugins/BField/ScalableBField.hpp>

#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/InterpolatedBFieldMap.hpp>
#include <Acts/MagneticField/SharedBField.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Helpers.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <Acts/Utilities/Units.hpp>
#include <Acts/Vertexing/FullBilloirVertexFitter.hpp>
#include <Acts/Vertexing/HelicalTrackLinearizer.hpp>
#include <Acts/Vertexing/ImpactPointEstimator.hpp>
#include <Acts/Vertexing/IterativeVertexFinder.hpp>
#include <Acts/Vertexing/LinearizedTrack.hpp>
#include <Acts/Vertexing/Vertex.hpp>
#include <Acts/Vertexing/VertexFinderConcept.hpp>
#include <Acts/Vertexing/VertexingOptions.hpp>
#include <Acts/Vertexing/ZScanVertexFinder.hpp>

#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>

#include <memory>

PHActsInitialVertexFinder::PHActsInitialVertexFinder(const std::string& name)
  : PHInitVertexing(name)
{}

int PHActsInitialVertexFinder::Setup(PHCompositeNode *topNode)
{
  if(getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;
  
  if(createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;
  
  m_seed = PHRandomSeed();
  m_random_number_generator.seed(m_seed);

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsInitialVertexFinder::Process(PHCompositeNode */*topNode*/)
{
  if(Verbosity() > 0)
    std::cout << "PHActsInitialVertexFinder processing event " 
	      << m_event << std::endl;

  if(m_trackMap->size() == 0)
    {
      std::cout << PHWHERE 
		<< "No silicon track seeds found. Can't run initial vertexing, setting dummy vertex of (0,0,0)" 
		<< std::endl;
      createDummyVertex(0,0,0);
    }
  else if(m_trackMap->size() == 1)
    {
      auto track = m_trackMap->get(0);
      createDummyVertex(track->get_x(),
			track->get_y(),
			track->get_z());
    }
  else
    {
      InitKeyMap keyMap;
      auto trackPointers = getTrackPointers(keyMap);
      
      auto vertices = findVertices(trackPointers);
      
      fillVertexMap(vertices, keyMap);
      
      /// Need to check that silicon stubs which may have been
      /// skipped over still have a vertex association
      checkTrackVertexAssociation();
      
      for(auto track : trackPointers)
	{
	  delete track;
	}
    }

 if(Verbosity() > 0)
    std::cout << "PHActsInitialVertexFinder processed event "
	      << m_event << std::endl;

  m_event++;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsInitialVertexFinder::ResetEvent(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsInitialVertexFinder::End(PHCompositeNode */*topNode*/)
{

  std::cout << "Acts IVF succeeded " << m_successFits 
	    << " out of " << m_totVertexFits << " total fits"
	    << std::endl;
  
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsInitialVertexFinder::checkTrackVertexAssociation()
{
  
  for(auto& [trackKey, track] : *m_trackMap)
    {
      /// If the track wasn't already given a vertex ID, it wasn't 
      /// included in the initial vertex finding due to Acts not liking
      /// tracks with large transverse position. So find the closest
      /// z vertex to it and assign it
      if(track->get_vertex_id() != UINT_MAX)
	continue;

      const auto trackZ = track->get_z();
      
      double closestVertZ = 9999;
      int vertId = -1;
      for(auto& [vertexKey, vertex] : *m_vertexMap)
	{
	  double dz = fabs(trackZ - vertex->get_z());

	  if(dz < closestVertZ) 
	    {
	      vertId = vertexKey;
	      closestVertZ = dz;
	    }

	}
      
      auto vertex = m_vertexMap->get(vertId);
      vertex->insert_track(trackKey);
      track->set_vertex_id(vertId);	
    }

}
void PHActsInitialVertexFinder::fillVertexMap(VertexVector& vertices,
					      InitKeyMap& keyMap)
{
  unsigned int vertexId = 0;
  if(vertices.size() != 0)
    m_vertexMap->clear();

  /// Create a fail safe for (e.g.) single particle events which 
  /// don't return a vertex
  if(vertices.size() == 0)
    {
      createDummyVertex(0,0,0);
      if(Verbosity() > 1)
	std::cout << "No vertices found. Adding a dummy vertex"
		  << std::endl;
      return;
    }
    
  for(auto vertex : vertices)
    {
      const auto &[chi2, ndf] = vertex.fitQuality();
      const auto numTracks = vertex.tracks().size();
      
      if(Verbosity() > 1)
	{
	  std::cout << "Found vertex at (" << vertex.position().x()
		    << ", " << vertex.position().y() << ", " 
		    << vertex.position().z() << ")" << std::endl;
	  std::cout << "Vertex has ntracks = " << numTracks
		    << " with chi2/ndf " << chi2 / ndf << std::endl;
	}

      /// Fill SvtxVertexMap
      #if __cplusplus < 201402L
      auto svtxVertex = boost::make_unique<SvtxVertex_v1>();
      #else
      auto svtxVertex = std::make_unique<SvtxVertex_v1>();
      #endif

      svtxVertex->set_x(vertex.position().x() / Acts::UnitConstants::cm);  
      svtxVertex->set_y(vertex.position().y() / Acts::UnitConstants::cm);
      svtxVertex->set_z(vertex.position().z() / Acts::UnitConstants::cm);
      for(int i = 0; i < 3; ++i) 
	for(int j = 0; j < 3; ++j)
	  svtxVertex->set_error(i, j,
				vertex.covariance()(i,j) 
				/ Acts::UnitConstants::cm2); 
	        
      svtxVertex->set_chisq(chi2);
      svtxVertex->set_ndof(ndf);
      svtxVertex->set_t0(vertex.time());
      svtxVertex->set_id(vertexId);
          
      for(const auto& track : vertex.tracks())
	{
	  const auto originalParams = track.originalParams;

	  const auto trackKey = keyMap.find(originalParams)->second;
	  svtxVertex->insert_track(trackKey);

	  /// Give the track the appropriate vertex id
	  const auto svtxTrack = m_trackMap->find(trackKey)->second;
	  
	  if(Verbosity() > 3)
	    {   
	      svtxTrack->identify();
	      std::cout << "Updating track key " << trackKey << " with vertex "
			<< vertexId << std::endl;
	    }

	  svtxTrack->set_vertex_id(vertexId);
	}

      m_vertexMap->insert(svtxVertex.release());

      ++vertexId;
    }
      
  return;
}

void PHActsInitialVertexFinder::createDummyVertex(const float x,
						  const float y,
						  const float z)
{

  /// If the Acts IVF finds 0 vertices or there weren't enough tracks
  /// for it to properly identify a good vertex. So just create
  /// a dummy vertex with large covariance for rest of track
  /// reconstruction to avoid seg faults

  #if __cplusplus < 201402L
  auto svtxVertex = boost::make_unique<SvtxVertex_v1>();
  #else
  auto svtxVertex = std::make_unique<SvtxVertex_v1>();
  #endif

  svtxVertex->set_x(x);  
  svtxVertex->set_y(y);
  svtxVertex->set_z(z);
  
  for(int i = 0; i < 3; ++i) 
    for(int j = 0; j < 3; ++j)
      {
	if( i == j)
	  svtxVertex->set_error(i, j, 10.); 
	else 
	  svtxVertex->set_error(i,j, 0);
      }

  float nan = NAN;
  svtxVertex->set_chisq(nan);
  svtxVertex->set_ndof(nan);
  svtxVertex->set_t0(nan);
  svtxVertex->set_id(0);

  m_vertexMap->insert(svtxVertex.release());
  std::cout << "Created dummy vertex at " << x << ", " << y << ", " << z
	    << std::endl;
  for(auto& [key, track] : *m_trackMap)
    track->set_vertex_id(0);

}
 
VertexVector PHActsInitialVertexFinder::findVertices(TrackParamVec& tracks)
{

  m_totVertexFits++;

  auto field = m_tGeometry->magField;

  /// Determine the input mag field type from the initial geometry
  /// and run the vertex finding with the determined mag field
  return std::visit([tracks, this](auto &inputField) {
      /// Setup aliases
      using InputMagneticField = 
	typename std::decay_t<decltype(inputField)>::element_type;
      using MagneticField = Acts::SharedBField<InputMagneticField>;
      
      using Stepper = Acts::EigenStepper<MagneticField>;
      using Propagator = Acts::Propagator<Stepper>;
      using TrackParameters = Acts::BoundTrackParameters;
      using Linearizer = Acts::HelicalTrackLinearizer<Propagator>;
      using VertexFitter = 
	Acts::FullBilloirVertexFitter<TrackParameters,Linearizer>;
      using ImpactPointEstimator = 
	Acts::ImpactPointEstimator<TrackParameters, Propagator>;
      using VertexSeeder = Acts::ZScanVertexFinder<VertexFitter>;
      using VertexFinder = 
	Acts::IterativeVertexFinder<VertexFitter, VertexSeeder>;
      using VertexFinderOptions = Acts::VertexingOptions<TrackParameters>;

      static_assert(Acts::VertexFinderConcept<VertexSeeder>,
		    "VertexSeeder does not fulfill vertex finder concept.");
      static_assert(Acts::VertexFinderConcept<VertexFinder>,
		    "VertexFinder does not fulfill vertex finder concept.");

      auto logLevel = Acts::Logging::FATAL;
      if(Verbosity() > 4)
	logLevel = Acts::Logging::VERBOSE;
      auto logger = Acts::getDefaultLogger("PHActsInitialVertexFinder", 
					   logLevel);
    
      MagneticField bField(inputField);
      auto propagator = std::make_shared<Propagator>(Stepper(bField));
      
      /// Setup vertex finder now
      typename VertexFitter::Config vertexFitterConfig;

      /// Vertex fitter seems to have no performance difference when
      /// iterating once vs. default of 5 times. Additionally, iterating
      /// more than once causes vertices with low numbers of tracks to
      /// fail fitting, causing an error to be thrown and 0 vertices 
      /// returned
      vertexFitterConfig.maxIterations = 1;

      VertexFitter vertexFitter(std::move(vertexFitterConfig));
      
      typename Linearizer::Config linearizerConfig(bField, propagator);
      Linearizer linearizer(std::move(linearizerConfig));
      
      typename ImpactPointEstimator::Config ipEstConfig(bField, propagator);
      ImpactPointEstimator ipEst(std::move(ipEstConfig));
      
      typename VertexSeeder::Config seederConfig(ipEst);

      /// Don't weight track contribution by pT, since the momentum
      /// resolution of the silicon seeds is poor
      if(m_disableWeights)
	seederConfig.disableAllWeights = true;
      VertexSeeder seeder(std::move(seederConfig));
      
      typename VertexFinder::Config finderConfig(std::move(vertexFitter), 
						 std::move(linearizer),
						 std::move(seeder), ipEst);
      finderConfig.maxVertices = m_maxVertices;
      finderConfig.reassignTracksAfterFirstFit = true;
      VertexFinder finder(finderConfig, std::move(logger));

      typename VertexFinder::State state(m_tGeometry->magFieldContext);
      VertexFinderOptions finderOptions(m_tGeometry->getGeoContext(),
					m_tGeometry->magFieldContext);
  
      auto result = finder.find(tracks, finderOptions, state);
    
      VertexVector vertexVector;

      if(result.ok())
	{
	  m_successFits++;

	  auto vertexCollection = *result;
	  
	  if(Verbosity() > 1)
	    {
	      std::cout << "Acts IVF found " << vertexCollection.size()
			<< " vertices in event" << std::endl;
	    }
	  
	  for(const auto& vertex : vertexCollection) 
	    {
	      vertexVector.push_back(vertex);
	    }
	}
      else
	{
	  if(Verbosity() > 0)
	    {
	      std::cout << "Acts initial vertex finder returned error: " 
			<< result.error().message() << std::endl;
	      std::cout << "Track positions IVF used are : " << std::endl;
	      for(const auto track : tracks)
		{
		  const auto position = track->position(m_tGeometry->getGeoContext());
		  std::cout << "(" << position(0) << ", " << position(1)
			    << ", " << position(2) << std::endl;
		}
	    }
	}
    
      return vertexVector;
      
    } /// end lambda
    , field
    ); /// end std::visit call

}

std::vector<SvtxTrack*> PHActsInitialVertexFinder::sortTracks()
{

  /// Implement a simple k-means clustering algorithm. Randomly select
  /// m_nCentroid track z PCAs (centroids), assign all tracks to 
  /// a centroid based on which they are closest to, and then iterate
  /// to update clusters and centroids

  std::vector<Acts::Vector3D> centroids(m_nCentroids);
  std::uniform_int_distribution<int> indices(0,m_trackMap->size() - 1);
  std::vector<int> usedIndices;

  /// Get the original centroids
  for(auto& centroid : centroids) 
    {
      auto index = indices(m_random_number_generator);
      for(const auto used : usedIndices)
	if(index == used)
	  index = indices(m_random_number_generator);

      usedIndices.push_back(index);
      
      centroid = Acts::Vector3D(m_trackMap->get(index)->get_x(),
				m_trackMap->get(index)->get_y(),
				m_trackMap->get(index)->get_z());
      
      if(Verbosity() > 3)
	{
	  std::cout << "Centroid is (" << centroid(0) 
		    << ", " << centroid(1) << ", " 
		    << centroid(2) << ")" << std::endl;
	}
    }
  
  /// This map holds the centroid index as the key and a
  /// vector of SvtxTracks that correspond to that centroid
  auto clusters = createCentroidMap(centroids);

  /// Take the map and identified centroids and remove tracks
  /// that aren't compatible
  auto sortedTracks = getIVFTracks(clusters, centroids);

  return sortedTracks;
}

std::vector<SvtxTrack*> PHActsInitialVertexFinder::getIVFTracks(
		        CentroidMap& clusters, 
			std::vector<Acts::Vector3D>& centroids)
{
  
  std::vector<SvtxTrack*> sortedTracks;

  if(Verbosity() > 2)
    {
      std::cout << "Final centroids are : " << std::endl;
      for(const auto& [centroidIndex, trackVec] : clusters)
	{
	  std::cout << "Centroid: " << centroids.at(centroidIndex).transpose()
		    << " has tracks " << std::endl;
	  for(const auto& track : trackVec)
	    {
	      std::cout << "(" << track->get_x() << ", "
			<< track->get_y() << ", " << track->get_z()
			<< ")" << std::endl;
	    }
	}
    }

  /// Note the centroid that has the most tracks
  int maxTrackCentroid = 0;
  std::vector<Acts::Vector3D> stddev(m_nCentroids);

  for(const auto& [centroidIndex, trackVec] : clusters)
    {
      Acts::Vector3D sum = Acts::Vector3D::Zero();
      if(trackVec.size() > maxTrackCentroid)
	{
	  maxTrackCentroid = trackVec.size();
	}

      for(const auto& track : trackVec)
	{
	  for(int i = 0; i < sum.rows(); i++)
	    {
	      sum(i) += pow(track->get_pos(i) - centroids.at(centroidIndex)(i), 2);
	    }
	}
      for(int i = 0; i < 3; i++)
	{ 
	  stddev.at(centroidIndex)(i) = sqrt(sum(i) / trackVec.size()); 
	}
    }
  
  for(const auto& [centroidIndex, trackVec] : clusters)
    {
      /// skip centroids that have a very small number of tracks
      /// compared to the largest centroid, as these are most likely
      /// composed of only a few (bad) stubs
      if(trackVec.size() < 0.2 * maxTrackCentroid)
	continue;

      /// Skip large transverse PCA centroids
      float centroidR = sqrt(pow(centroids.at(centroidIndex)(0), 2) +
			     pow(centroids.at(centroidIndex)(1), 2));
    
      if(Verbosity() > 2)
	{
	  std::cout << "Checking to add tracks from centroid " 
		    << centroids.at(centroidIndex).transpose() << std::endl;
	}

      if(centroidR > m_pcaCut)
	{
	  continue;
	}
    
      for(const auto& track : trackVec)
	{
	  Acts::Vector3D pulls = Acts::Vector3D::Zero();
	  for(int i = 0; i < 3; i++)
	    {
	      pulls(i) = fabs(track->get_pos(i) - centroids.at(centroidIndex)(i)) / stddev.at(centroidIndex)(i);
	    }
	  
	  if(Verbosity() > 3)
	    {
	      std::cout << "Track pos is (" << track->get_x() << ", " 
			<< track->get_y() << ", " << track->get_z() 
			<< ") and pull is " << pulls.transpose() << std::endl;
	    }
	 
	  if ((pulls(0) < 2 and pulls(1) < 2 and pulls(2) < 2))
	    {
	      sortedTracks.push_back(track);
	    }
	  else
	    {
	      if(m_removeSeeds)
		{
		  m_trackMap->erase(track->get_id());
		}

	      if(Verbosity() > 3)
		{
		  std::cout << "Not adding track with pos (" << track->get_x()
			    << ", " << track->get_y() << ", " << track->get_z() 
			    << ") as it is incompatible with centroid " 
			    << centroids.at(centroidIndex).transpose() 
			    << " with std dev " 
			    << stddev.at(centroidIndex).transpose() << std::endl;
		}
	    }
	}
    }

  return sortedTracks;

}

CentroidMap PHActsInitialVertexFinder::createCentroidMap(std::vector<Acts::Vector3D>& centroids)
{
  CentroidMap clusters;
  
  for(int niter = 0; niter < m_nIterations; niter++)
    {
      /// reset the centroid-track map
      clusters.clear();
      for(unsigned int i = 0; i<m_nCentroids; i++)
	{
	  std::vector<SvtxTrack*> vec;
	  clusters.insert(std::make_pair(i, vec));
	}  
      
      if(Verbosity() > 3)
	{
	  for(int i =0; i< m_nCentroids; i++)
	    std::cout << "Starting centroid is : " 
		      << centroids.at(i).transpose() << std::endl;
	}
      for(const auto& [key, track] : *m_trackMap)
	{
	  Acts::Vector3D trackPos(track->get_x(),
				  track->get_y(),
				  track->get_z());
	      
	  double minDist = std::numeric_limits<double>::max();
	  unsigned int centKey = std::numeric_limits<unsigned int>::max();
	  for(int i = 0; i < centroids.size(); i++)
	    {
	      double dist = sqrt(pow(trackPos(0) - centroids.at(i)(0), 2) +
				 pow(trackPos(1) - centroids.at(i)(1), 2) +
				 pow(trackPos(2) - centroids.at(i)(2), 2));
	      if(dist < minDist)
		{
		  minDist = dist;
		  centKey = i;
		  if(Verbosity() > 3)
		    {
		      std::cout << "mindist and centkey are " 
				<< minDist << ", " << centKey 
				<< std::endl;
		    }
		}
	    }
	  
	  /// Add this track to the map that associates centroids with tracks
	  if(Verbosity() > 3)
	    {
	      std::cout << "adding track with " << trackPos.transpose() 
			<< " to centroid " 
			<< centroids.at(centKey).transpose() << std::endl;
	    }

	  clusters.find(centKey)->second.push_back(track);	  
	}
      
      /// Update pos centroids
      std::vector<Acts::Vector3D> newCentroids(m_nCentroids);
      for(auto& centroid : newCentroids)
	centroid = Acts::Vector3D(0,0,0);

      for(const auto& [centroidVal, trackVec] : clusters)
	{
	  for(const auto& track : trackVec)
	    {
	      for(int i = 0; i < 3; i++)
		newCentroids.at(centroidVal)(i) += track->get_pos(i);
	    }

	  /// Sets the centroid as the average value
	  centroids.at(centroidVal) = 
	    newCentroids.at(centroidVal) / trackVec.size();
	}
      
      if(Verbosity() > 3)
	{
	  for(int i = 0; i < m_nCentroids; i++)
	    std::cout << "new centroids " << centroids.at(i).transpose() 
		      << std::endl;
   
	  for(const auto& [centKey, trackVec] : clusters)
	    {
	      std::cout << "cent key : " << centKey << " has " << trackVec.size() 
			<< "tracks" << std::endl;
	      for(const auto track : trackVec) 
		{
		  std::cout << "track id : " << track->get_id() 
			    << " with pos (" << track->get_x() << ", "
			    << track->get_y() << ", " <<  track->get_z()
			    << ")" << std::endl;
		  
		}
	    }
	}
    }

  return clusters;

}

TrackParamVec PHActsInitialVertexFinder::getTrackPointers(InitKeyMap& keyMap)
{
  TrackParamVec tracks;

  /// If there are fewer tracks than centroids, just one with 1 centroid
  /// Otherwise algorithm does not converge. nCentroids should only be
  /// a handful, so this only affects small nTrack events.
  if(m_trackMap->size() < m_nCentroids) 
    {
      m_nCentroids = 1;
    }

  std::vector<SvtxTrack*> sortedTracks;
  if(m_svtxTrackMapName.find("Silicon") != std::string::npos)
    {
      sortedTracks = sortTracks();
    }
  else
    {
      for(const auto& [key, track] : *m_trackMap)
	sortedTracks.push_back(track);
    }

  for(const auto& track : sortedTracks)
    {
      if(Verbosity() > 3)
	{
	  std::cout << "Adding track seed to vertex finder " 
		    << std::endl;
	  track->identify();
	}
      
      /// Only vertex with stubs that have five clusters
      if(m_svtxTrackMapName.find("Silicon") != std::string::npos)
	{
	  if(track->size_cluster_keys() < 5)
	    {
	      continue;
	    }
	}
    
      const Acts::Vector4D stubVec(
                  track->get_x() * Acts::UnitConstants::cm,
		  track->get_y() * Acts::UnitConstants::cm,
		  track->get_z() * Acts::UnitConstants::cm,
		  10 * Acts::UnitConstants::ns);
     
      const Acts::Vector3D stubMom(track->get_px(),
				   track->get_py(),
				   track->get_pz());
      int trackQ = track->get_charge() * Acts::UnitConstants::e;
      
      /// The vertexing has an expectation about the field direction
      /// and associated charge. For the 3D map, to work with the 
      /// required swap of the field for the silicon seeds, we need
      /// to flip the charge
      if(m_magField.find("3d") != std::string::npos)
	trackQ *= -1;

      const double p = track->get_p();
      
      /// Make a dummy covariance matrix for Acts that corresponds
      /// to the resolutions of the silicon seeds
      Acts::BoundSymMatrix cov;
      if(m_resetTrackCovariance)
	cov << 50 * Acts::UnitConstants::um, 0., 0., 0., 0., 0.,
	       0., 30 * Acts::UnitConstants::um, 0., 0., 0., 0.,
	       0., 0., 0.005, 0., 0., 0.,
	       0., 0., 0., 0.001, 0., 0.,
	       0., 0., 0., 0., 0.3 , 0.,
	       0., 0., 0., 0., 0., 1.;
      
      else 
	{
	  ActsTransformations transform;
	  transform.setVerbosity(Verbosity());
	  cov = transform.rotateSvtxTrackCovToActs(track,
						   m_tGeometry->getGeoContext());
	}

      /// Make a dummy perigee surface to bound the track to
      auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(
				    Acts::Vector3D(stubVec(0),
						   stubVec(1),
						   stubVec(2)));
								     

      const auto param = new Acts::BoundTrackParameters(
			           perigee,
				   m_tGeometry->getGeoContext(),
				   stubVec, stubMom,
				   p, trackQ, cov);

      tracks.push_back(param);
      keyMap.insert(std::make_pair(param, track->get_id()));
    }

  return tracks;
}

int PHActsInitialVertexFinder::getNodes(PHCompositeNode *topNode)
{

  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_svtxTrackMapName.c_str());
  if(!m_trackMap)
    {
      std::cout << PHWHERE << "No " << m_svtxTrackMapName.c_str() 
		<< " on node tree, bailing."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

   m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode, 
							 "ActsTrackingGeometry");
  if(!m_tGeometry)
    {
      std::cout << PHWHERE << "ActsTrackingGeometry not on node tree. Exiting"
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsInitialVertexFinder::createNodes(PHCompositeNode *topNode)
{
  
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  /// Check that it is there
  if (!dstNode)
  {
    std::cerr << "DST Node missing, quitting" << std::endl;
    throw std::runtime_error("failed to find DST node in PHActsInitialVertexFinder::createNodes");
  }

  /// Get the tracking subnode
  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode,
						  m_svtxVertexMapName.c_str());
  
  if(!m_vertexMap)
    {
      m_vertexMap = new SvtxVertexMap_v1;
      PHIODataNode<PHObject>* vertexNode = new PHIODataNode<PHObject>( 
		   m_vertexMap, m_svtxVertexMapName.c_str(),"PHObject");

      svtxNode->addNode(vertexNode);

    }


  return Fun4AllReturnCodes::EVENT_OK;
}
