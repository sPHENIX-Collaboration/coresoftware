#include "PHActsVertexFinder.h"
#include "ActsTransformations.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#if __cplusplus < 201402L
#include <boost/make_unique.hpp>
#endif

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertex_v1.h>
#include <trackbase_historic/SvtxVertexMap_v1.h>

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
#include <iostream>

PHActsVertexFinder::PHActsVertexFinder(const std::string &name)
  : PHInitVertexing(name)
  , m_actsFitResults(nullptr)
  , m_actsVertexMap(nullptr)
{
}

int PHActsVertexFinder::Setup(PHCompositeNode *topNode)
{
  int ret = createNodes(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK)
    return ret;
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsVertexFinder::Process(PHCompositeNode *topNode)
{
  if(Verbosity() > 0)
    {
      std::cout << "Starting event " << m_event << " in PHActsVertexFinder"
		<< std::endl;
    }
  
  int ret = getNodes(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK)
    return ret;

  /// Create a map that correlates the track momentum to the track key
  KeyMap keyMap;

  /// Get the list of tracks in Acts form
  auto trackPointers = getTracks(keyMap);

  auto vertices = findVertices(trackPointers);

  fillVertexMap(vertices, keyMap);
  
  /// Clean up the track pointer vector memory
  for(auto track : trackPointers)
    {
      delete track;
    }

  if(Verbosity() > 0)
    std::cout << "Finished PHActsVertexFinder::process_event" << std::endl;

  m_event++;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsVertexFinder::ResetEvent(PHCompositeNode *topNode)
{
  m_actsVertexMap->clear();
  m_svtxVertexMap->clear();

  return Fun4AllReturnCodes::EVENT_OK;

}

int PHActsVertexFinder::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

TrackPtrVector PHActsVertexFinder::getTracks(KeyMap& keyMap)
{
  std::vector<const Acts::BoundTrackParameters*> trackPtrs;

  for(const auto &[key, traj] : *m_actsFitResults)
  {
    const auto &[trackTips, mj] = traj.trajectory();
    
    for(const size_t &trackTip : trackTips)
      {
	if(traj.hasTrackParameters(trackTip))
	  {
	    const auto param = new Acts::BoundTrackParameters(traj.trackParameters(trackTip));
	    keyMap.insert(std::make_pair(param, key));
	    trackPtrs.push_back(param);
	  }
      }
  }
  
  if(Verbosity() > 3)
    {
      std::cout << "Finding vertices for the following number of tracks "
		<< trackPtrs.size()
		<< std::endl;
     
      for(const auto param : trackPtrs)
	{
	  std::cout << "Track position: (" 
		    << param->position(m_tGeometry->geoContext)(0)
		    <<", " << param->position(m_tGeometry->geoContext)(1) << ", "
		    << param->position(m_tGeometry->geoContext)(2) << ")" 
		    << std::endl;
	}
    }

  return trackPtrs;

}

VertexVector PHActsVertexFinder::findVertices(TrackPtrVector& tracks)
{
  /// Determine the input mag field type from the initial geometry
  /// and run the vertex finding with the determined mag field

  return std::visit([tracks, this](auto &inputField) {
      /// Setup aliases
      using InputMagneticField = 
	typename std::decay_t<decltype(inputField)>::element_type;
      using MagneticField = Acts::SharedBField<InputMagneticField>;
      
      using Stepper = Acts::EigenStepper<MagneticField>;
      using Propagator = Acts::Propagator<Stepper>;
      using PropagatorOptions = Acts::PropagatorOptions<>;
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
      auto logger = Acts::getDefaultLogger("PHActsVertexFinder", logLevel);

      MagneticField bField(inputField);
      auto propagator = std::make_shared<Propagator>(Stepper(bField));
      
      /// Setup vertex finder now
      typename VertexFitter::Config vertexFitterConfig;
      VertexFitter vertexFitter(std::move(vertexFitterConfig));
      
      typename Linearizer::Config linearizerConfig(bField, propagator);
      Linearizer linearizer(std::move(linearizerConfig));
      
      typename ImpactPointEstimator::Config ipEstConfig(bField, propagator);
      ImpactPointEstimator ipEst(std::move(ipEstConfig));
      
      typename VertexSeeder::Config seederConfig(ipEst);
      VertexSeeder seeder(std::move(seederConfig));
      
      typename VertexFinder::Config finderConfig(std::move(vertexFitter), 
						 std::move(linearizer),
						 std::move(seeder), ipEst);
      finderConfig.maxVertices = m_maxVertices;
      finderConfig.reassignTracksAfterFirstFit = true;
      VertexFinder finder(finderConfig, std::move(logger));
      
      typename VertexFinder::State state(m_tGeometry->magFieldContext);
      VertexFinderOptions finderOptions(m_tGeometry->geoContext,
					m_tGeometry->magFieldContext);
      
      auto result = finder.find(tracks, finderOptions, state);
    
      VertexVector vertexVector;

      if(result.ok())
	{
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
	  if(Verbosity() > 1)
	    {
	      std::cout << "Acts vertex finder returned error: " 
			<< result.error().message() << std::endl;
	    }	  
	}

      return vertexVector;
      
    } /// end lambda
    , m_tGeometry->magField
    ); /// end std::visit call
}




void PHActsVertexFinder::fillVertexMap(VertexVector& vertices,
				       KeyMap& keyMap)
{
  unsigned int key = 0;
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

      /// Make some basic QA cuts on the vertices 
      if(numTracks < 3)
	continue;

      /// Fill Acts vertex map
      auto pair = std::make_pair(key, vertex);
      m_actsVertexMap->insert(pair);

      /// Fill SvtxVertexMap
      #if __cplusplus < 201402L
      auto svtxVertex = boost::make_unique<SvtxVertex_v1>();
      #else
      auto svtxVertex = std::make_unique<SvtxVertex_v1>();
      #endif

      const auto vertexX = vertex.position().x() / Acts::UnitConstants::cm;
      const auto vertexY = vertex.position().y() / Acts::UnitConstants::cm;
      const auto vertexZ = vertex.position().z() / Acts::UnitConstants::cm;
      
      svtxVertex->set_x(vertexX);  
      svtxVertex->set_y(vertexY);
      svtxVertex->set_z(vertexZ);
      for(int i = 0; i < 3; ++i) 
	{
	  for(int j = 0; j < 3; ++j)
	    {
	      svtxVertex->set_error(i, j,
	       vertex.covariance()(i,j) / Acts::UnitConstants::cm2); 
	    }
	}

      for(const auto track : vertex.tracks())
	{
	  const auto originalParams = track.originalParams;
	  const auto trackKey = keyMap.find(originalParams)->second;
	  svtxVertex->insert_track(trackKey);

	  updateTrackDCA(trackKey, Acts::Vector3D(vertexX,
						  vertexY,
						  vertexZ));
	}

      svtxVertex->set_chisq(chi2);
      svtxVertex->set_ndof(ndf);
      svtxVertex->set_t0(vertex.time());
      svtxVertex->set_id(key);

      m_svtxVertexMap->insert(svtxVertex.release());

      ++key;
    }
      
  return;
}

void PHActsVertexFinder::updateTrackDCA(const unsigned int trackKey,
					const Acts::Vector3D vertex)
{
  
  auto svtxTrack = m_svtxTrackMap->find(trackKey)->second;
  
  Acts::Vector3D pos(svtxTrack->get_x(),
		     svtxTrack->get_y(),
		     svtxTrack->get_z());
  Acts::Vector3D mom(svtxTrack->get_px(),
		     svtxTrack->get_py(),
		     svtxTrack->get_pz());
  
  pos -= vertex;

  Acts::ActsSymMatrixD<3> posCov;
  for(int i = 0; i < 3; ++i)
    {
      for(int j = 0; j < 3; ++j)
	{
	  posCov(i, j) = svtxTrack->get_error(i, j);
	} 
    }
  
  Acts::Vector3D r = mom.cross(Acts::Vector3D(0.,0.,1.));
  float phi = atan2(r(1), r(0));

  Acts::RotationMatrix3D rot;
  Acts::RotationMatrix3D rot_T;
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

  Acts::Vector3D pos_R = rot * pos;
  Acts::ActsSymMatrixD<3> rotCov = rot * posCov * rot_T;

  const auto dca3Dxy = pos_R(0);
  const auto dca3Dz = pos_R(2);
  const auto dca3DxyCov = rotCov(0,0);
  const auto dca3DzCov = rotCov(2,2);

  svtxTrack->set_dca3d_xy(dca3Dxy);
  svtxTrack->set_dca3d_z(dca3Dz);
  svtxTrack->set_dca3d_xy_error(sqrt(dca3DxyCov));
  svtxTrack->set_dca3d_z_error(sqrt(dca3DzCov));

}

int PHActsVertexFinder::createNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
    {
      std::cerr << "DST node is missing, quitting" << std::endl;
      throw std::runtime_error("Failed to find DST node in PHActsTracks::createNodes");
    }
  
  PHCompositeNode *svtxNode = 
    dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "SVTX"));
  
  if (!svtxNode)
    {
      svtxNode = new PHCompositeNode("SVTX");
      dstNode->addNode(svtxNode);
    }
 
  m_actsVertexMap = 
    findNode::getClass<VertexMap>(topNode, "ActsVertexMap");
  if(!m_actsVertexMap)
    {
      m_actsVertexMap = new VertexMap;
      
      PHDataNode<VertexMap> *node = 
	new PHDataNode<VertexMap>(m_actsVertexMap,
				  "ActsVertexMap");
   
      svtxNode->addNode(node);
    }

  m_svtxVertexMap = 
    findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMapActs");
  if(!m_svtxVertexMap)
    {
      m_svtxVertexMap = new SvtxVertexMap_v1;
      PHIODataNode<PHObject> *node = 
	new PHIODataNode<PHObject>(m_svtxVertexMap,
				   "SvtxVertexMapActs", "PHObject");
      svtxNode->addNode(node);
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsVertexFinder::getNodes(PHCompositeNode *topNode)
{
  
  m_actsFitResults = findNode::getClass<std::map<const unsigned int, Trajectory>>
    (topNode, "ActsFitResults");
  if(!m_actsFitResults)
    {
      std::cout << PHWHERE << "Acts Trajectories not found on node tree, exiting."
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

  m_svtxTrackMap = findNode::getClass<SvtxTrackMap>(topNode,
						    "SvtxTrackMap");
  if(!m_svtxTrackMap)
    {
      std::cout << PHWHERE << "No SvtxTrackMap on node tree, exiting."
		<< std::endl;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}
