#include "PHActsVertexFitter.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertexMap_v1.h>
#include <trackbase_historic/SvtxVertex_v1.h>

/// Tracking includes
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

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
#include <Acts/Vertexing/FullBilloirVertexFitter.hpp>
#include <Acts/Vertexing/HelicalTrackLinearizer.hpp>
#include <Acts/Vertexing/LinearizedTrack.hpp>
#include <Acts/Vertexing/VertexingOptions.hpp>

#include <iostream>

PHActsVertexFitter::PHActsVertexFitter(const std::string &name)
  : SubsysReco(name)
  , m_actsFitResults(nullptr)
  , m_tGeometry(nullptr)
{
}

int PHActsVertexFitter::Init(PHCompositeNode */*topNode*/)
{
  if (Verbosity() > 1)
    std::cout << "PHActsVertexFitter::Init" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsVertexFitter::End(PHCompositeNode */*topNode*/)
{
  if (Verbosity() > 1)
    std::cout << "PHActsVertexFitter::End " << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsVertexFitter::ResetEvent(PHCompositeNode */*topNode*/)
{
  m_actsVertexMap->clear();

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsVertexFitter::InitRun(PHCompositeNode *topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTRUN;

  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTRUN;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsVertexFitter::process_event(PHCompositeNode */*topNode*/)
{
  auto logLevel = Acts::Logging::FATAL;
  if (Verbosity() > 0)
  {
    std::cout << "Beginning PHActsVertexFitter::process_event number "
              << m_event << std::endl;
    if (Verbosity() > 5)
      logLevel = Acts::Logging::VERBOSE;
  }

  const auto vertexTrackMap = getTracks();

  for (const auto &[vertexId, trackVec] : vertexTrackMap)
  {
    const auto vertex = fitVertex(trackVec, logLevel);

    createActsSvtxVertex(vertexId, vertex);

    if (m_updateSvtxVertexMap)
      updateSvtxVertex(vertexId, vertex);
  }

  if (Verbosity() > 1)
    std::cout << "Finished PHActsVertexFitter::process_event"
              << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsVertexFitter::updateSvtxVertex(const unsigned int vertexId,
                                          ActsVertex vertex)
{
  auto svtxVertex = m_vertexMap->get(vertexId);

  if (Verbosity() > 1)
    std::cout << "Updating SvtxVertex id " << vertexId << std::endl;

  svtxVertex->set_x(vertex.position().x() / Acts::UnitConstants::cm);
  svtxVertex->set_y(vertex.position().y() / Acts::UnitConstants::cm);
  svtxVertex->set_z(vertex.position().z() / Acts::UnitConstants::cm);

  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      svtxVertex->set_error(i, j,
                            vertex.covariance()(i, j) / Acts::UnitConstants::cm2);
    }
  }

  const auto &[chi2, ndf] = vertex.fitQuality();
  svtxVertex->set_ndof(ndf);
  svtxVertex->set_chisq(chi2);
  svtxVertex->set_t0(vertex.time());
}

void PHActsVertexFitter::createActsSvtxVertex(const unsigned int vertexId,
                                              ActsVertex vertex)
{
#if __cplusplus < 201402L
  auto svtxVertex = boost::make_unique<SvtxVertex_v1>();
#else
  auto svtxVertex = std::make_unique<SvtxVertex_v1>();
#endif

  if (Verbosity() > 2)
  {
    std::cout << "Creating vertex for id " << vertexId
              << std::endl;
  }

  svtxVertex->set_x(vertex.position().x() / Acts::UnitConstants::cm);
  svtxVertex->set_y(vertex.position().y() / Acts::UnitConstants::cm);
  svtxVertex->set_z(vertex.position().z() / Acts::UnitConstants::cm);

  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      svtxVertex->set_error(i, j,
                            vertex.covariance()(i, j) / Acts::UnitConstants::cm2);
    }
  }

  const auto &[chi2, ndf] = vertex.fitQuality();
  svtxVertex->set_ndof(ndf);
  svtxVertex->set_chisq(chi2);
  svtxVertex->set_t0(vertex.time());
  svtxVertex->set_id(vertexId);

  m_actsVertexMap->insert(svtxVertex.release());
}

ActsVertex PHActsVertexFitter::fitVertex(BoundTrackParamVec tracks, Acts::Logging::Level logLevel) const
{
  /// Determine the input mag field type from the initial
  /// geometry created in MakeActsGeometry
  return std::visit([tracks, logLevel, this](auto &inputField) {
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
        Acts::FullBilloirVertexFitter<TrackParameters, Linearizer>;
    using VertexFitterOptions = Acts::VertexingOptions<TrackParameters>;

    auto logger = Acts::getDefaultLogger("PHActsVertexFitter",
                                         logLevel);

    /// Create necessary templated inputs for Acts vertex fitter
    MagneticField bField(inputField);
    auto propagator = std::make_shared<Propagator>(Stepper(bField));
    PropagatorOptions propagatorOpts(m_tGeometry->getGeoContext(),
                                     m_tGeometry->magFieldContext,
                                     Acts::LoggerWrapper(*logger));

    typename VertexFitter::Config vertexFitterCfg;
    VertexFitter fitter(vertexFitterCfg);
    typename VertexFitter::State state(m_tGeometry->magFieldContext);

    typename Linearizer::Config linConfig(bField, propagator);
    Linearizer linearizer(linConfig);

    /// Can add a vertex fitting constraint as an option, if desired
    VertexFitterOptions vfOptions(m_tGeometry->getGeoContext(),
                                  m_tGeometry->magFieldContext);

    /// Call the fitter and get the result
    auto fitRes = fitter.fit(tracks, linearizer,
                             vfOptions, state);

    Acts::Vertex<TrackParameters> fittedVertex;

    if (fitRes.ok())
    {
      fittedVertex = *fitRes;
      if (Verbosity() > 3)
      {
        std::cout << "Fitted vertex position "
                  << fittedVertex.position().x()
                  << ", "
                  << fittedVertex.position().y()
                  << ", "
                  << fittedVertex.position().z()
                  << std::endl;
      }
    }
    else
    {
      if (Verbosity() > 3)
      {
        std::cout << "Acts vertex fit error: "
                  << fitRes.error().message()
                  << std::endl;
      }
    }

    return fittedVertex;
  },
                    m_tGeometry->magField);  /// end std::visit call
}

VertexTrackMap PHActsVertexFitter::getTracks()
{
  VertexTrackMap trackPtrs;

  for (const auto &[key, track] : *m_trackMap)
  {
    const unsigned int vertexId = track->get_vertex_id();
    const auto trackParam = makeTrackParam(track);
    auto trackVecPos = trackPtrs.find(vertexId);

    if (trackVecPos == trackPtrs.end())
    {
      BoundTrackParamVec trackVec;

      trackVec.push_back(trackParam);
      auto pair = std::make_pair(vertexId, trackVec);
      trackPtrs.insert(pair);
    }
    else
    {
      trackVecPos->second.push_back(trackParam);
    }
  }

  if (Verbosity() > 3)
  {
    for (const auto &[vertexId, trackVec] : trackPtrs)
    {
      std::cout << "Fitting vertexId : " << vertexId
                << " with the following number of tracks "
                << trackVec.size()
                << std::endl;

      for (const auto param : trackVec)
      {
        std::cout << "Track position: ("
                  << param->position(m_tGeometry->getGeoContext())(0)
                  << ", " << param->position(m_tGeometry->getGeoContext())(1) << ", "
                  << param->position(m_tGeometry->getGeoContext())(2) << ")"
                  << std::endl;
      }
    }
  }

  return trackPtrs;
}

const Acts::BoundTrackParameters *PHActsVertexFitter::makeTrackParam(const SvtxTrack *track) const
{
  const Acts::Vector4D trackPos(
      track->get_x() * Acts::UnitConstants::cm,
      track->get_y() * Acts::UnitConstants::cm,
      track->get_z() * Acts::UnitConstants::cm,
      10 * Acts::UnitConstants::ns);

  const Acts::Vector3D trackMom(track->get_px(),
                                track->get_py(),
                                track->get_pz());

  const int trackQ = track->get_charge() * Acts::UnitConstants::e;
  const double p = track->get_p();
  Acts::BoundSymMatrix cov;

  for (int i = 0; i < 6; ++i)
  {
    for (int j = 0; j < 6; ++j)
    {
      cov(i, j) = track->get_error(i, j);
    }
  }

  auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(
      Acts::Vector3D(track->get_x() * Acts::UnitConstants::cm,
                     track->get_y() * Acts::UnitConstants::cm,
                     track->get_z() * Acts::UnitConstants::cm));

  const auto param = new Acts::BoundTrackParameters(
      perigee, m_tGeometry->getGeoContext(),
      trackPos, trackMom, p, trackQ, cov);

  return param;
}

int PHActsVertexFitter::getNodes(PHCompositeNode *topNode)
{
  m_actsFitResults = findNode::getClass<std::map<const unsigned int, Trajectory>>(topNode, "ActsFitResults");
  if (!m_actsFitResults)
  {
    std::cout << PHWHERE << "Acts Trajectories not found on node tree, exiting."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode,
                                                "SvtxTrackMap");
  if (!m_trackMap)
  {
    std::cout << PHWHERE << "No SvtxTrackMap on node tree. Bailing."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode,
                                                  "SvtxVertexMap");
  if (!m_vertexMap)
  {
    std::cout << PHWHERE << "No SvtxVertexMap on node tree, bailing"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode, "ActsTrackingGeometry");
  if (!m_tGeometry)
  {
    std::cout << PHWHERE << "ActsTrackingGeometry not on node tree. Exiting"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsVertexFitter::createNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

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
      findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMapActs");
  if (!m_actsVertexMap)
  {
    m_actsVertexMap = new SvtxVertexMap_v1;
    PHIODataNode<PHObject> *node =
        new PHIODataNode<PHObject>(m_actsVertexMap,
                                   "SvtxVertexMapActs", "PHObject");
    svtxNode->addNode(node);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
