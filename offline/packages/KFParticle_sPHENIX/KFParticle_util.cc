#include "KFParticle_util.h"

#include "KFParticle_Config.h"

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackAnalysisUtils.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrClusterContainer.h>

#include <phool/getClass.h>

#include <KFParticle.h>

#include <TDatabasePDG.h>
#include <TMatrixD.h>
#include <TMatrixDfwd.h>  // for TMatrixD
#include <TMatrixT.h>     // for TMatrixT, operator*
#include <TParticlePDG.h>

#include <cmath>     // for sqrt
#include <cstdlib>   // for abs
#include <iostream>

namespace
{
  // ANSI colour codes for printSelectionCheck
  const std::string accept_colour = "32";
  const std::string reject_colour = "31";
}  // namespace

SvtxTrack *KFParticle_util::getTrack(unsigned int track_id, const SvtxTrackMap *trackmap)
{
  // Track ids are the keys of the map, so look the track up directly instead of scanning the whole map
  auto iter = trackmap->find(track_id);
  if (iter == trackmap->end())
  {
    return nullptr;
  }

  return iter->second;
}

const Vertex *KFParticle_util::findDstVertex(const DstVertices &vertices, unsigned int id)
{
  auto iter = vertices.find(id);
  return iter == vertices.end() ? nullptr : iter->second;
}

bool KFParticle_util::isMbdVertex(const Vertex *vertex)
{
  return dynamic_cast<const MbdVertex *>(vertex) != nullptr;
}

KFParticle KFParticle_util::makeVertex(const Vertex *vertex)
{
  const bool isMbd = isMbdVertex(vertex);

  float vtxX = isMbd ? 0 : vertex->get_x();
  float vtxY = isMbd ? 0 : vertex->get_y();
  float vtxZ = vertex->get_z();

  float f_vertexParameters[6] = {vtxX, vtxY, vtxZ, 0, 0, 0};

  float f_vertexCovariance[21] = {0};
  unsigned int iterate = 0;
  if (isMbd)
  {
    f_vertexCovariance[5] = vertex->get_z_err();
  }
  else
  {
    for (unsigned int i = 0; i < 3; ++i)
    {
      for (unsigned int j = 0; j <= i; ++j)
      {
        f_vertexCovariance[iterate] = vertex->get_error(i, j);
        ++iterate;
      }
    }
  }

  KFParticle kfp_vertex;
  kfp_vertex.Create(f_vertexParameters, f_vertexCovariance, 0, -1);
  kfp_vertex.NDF() = isMbd ? 0 : vertex->get_ndof();
  kfp_vertex.Chi2() = isMbd ? 0 : vertex->get_chisq();

  return kfp_vertex;
}

std::vector<KFParticle> KFParticle_util::makeAllPrimaryVertices(const DstVertices &dstVertices, const GlobalVertexMap *globalVertexMap,
                                                                const KFParticle_Config &config)
{
  std::vector<KFParticle> primaryVertices;

  if (config.dontUseGlobalVertex)
  {
    for (const auto &[id, vertex] : dstVertices)
    {
      primaryVertices.push_back(makeVertex(vertex));
      primaryVertices.back().SetId(id);
    }

    return primaryVertices;
  }

  if (!globalVertexMap)
  {
    std::cout << "Can't continue in KFParticle_util::makeAllPrimaryVertices" << std::endl;
    return primaryVertices;
  }

  for (const auto &entry : *globalVertexMap)
  {
    const GlobalVertex *globalVertex = entry.second;

    GlobalVertex::VTXTYPE whichVtx = config.useMbdVertex ? GlobalVertex::MBD : GlobalVertex::SVTX;

    auto svtxiter = globalVertex->find_vertexes(whichVtx);
    // check that it contains a track vertex
    if (svtxiter == globalVertex->end_vertexes())
    {
      continue;
    }

    const auto &svtxvertexvector = svtxiter->second;

    for (const auto &vertex : svtxvertexvector)
    {
      primaryVertices.push_back(makeVertex(KFParticle_util::findDstVertex(dstVertices, vertex->get_id())));
      primaryVertices.back().SetId(globalVertex->get_id());
    }
  }

  return primaryVertices;
}

std::vector<KFParticle> KFParticle_util::makeAllPrimaryVertices(PHCompositeNode *topNode, const KFParticle_Config &config)
{
  const DstVertices vertexById = config.useMbdVertex
                                     ? dstVertices(findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap"))
                                     : dstVertices(findNode::getClass<SvtxVertexMap>(topNode, config.vertexMapNodeName));

  const GlobalVertexMap *globalVertexMap = config.dontUseGlobalVertex ? nullptr : findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");

  return makeAllPrimaryVertices(vertexById, globalVertexMap, config);
}

float KFParticle_util::get_dEdx(PHCompositeNode *topNode, const KFParticle &daughter, const std::string &trackMapName)
{
  SvtxTrackMap *trackmap = findNode::getClass<SvtxTrackMap>(topNode, trackMapName);
  TrkrClusterContainer *clusterMap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!clusterMap)
  {
    clusterMap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER_SEED");
  }
  auto *geomContainer = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  auto *geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!clusterMap || !geomContainer || !geometry)
  {
    return -1.0;
  }

  SvtxTrack *daughter_track = KFParticle_util::getTrack(daughter.Id(), trackmap);
  if (!daughter_track)
  {
    // No backing SvtxTrack -- e.g. daughter is itself a built intermediate.
    // dE/dx PID is only meaningful for a real reconstructed track.
    return -1.0;
  }

  bool isTrack_v5 = daughter_track->InheritsFrom("SvtxTrack_v5");
  TrackSeedContainer *tpcSeedContainer = isTrack_v5 ? findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer") : nullptr;

  TrackSeed *tpcseed = isTrack_v5 ? KFParticle_util::findSeed(tpcSeedContainer, daughter_track->get_tpc_seed_id()) : daughter_track->get_tpc_seed();
  if (!tpcseed)
  {
    return -1;
  }

  float layerThicknesses[4] = {0.0, 0.0, 0.0, 0.0};
  // These are randomly chosen layer thicknesses for the TPC, to get the
  // correct region thicknesses in an easy to pass way to the helper fxn
  layerThicknesses[0] = geomContainer->GetLayerCellGeom(7)->get_thickness();
  layerThicknesses[1] = geomContainer->GetLayerCellGeom(8)->get_thickness();
  layerThicknesses[2] = geomContainer->GetLayerCellGeom(27)->get_thickness();
  layerThicknesses[3] = geomContainer->GetLayerCellGeom(50)->get_thickness();

  return TrackAnalysisUtils::calc_dedx(tpcseed, clusterMap, geometry, layerThicknesses);
}

int KFParticle_util::getPVID(PHCompositeNode *topNode, const KFParticle &kfpvertex,
                             bool dontUseGlobalVertex, bool useMbdVertex, const std::string &vtxMapName)
{
  if (dontUseGlobalVertex)
  {
    if (useMbdVertex)
    {
      auto *dst_mbdvertexmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
      if (dst_mbdvertexmap)
      {
        MbdVertex *m_dst_vertex = dst_mbdvertexmap->get(kfpvertex.Id());
        return m_dst_vertex->get_beam_crossing();
      }

      std::cout << "KFParticle vertex matching: " << vtxMapName << " does not exist" << std::endl;
    }
    else
    {
      auto *dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, vtxMapName);
      if (dst_vertexmap)
      {
        SvtxVertex *m_dst_vertex = dst_vertexmap->get(kfpvertex.Id());
        return m_dst_vertex->get_beam_crossing();
      }

      std::cout << "KFParticle vertex matching: " << vtxMapName << " does not exist" << std::endl;
    }
  }
  else
  {
    auto *globalvertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
    if (!globalvertexmap)
    {
      return -100;
    }

    GlobalVertex *gvertex = globalvertexmap->get(kfpvertex.Id());
    return gvertex->get_beam_crossing();
  }

  return -100;
}

void KFParticle_util::identify(const KFParticle &particle)
{
  std::cout << "Track ID: " << particle.Id() << std::endl;
  std::cout << "PDG ID: " << particle.GetPDG() << ", charge: " << (int) particle.GetQ() << ", mass: " << particle.GetMass() << " GeV" << std::endl;
  std::cout << "(px,py,pz) = (" << particle.GetPx() << " +/- " << std::sqrt(particle.GetCovariance(3, 3)) << ", ";
  std::cout << particle.GetPy() << " +/- " << std::sqrt(particle.GetCovariance(4, 4)) << ", ";
  std::cout << particle.GetPz() << " +/- " << std::sqrt(particle.GetCovariance(5, 5)) << ") GeV" << std::endl;
  std::cout << "(x,y,z) = (" << particle.GetX() << " +/- " << std::sqrt(particle.GetCovariance(0, 0)) << ", ";
  std::cout << particle.GetY() << " +/- " << std::sqrt(particle.GetCovariance(1, 1)) << ", ";
  std::cout << particle.GetZ() << " +/- " << std::sqrt(particle.GetCovariance(2, 2)) << ") cm\n" << std::endl;
}

float KFParticle_util::eventDIRA(const KFParticle &particle, const KFParticle &vertex, bool do3D)
{
  const int nDimensions = do3D ? 3 : 2;
  TMatrixD flightVector(nDimensions, 1);
  TMatrixD momVector(nDimensions, 1);
  flightVector(0, 0) = particle.GetX() - vertex.GetX();
  flightVector(1, 0) = particle.GetY() - vertex.GetY();

  momVector(0, 0) = particle.GetPx();
  momVector(1, 0) = particle.GetPy();

  if (do3D)
  {
    flightVector(2, 0) = particle.GetZ() - vertex.GetZ();
    momVector(2, 0) = particle.GetPz();
  }

  TMatrixD momDotFD(1, 1);  // Calculate momentum dot flight distance
  momDotFD = TMatrixD(momVector, TMatrixD::kTransposeMult, flightVector);
  float f_momDotFD = momDotFD(0, 0);

  TMatrixD sizeOfMom(1, 1);  // Calculates the size of the momentum vector
  sizeOfMom = TMatrixD(momVector, TMatrixD::kTransposeMult, momVector);
  float f_sizeOfMom = sqrt(sizeOfMom(0, 0));

  TMatrixD sizeOfFD(1, 1);  // Calculates the size of the flight distance vector
  sizeOfFD = TMatrixD(flightVector, TMatrixD::kTransposeMult, flightVector);
  float f_sizeOfFD = sqrt(sizeOfFD(0, 0));

  return f_momDotFD / (f_sizeOfMom * f_sizeOfFD);
}

float KFParticle_util::flightDistanceChi2(const KFParticle &particle, const KFParticle &vertex)
{
  TMatrixD flightVector(3, 1);
  TMatrixD flightDistanceCovariance(3, 3);

  const KFParticle &kfp_vertex(vertex);

  flightVector(0, 0) = particle.GetX() - kfp_vertex.GetX();
  flightVector(1, 0) = particle.GetY() - kfp_vertex.GetY();
  flightVector(2, 0) = particle.GetZ() - kfp_vertex.GetZ();

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      flightDistanceCovariance(i, j) = particle.GetCovariance(i, j) + kfp_vertex.GetCovariance(i, j);
    }
  }

  TMatrixD anInverseMatrix(3, 3);
  anInverseMatrix = flightDistanceCovariance.Invert();
  TMatrixD m_chi2Value(1, 1);
  m_chi2Value = TMatrixD(flightVector, TMatrixD::kTransposeMult, anInverseMatrix * flightVector);
  return m_chi2Value(0, 0);
}

float KFParticle_util::calculateEllipsoidVolume(const KFParticle &particle)
{
  TMatrixD cov_matrix(3, 3);

  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      cov_matrix(i, j) = particle.GetCovariance(i, j);
    }
  }

  float volume;
  if (cov_matrix(0, 0) * cov_matrix(1, 1) * cov_matrix(2, 2) == 0)
  {
    volume = 0;
  }
  else
  {
    volume = (4. / 3.) * M_PI * sqrt((std::abs(cov_matrix.Determinant())));  // The covariance matrix is error-squared
  }

  return volume;
}

int KFParticle_util::getParticleID(const std::string &particle)
{
  return TDatabasePDG::Instance()->GetParticle(particle.c_str())->PdgCode();
}

float KFParticle_util::getParticleMass(const std::string &particle)
{
  return TDatabasePDG::Instance()->GetParticle(particle.c_str())->Mass();
}

float KFParticle_util::getParticleMass(const int PDGID)
{
  return TDatabasePDG::Instance()->GetParticle(PDGID)->Mass();
}

bool KFParticle_util::findParticle(const std::string &particle)
{
  bool particleFound = true;
  if (!TDatabasePDG::Instance()->GetParticle(particle.c_str()))
  {
    std::cout << "The particle, " << particle << " is not recognized" << std::endl;
    std::cout << "Check TDatabasePDG for a list of available particles" << std::endl;
    particleFound = false;
  }

  return particleFound;
}

//Taken from PHTrackCleaner.cc but this was declared private there
TrackSeed *KFParticle_util::findSeed(const TrackSeedContainer *container, unsigned short id)
{
  // find() does no bounds checking, so guard against unset (USHRT_MAX) or stale ids
  if (!container || id >= container->size())
  {
    return nullptr;
  }
  return *container->find(id);
}

bool KFParticle_util::isInRange(float min, float value, float max)
{
  return min <= value && value <= max;
}

void KFParticle_util::printSelectionCheck(const std::string &parameter, float min, float val, float max)
{
  std::string trailer = "the " + parameter + " requirement\033[0m";
  std::string passOrFail = isInRange(min, val, max) ? "\033[1;" + accept_colour + "mPassed " + trailer
                                                    : "\033[1;" + reject_colour + "mFailed " + trailer;
  std::cout << passOrFail << ". Lower bound = " << min << ", measured value =  " << val << ", upper bound = " << max << std::endl; 
}

void KFParticle_util::printSelectionCheck(const std::string &start, const std::string &accept, const std::string &reject, const std::string &end, bool equality)
{
  std::string decision = equality ? accept : reject;
  std::string colour = equality ? accept_colour : reject_colour;
  std::string spacing = start.empty() ? "" : " ";
  std::cout << "\033[1;" << colour << "m" << start << spacing << decision << " " << end << "\033[0m" << std::endl;
}

void KFParticle_util::printSelectionCheck(const std::string &info, unsigned int value)
{
  std::string colour = value > 0 ? accept_colour : reject_colour;
  std::cout << info << " = \033[1;" << colour << "m" << value << "\033[0m" << std::endl;
}
