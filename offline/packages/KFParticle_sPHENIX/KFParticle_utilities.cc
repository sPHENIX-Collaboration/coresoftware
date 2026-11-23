#include "KFParticle_utilities.h"

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <phool/getClass.h>

#include <KFParticle.h>

#include <TDatabasePDG.h>
#include <TMatrixD.h>
#include <TMatrixDfwd.h>  // for TMatrixD
#include <TMatrixT.h>     // for TMatrixT, operator*
#include <TParticlePDG.h>

#include <cmath>     // for sqrt, pow
#include <cstdlib>   // for abs
#include <iostream>

SvtxTrack *KFParticle_utilities::getTrack(unsigned int track_id, SvtxTrackMap *trackmap)
{
  // Track ids are the keys of the map, so look the track up directly instead of scanning the whole map
  auto iter = trackmap->find(track_id);
  if (iter == trackmap->end())
  {
    return nullptr;
  }

  return iter->second;
}

int KFParticle_utilities::getPVID(PHCompositeNode *topNode, const KFParticle &kfpvertex,
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

void KFParticle_utilities::identify(const KFParticle &particle)
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

double KFParticle_utilities::calculateEllipsoidRadius(int posOrNeg, double sigma_ii, double sigma_jj, double sigma_ij)
{  // Note - Only works for a 2D ellipsoid OR rotated nD ellipsoid to avoid projections
  if (std::abs(posOrNeg) != 1)
  {
    std::cout << "You have set posOrNeg to " << posOrNeg << ". This value must be  +/- 1! Skipping" << std::endl;
    return 0;
  }

  double r_ij = sqrt((sigma_ii + sigma_jj) / 2 + posOrNeg * (sqrt(pow((sigma_ii - sigma_jj) / 2, 2) + pow(sigma_ij, 2))));

  return r_ij;
}

float KFParticle_utilities::eventDIRA(const KFParticle &particle, const KFParticle &vertex, bool do3D)
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

float KFParticle_utilities::flightDistanceChi2(const KFParticle &particle, const KFParticle &vertex)
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

float KFParticle_utilities::calculateEllipsoidVolume(const KFParticle &particle)
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

int KFParticle_utilities::getParticleID(const std::string &particle)
{
  return TDatabasePDG::Instance()->GetParticle(particle.c_str())->PdgCode();
}

float KFParticle_utilities::getParticleMass(const std::string &particle)
{
  return TDatabasePDG::Instance()->GetParticle(particle.c_str())->Mass();
}

float KFParticle_utilities::getParticleMass(const int PDGID)
{
  return TDatabasePDG::Instance()->GetParticle(PDGID)->Mass();
}

//Taken from PHTrackCleaner.cc but this was declared private there
TrackSeed *KFParticle_utilities::findSeed(TrackSeedContainer *container, unsigned short id)
{
  // find() does no bounds checking, so guard against unset (USHRT_MAX) or stale ids
  if (!container || id >= container->size())
  {
    return nullptr;
  }
  return *container->find(id);
}
