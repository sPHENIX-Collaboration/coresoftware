#ifndef KFPARTICLESPHENIX_KFPARTICLEUTILITIES_H
#define KFPARTICLESPHENIX_KFPARTICLEUTILITIES_H

#include <string>

class KFParticle;
class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class TrackSeed;
class TrackSeedContainer;

// Stateless helpers shared across the package
class KFParticle_utilities
{
 public:
  KFParticle_utilities() = delete;

  // Size of the per-daughter branch arrays
  static constexpr int max_tracks{20};

  // Speed of light in cm/ps: KFParticle decay times come out as lengths (c*tau, cm), divide by this to get ps
  static constexpr float speedOfLight{2.99792458e-2};

  static SvtxTrack *getTrack(unsigned int track_id, SvtxTrackMap *trackmap);

  // Beam crossing of the vertex the candidate was built with; -100 if it can't be found
  static int getPVID(PHCompositeNode *topNode, const KFParticle &kfpvertex,
                     bool dontUseGlobalVertex, bool useMbdVertex, const std::string &vtxMapName);

  // Prints the ID, PDG ID, charge, mass, momentum and position (with errors) of a particle
  static void identify(const KFParticle &particle);

  // Note - Only works for a 2D ellipsoid OR rotated nD ellipsoid to avoid projections. No callers
  static double calculateEllipsoidRadius(int posOrNeg, double sigma_ii, double sigma_jj, double sigma_ij);

  // Cosine of the angle between the flight direction (vertex -> particle) and the momentum
  static float eventDIRA(const KFParticle &particle, const KFParticle &vertex, bool do3D = true);
  static float flightDistanceChi2(const KFParticle &particle, const KFParticle &vertex);
  static float calculateEllipsoidVolume(const KFParticle &particle);

  // TDatabasePDG lookups; the particle must exist (no null check)
  static int getParticleID(const std::string &particle);
  static float getParticleMass(const std::string &particle);
  static float getParticleMass(int PDGID);

  // False (with a message) if the name is not in TDatabasePDG
  static bool findParticle(const std::string &particle);

  // Seed with this index, or nullptr for an unset or out-of-range id
  static TrackSeed *findSeed(TrackSeedContainer *container, unsigned short id);
};

#endif  // KFPARTICLESPHENIX_KFPARTICLEUTILITIES_H
