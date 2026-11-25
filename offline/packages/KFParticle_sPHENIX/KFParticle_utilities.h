#ifndef KFPARTICLESPHENIX_KFPARTICLEUTILITIES_H
#define KFPARTICLESPHENIX_KFPARTICLEUTILITIES_H

#include <algorithm>  // for remove, stable_sort
#include <cstddef>    // for size_t
#include <numeric>    // for iota
#include <string>
#include <utility>    // for move
#include <vector>

class KFParticle;
class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class TrackSeed;
class TrackSeedContainer;

// Stateless helpers shared across the package
namespace KFParticle_utilities
{
  // Size of the per-daughter branch arrays
  inline constexpr int max_tracks{20};

  // Speed of light in cm/ps: KFParticle decay times come out as lengths (c*tau, cm), divide by this to get ps
  inline constexpr float speedOfLight{2.99792458e-2};

  SvtxTrack *getTrack(unsigned int track_id, SvtxTrackMap *trackmap);

  // Beam crossing of the vertex the candidate was built with; -100 if it can't be found
  int getPVID(PHCompositeNode *topNode, const KFParticle &kfpvertex,
              bool dontUseGlobalVertex, bool useMbdVertex, const std::string &vtxMapName);

  // Prints the ID, PDG ID, charge, mass, momentum and position (with errors) of a particle
  void identify(const KFParticle &particle);

  // Note - Only works for a 2D ellipsoid OR rotated nD ellipsoid to avoid projections. No callers
  double calculateEllipsoidRadius(int posOrNeg, double sigma_ii, double sigma_jj, double sigma_ij);

  // Cosine of the angle between the flight direction (vertex -> particle) and the momentum
  float eventDIRA(const KFParticle &particle, const KFParticle &vertex, bool do3D = true);
  float flightDistanceChi2(const KFParticle &particle, const KFParticle &vertex);
  float calculateEllipsoidVolume(const KFParticle &particle);

  // TDatabasePDG lookups; the particle must exist (no null check)
  int getParticleID(const std::string &particle);
  float getParticleMass(const std::string &particle);
  float getParticleMass(int PDGID);

  // False (with a message) if the name is not in TDatabasePDG
  bool findParticle(const std::string &particle);

  // Seed with this index, or nullptr for an unset or out-of-range id
  TrackSeed *findSeed(TrackSeedContainer *container, unsigned short id);

  // min <= value <= max
  bool isInRange(float min, float value, float max);

  // Removes duplicates keeping the first occurrence of each element and the original order, exactly like the
  // previous quadratic implementation, but in O(n log n) for long vectors. Short vectors use the simple loop.
  template <typename T>
  void removeDuplicates(std::vector<T> &v)
  {
    const std::size_t n = v.size();
    if (n < 2)
    {
      return;
    }
    if (n <= 16)
    {
      auto end = v.end();
      for (auto it = v.begin(); it != end; ++it)
      {
        end = std::remove(it + 1, end, *it);
      }
      v.erase(end, v.end());
      return;
    }

    std::vector<std::size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(), [&v](std::size_t a, std::size_t b) { return v[a] < v[b]; });

    // In each run of equal elements the first entry of the stable sort has the smallest original position
    std::vector<char> keep(n, 0);
    for (std::size_t i = 0; i < n; ++i)
    {
      if (i == 0 || !(v[idx[i]] == v[idx[i - 1]]))
      {
        keep[idx[i]] = 1;
      }
    }

    std::size_t out = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
      if (keep[i])
      {
        if (out != i)
        {
          v[out] = std::move(v[i]);
        }
        ++out;
      }
    }
    v.erase(v.begin() + out, v.end());
  }

  // Coloured pass/fail lines for the verbose selection printout
  void printSelectionCheck(const std::string &parameter, float min, float val, float max);
  void printSelectionCheck(const std::string &start, const std::string &accept, const std::string &reject, const std::string &end, bool equality);
  void printSelectionCheck(const std::string &info, unsigned int value);
}  // namespace KFParticle_utilities

#endif  // KFPARTICLESPHENIX_KFPARTICLEUTILITIES_H
