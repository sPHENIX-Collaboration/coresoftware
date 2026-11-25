#ifndef KFPARTICLESPHENIX_KFPARTICLEUTIL_H
#define KFPARTICLESPHENIX_KFPARTICLEUTIL_H

#include <algorithm>  // for remove, stable_sort
#include <cstddef>    // for size_t
#include <map>
#include <numeric>    // for iota
#include <string>
#include <utility>    // for move
#include <vector>

class KFParticle;
struct KFParticle_Config;
class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class TrackSeed;
class TrackSeedContainer;
class Vertex;

// Stateless helpers shared across the package
namespace KFParticle_util
{
  // Constants

  // Size of the per-daughter branch arrays
  inline constexpr int maxTracks{20};

  // Speed of light in cm/ps: KFParticle decay times come out as lengths (c*tau, cm), divide by this to get ps
  inline constexpr float speedOfLight{2.99792458e-2};

  // Particle data (TDatabasePDG)

  // The particle must exist (no null check)
  int getParticleID(const std::string &particle);
  float getParticleMass(const std::string &particle);
  float getParticleMass(int PDGID);

  // False (with a message) if the name is not in TDatabasePDG
  bool findParticle(const std::string &particle);

  // Kinematics

  // Cosine of the angle between the flight direction (vertex -> particle) and the momentum
  float eventDIRA(const KFParticle &particle, const KFParticle &vertex, bool do3D = true);
  float flightDistanceChi2(const KFParticle &particle, const KFParticle &vertex);
  float calculateEllipsoidVolume(const KFParticle &particle);

  // Node and DST access

  // Track with this id, or nullptr if it is not in the map
  SvtxTrack *getTrack(unsigned int track_id, const SvtxTrackMap *trackmap);

  // Seed with this index, or nullptr for an unset or out-of-range id
  TrackSeed *findSeed(const TrackSeedContainer *container, unsigned short id);

  // TPC dE/dx of the daughter's track (looked up in the track map trackMapName); -1 if the track, its TPC seed
  // or one of the cluster/geometry nodes is missing
  float get_dEdx(PHCompositeNode *topNode, const KFParticle &daughter, const std::string &trackMapName);

  // MBD or SVTX vertices by id. MbdVertexMap and SvtxVertexMap have no common base class, but both hold
  // id -> vertex; copying the one in use into this map lets the rest of the code treat both alike.
  using DstVertices = std::map<unsigned int, const Vertex *>;

  // Copy of an MbdVertexMap or SvtxVertexMap; empty for a null map
  template <class VertexMap>
  DstVertices dstVertices(const VertexMap *vertexMap)
  {
    DstVertices vertices;
    if (vertexMap)
    {
      for (const auto &[id, vertex] : *vertexMap)
      {
        vertices.emplace(id, vertex);
      }
    }
    return vertices;
  }

  // Vertex with this id, or nullptr if it is not there
  const Vertex *findDstVertex(const DstVertices &vertices, unsigned int id);

  // The MBD measures only z. Its vertices have no x or y (MbdVertexv1 returns NaN, v2 and v3 return 0), no covariance
  // apart from the z error, and no chi2 or NDF, so they are read differently from SVTX vertices.
  bool isMbdVertex(const Vertex *vertex);

  // The vertex as a KFParticle (MBD: x = y = 0, only the z error, NDF = chi2 = 0)
  KFParticle makeVertex(const Vertex *vertex);

  // Every primary vertex of the event, from the MBD or SVTX vertices (config.useMbdVertex) through the
  // GlobalVertexMap (id = the global vertex id) or, with config.dontUseGlobalVertex, from the vertex map
  // directly (id = the map's id)
  std::vector<KFParticle> makeAllPrimaryVertices(PHCompositeNode *topNode, const KFParticle_Config &config);

  // Beam crossing of the vertex the candidate was built with; -100 if it can't be found
  int getPVID(PHCompositeNode *topNode, const KFParticle &kfpvertex,
              bool dontUseGlobalVertex, bool useMbdVertex, const std::string &vtxMapName);

  // Generic helpers

  // min <= value <= max
  bool isInRange(float min, float value, float max);

  // Removes duplicates in place, keeping the first occurrence of each element and the original order.
  // O(n log n) for long vectors; short vectors (up to 16 elements) use the simple quadratic loop.
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

  // Console output

  // Prints the ID, PDG ID, charge, mass, momentum and position (with errors) of a particle
  void identify(const KFParticle &particle);

  // Coloured pass/fail lines for the verbose selection printout
  void printSelectionCheck(const std::string &parameter, float min, float val, float max);
  void printSelectionCheck(const std::string &start, const std::string &accept, const std::string &reject, const std::string &end, bool equality);
  void printSelectionCheck(const std::string &info, unsigned int value);
}  // namespace KFParticle_util

#endif  // KFPARTICLESPHENIX_KFPARTICLEUTIL_H
