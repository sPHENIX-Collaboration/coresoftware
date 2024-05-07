#ifndef TRACKRECO_PHACTSINITIALVERTEXFINDER_H
#define TRACKRECO_PHACTSINITIALVERTEXFINDER_H

#include "PHInitVertexing.h"
#include <trackbase/ActsTrackingGeometry.h>

#include <trackbase/TrkrDefs.h>

#include <Acts/Utilities/Result.hpp>
#include <Acts/Vertexing/Vertex.hpp>

#include <ActsExamples/EventData/Track.hpp>
#include <ActsExamples/EventData/TrkrClusterMultiTrajectory.hpp>

#include <random>

class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class SvtxVertexMap;
class SvtxVertex;

using VertexVector = std::vector<Acts::Vertex<Acts::BoundTrackParameters>>;

using TrackParamVec = std::vector<const Acts::BoundTrackParameters*>;

using InitKeyMap = std::map<const Acts::BoundTrackParameters*, const unsigned int>;

using CentroidMap = std::map<unsigned int, std::vector<SvtxTrack*>>;

class PHActsInitialVertexFinder: public PHInitVertexing
{
 public: 
  PHActsInitialVertexFinder(const std::string& name="PHActsInitialVertexFinder");
  ~PHActsInitialVertexFinder() override {}

  void setMaxVertices(const int maxVertices)
  { m_maxVertices = maxVertices;}

  void setSvtxVertexMapName(const std::string& name)
  { m_svtxVertexMapName = name; }
  
  void setSvtxTrackMapName(const std::string& name)
  { m_svtxTrackMapName = name; }

  void disablePtWeights(const bool weight)
  { m_disableWeights = weight; }

  void resetTrackCovariance(const bool initial)
  { m_resetTrackCovariance = initial; }

  void setCentroids(const int centroids)
  { m_nCentroids = centroids;}

  void setIterations(const int iterations)
  {m_nIterations = iterations;}

  void removeSiliconSeeds(const bool removeSeeds)
  {m_removeSeeds = removeSeeds;}
  void magFieldName(const std::string& magfield)
  {m_magField = magfield;}

  void setPCACut(const float pcaCut)
  {m_pcaCut = pcaCut;}

 protected:
  int Setup(PHCompositeNode *topNode) override;
  int Process(PHCompositeNode *topNode) override;
  int ResetEvent(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 private:
  int getNodes(PHCompositeNode *topNode);
  int createNodes(PHCompositeNode *topNode);
  
  /// Gets silicon stubs to send to Acts IVF
  TrackParamVec getTrackPointers(InitKeyMap& keyMap);

  /// Calls Acts IVF
  VertexVector findVertices(TrackParamVec& tracks);

  /// Creates SvtxVertexMap
  void fillVertexMap(VertexVector& vertices, InitKeyMap& keyMap);

  /// Makes a dummy vertex only if Acts returns 0 vertices or there 
  /// are either 0 or 1 silicon seeds
  void createDummyVertex(const float x, const float y, const float z);

  /// Assigns silicon seed a vertex ID if it was left out of Acts IVF
  void checkTrackVertexAssociation();
  
  /// Implements a k-means cluster algorithm to identify bad seeds
  /// to remove from Acts initial vertexing
  std::vector<SvtxTrack*> sortTracks();
  
  /// Helper functions for the k-means cluster algorithm
  CentroidMap createCentroidMap(std::vector<Acts::Vector3D>& centroids);

  std::vector<SvtxTrack*> getIVFTracks(CentroidMap& clusters, 
				       std::vector<Acts::Vector3D>& centroids);
  
  /// Number of centroids for k-means clustering algorithm
  int m_nCentroids = 5;
  /// Number of times to iterate for clusters to converge
  int m_nIterations = 15;
  /// Max number of vertices allowed by the Acts IVF
  int m_maxVertices = 5;
  /// Maximum centroid transverse PCA cut
  float m_pcaCut = 0.03; // cm
  /// Event num
  int m_event = 0;
  /// Diagnostic vertex numbers
  unsigned int m_totVertexFits = 0;
  unsigned int m_successFits = 0;

  unsigned int m_seed = 0;
  std::mt19937 m_random_number_generator;

  std::string m_svtxTrackMapName = "SvtxSiliconTrackMap";
  std::string m_svtxVertexMapName = "SvtxVertexMap";
  std::string m_magField = "";
  bool m_resetTrackCovariance = true;
  bool m_disableWeights = true;
  bool m_removeSeeds = false;

  SvtxTrackMap *m_trackMap = nullptr;
  SvtxVertexMap *m_vertexMap = nullptr;
  ActsTrackingGeometry *m_tGeometry = nullptr;

};


#endif
