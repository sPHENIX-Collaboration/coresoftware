#ifndef ACTSTRACK_H
#define ACTSTRACK_H

#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/EventData/MultiTrajectory.hpp>

#include <ActsExamples/EventData/Track.hpp>
#include <ActsExamples/EventData/TrkrClusterSourceLink.hpp>
#include <ActsExamples/EventData/TrkrClusterMultiTrajectory.hpp>

using SourceLink = ActsExamples::TrkrClusterSourceLink;
using Trajectory = ActsExamples::TrkrClusterMultiTrajectory;

/**
 * A class that contains an Acts track seed and the corresponding source links.
 * Also contains a vertex id corresponding to the original track seed
 * vertex id.
 * Used by a variety of the Acts track fitting classes.
 */
class ActsTrack
{

 public:
  ActsTrack();

  ActsTrack(ActsExamples::TrackParameters params, 
	    const std::vector<SourceLink> &links,
	    Acts::Vector3D vertex)
    : m_trackParams(params)
    , m_sourceLinks(links)
    , m_vertex(vertex)
    {}

  ~ActsTrack(){}

  ActsExamples::TrackParameters getTrackParams()
    { return m_trackParams; }
  void setTrackParams(ActsExamples::TrackParameters params) 
    { m_trackParams = params; }

  std::vector<SourceLink> getSourceLinks(){ return m_sourceLinks; }
  void setSourceLinks(const std::vector<SourceLink> &srcLinks)
    { m_sourceLinks = srcLinks; }
  
  Acts::Vector3D  getVertex() { return m_vertex;}
  void setVertex(Acts::Vector3D vertex){m_vertex = vertex;}

 private:
  /// Initial track seed parameters
  ActsExamples::TrackParameters m_trackParams;

  /// SourceLinks associated to this track
  std::vector<SourceLink> m_sourceLinks;

  /// Initial x,y,z vertex estimate
  Acts::Vector3D m_vertex;

};

#endif
