#ifndef ACTSTRACK_H
#define ACTSTRACK_H

#include <Acts/EventData/TrackParameters.hpp>

#include <ACTFW/EventData/TrkrClusterSourceLink.hpp>

using SourceLink = FW::Data::TrkrClusterSourceLink;


/**
 * Aclass that contains an Acts track seed and the corresponding source links.
 * Also contains a vertex id corresponding to the original track seed
 * vertex id.
 * Used by a variety of the Acts track fitting classes.
 */
class ActsTrack
{

 public:
  ActsTrack(FW::TrackParameters params, 
	    const std::vector<SourceLink> &links,
	    unsigned int vertexId)
    : m_trackParams(params)
    , m_sourceLinks(links)
    , m_vertexId(vertexId)
    {}

  ~ActsTrack(){}

  FW::TrackParameters getTrackParams(){ return m_trackParams; }
  void setTrackParams(FW::TrackParameters params) { m_trackParams = params; }

  std::vector<SourceLink> getSourceLinks(){ return m_sourceLinks; }
  void setSourceLinks(std::vector<SourceLink> srcLinks)
  { m_sourceLinks = srcLinks; }
  
  unsigned int getVertexId() { return m_vertexId;}
  void setVertexId(unsigned int id){m_vertexId = id;}

 private:
  FW::TrackParameters m_trackParams;
  std::vector<SourceLink> m_sourceLinks;
  unsigned int m_vertexId;

};

#endif
