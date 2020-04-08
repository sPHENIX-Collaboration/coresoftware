#ifndef ACTSTRACK_H
#define ACTSTRACK_H

#include <Acts/EventData/TrackParameters.hpp>

#include <ACTFW/EventData/TrkrClusterSourceLink.hpp>

using SourceLink = FW::Data::TrkrClusterSourceLink;


/**
 * A struct that contains an Acts track seed and the corresponding source links.
 * Used by a variety of the Acts track fitting classes.
 */
class ActsTrack
{

 public:
  ActsTrack(FW::TrackParameters params, const std::vector<SourceLink> &links)
    : m_trackParams(params)
    , m_sourceLinks(links)
    {}

  ~ActsTrack(){}

  FW::TrackParameters getTrackParams(){ return m_trackParams; }
  void setTrackParams(FW::TrackParameters params) { m_trackParams = params; }

  std::vector<SourceLink> getSourceLinks(){ return m_sourceLinks; }
  void setSourceLinks(std::vector<SourceLink> srcLinks)
  { m_sourceLinks = srcLinks; }
  

 private:
  
  FW::TrackParameters m_trackParams;
  std::vector<SourceLink> m_sourceLinks;

};

#endif
