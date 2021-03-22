/*!
 *  \file       PHTpcEventExporter.h
 *  \brief      
 *  \author     Dmitry Arkhipkin <arkhipkin@bnl.gov>
 */

#ifndef PHTPCEVENTEXPORTER_H_
#define PHTPCEVENTEXPORTER_H_

#include <string>  // for string
#include <vector>  // for vector

class TrkrClusterContainer;
class TrkrHitSetContainer;

namespace PHGenFit2
{
  class Track;
}
namespace kdfinder
{
  template <class T>
  class TrackCandidate;
}

/// \class PHTpcEventExporter
///
/// \brief
///
class PHTpcEventExporter
{
 public:
  PHTpcEventExporter();
  virtual ~PHTpcEventExporter() {}

  void exportEvent(TrkrClusterContainer* cluster_map, TrkrHitSetContainer* hitsets, std::vector<kdfinder::TrackCandidate<double>*> candidates,
                   double B, const std::string& filename);

  void exportEvent(TrkrClusterContainer* cluster_map, TrkrHitSetContainer* hitsets,std::vector<PHGenFit2::Track*> gtracks,
                   double B, const std::string& filename);

  void exportEvent(std::vector<PHGenFit2::Track*> gtracks, double B, const std::string& filename);

 protected:
 private:
};

#endif /* PHTPCEVENTEXPORTER_H_ */
