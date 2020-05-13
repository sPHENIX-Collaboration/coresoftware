/*!
 *  \file       PHTpcVertexFinder.h
 *  \brief      
 *  \author     Dmitry Arkhipkin <arkhipkin@gmail.com>
 */

#ifndef PHTPCVERTEXFINDER_H_
#define PHTPCVERTEXFINDER_H_

#include <vector>

namespace PHGenFit2 { class Track; }
namespace genfit { class GFRaveVertex; }
namespace genfit { class GFRaveVertexFactory; }

/// \class PHTpcVertexFinder
///
/// \brief
///
class PHTpcVertexFinder
{
 public:
  PHTpcVertexFinder();
  ~PHTpcVertexFinder();

  std::vector<genfit::GFRaveVertex*> findVertices(std::vector<PHGenFit2::Track*>& gtracks);

 protected:
  genfit::GFRaveVertexFactory* _vertex_finder;

 private:
};

#endif /* PHTPCVERTEXFINDER_H_ */
