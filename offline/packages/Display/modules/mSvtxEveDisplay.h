/*
        \file mSvtxEveDisplay.h
        \author Sookhyun Lee
        \brief reconstructed charged tracks and their clusters
        \version $Revision: 1.1 $
        \date    $Date: 07/26/2016
*/

#ifndef __MSVTXEVEDISPLAY_H__
#define __MSVTXEVEDISPLAY_H__

#include <mPHEveModuleBase.h>
#include <boost/shared_ptr.hpp>
#include <set>

class TEveManager;
class TEveTrackPropagator;
class TEveElementList;
class TEveTrackList;
class PHEveDisplay;
class TH2F;

class SvtxVertexMap;
class SvtxTrackMap;
class SvtxClusterMap;

class mSvtxEveDisplay : public mPHEveModuleBase
{
 public:
  mSvtxEveDisplay(boost::shared_ptr<PHEveDisplay>);
  ~mSvtxEveDisplay();
  
  void init(PHCompositeNode* topNode);
  void init_run(PHCompositeNode* topNode);
  bool event(PHCompositeNode* topNode);
  void end(PHCompositeNode* topNode);
  void draw_event();
  
  void create_nodes(PHCompositeNode* topNode);
  void draw_tracks();
  bool pid_cut(int pid);
  void clear();
  
private:

  boost::shared_ptr<PHEveDisplay> _evedisp;
  SvtxVertexMap *_vertexmap;
  SvtxTrackMap *_trackmap;
  SvtxClusterMap *_clustermap;

  TEveTrackPropagator* _prop;
  TEveTrackList* _svtx_tracks;
  
};

#endif // __MSVTXEVEDISPLAY_H__
