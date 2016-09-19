/*!
        \file mG4EveDisplay.h
        \author Sookhyun Lee
        \brief true tracks and true jets from truthinfo container
        \version $Revision: 1.1 $
        \date    $Date: 07/26/2016
*/

#ifndef __MG4EVEDISPLAY_H__
#define __MG4EVEDISPLAY_H__

#include <mPHEveModuleBase.h>
#include <boost/shared_ptr.hpp>
#include <set>

class TEveManager;
class TEveTrackPropagator;
class TEveElementList;
class TEveTrackList;
class PHEveDisplay;
class TH2F;

class PHG4TruthInfoContainer;
class SvtxEvalStack;
class PHG4VtxPoint;
class JetMap;

class mG4EveDisplay : public mPHEveModuleBase
{

 public:

  mG4EveDisplay(boost::shared_ptr<PHEveDisplay>);
  ~mG4EveDisplay();
  
  void init(PHCompositeNode* topNode);
  void init_run(PHCompositeNode* topNode);
  bool event(PHCompositeNode* topNode);
  void end(PHCompositeNode* topNode);
  
  void create_nodes(PHCompositeNode* topNode);
  void draw_ideal_tracks();
  void draw_g4_tracks();
  void draw_jets();
  void clear();

  
 private:

  boost::shared_ptr<PHEveDisplay> _evedisp;
  PHG4TruthInfoContainer* _truth;
  SvtxEvalStack* _svtxevalstack;
  JetMap *_jetmap;

  TEveTrackPropagator* _prop;
  TEveTrackList* _true_tracks;
  TEveElementList* _true_jets;

  float radius;
  float length;

  int verbosity;  

};

#endif // __MG4EVEDISPLAY_H__
