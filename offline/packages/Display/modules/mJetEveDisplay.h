/*
        \file mJetEveDisplay.h
        \author Sookhyun Lee
        \brief reconstructed jets
        \version $Revision: 1.1 $
        \date    $Date: 07/26/2016
*/

#ifndef __MJETEVEDISPLAY_H__
#define __MJETEVEDISPLAY_H__

#include <mPHEveModuleBase.h>
#include <boost/shared_ptr.hpp>
#include <set>

class TEveManager;
class TEveElementList;
class PHEveDisplay;

class SvtxVertexMap;
class JetMap;

class mJetEveDisplay : public mPHEveModuleBase
{
 public:
  mJetEveDisplay(boost::shared_ptr<PHEveDisplay>);
  ~mJetEveDisplay();
  
  void init(PHCompositeNode* topNode);
  void init_run(PHCompositeNode* topNode);
  bool event(PHCompositeNode* topNode);
  void end(PHCompositeNode* topNode);
  void draw_event();
  
  void create_nodes(PHCompositeNode* topNode);
  void draw_jets();
  void clear();


private:
  boost::shared_ptr<PHEveDisplay> _evedisp;

  SvtxVertexMap* _vtxmap;
  JetMap* _jetmap;
  TEveElementList * _reco_jets;
  float radius;
  float length;
};

#endif // __MJETSEVEDISPLAY_H__
