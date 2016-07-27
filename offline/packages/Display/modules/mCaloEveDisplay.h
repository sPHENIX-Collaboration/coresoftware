/*!
        \file mCaloEveDisplay.h
        \author Sookhyun Lee
        \brief reconstructed energy clusters from cemc/hcalin/hcalout
        \version $Revision: 1.1 $
        \date    $Date: 07/26/2016
*/

#ifndef __MCALOEVEDISPLAY_H__
#define __MCALOEVEDISPLAY_H__

#include <mPHEveModuleBase.h>
#include <boost/shared_ptr.hpp>
#include <set>

class TEveManager;
class TEveTrackPropagator;
class TEveElementList;
class TEveBoxSet;
class TEveRGBAPalette;

class PHEveDisplay;
class TH2F;

class RawClusterContainer;
class RawTowerContainer;
class RawTowerGeomContainer;

class mCaloEveDisplay : public mPHEveModuleBase
{
 public:
  mCaloEveDisplay(boost::shared_ptr<PHEveDisplay>);
  ~mCaloEveDisplay();
  
  void init(PHCompositeNode* topNode);
  void init_run(PHCompositeNode* topNode);
  bool event(PHCompositeNode* topNode);
  void end(PHCompositeNode* topNode);
  void draw_event();
  
  void create_nodes(PHCompositeNode* topNode);
  void draw_clusters(bool _is_cemc, bool _is_hcalin, bool _is_hcalout);
  void clear();

private:
  boost::shared_ptr<PHEveDisplay> _evedisp;

  TEveElementList* _cemc_list;
  TEveElementList* _hcalin_list;
  TEveElementList* _hcalout_list;
  
  TEveBoxSet* _cemc_boxset;
  TEveBoxSet* _hcalin_boxset;
  TEveBoxSet* _hcalout_boxset; 
  TEveRGBAPalette *_pal;

  RawClusterContainer *_cemc_clusters;
  RawClusterContainer *_hcalin_clusters;
  RawClusterContainer *_hcalout_clusters;
  RawTowerContainer *_cemc_towers;
  RawTowerContainer *_hcalin_towers;
  RawTowerContainer *_hcalout_towers;
  RawTowerGeomContainer *_cemc_towergeo;
  RawTowerGeomContainer *_hcalin_towergeo;
  RawTowerGeomContainer *_hcalout_towergeo;


};

#endif // __MCALOEVEDISPLAY_H__
