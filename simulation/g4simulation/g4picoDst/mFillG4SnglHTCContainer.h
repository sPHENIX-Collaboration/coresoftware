/*!
 *** \\ class mFillG4SnglHTCContainer for filling G4 hits, towers, clusters into a root file
 *** \\ author: Dr. Liang Xue from Georgia State University
*/

#ifndef MFILLG4SNGLHTCCONTAINER_H__
#define MFILLG4SNGLHTCCONTAINER_H__

#include <fun4all/SubsysReco.h>
#include <map>
#include <set>
#include <string>



// Forward declerations
class Fun4AllHistoManager;
class PHCompositeNode;
class PHG4HitContainer;
class RawTowerContainer;
class RawTowerGeomContainer;
class RawClusterContainer;

class mFillG4SnglHTCContainer: public SubsysReco
{
 public:

  //! constructor
  mFillG4SnglHTCContainer( const std::string &name = "mFillG4SnglHTCContainer");

  //! destructor
  virtual ~mFillG4SnglHTCContainer();

  //! full initialization
  int Init(PHCompositeNode *);

  //! event processing method
  int process_event(PHCompositeNode *);

  //! G4 hits processing method
  int process_hit(int detid, PHG4HitContainer *hits, G4SnglHit snglhit, G4SnglHTCContainer *snglhits);

  //! G4 towers processing method
  int process_twr(int detid, RawTowerContainer *twrs, RawTowerGeomContainer *twrgeom, G4SnglTower sngltwr, G4SnglHTCContainer *sngltwrs);

  //! G4 clusters processing method
  int process_clr(int detid, RawClusterContainer *clrs, G4SnglCluster snglclr, G4SnglHTCContainer *sngltcls);

  //! end of run method
  int End(PHCompositeNode *);

  void CreateNode(const std::string &name) { _nodename = name; }

  void AddNode(const std::string &name, const int detid=0);

  void Set_MakeHit(bool a=true) { make_hit = a; }
  void Set_MakeTower(bool a=true) { make_twr = a; }
  void Set_MakeCluster(bool a=true) { make_clr = a; }

 public:
 bool make_hit;
 bool make_twr;
 bool make_clr;
 int nevents;
 std::string _nodename;
 std::set<std::string> _node_postfix;
 std::map<std::string, int> _detid;

};

#endif 
