// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHSILICONTPCTRACKMATCHING_H
#define PHSILICONTPCTRACKMATCHING_H

#include <fun4all/SubsysReco.h>
#include <trackbase/ActsGeometry.h>
#include <phparameter/PHParameterInterface.h>

#include <string>
#include <map>

class PHCompositeNode;
class TrackSeedContainer;
class TrackSeed;
class TrkrClusterContainer;
class TF1;
class TrkrClusterCrossingAssoc;

class PHSiliconTpcTrackMatching : public SubsysReco, public PHParameterInterface
{
 public:

  PHSiliconTpcTrackMatching(const std::string &name = "PHSiliconTpcTrackMatching");

  ~PHSiliconTpcTrackMatching() override;

 void SetDefaultParameters() override;

  void set_phi_search_window(const double win){_phi_search_win = win;}
  void set_eta_search_window(const double win){_eta_search_win = win;}
  void set_x_search_window(const double win){_x_search_win = win;}
  void set_y_search_window(const double win){_y_search_win = win;}
  void set_z_search_window(const double win){_z_search_win = win;}

  void set_test_windows_printout(const bool test){_test_windows = test ;}
  void set_pp_mode(const bool flag){_pp_mode = flag ;}
  void set_use_intt_time(const bool flag){_use_intt_time = flag ;}

  int InitRun(PHCompositeNode* topNode) override;

  int process_event(PHCompositeNode*) override;

  int End(PHCompositeNode*) override;

  void set_silicon_track_map_name(const std::string &map_name) { _silicon_track_map_name = map_name; }
  void set_track_map_name(const std::string &map_name) { _track_map_name = map_name; }
  void SetIteration(int iter){_n_iteration = iter;}
 private:

  int GetNodes(PHCompositeNode* topNode);

  void findEtaPhiMatches( std::set<unsigned int> &tpc_matched_set,
                            std::set<unsigned int> &tpc_unmatched_set,
			    std::multimap<unsigned int, unsigned int> &tpc_matches );
  std::vector<short int> getInttCrossings(TrackSeed *si_track);
   void checkCrossingMatches( std::multimap<unsigned int, unsigned int> &tpc_matches);
   short int getCrossingIntt(TrackSeed *_tracklet_si);

   //   void checkCrossingMatches( std::multimap<short int, std::pair<unsigned int, unsigned int>> &crossing_matches,  std::map<unsigned int, short int> &tpc_crossing_map );
  //double getBunchCrossing(unsigned int trid, double z_mismatch);
  //double getMedian(std::vector<double> &v);
  //void addSiliconClusters( std::multimap<short int, std::pair<unsigned int, unsigned int>> &crossing_matches);
  //void addSiliconClusters(  std::multimap<unsigned int, unsigned int> &tpc_matches);
  //void tagInTimeTracks(  std::multimap<unsigned int, unsigned int> &tpc_matches,
  //			 std::multimap<int, std::pair<unsigned int, unsigned int>> &crossing_matches,
  //			 std::map<unsigned int, int> &tpc_crossing_map );
  //void tagMatchCrossing( std::multimap<unsigned int, unsigned int> &tpc_matches,
  //			 std::multimap<short int, std::pair<unsigned int, unsigned int>> &crossing_matches,
  //			 std::map<unsigned int, short int> &tpc_crossing_map );
  //   void copySiliconClustersToCorrectedMap( );
  //void correctTpcClusterZIntt(  std::map<unsigned int, short int> &tpc_crossing_map );
   //void getMatchCrossingIntt(  
   //			       std::multimap<unsigned int, unsigned int> &tpc_matches,
   //			       std::multimap<short int, std::pair<unsigned int, unsigned int>> &crossing_matches,
   //			       std::map<unsigned int, short int> &tpc_crossing_map );
   //  void addTrackBunchCrossing(std::multimap<unsigned int, unsigned int> &tpc_matches);	  
    //  void addTrackBunchCrossing( std::map<unsigned int, short int> &tpc_crossing_map);	  
  //  void addTrackBunchCrossing(
   //						   std::map<unsigned int, short int> &vertex_crossings_map,
   //						   std::multimap<unsigned int, std::pair<unsigned int, unsigned int>>  &vertex_map);	  
  //  void cleanVertexMap( std::map<unsigned int, short int> &vertex_crossings_map,
  //		       std::multimap<unsigned int, std::pair<unsigned int, unsigned int>>  &vertex_map,
  //		       std::map<unsigned int, short int> &tpc_crossing_map );
  // void getCrossingNumber( std::vector<double> &vertex_list,
  //			    std::multimap<unsigned int, std::pair<unsigned int, unsigned int>>  &vertex_map, 
  //			    std::map<unsigned int, short int> &vertex_crossings_map);
  //void getSiVertexList( std::multimap<double, std::pair<unsigned int, unsigned int>> &si_sorted_map,
  //			  std::vector<double> &vertex_list,
  //			  std::multimap<unsigned int, std::pair<unsigned int, unsigned int>>  &vertex_map);
  //  void addSiliconClusters(  std::multimap<unsigned int, std::pair<unsigned int, unsigned int>> &vertex_map);
   //  void correctTpcClusterZ( std::map<unsigned int, double> &vertex_crossings_map,
   //			     std::multimap<unsigned int, std::pair<unsigned int, unsigned int>>  &vertex_map );

  // default values, can be replaced from the macro
  double _phi_search_win = 0.01;
  double _eta_search_win = 0.004;
  double _x_search_win = 0.3;
  double _y_search_win = 0.3;
  double _z_search_win = 0.4;
  
  TrackSeedContainer *_svtx_seed_map{nullptr};
  TrackSeedContainer *_track_map{nullptr};
  TrackSeedContainer *_track_map_silicon{nullptr};
  TrackSeed *_tracklet_tpc{nullptr};
  TrackSeed *_tracklet_si{nullptr};
  TrkrClusterContainer *_cluster_map{nullptr};
  ActsGeometry *_tGeometry{nullptr};
  TrkrClusterCrossingAssoc *_cluster_crossing_map{nullptr};

  std::map<unsigned int, double> _z_mismatch_map;

  double _collision_rate = 50e3;  // input rate for phi correction
  double _reference_collision_rate = 50e3;  // reference rate for phi correction
  double _si_vertex_dzmax = 0.25;  // mm
  double crossing_period = 106.0;  // ns

  bool _test_windows = false;
  bool _pp_mode = false;
  bool _use_intt_time = false;

  int _n_iteration = 0;
  std::string _track_map_name = "TpcTrackSeedContainer";
  std::string _silicon_track_map_name = "SiliconTrackSeedContainer";
};

#endif // PHTRUTHSILICONASSOCIATION_H
