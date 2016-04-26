#ifndef __PHG4TrackKalmanFitter_H__
#define __PHG4TrackKalmanFitter_H__

#include <fun4all/SubsysReco.h>
#include <string>

//Forward declerations
class PHCompositeNode;
class PHG4TruthInfoContainer;
class SvtxClusterMap;
class SvtxEvalStack;
class TFile;
class TTree;


//Brief: basic ntuple and histogram creation for sim evaluation
class PHG4TrackKalmanFitter: public SubsysReco
{
 public: 
  //Default constructor
  PHG4TrackKalmanFitter(const std::string &name="PHG4TrackKalmanFitter");

  //Initialization, called for initialization
  int Init(PHCompositeNode *);

  //Process Event, called for each event
  int process_event(PHCompositeNode *);

  //End, write and close files
  int End(PHCompositeNode *);

  //Change output filename
  void set_filename(const char* file)
  { if(file) _outfile = file; }

  //Flags of different kinds of outputs
  enum Flag
  {
    //all disabled
    NONE = 0,
  };

  //Set the flag
  //Flags should be set like set_flag(PHG4TrackKalmanFitter::TRUTH, true) from macro
  void set_flag(const Flag& flag, const bool& value)
  {
   if(value) _flags |= flag;
   else _flags &= (~flag);
  }

  //User modules
  void fill_tree(PHCompositeNode*);
  void reset_variables();

 private:
  //output filename
  std::string _outfile;
   
  //Event counter
  int _event;

  //Get all the nodes
  void GetNodes(PHCompositeNode *);
  
  //flags
  unsigned int _flags;

  //TTrees
  TTree* _tracks;
  int event;
  //-- truth
  int gtrackID;
  int gflavor;
  float gpx;
  float gpy;
  float gpz;
  float gvx;
  float gvy;
  float gvz;
  //-- reco
  int trackID;
  int charge;
  int nhits;
  float px;
  float py;
  float pz;
  float dca2d;
  //-- clusters
  int clusterID[7];
  int layer[7];
  float x[7];
  float y[7];
  float z[7];
  float size_dphi[7];
  float size_dz[7];

  //Node pointers
  PHG4TruthInfoContainer* _truth_container;
  SvtxClusterMap* _clustermap;

  // eval stack
  SvtxEvalStack* _svtxevalstack;

};

#endif //* __PHG4TrackKalmanFitter_H__ *//
