// Tell emacs that this is a C++ source
//  -*- C++ -*-.
/*!
 *  \file		PHG4TrackFastSimEval.h
 *  \brief		Evaluation module for PHG4TrackFastSim output
 *  \details	input: PHG4TruthInfoContainer, SvtxTrackMap with SvtxTrack_FastSim inside
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef G4TRACKFASTSIM_PHG4TRACKFASTSIMEVAL_H
#define G4TRACKFASTSIM_PHG4TRACKFASTSIMEVAL_H

#include <fun4all/SubsysReco.h>
#include <string>

//Forward declerations
class PHCompositeNode;
class PHG4TruthInfoContainer;
class SvtxTrackMap;
class SvtxVertexMap;
class TTree;
class TH2D;

//Brief: basic ntuple and histogram creation for sim evaluation
class PHG4TrackFastSimEval : public SubsysReco
{
 public:
  //Default constructor
  PHG4TrackFastSimEval(const std::string& name = "PHG4TrackFastSimEval",
                       const std::string& filename = "g4eval.root",
                       const std::string& trackmapname = "SvtxTrackMap");

  //Initialization, called for initialization
  int Init(PHCompositeNode*);

  //Process Event, called for each event
  int process_event(PHCompositeNode*);

  //End, write and close files
  int End(PHCompositeNode*);

  //Change output filename
  void set_filename(const char* file)
  {
    if (file) _outfile_name = file;
  }

  //Flags of different kinds of outputs
  enum Flag
  {
    //all disabled
    NONE = 0,
  };

  //Set the flag
  //Flags should be set like set_flag(PHG4TrackFastSimEval::TRUTH, true) from macro
  void set_flag(const Flag& flag, const bool& value)
  {
    if (value)
      _flags |= flag;
    else
      _flags &= (~flag);
  }

  //User modules
  void reset_variables();

 private:
  void fill_track_tree(PHCompositeNode*);
  void fill_vertex_tree(PHCompositeNode*);

  //output filename
  std::string _outfile_name;

  //name of SvtxTrackMap collection
  std::string _trackmapname;

  //Event counter
  int _event;

  //Get all the nodes
  int GetNodes(PHCompositeNode*);

  //flags
  unsigned int _flags;

  //TTrees
  TTree* _eval_tree_tracks;
  TTree* _eval_tree_vertex;
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
  float gvt;

  //-- reco
  int trackID;
  int charge;
  int nhits;
  float px;
  float py;
  float pz;
  float pcax;
  float pcay;
  float pcaz;
  float dca2d;

  //vertex
  float vx;
  float vy;
  float vz;
  float deltavx;
  float deltavy;
  float deltavz;
  int ntracks;
  int n_from_truth;

  //Histos
  TH2D* _h2d_Delta_mom_vs_truth_mom;
  TH2D* _h2d_Delta_mom_vs_truth_eta;

  //Node pointers
  PHG4TruthInfoContainer* _truth_container;
  SvtxTrackMap* _trackmap;
  SvtxVertexMap* _vertexmap;
};

#endif  //* __PHG4TrackFastSimEval_H__ *//
