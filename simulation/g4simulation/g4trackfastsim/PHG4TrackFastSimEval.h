/*!
 *  \file		PHG4TrackFastSimEval.h
 *  \brief		Evaluation module for PHG4TrackFastSim output
 *  \details	input: PHG4TruthInfoContainer, SvtxTrackMap with SvtxTrack_FastSim inside
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef __PHG4TrackFastSimEval_H__
#define __PHG4TrackFastSimEval_H__

#include <fun4all/SubsysReco.h>
#include <string>

//Forward declerations
class PHCompositeNode;
class PHG4TruthInfoContainer;
class SvtxClusterMap;
class SvtxTrackMap;
class TFile;
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
  void fill_tree(PHCompositeNode*);
  void reset_variables();

 private:
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

  //Histos
  TH2D* _h2d_Delta_mom_vs_truth_mom;
  TH2D* _h2d_Delta_mom_vs_truth_eta;

  //Node pointers
  PHG4TruthInfoContainer* _truth_container;
  SvtxTrackMap* _trackmap;
};

#endif  //* __PHG4TrackFastSimEval_H__ *//
