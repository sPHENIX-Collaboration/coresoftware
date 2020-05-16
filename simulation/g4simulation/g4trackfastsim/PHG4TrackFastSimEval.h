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

#include <map>
#include <string>

//Forward declarations
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
  void set_filename(const std::string &file)
  {
    _outfile_name = file;
  }

  //User modules
  void reset_variables();

  void AddProjection(const std::string& name);

 private:
  void fill_track_tree(PHCompositeNode*);
  void fill_vertex_tree(PHCompositeNode*);

  //Get all the nodes
  int GetNodes(PHCompositeNode*);

  //Node pointers
  PHG4TruthInfoContainer* m_TruthInfoContainer;
  SvtxTrackMap* m_TrackMap;
  SvtxVertexMap* m_VertexMap;

  //TTrees
  TTree* m_TracksEvalTree;
  TTree* m_VertexEvalTree;

  //Histos
  TH2D* m_H2D_DeltaMomVsTruthMom;
  TH2D* m_H2D_DeltaMomVsTruthEta;

  //Event counter
  int m_EventCounter;

// TTree variables
  int m_TTree_Event;
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

  static const int nproj = 3;
  // projections hits/mom
  float proj[3][nproj];
  float proj_p[3][nproj];
  // hits/mom at reference
  float ref[3][nproj];
  float ref_p[3][nproj];

  //vertex
  float vx;
  float vy;
  float vz;
  float deltavx;
  float deltavy;
  float deltavz;
  int ntracks;
  int n_from_truth;

  //output filename
  std::string _outfile_name;

  //name of SvtxTrackMap collection
  std::string _trackmapname;

// names and index of projections
  std::map<std::string, int> m_ProjectionNameMap;
};

#endif  //* G4TRACKFASTSIM_PHG4TRACKFASTSIMEVAL_H *//
