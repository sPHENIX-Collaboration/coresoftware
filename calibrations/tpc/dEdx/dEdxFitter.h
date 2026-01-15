#ifndef __DEDXFITTER_H__
#define __DEDXFITTER_H__

#include <fun4all/SubsysReco.h>
#include <string>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/ActsGeometry.h>
#include <g4detectors/PHG4TpcGeomContainer.h>
#include <globalvertex/SvtxVertexMap.h>

#include "GlobaldEdxFitter.h"

//Forward declerations
class PHCompositeNode;
class TFile;

// dEdx fit analysis module
class dEdxFitter: public SubsysReco
{
 public: 
  //Default constructor
  dEdxFitter(const std::string &name="dEdxFitter");

  //Initialization, called for initialization
  int InitRun(PHCompositeNode *);

  //Process Event, called for each event
  int process_event(PHCompositeNode *);

  //End, write and close files
  int End(PHCompositeNode *);

  //Change output filename
  void set_filename(const char* file)
  { if(file) _outfile = file; }

  void set_nmaps_cut(int nmaps)
  { nmaps_cut = nmaps; }

  void set_nintt_cut(int nintt)
  { nintt_cut = nintt; }

  void set_ntpc_cut(int ntpc)
  { ntpc_cut = ntpc; }

  void set_eta_cut(float eta)
  { eta_cut = eta; }

  void set_dcaxy_cut(float dcaxy)
  { dcaxy_cut = dcaxy; }

  void set_ntracks_to_fit(size_t ntrk)
  { ntracks_to_fit = ntrk; }

 private:
  //output filename
  std::string _outfile = "dedx_outfile.root";
  size_t _event;
   
  SvtxTrackMap* _trackmap = nullptr;
  TrkrClusterContainer* _clustermap = nullptr;
  ActsGeometry* _geometry = nullptr;
  PHG4TpcGeomContainer* _tpcgeom = nullptr;
  SvtxVertexMap* _vertexmap = nullptr;
  
  //Get all the nodes
  void GetNodes(PHCompositeNode *);

  void process_tracks(PHCompositeNode *);

  int nmaps_cut = 1;
  int nintt_cut = 1;
  int ntpc_cut = 30;
  float eta_cut = 1.;
  float dcaxy_cut = 0.5;

  size_t ntracks_to_fit = 40000;
  std::vector<double> minima;
  std::unique_ptr<GlobaldEdxFitter> fitter;

  std::tuple<int,int,int> get_nclus(SvtxTrack* track);
  double get_dedx(SvtxTrack* track);
  double get_dcaxy(SvtxTrack* track);

};

#endif //* __DEDXFITTER_H__ *//
