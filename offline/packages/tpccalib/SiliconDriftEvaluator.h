#ifndef G4EVAL_SiliconDriftEvaluator_H
#define G4EVAL_SiliconDriftEvaluator_H

/*
 * Bade Sayki June 16th, 2026 -- LANL
 * This module is created to calibrate the drift velocity in the TPC by projecting the silicon seeds and the TPC seeds to the beam axis and calculating the z residuals.
 * This is heavily inspired by and distilled from Dr. Hugo Pereira Da Costa's TrackingEvaluator_hp module. It is meant to be a more lightweight and specialized version.
 */

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>

#include <string>
#include <vector>

class SvtxTrackMap;
class TH2F;
class TH3F;

class SiliconDriftEvaluator : public SubsysReco
{
  public:

  //! constructor
  SiliconDriftEvaluator( const std::string& = "SiliconDriftEvaluator" );

  //! global initialization
  virtual int Init(PHCompositeNode*);

  //! run initialization
  virtual int InitRun(PHCompositeNode*);

  //! event processing
  virtual int process_event(PHCompositeNode*);

  //! end of processing
  virtual int End(PHCompositeNode*);

  // track information stored in the Container
  class TrackStruct
  {
    public:

    using List = std::vector<TrackStruct>;

    //cluster counts
    unsigned int _nclusters_mvtx = 0;
    unsigned int _nclusters_intt = 0;
    unsigned int _nclusters_tpc  = 0;


    //tpc seed kinematics
    float _pt  = 0;
    float _eta = 0;
    float _phi = 0;

    //seed z positions

    // z position of the TPC seed at the beamline
    float _z_tpc = 0;

    //! z position of the silicon seed at the beamline
    float _z_si  = 0;


    //! beam-bunch crossing number
    short int _crossing = 0;

    //! crossing-corrected dz = z_tpc_corr - z_si (cm)
    float _dz = 0;

  };

  //! track container stored on the node tree
  class Container : public PHObject
  {
    public:

    //! constructor
    explicit Container() = default;

    //! copy constructor
    explicit Container( const Container& ) = delete;

    //! assignment operator
    Container& operator=( const Container& ) = delete;

    //! reset
    void Reset() override
    { _tracks.clear(); }

    //!@name accessors
    //@{

    const TrackStruct::List& tracks() const
    { return _tracks; }


    // modifiers

    void addTrack( const TrackStruct& track )
    { _tracks.push_back( track ); }

    void clearTracks()
    { _tracks.clear(); }


    private:

    //! tracks array
    TrackStruct::List _tracks;

    ClassDefOverride(Container, 1)

  };

  //! track map name
  void set_trackmapname( const std::string& value )
  { m_trackmapname = value; }

  // initial drift velocity (cm/ns); used as starting point for the fit and for the crossing correction
  void set_drift_velocity( double value )
  { m_drift_velocity = value; }

  // bunch-crossing interval in ns (default: 106.65237 ns)
  void set_crossing_interval( double value )
  { m_crossing_interval = value; }

  // minimum pT cut on tracks (GeV)
  void set_min_pt( double value )
  { m_min_pt = value; }

  // minimum number of TPC clusters required
  void set_min_nclusters_tpc( unsigned int value )
  { m_min_nclusters_tpc = value; }

  // minimum number of MVTX clusters required
  void set_min_nclusters_mvtx( unsigned int value )
  { m_min_nclusters_mvtx = value; }

  // minimum number of INTT clusters required
  void set_min_nclusters_intt( unsigned int value )
  { m_min_nclusters_intt = value; }

  // maximum abs(eta) of TPC seed accepted
  void set_max_eta( double value )
  { m_max_eta = value; }

  // half-range of the z_si histogram axis (cm)
  void set_max_z( double value )
  { m_max_z = value; }

  // half-range of the dz histogram axis (cm)
  void set_max_dz( double value )
  { m_max_dz = value; }

  //! minimum entries per z slice required by FitSlicesY
  void set_min_slice_entries( int value )
  { m_min_slice_entries = value; }

  // output drift plot filename (.png)
  void set_plot_filename( const std::string& value )
  { m_plot_filename = value; }

  private:

  //! load nodes
  int load_nodes( PHCompositeNode* );

  //! evaluate tracks
  void evaluate_tracks();

  //! evaluation node
  Container* m_container = nullptr;

  //! track map
  SvtxTrackMap* m_track_map = nullptr;

  //! 3D accumulator histogram: x = eta bin [2], y = z_si, z = dz
  TH3F* m_hist3D = nullptr;

  //! track map name
  std::string m_trackmapname = "SvtxTrackMap";

  //! initial drift velocity (cm/ns)
  double m_drift_velocity = 0.00749;

  //! bunch-crossing interval (ns)
  double m_crossing_interval = 106.65237;

  // track selection cuts
  
  double       m_min_pt            = 0.5;
  unsigned int m_min_nclusters_tpc  = 20;
  unsigned int m_min_nclusters_mvtx = 3;
  unsigned int m_min_nclusters_intt = 2;
  double       m_max_eta            = 0.9;


  //! histogram range
  double m_max_z  = 20.0;
  double m_max_dz = 10.0;

  //! minimum entries per z slice for FitSlicesY
  int m_min_slice_entries = 10;

  //! output QA plot filename
  std::string m_plot_filename = "silicon_drift_calib_QA.png";

};

#endif  // G4EVAL_SiliconDriftEvaluator_H
