#ifndef G4EVAL_MicromegasClusterEvaluator_hp_H
#define G4EVAL_MicromegasClusterEvaluator_hp_H

#include <fun4all/SubsysReco.h>
#include <gsl/gsl_rng.h>
#include <micromegas/MicromegasCalibrationData.h>
#include <micromegas/MicromegasDefs.h>
#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>

#include <TVector3.h>

#include <map>
#include <memory>

class ActsGeometry;
class PHG4CylinderGeomContainer;
class TrkrClusterHitAssoc;
class TrkrClusterContainer;
class TrkrHitSetContainer;

class MicromegasClusterEvaluator_hp : public SubsysReco
{
  public:

  /// constructor
  MicromegasClusterEvaluator_hp( const std::string& = "MicromegasClusterEvaluator_hp" );

  /// global initialization
  virtual int Init(PHCompositeNode*);

  /// run initialization
  virtual int InitRun(PHCompositeNode*);

  /// event processing
  virtual int process_event(PHCompositeNode*);

  /// end of processing
  virtual int End(PHCompositeNode*);

  /// set default pedestal
  void set_default_pedestal( double value )
  { m_default_pedestal = value; }

  /// set whether default pedestal is used or not
  void set_use_default_pedestal( bool value )
  { m_use_default_pedestal = value; }

  /// calibration file
  void set_calibration_file( const std::string& value )
  { m_calibration_filename = value; }

  class Cluster
  {
    public:

    unsigned short layer = 0;

    unsigned short tile = 0;

    unsigned short size = 0;

    double charge = 0;

    int strip = 0;

    int region = 0;

    double x_local = 0;
    double y_local = 0;

    double x = 0;
    double y = 0;
    double z = 0;

    bool is_signal = true;

    using List = std::vector<Cluster>;
  };

  class Container: public PHObject
  {

    public:

    void Reset();

    // number of signal waveforms
    unsigned short n_waveforms_signal = 0;

    // number of clusters
    unsigned short n_clusters = 0;

    // clusters
    Cluster::List clusters;

    // number of clusters per detector
    std::vector<unsigned short> n_detector_clusters;

    // number of clusters per region
    std::vector<unsigned short> n_region_clusters;

    // minimum cluster charge per detector
    std::vector<unsigned short> min_cluster_size;

    // minimum cluster charge per detector
    std::vector<double> min_cluster_charge;

    // strip of the first cludster, per detector
    std::vector<int> first_cluster_strip;

    // number of clusters
    unsigned short n_phi_clusters = 0;

    // number of clusters
    unsigned short n_z_clusters = 0;

    ClassDef(Container,1)

  };

  private:

  /// load nodes
  int load_nodes( PHCompositeNode* );

  /// evaluate hits
  void evaluate_clusters();

  /// cluster array
  Container* m_container = nullptr;

  /// Acts tracking geometry for surface lookup
  ActsGeometry *m_tGeometry = nullptr;

  /// gemometry
  PHG4CylinderGeomContainer* m_geonode = nullptr;

  //! hits
  TrkrHitSetContainer* m_hitsetcontainer = nullptr;

  //! clusters
  TrkrClusterContainer* m_cluster_map = nullptr;

  //! cluster to hit association
  TrkrClusterHitAssoc* m_cluster_hit_map = nullptr;

  //!@name calibration filename
  //@{

  /// if true, use default pedestal to get hit charge. Relies on calibration data otherwise
  bool m_use_default_pedestal = true;

  /// default pedestal
  double m_default_pedestal = 74.6;

  /// calibration filename
  std::string m_calibration_filename;

  /// calibration data
  MicromegasCalibrationData m_calibration_data;

  //@}

};

#endif  // G4EVAL_MicromegasClusterEvaluator_hp_H
