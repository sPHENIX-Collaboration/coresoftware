#ifndef TPC_TPCSPACECHARGECORRECTION_H
#define TPC_TPCSPACECHARGECORRECTION_H

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>
#include <trackbase/TrkrDefs.h>

#include <array>

class TrkrCluster;
class TrkrClusterContainer;
class TFile;
class TH3;

class TpcSpaceChargeCorrection : public SubsysReco
{
  public:

  //! constructor
  TpcSpaceChargeCorrection(  const std::string& = "TpcSpaceChargeCorrection" );

  //! global initialization
  virtual int InitRun(PHCompositeNode*);

  //! event processing
  virtual int process_event(PHCompositeNode*);

  //! distortion filename
  void set_distortion_filename( const std::string& value )
  { m_distortion_filename = value; }

  enum CoordMask
  {
    COORD_PHI = 1<<0,
    COORD_R = 1<<1,
    COORD_Z = 1<<2
  };

  // use corrections on specific set of coordinates
  void set_coordinates( unsigned int value )
  { m_coordinates = value; }

  private:

  //! load nodes
  int load_nodes( PHCompositeNode* );

  //! transform clusters
  void transform_clusters();

  //! transform clusters
  void transform_cluster( TrkrCluster* );

  //! cluster container
  TrkrClusterContainer* m_cluster_map = nullptr;

  //! space charge correction file name
  std::string m_distortion_filename;
  TFile *m_distortion_tfile = nullptr;

  //!@name space charge distortion histograms
  //@{
  std::array<TH3*,2> m_hDRint = {{nullptr, nullptr}};
  std::array<TH3*,2> m_hDPint = {{nullptr, nullptr}};
  std::array<TH3*,2> m_hDZint = {{nullptr, nullptr}};
  //@}

  //* coordinated for which corrections are applied
  unsigned int m_coordinates = COORD_PHI|COORD_R|COORD_Z;

};

#endif
