#ifndef TPC_TpcSpaceChargeCorrection_hp_H
#define G4EVAL_TRACKINGEVALUATOR_HP_H

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>
#include <trackbase/TrkrDefs.h>

class TrkrCluster;
class TrkrClusterContainer;
class TFile;
class TH3;

class TpcSpaceChargeCorrection_hp : public SubsysReco
{
  public:

  //! constructor
  TpcSpaceChargeCorrection_hp(  const std::string& = "TpcSpaceChargeCorrection_hp" );

  //! global initialization
  virtual int InitRun(PHCompositeNode*);

  //! event processing
  virtual int process_event(PHCompositeNode*);

  //! distortion filename
  void set_distortion_filename( const std::string& value )
  { m_distortion_filename = value; }

  private:

  //! load nodes
  int load_nodes( PHCompositeNode* );

  //! transform clusters
  void transform_clusters();

  //! transform clusters
  void transform_cluster( TrkrCluster* );

  //! cluster container
  TrkrClusterContainer* _cluster_map = nullptr;

  //! space charge correction file name
  std::string m_distortion_filename;
  TFile *m_distortion_tfile = nullptr;

  //!@name space charge distortion histograms
  //@{
  TH3 *hDRint = nullptr;
  TH3 *hDPint = nullptr;
  TH3 *hDZint = nullptr;
  //@}

  /*! \brief
   true if the maps contain the full z range
   assume it only contains positive z otherwise
  */
  bool m_fullzrange = true;

};

#endif
