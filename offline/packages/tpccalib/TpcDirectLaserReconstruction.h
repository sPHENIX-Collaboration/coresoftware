#ifndef TPCCALIB_TPCDIRECTLASERRECONSTRUCTION_H
#define TPCCALIB_TPCDIRECTLASERRECONSTRUCTION_H

/**
 * \file TpcDirectLaserReconstruction.h
 * \brief performs the reconstruction of TPC direct laser tracks
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <fun4all/SubsysReco.h>
#include <phparameter/PHParameterInterface.h>

class SvtxTrack;
class SvtxTrackMap;
class TrkrClusterContainer;
class TrkrHitSetContainer;

class TpcDirectLaserReconstruction: public SubsysReco, public PHParameterInterface
{
  
  public:
  
  /// constructor
  TpcDirectLaserReconstruction( const std::string& = "TpcDirectLaserReconstruction" );
 
  /// global initialization
  int Init(PHCompositeNode*) override;

  /// run initialization
  int InitRun(PHCompositeNode*) override;

  /// event processing
  int process_event(PHCompositeNode*) override;

  /// end of processing
  int End(PHCompositeNode*) override;

  /// parameters
  void SetDefaultParameters() override;
   
  private:

  /// load nodes
  int load_nodes( PHCompositeNode* );

  /// process tracks
  void process_tracks();

  /// process track
  void process_track( SvtxTrack* );
    
  ///@name nodes
  //@{
  TrkrHitSetContainer* m_hitsetcontainer = nullptr;
  SvtxTrackMap* m_track_map = nullptr;
  TrkrClusterContainer* m_cluster_map = nullptr;
  //@}
  
};

#endif
