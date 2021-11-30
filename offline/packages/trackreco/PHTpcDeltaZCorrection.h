// Tell emacs that this is a C++ source
//  -*- C++ -*-.

/*!
 *  \file PHTpcDeltaZCorrection.h
 *  \brief Navigate along a given track and correct TPC cluster z position to account for particle finite velocity in the TPC
 *  \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#ifndef TRACKRECO_PHTPCDELTAZCORRECTION_H
#define TRACKRECO_PHTPCDELTAZCORRECTION_H

#include <fun4all/SubsysReco.h>
#include <phparameter/PHParameterInterface.h>
#include <trackbase_historic/ActsTransformations.h>
#include <trackbase/TrkrDefs.h>

class SvtxTrackMap;
class TrkrClusterContainer;

class PHTpcDeltaZCorrection : public SubsysReco, public PHParameterInterface
{
 public:

  /// constructor
  PHTpcDeltaZCorrection(const std::string &name = "PHTpcDeltaZCorrection");

  /// destructor
  ~PHTpcDeltaZCorrection() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  void SetDefaultParameters() override;

  private:

  /// load nodes
  int load_nodes( PHCompositeNode* );

  /// process tracks
  void process_tracks();

  /// process track
  void process_track( unsigned int, SvtxTrack* );

  /// Acts surface maps for surface lookup
  ActsSurfaceMaps *m_surfmaps = nullptr;

  /// Acts tracking geometry for surface lookup
  ActsTrackingGeometry *m_tGeometry = nullptr;

  /// acts transformation
  ActsTransformations m_transformer;

  /// track map
  SvtxTrackMap *m_track_map = nullptr;

  /// cluster map
  TrkrClusterContainer *m_cluster_map = nullptr;

  /// list of corrected cluster keys
  /** needed to prevent clusters to be corrected twice, when same cluster is used for two different tracks */
  std::set<TrkrDefs::cluskey> m_corrected_clusters;

  /// constant magnetic field
  /** it is used to get helix trajectory from momentum at origin */
  double m_bz_const = 1.4; // Tesla

  /// electron drift velocity in gas
  double m_drift_velocity = NAN;

};

#endif // PHTpcDeltaZCorrection_H
