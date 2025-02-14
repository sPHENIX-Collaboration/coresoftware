// Tell emacs that this is a C++ source
//  -*- C++ -*-.

/*!
 *  \file TPCClusterTzeroCorrection.h
 *  \brief Correct TPC cluster time for TPC ADC time zero
 *  \author Tony Frawley <afrawley@fsu.edu>
 */

#ifndef TRACKRECO_TPCCLUSTERTZEROCORRECTION_H
#define TRACKRECO_TPCCLUSTERTZEROCORRECTION_H

#include <fun4all/SubsysReco.h>
#include <phparameter/PHParameterInterface.h>
#include <trackbase/TrkrDefs.h>

class TrkrClusterContainer;

class TpcClusterTzeroCorrection : public SubsysReco, public PHParameterInterface
{
 public:

  /// constructor
  TpcClusterTzeroCorrection(const std::string &name = "TpcClusterTzeroCorrection");

  /// destructor
  ~TpcClusterTzeroCorrection() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  void SetDefaultParameters() override;
  void setTrkrClusterContainerName(std::string &name){ m_clusterContainerName = name; }

  void setTpcTzeroCorrection(const float tz) {m_tzero = tz  ;}

  private:

  /// load nodes
  int load_nodes( PHCompositeNode* );

  /// process clusters
  void process_clusters();

  /// cluster map
  TrkrClusterContainer *m_cluster_map = nullptr;

  //cluster container name
  std::string m_clusterContainerName = "TRKR_CLUSTER";

  float m_tzero = 0.0;
};

#endif // TpcClusterTzeroCorrection_H
