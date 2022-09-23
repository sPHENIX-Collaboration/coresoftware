// Tell emacs that this is a C++ source
//  -*- C++ -*-.

/*!
 *  \file		  TpcClusterCleaner.h
 *  \brief		 Class for removing bad TPC clusters
 *  \author	 Tony Frawley <afrawley@fsu.edu>
 */

#ifndef TPCCLUSTERCLEANER_H
#define TPCCLUSTERCLEANER_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;
class TrkrCluster;
class TrkrClusterContainer;

class TpcClusterCleaner : public SubsysReco
{
 public:

  TpcClusterCleaner(const std::string &name = "TpcClusterCleaner");

  ~TpcClusterCleaner() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void set_rphi_error_low_cut(double cut){_rphi_error_low_cut = cut;}
  void set_rphi_error_high_cut(double cut){_rphi_error_high_cut = cut;}

  void set_new_rphi_error(double err){_new_rphi_error = err;}
  void set_new_z_error(double err){_new_z_error = err;}
  void set_cluster_version(int value) { m_cluster_version = value; }

 private:

  void rotate_error(double erphi, double ez, double phi, double error[][3]);

  int GetNodes(PHCompositeNode* topNode);
  TrkrClusterContainer *_cluster_map = nullptr;

  double _rphi_error_low_cut = 0.01;
  double _rphi_error_high_cut = 0.1;  // made large enough to not matter for now

  double _new_rphi_error = 0.05;
  double _new_z_error = 0.1;
  int m_cluster_version = 4;
};

#endif // TPCCLUSTERCLEANER_H
