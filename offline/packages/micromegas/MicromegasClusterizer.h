#ifndef MICROMEGAS_MICROMEGASCLUSTERIZER_H
#define MICROMEGAS_MICROMEGASCLUSTERIZER_H

/*!
 * \file MicromegasClusterizer.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <fun4all/SubsysReco.h>

#include <string>        

class PHCompositeNode;

//! micromegas clusterizer
class MicromegasClusterizer : public SubsysReco
{
 public:

  //! constructor
  MicromegasClusterizer( const std::string &name = "MicromegasClusterizer" );

  //! run initialization
  int InitRun(PHCompositeNode*) override;

  //! event processing
  int process_event(PHCompositeNode*) override;

  //! cluster version
  void set_cluster_version(int value) { m_cluster_version = value; }

  //! read raw data 
  /** not implemented for now */
  void set_read_raw(bool read_raw){ do_read_raw = read_raw;}

  private:

  bool do_read_raw = false;
  int m_cluster_version = 4;
};

#endif
