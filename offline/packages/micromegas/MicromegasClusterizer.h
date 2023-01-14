// Tell emacs that this is a C++ source
// -*- C++ -*-.

#ifndef MICROMEGAS_MICROMEGASCLUSTERIZER_H
#define MICROMEGAS_MICROMEGASCLUSTERIZER_H

/*!
 * \file MicromegasClusterizer.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <fun4all/SubsysReco.h>

#include <string>                // for string

class PHCompositeNode;

//! micromegas clusterizer
class MicromegasClusterizer : public SubsysReco
{
 public:

  //! constructor
  MicromegasClusterizer(
    const std::string &name = "MicromegasClusterizer",
    const std::string &detector = "MICROMEGAS"
    );

  //! run initialization
  int InitRun(PHCompositeNode*) override;

  //! event processing
  int process_event(PHCompositeNode*) override;

  void set_cluster_version(int value) { m_cluster_version = value; }
  void set_read_raw(bool read_raw){ do_read_raw = read_raw;}

  private:

  //! detector name
  std::string m_detector;
  bool do_read_raw = false;
  int m_cluster_version = 4;
};

#endif
