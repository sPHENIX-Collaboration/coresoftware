// Tell emacs that this is a C++ source
// -*- C++ -*-.

#ifndef MICROMEGAS_MICROMEGASCLUSTERIZER_H
#define MICROMEGAS_MICROMEGASCLUSTERIZER_H

/*!
 * \file MicromegasClusterizer.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase/ActsSurfaceMaps.h>

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

  private:

  //! detector name
  std::string m_detector;

  ActsTrackingGeometry *m_tGeometry = nullptr;

  Surface getMmSurfaceFromCoords(PHCompositeNode *topNode, 
				 TrkrDefs::hitsetkey hitsetkey, 
				 Acts::Vector3D world);

};

#endif
