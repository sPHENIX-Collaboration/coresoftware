/*!
 * \file DSTContainerTcl.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 * \modified Christof Roland <cer@mit.edu>
 */

#include "DSTContainerTcl.h"

//_____________________________________________________________________
void DSTContainerTcl::Reset()
{
  m_events.clear();
  m_clusters.clear();
  m_tracks.clear();
}
