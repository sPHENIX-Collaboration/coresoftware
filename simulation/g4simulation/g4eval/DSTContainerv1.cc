/*!
 * \file DSTContainerv1.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 * \modified Christof Roland <cer@mit.edu>
 */

#include "DSTContainerv1.h"

//_____________________________________________________________________
void DSTContainerv1::Reset()
{
  m_events.clear();
  m_clusters.clear();
  m_tracks.clear();
}
