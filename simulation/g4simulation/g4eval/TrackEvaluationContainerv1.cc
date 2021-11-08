/*!
 * \file TrackEvaluationContainerv1.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TrackEvaluationContainerv1.h"

//_____________________________________________________________________
void TrackEvaluationContainerv1::Reset()
{
  m_events.clear();
  m_clusters.clear();
  m_tracks.clear();
}
