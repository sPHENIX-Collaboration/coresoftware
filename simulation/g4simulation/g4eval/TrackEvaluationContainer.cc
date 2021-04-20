/*!
 * \file PHG4MicromegasDigitizer.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TrackEvaluationContainer.h"

//_____________________________________________________________________
void TrackEvaluationContainer::Reset()
{
  m_events.clear();
  m_clusters.clear();
  m_tracks.clear();
}
