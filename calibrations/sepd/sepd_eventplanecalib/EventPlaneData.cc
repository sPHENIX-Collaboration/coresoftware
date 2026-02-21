#include "EventPlaneData.h"

#include <cmath>

EventPlaneData::EventPlaneData()
{
  sepd_charge.fill(0);
}

void EventPlaneData::Reset()
{
  event_id = 0;
  event_zvertex = std::numeric_limits<double>::quiet_NaN();
  event_centrality = std::numeric_limits<double>::quiet_NaN();
  sepd_totalcharge = std::numeric_limits<double>::quiet_NaN();
  sepd_charge.fill(0);
}

void EventPlaneData::identify(std::ostream& os) const
{
  os << "--- EventPlaneData Identify ---" << std::endl;
  os << "Event ID: " << event_id << std::endl;
  os << "Z-Vertex: " << event_zvertex << std::endl;
  os << "Centrality: " << event_centrality << std::endl;
  os << "sEPD Total Charge: " << sepd_totalcharge << std::endl;
  os << "-------------------------------" << std::endl;
}

int EventPlaneData::isValid() const
{
  // An object is considered invalid if the Z-vertex is still NaN
  // or if the event ID hasn't been set (remains 0).
  if (std::isnan(event_zvertex))
  {
    return 0;
  }

  if (event_id == 0)
  {
    return 0;
  }

  return 1;
}
