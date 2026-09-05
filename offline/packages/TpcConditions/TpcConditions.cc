#include "TpcConditions.h"

TpcConditions::TpcConditions()
  : m_LoadCurrent(0.0)
  , m_LoadNorth  (0.0)
  , m_LoadSouth  (0.0)
  , m_LoadSR1    (0.0)
  , m_LoadSR2    (0.0)
  , m_LoadSR3    (0.0)
  , m_LoadNR1    (0.0)
  , m_LoadNR2    (0.0)
  , m_LoadNR3    (0.0)
  , m_Temperature(0.0)
  , m_Pressure(0.0)
  , m_FieldOK(false)
  , m_GainOK(false)
{
}
