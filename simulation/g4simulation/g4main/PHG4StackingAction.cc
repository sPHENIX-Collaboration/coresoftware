#include "PHG4StackingAction.h"

#include <string>

PHG4StackingAction::PHG4StackingAction(const std::string& name, const int i)
  : m_Verbosity(i)
  , m_Name(name)
{
}
