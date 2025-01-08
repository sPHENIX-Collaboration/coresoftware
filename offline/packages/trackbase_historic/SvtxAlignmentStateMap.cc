#include "SvtxAlignmentStateMap.h"

namespace
{
  SvtxAlignmentStateMap::StateMap map;
}

SvtxAlignmentStateMap::ConstIter SvtxAlignmentStateMap::begin() const
{
  return map.end();
}

SvtxAlignmentStateMap::ConstIter SvtxAlignmentStateMap::find(unsigned int /*unused*/) const
{
  return map.end();
}

SvtxAlignmentStateMap::ConstIter SvtxAlignmentStateMap::end() const
{
  return map.end();
}

SvtxAlignmentStateMap::Iter SvtxAlignmentStateMap::begin()
{
  return map.end();
}

SvtxAlignmentStateMap::Iter SvtxAlignmentStateMap::find(unsigned int /*unused*/)
{
  return map.end();
}

SvtxAlignmentStateMap::Iter SvtxAlignmentStateMap::end()
{
  return map.end();
}
