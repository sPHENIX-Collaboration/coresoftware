#include "Eventplaneinfo.h"

std::map<Eventplaneinfo::EPTYPE, unsigned int> DummyEventplaneinfo;

Eventplaneinfo::ConstEpIter Eventplaneinfo::begin_ep_ids() const
{
  return DummyEventplaneinfo.end();
}

Eventplaneinfo::ConstEpIter Eventplaneinfo::find_ep_ids(EPTYPE /*type*/) const
{
  return DummyEventplaneinfo.end();
}

Eventplaneinfo::ConstEpIter Eventplaneinfo::end_ep_ids() const
{
  return DummyEventplaneinfo.end();
}

Eventplaneinfo::EpIter Eventplaneinfo::begin_ep_ids()
{
  return DummyEventplaneinfo.end();
}

Eventplaneinfo::EpIter Eventplaneinfo::find_ep_ids(EPTYPE /*type*/)
{
  return DummyEventplaneinfo.end();
}

Eventplaneinfo::EpIter Eventplaneinfo::end_ep_ids()
{
  return DummyEventplaneinfo.end();
}
