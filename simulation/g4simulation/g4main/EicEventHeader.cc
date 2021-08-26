#include "EicEventHeader.h"

#include <TSystem.h>

#include <cassert>
#include <cstdlib>

using namespace std;

void EicEventHeader::CopyFrom(const PHObject *phobj)
{
  const EicEventHeader *evthead = dynamic_cast<const EicEventHeader *>(phobj);
  assert(evthead);
  // This is a generic copy of ALL properties an eic event header has
  // do not add explicit copies, they will be added to
  // the new eic event header with their default value increasing memory use
  for (unsigned char ic = 0; ic < UCHAR_MAX; ic++)
  {
    PROPERTY prop_id = static_cast<EicEventHeader::PROPERTY>(ic);
    if (evthead->has_property(prop_id))
    {
      set_property_nocheck(prop_id, evthead->get_property_nocheck(prop_id));
    }
  }
}

EicEventHeader::~EicEventHeader()
{
  return;
}

void EicEventHeader::Reset()
{
  cout << "Reset not implemented by daughter class" << endl;
  return;
}

std::pair<const std::string, EicEventHeader::PROPERTY_TYPE>
EicEventHeader::get_property_info(const PROPERTY prop_id)
{
  switch (prop_id)
  {
  case prop_eventgen:
    return make_pair("Event Generator", EicEventHeader::type_int);
  case prop_milou_weight:
    return make_pair("Milou weight", EicEventHeader::type_float);
  case prop_milou_truex:
    return make_pair("Milou True X", EicEventHeader::type_float);
  case prop_milou_trueq2:
    return make_pair("Milou True Q2", EicEventHeader::type_float);
  case prop_demp_weight:
    return make_pair("DEMP weight", EicEventHeader::type_float);

  default:
    cout << "EicEventHeader::get_property_info - Fatal Error - unknown index " << prop_id << endl;
    gSystem->Exit(1);
    exit(1);
  }
}

bool EicEventHeader::check_property(const PROPERTY prop_id, const PROPERTY_TYPE prop_type)
{
  pair<const string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
  if (property_info.second != prop_type)
  {
    return false;
  }
  return true;
}

string
EicEventHeader::get_property_type(const PROPERTY_TYPE prop_type)
{
  switch (prop_type)
  {
  case type_int:
    return "int";
  case type_uint:
    return "unsigned int";
  case type_float:
    return "float";
  default:
    return "unkown";
  }
}

void EicEventHeader::identify(ostream &os) const
{
  os << "Class " << this->ClassName() << endl;
  os << "Event Generator: ";
  switch (get_eventgenerator_type())
  {
  case EvtGen::Milou:
    os << "Milou";
    break;
  case EvtGen::DEMP:
    os << "DEMP";
    break;
  default:
    os << "Unknown";
    break;
  }
  os << endl;
}
