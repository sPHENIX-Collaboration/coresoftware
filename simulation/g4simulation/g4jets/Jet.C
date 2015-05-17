#include "Jet.h"

#include <cmath>
#include <cstdlib>

using namespace std;

ClassImp(Jet);

Jet::typ_comp_ids Jet::_dummy_ids;

Jet::Jet()
{
}

void
Jet::identify(ostream& os) const
{
  os << "---Jet-----------------------" << endl;
  return;
}

void
Jet::Reset()
{
}

int
Jet::isValid() const
{
  return 0;
}

unsigned int
Jet::get_id() const
{
  return 0xFFFFFFFF;
}
void
Jet::set_id(unsigned int id)
{
}

float
Jet::get_px() const
{
  return NAN;
}
void
Jet::set_px(float px)
{
}

float
Jet::get_py() const
{
  return NAN;
}
void
Jet::set_py(float py)
{
}

float
Jet::get_pz() const
{
  return NAN;
}
void
Jet::set_pz(float pz)
{
}

float
Jet::get_e() const
{
  return NAN;
}
void
Jet::set_e(float e)
{
}

Jet::ConstIter
Jet::begin_comp() const
{

  return _dummy_ids.begin();
}

Jet::ConstIter
Jet::lower_bound_comp(SRC source) const
{
  return _dummy_ids.lower_bound(source);
}
Jet::ConstIter
Jet::upper_bound_comp(SRC source) const
{
  return _dummy_ids.upper_bound(source);
}

Jet::ConstIter
Jet::find(SRC source) const
{
  return _dummy_ids.find(source);
}

Jet::ConstIter
Jet::end_comp() const
{
  return _dummy_ids.end();
}

Jet::Iter
Jet::begin_comp()
{
  cout << __PRETTY_FUNCTION__ << " is not implemented!" << endl;
  exit(7391);
  return _dummy_ids.begin();
}

Jet::Iter
Jet::lower_bound_comp(SRC source)
{
  cout << __PRETTY_FUNCTION__ << " is not implemented!" << endl;
  exit(7391);
  return _dummy_ids.lower_bound(source);
}

Jet::Iter
Jet::upper_bound_comp(SRC source)
{
  cout << __PRETTY_FUNCTION__ << " is not implemented!" << endl;
  exit(7391);
  return _dummy_ids.upper_bound(source);
}

Jet::Iter
Jet::find(SRC source)
{
  cout << __PRETTY_FUNCTION__ << " is not implemented!" << endl;
  exit(7391);
  return _dummy_ids.find(source);
}

Jet::Iter
Jet::end_comp()
{
  cout << __PRETTY_FUNCTION__ << " is not implemented!" << endl;
  exit(7391);
  return _dummy_ids.end();
}


float
Jet::
get_property(PROPERTY prop_id) const
{
  return NAN;
}

void
Jet::
set_property(PROPERTY prop_id, float value)
{
  return;
}
