
#include "MicromegasGeometryContainer.h"

//___________________________________________________________
void MicromegasGeometryContainer::identify( std::ostream& out ) const
{ out << "MicromegasGeometryContainer" << std::endl; }

//___________________________________________________________
TVector3 MicromegasGeometryContainer::get_strip_begin( unsigned int layer, unsigned int tile, unsigned int strip ) const
{ return m_strip_begin.at({.layer=layer, .tile=tile, .strip=strip}); }

//___________________________________________________________
TVector3 MicromegasGeometryContainer::get_strip_end( unsigned int layer, unsigned int tile, unsigned int strip ) const
{ return m_strip_end.at({.layer=layer, .tile=tile, .strip=strip}); }

//___________________________________________________________
void MicromegasGeometryContainer::Reset()
{
  m_strip_begin.clear();
  m_strip_end.clear();
}

//___________________________________________________________
void MicromegasGeometryContainer::add_strip( 
  unsigned int layer, unsigned int tile, unsigned int strip, 
  const TVector3& begin, const TVector3& end )
{
  const StripId strip_id{.layer=layer, .tile=tile, .strip=strip};
  m_strip_begin[strip_id] = begin;
  m_strip_end[strip_id] = end;
}
