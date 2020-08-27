#include "ParticleFlowElementContainer.h"

#include "ParticleFlowElement.h"

#include <cstdlib>
#include <iostream>
#include <utility>


ParticleFlowElementContainer::ConstRange
ParticleFlowElementContainer::getParticleFlowElements(void) const
{
  return make_pair( _pflowElementMap.begin(), _pflowElementMap.end() );
}

ParticleFlowElementContainer::Range
ParticleFlowElementContainer::getParticleFlowElements(void)
{
  return make_pair( _pflowElementMap.begin(), _pflowElementMap.end() );
}

void
ParticleFlowElementContainer::AddParticleFlowElement( int index , ParticleFlowElement *pflowElement )
{

  _pflowElementMap[index] = pflowElement;

}

ParticleFlowElement *
ParticleFlowElementContainer::getParticleFlowElement( int index ) 
{
  ConstIterator it = _pflowElementMap.find( index );
  if (it != _pflowElementMap.end())
  {
    return it->second;
  }
  return NULL;
}

const ParticleFlowElement *
ParticleFlowElementContainer::getParticleFlowElement( int index ) const
{
  ConstIterator it = _pflowElementMap.find( index );
  if (it != _pflowElementMap.end())
  {
    return it->second;
  }
  return NULL;
}

int ParticleFlowElementContainer::isValid() const
{
  return (!_pflowElementMap.empty());
}

void ParticleFlowElementContainer::Reset()
{
  while (_pflowElementMap.begin() != _pflowElementMap.end())
  {
    delete _pflowElementMap.begin()->second;
    _pflowElementMap.erase(_pflowElementMap.begin());
  }
}

void ParticleFlowElementContainer::identify(std::ostream &os) const
{
  os << "ParticleFlowElementContainer, number of elements: " << size() << std::endl;
}

