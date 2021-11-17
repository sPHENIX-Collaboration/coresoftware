#ifndef PARTICLEFLOW_PARTICLEFLOWELEMENTCONTAINER_H
#define PARTICLEFLOW_PARTICLEFLOWELEMENTCONTAINER_H

//===========================================================
/// \file ParticleFlowElementContainer.h
/// \brief Simple container for particle flow elements
/// \author Dennis V. Perepelitsa
//===========================================================

#include <phool/PHObject.h>

#include <iostream>
#include <map>

class ParticleFlowElement;

class ParticleFlowElementContainer : public PHObject
{
 public:
  typedef std::map<int, ParticleFlowElement *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;
  
  ParticleFlowElementContainer( )
    {
      
    }
  
  ~ParticleFlowElementContainer() override {}

  void Reset() override;
  int isValid() const override;
  void identify(std::ostream &os = std::cout) const override;

  void AddParticleFlowElement(int index, ParticleFlowElement *pflowElement);
  ParticleFlowElement *getParticleFlowElement(int index);
  const ParticleFlowElement *getParticleFlowElement(int index) const;

  //! return all elements
  ConstRange getParticleFlowElements(void) const;
  Range getParticleFlowElements(void);

  unsigned int size() const { return _pflowElementMap.size(); }

 protected:
  Map _pflowElementMap;

  ClassDefOverride(ParticleFlowElementContainer, 1)
};

#endif
