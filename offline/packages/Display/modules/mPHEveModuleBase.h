/*
        \file mPHEveModuleBase.h
        \author Sookhyun Lee
        \brief abstract interface for PHEve framework modules
        \version $Revision: 1.1 $
        \date    $Date: 07/26/2016
*/

#ifndef __MPHEVEMODULEBASE_H__
#define __MPHEVEMODULEBASE_H__

#include <boost/bind.hpp>

#include <phool/phool.h>
#include <map>

#include <TEveManager.h>

using boost::bind;

class PHCompositeNode;
class TEveElement;

class mPHEveModuleBase 
{
 public:
  virtual ~mPHEveModuleBase() {}
  virtual void init(PHCompositeNode* top_node) = 0;
  virtual void init_run(PHCompositeNode* top_node) = 0;
  virtual bool event(PHCompositeNode* top_node) = 0;
  virtual void draw_event(){};
  void add_elements()
  {
    std::for_each(_elmt_buffer.begin(),
		  _elmt_buffer.end(),
		  bind(&TEveManager::AddElement,
		       _evemanager,
		       bind(&mvt::first, _1),
		       bind(&mvt::second, _1)));    
  };
  void buffer_element(TEveElement* el, 
		     TEveElement* parent_elmt)
  {
    _elmt_buffer.insert(std::make_pair(el, parent_elmt));
  }

  void clear_element_buffer()
  {
    _elmt_buffer.clear();
  }

  virtual void clear() = 0;
  
  virtual void end(PHCompositeNode* top_node) = 0;
  
 protected:

  TEveManager* _evemanager;

  typedef std::multimap<TEveElement*, TEveElement*> elmtmap;
  elmtmap _elmt_buffer;
  typedef elmtmap::iterator::value_type mvt;

};

#endif // __MPHEVEMODULEBASE_H__
