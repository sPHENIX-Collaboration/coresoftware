#ifndef __SUBSYSRECOSTACK_H__
#define __SUBSYSRECOSTACK_H__

#include <iostream>
#include <string>
#include <list>

#include <SubsysReco.h>



class PHCompositeNode;


/**
 * a SubsysReco that can stack other modules.  use this module to group
 * together a set of SubsysReco modules.  an alternate topnode to be processed 
 * can be specified to the constructor (for example when you want to process 
 * an alternate topnode during embedding).
 *
 * modules put in this container are adopted by the container, that is,
 * they will be destructed when the container gets destructed.
 *
 */
class SubsysRecoStack: public SubsysReco, public std::list< SubsysReco * > {
public:
  SubsysRecoStack(const char* name = "NONAME", PHCompositeNode * chroot = 0);
  ~SubsysRecoStack();

  // ROOT still does not like list::push_back(), so here is a wrapper
  void x_push_back(SubsysReco * r){ push_back( r ); }

  virtual int Init(PHCompositeNode *topNode);
  virtual int InitRun(PHCompositeNode *topNode);
  virtual int process_event(PHCompositeNode *topNode);  
  virtual int Reset(PHCompositeNode *topNode);
  virtual int ResetEvent(PHCompositeNode *topNode);  
  virtual int End(PHCompositeNode *topNode);  
  virtual int EndRun(const int runnumber);


  void Print(const std::string &what = "ALL") const;
  void Print2(const std::string &what = "ALL", std::string prefix = "") const;
  void Verbosity(const int ival);


protected:
  PHCompositeNode * getroot(PHCompositeNode * root){ return chroot ? chroot : root; }

  PHCompositeNode * chroot;
};


#endif
