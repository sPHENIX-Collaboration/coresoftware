#ifndef __PHNODE_H__
#define __PHNODE_H__

//  Declaration of class PHNode
//  Purpose: abstract base class for all node classes

#include "phool.h"
#include "PHString.h"

#include <iosfwd>

class PHIOManager;

class PHNode 
{ 
public: 

  // Note that the constructor makes a node transient by default.
  PHNode(const PHString&); 
  PHNode(const PHString&, const PHString&);
  virtual ~PHNode(); 

public:

  PHNode* getParent() const { return parent; }
  
  PHBoolean isPersistent() const { return persistent; }
  void makePersistent() { persistent = True;}
  
  const PHString& getObjectType() const { return objecttype; }
  const PHString& getType() const { return type; }
  const PHString& getName() const { return name; }
  void setName(const PHString& n);
 
  void setParent(PHNode *p) { parent = p; }
  void setObjectType(const PHString& type) {objecttype = type;} 
  virtual void prune() = 0;
  virtual void print(const PHString&) = 0;
  virtual void forgetMe(PHNode*) = 0;
  virtual PHBoolean write(PHIOManager *, const PHString& = "") = 0;

  virtual void setResetFlag(const int val);
  virtual PHBoolean getResetFlag() const;
  void makeTransient()  { persistent = False;}
  
protected:
  
  PHNode();
  PHNode(const PHNode&); // implement invalid copy ctor
  PHNode & operator=(const PHNode&);
  
protected:

  PHNode*   parent;
  PHBoolean persistent;
  PHString  type;
  PHString  objecttype;
  PHString  name;
  PHBoolean reset_able;
};

std::ostream & operator << (std::ostream &, const PHNode &);

#endif /* __PHNODE_H__ */
