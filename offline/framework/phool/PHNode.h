#ifndef __PHNODE_H__
#define __PHNODE_H__

//  Declaration of class PHNode
//  Purpose: abstract base class for all node classes

#include "phool.h"

#include <iosfwd>
#include <string>

class PHIOManager;

class PHNode 
{ 
public: 

  // Note that the constructor makes a node transient by default.
  PHNode(const std::string &); 
  PHNode(const std::string &, const std::string &);
  virtual ~PHNode(); 

public:

  PHNode* getParent() const { return parent; }
  
  PHBoolean isPersistent() const { return persistent; }
  void makePersistent() { persistent = True;}
  
  const std::string getObjectType() const { return objecttype; }
  const std::string getType() const { return type; }
  const std::string getName() const { return name; }

  void setParent(PHNode *p) { parent = p; }
  void setName(const std::string &n) {name = n;}
  void setObjectType(const std::string &type) {objecttype = type;} 
  virtual void prune() = 0;
  virtual void print(const std::string &) = 0;
  virtual void forgetMe(PHNode*) = 0;
  virtual PHBoolean write(PHIOManager *, const std::string& = "") = 0;

  virtual void setResetFlag(const int val);
  virtual PHBoolean getResetFlag() const;
  void makeTransient()  { persistent = False;}
  
protected:
  
  PHNode();
  PHNode(const PHNode&); // implement invalid copy ctor
  PHNode & operator=(const PHNode&);
  
  PHNode*   parent;
  bool persistent;
  std::string  type;
  std::string  objecttype;
  std::string  name;
  bool reset_able;
};

std::ostream & operator << (std::ostream &, const PHNode &);

#endif /* __PHNODE_H__ */
