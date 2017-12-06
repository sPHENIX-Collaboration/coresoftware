#ifndef PHNODE_H__
#define PHNODE_H__

//  Declaration of class PHNode
//  Purpose: abstract base class for all node classes

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
  PHNode *getParent() const { return parent; }
  bool isPersistent() const { return persistent; }
  void makePersistent() { persistent = true; }
  const std::string getObjectType() const { return objecttype; }
  const std::string getType() const { return type; }
  const std::string getName() const { return name; }
  const std::string getClass() const { return objectclass; }
  void setParent(PHNode *p) { parent = p; }
  void setName(const std::string &n) { name = n; }
  void setObjectType(const std::string &type) { objecttype = type; }
  virtual void prune() = 0;
  virtual void print(const std::string &) = 0;
  virtual void forgetMe(PHNode *) = 0;
  virtual bool write(PHIOManager *, const std::string & = "") = 0;

  virtual void setResetFlag(const bool b) { reset_able = b; }
  virtual bool getResetFlag() const { return reset_able; }
  void makeTransient() { persistent = false; }
 protected:
  PHNode();
  PHNode(const PHNode &);  // implement invalid copy ctor
  PHNode &operator=(const PHNode &);

  PHNode *parent;
  bool persistent;
  std::string type;
  std::string objecttype;
  std::string name;
  std::string objectclass;
  bool reset_able;
};

std::ostream &operator<<(std::ostream &, const PHNode &);

#endif /* __PHNODE_H__ */
