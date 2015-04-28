#ifndef __FUN4ALLBASE_H__
#define __FUN4ALLBASE_H__

#include <string>

/** Base class for all Fun4All Classes
 *
 *  It implements the Name, the Verbosity and the print method  
 */

class Fun4AllBase
{
 public:


  /** dtor. 
      Does nothing as this is a base class only.
  */
  virtual ~Fun4AllBase();

  /// Returns the name of this module.
  virtual const char *Name() const {return ThisName.c_str();}

  /// Sets the name of this module.
  virtual void Name(const char *name) {ThisName = name;}

  /** Print out some info about this module. 
      @param what can be used to specify what to print exactly.
  */
  virtual void Print(const std::string &what = "ALL") const;

  /// Sets the verbosity of this module (0 by default=quiet).
  virtual void Verbosity(const int ival) {verbosity = ival;}

  /// Gets the verbosity of this module.
  virtual int Verbosity() const {return verbosity;}

 protected:

  /** ctor.
  */
  Fun4AllBase(const std::string &name = "NONAME");

  std::string ThisName;

  /// The verbosity level. 0 means not verbose at all.
  int verbosity;
};

#endif /* __FUN4ALLBASE_H__ */

