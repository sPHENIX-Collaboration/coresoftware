#ifndef __FUN4ALLBASE_H__
#define __FUN4ALLBASE_H__

#include <string>
#include <climits>       // std::numeric_limits

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
  virtual const std::string Name() const {return ThisName;}

  /// Sets the name of this module.
  virtual void Name(const std::string &name) {ThisName = name;}

  /** Print out some info about this module. 
      @param what can be used to specify what to print exactly.
  */
  virtual void Print(const std::string &what = "ALL") const;

  enum enu_Verbosity{

    //! Quiet mode. Only output critical messages. Intended for batch production mode.
    VERBOSITY_QUIET = 0,

    //! Output some useful messages during manual command line running
    VERBOSITY_SOME = 1,

    //! Output more messages
    VERBOSITY_MORE = 2,

    //! Output even more messages
    VERBOSITY_EVEN_MORE = 3,

    //! Output a lot of messages
    VERBOSITY_A_LOT = 4,

    // ... use your imagination ...

    //! Show all messages. Useful for step-by-step debugging
    VERBOSITY_MAX =  INT_MAX - 10

  } ;

  /// Sets the verbosity of this module (0 by default=quiet).
  virtual void Verbosity(const int ival) {verbosity = static_cast<enu_Verbosity>(ival);}

  /// Sets the verbosity of this module (0 by default=quiet).
  virtual void Verbosity(enu_Verbosity ival) {verbosity = ival;}

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

