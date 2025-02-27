// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALL_FUN4ALLBASE_H
#define FUN4ALL_FUN4ALLBASE_H

#include <cstdint>
#include <limits>  // std::numeric_limits
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
  virtual const std::string Name() const { return m_ThisName; }

  /// Sets the name of this module.
  virtual void Name(const std::string &name) { m_ThisName = name; }

  /** Print out some info about this module.
      @param what can be used to specify what to print exactly.
  */
  virtual void Print(const std::string &what = "ALL") const;

  enum enu_Verbosity
  {

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
    VERBOSITY_MAX = std::numeric_limits<int>::max() - 10

  };

  /// Sets the verbosity of this module (0 by default=quiet).
  virtual void Verbosity(const uint64_t ival) { m_Verbosity = ival; }

  /// Sets the verbosity of this module (0 by default=quiet).
  virtual void Verbosity(enu_Verbosity ival) { m_Verbosity = ival; }

  /// Gets the verbosity of this module.
  virtual uint64_t Verbosity() const { return m_Verbosity; }

  /// Gets the downscale of printouts for this module.
  virtual uint32_t VerbosityDownscale() const { return m_VerbosityDownscale; }

  /// Sets the downscale of printouts for this module.
  /// Use (cnt%VerbosityDownscale()) in the code to apply
  virtual void VerbosityDownscale(uint32_t ival) { m_VerbosityDownscale = std::max(1U, ival); }

 protected:
  /** ctor.
   */
  Fun4AllBase(const std::string &name = "NONAME");

 private:
  std::string m_ThisName;

  /// The verbosity level. 0 means not verbose at all.
  uint64_t m_Verbosity{VERBOSITY_QUIET};
  /// The frequency of printouts
  uint32_t m_VerbosityDownscale{1};
};

#endif
