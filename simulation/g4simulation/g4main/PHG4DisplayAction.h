// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4DISPLAYACTION_H
#define G4MAIN_PHG4DISPLAYACTION_H

#include <string>

class G4VPhysicalVolume;

class PHG4DisplayAction
{
 public:
  //! constructor
  // delete default ctor, nobody should use it
  PHG4DisplayAction() = delete;
  // this is the ctor we use
  PHG4DisplayAction(const std::string &name)
    : m_Detector(name)
  {
  }

  //! destructor
  virtual ~PHG4DisplayAction() {}

  //! ApplyDisplayAction method
  /**
   pure virtual - has to be implemented by derived class
   creates and set VisAttributes for volumes
   @param[in] physvol starting volume in hierarchy (typically world volume)
  */

  virtual void ApplyDisplayAction(G4VPhysicalVolume *physvol) = 0;

  virtual void SetName(const std::string &name) { m_Detector = name; }

  virtual std::string GetName() const { return m_Detector; }

  virtual void Print(const std::string &/*what*/="ALL") {}

  enum CheckReturnCodes
  {
    ABORT = -1,
    FAILED = 0,
    ACCEPT = 1
  };

 protected:
  //! find FindVolume method
  /*
 * @param[in] starting volume
 */
  int FindVolumes(G4VPhysicalVolume *physvol);

  //! find CheckVolume method
  /*
 * @param[in] physical volume to be checked
 */
  virtual int CheckVolume(G4VPhysicalVolume */*physvol*/) { return 0; }

  //! ApplyVisAttributes method
  /**
 *@param[in] physvol selected physical volume
 */
  virtual void ApplyVisAttributes(G4VPhysicalVolume */*physvol*/) { return; }

 private:
  std::string m_Detector;
};

#endif  // G4MAIN_PHG4DISPLAYACTION_H
