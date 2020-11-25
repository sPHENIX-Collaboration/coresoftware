// $Id: $

/*!
 * \file PHG4TPCDistortion.h
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>, Henry Klest <henry.klest@stonybrook.edu>, Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef G4TPC_PHG4TPCDISTORTION_H
#define G4TPC_PHG4TPCDISTORTION_H

#include <memory>

class TFile;
class TH3;
class TTree;

//! handle distortions (static and time-ordered)
class PHG4TpcDistortion
{
  public:

  //! constructor
  explicit PHG4TpcDistortion(bool do_time_ordered_distortion = false, bool do_static_distortion = false);

  //! destructor
  ~PHG4TpcDistortion() = default;

  //! x distortion for a given truth location of the primary ionization
  double get_x_distortion(double x, double y, double z)
  { return get_distortion(hDXint, TimehDX, x, y, z ); }

  //! y distortion for a given truth location of the primary ionization
  double get_y_distortion(double x, double y, double z)
  { return get_distortion(hDYint, TimehDY, x, y, z ); }

  //! z distortion for a given truth location of the primary ionization
  double get_z_distortion(double x, double y, double z)
  { return get_distortion(hDZint, TimehDZ, x, y, z ); }

  //! get relevant histogram from time ordered TTrees
  void load_event(int event_num);

  //! Sets the verbosity of this module (0 by default=quiet).
  virtual void Verbosity(const int ival)
  { verbosity = ival; }

  //! Gets the verbosity of this module.
  int Verbosity() const
  { return verbosity; }

  private:

  static double get_distortion( TH3* hstatic, TH3* htimeOrdered, double x, double y, double z );

  //! The verbosity level. 0 means not verbose at all.
  int verbosity = 0;

  //!@name static histograms
  //@{
  std::unique_ptr<TFile> m_staticTFile;
  TH3 *hDXint = nullptr;
  TH3 *hDYint = nullptr;
  TH3 *hDZint = nullptr;
  //@}

  //!@name time ordered histograms
  //@{
  std::unique_ptr<TFile> m_timeOrderedTFile;
  TTree *TimeTree = nullptr;
  TH3 *TimehDX = nullptr;
  TH3 *TimehDY = nullptr;
  TH3 *TimehDZ = nullptr;
  //@}


};

#endif /* G4TPC_PHG4TPCDISTORTION_H */
