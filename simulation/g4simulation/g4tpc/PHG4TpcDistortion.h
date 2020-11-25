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
  explicit PHG4TpcDistortion() = default;

  //! destructor
  ~PHG4TpcDistortion() = default;

  //!@name accessors
  //@{

  //! x distortion for a given truth location of the primary ionization
  double get_x_distortion(double x, double y, double z);

  //! y distortion for a given truth location of the primary ionization
  double get_y_distortion(double x, double y, double z);

  //! z distortion for a given truth location of the primary ionization
  double get_z_distortion(double x, double y, double z);

  //! Gets the verbosity of this module.
  int Verbosity() const
  { return verbosity; }

  //@}

  //!@name modifiers
  //@{

  //! static distortion filename
  void set_static_distortion_filename( const std::string& value )
  { m_static_distortion_filename = value; }

  //! time ordered distortion filename
  void set_time_ordered_distortion_filename( const std::string& value )
  { m_time_ordered_distortion_filename = value; }

  //! initialize
  void Init( bool do_time_ordered_distortion, bool do_static_distortion );

  //! get relevant histogram from time ordered TTrees
  void load_event(int event_num);

  //! Sets the verbosity of this module (0 by default=quiet).
  virtual void Verbosity(const int ival)
  { verbosity = ival; }

  //@}

  private:

  //! The verbosity level. 0 means not verbose at all.
  int verbosity = 0;

  //!@name static histograms
  //@{
  std::string m_static_distortion_filename = "$CALIBRATIONROOT/TPC/DistortionMaps/fluct_average.rev3.1side.3d.file0.h_negz.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root";
  std::unique_ptr<TFile> m_static_tfile;
  TH3 *hDXint = nullptr;
  TH3 *hDYint = nullptr;
  TH3 *hDZint = nullptr;
  //@}

  //!@name time ordered histograms
  //@{
  std::string m_time_ordered_distortion_filename = "/gpfs/mnt/gpfs02/sphenix/user/klest/TimeOrderedDistortions.root";
  std::unique_ptr<TFile> m_time_ordered_tfile;
  TTree *TimeTree = nullptr;
  TH3 *TimehDX = nullptr;
  TH3 *TimehDY = nullptr;
  TH3 *TimehDZ = nullptr;
  //@}


};

#endif /* G4TPC_PHG4TPCDISTORTION_H */
