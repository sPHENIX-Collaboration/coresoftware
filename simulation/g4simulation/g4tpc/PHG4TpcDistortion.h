// $Id: $

/*!
 * \file PHG4TPCDistortion.h
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>, Henry Klest <henry.klest@stonybrook.edu>, Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>, Ross Corliss <ross.corliss@stonybrook.edu>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef G4TPC_PHG4TPCDISTORTION_H
#define G4TPC_PHG4TPCDISTORTION_H

#include <memory>
#include <string>

class TFile;
class TH3;
class TTree;

//! handle distortions (static and time-ordered)
class PHG4TpcDistortion
{
 public:
  //! constructor
  explicit PHG4TpcDistortion() = default;

  //!@name accessors
  //@{

  //! x distortion for a given truth location of the primary ionization
  double get_x_distortion_cartesian(double x, double y, double z) const;

  //! y distortion for a given truth location of the primary ionization
  double get_y_distortion_cartesian(double x, double y, double z) const;

  //! z distortion for a given truth location of the primary ionization
  double get_z_distortion_cartesian(double x, double y, double z) const;

  //! radial distortion for a given cylindrical truth location of the primary ionization
  double get_r_distortion(double r, double phi, double z) const;

  //! R*phi (hence the unitful phi-hat) distortion for a given cylindrical truth location of the primary ionization
  double get_rphi_distortion(double r, double phi, double z) const;

  //! z distortion for a given cylindrical truth location of the primary ionization
  double get_z_distortion(double r, double phi, double z) const;
  
  //The ReachesReadout serves as a fourth axis in the distortion histogram
  double get_reaches_readout(double r, double phi, double z) const;

  //! Gets the verbosity of this module.
  int Verbosity() const
  {
    return verbosity;
  }

  //@}

  //!@name modifiers
  //@{

  //! enable distortions from a single fixed source
  void set_do_static_distortions(bool value)
  {
    m_do_static_distortions = value;
  }

  //! set the filename for the single fixed distortion
  void set_static_distortion_filename(const std::string &value)
  {
    m_static_distortion_filename = value;
  }

  //! enable time ordered distortions
  void set_do_time_ordered_distortions(bool value)
  {
    m_do_time_ordered_distortions = value;
  }

  //! time ordered distortion filename
  void set_time_ordered_distortion_filename(const std::string &value)
  {
    m_time_ordered_distortion_filename = value;
  }

  void set_read_phi_as_radians(bool flag=true){
    m_phi_hist_in_radians=flag;
  }

  //! initialize
  void Init();

  //! get relevant histogram from time ordered TTrees
  void load_event(int event_num);

  //! Sets the verbosity of this module (0 by default=quiet).
  void Verbosity(const int ival)
  {
    verbosity = ival;
  }

  //@}

 private:
  //! get distortion for a set of histogram and an input momentum distribution
  double get_distortion(char axis, double r, double phi, double z) const;

  //! The verbosity level. 0 means not verbose at all.
  int verbosity = 0;

  //! Flag controls whether to assume the phi hist units are radians or cm.
  bool m_phi_hist_in_radians=true;

  
  //!@name static histograms
  //@{
  bool m_do_static_distortions = false;
  std::string m_static_distortion_filename;
  std::unique_ptr<TFile> m_static_tfile;
  TH3 *hDRint[2] = {nullptr, nullptr};
  TH3 *hDPint[2] = {nullptr, nullptr};
  TH3 *hDZint[2] = {nullptr, nullptr};
  TH3 *hReach[2] = {nullptr, nullptr};
  //@}

  //!@name time ordered histograms
  //@{
  bool m_do_time_ordered_distortions = false;
  std::string m_time_ordered_distortion_filename;
  std::unique_ptr<TFile> m_time_ordered_tfile;
  TTree *TimeTree = nullptr;
  TH3 *TimehDR[2] = {nullptr, nullptr};
  TH3 *TimehDP[2] = {nullptr, nullptr};
  TH3 *TimehDZ[2] = {nullptr, nullptr};
  TH3 *TimehRR[2] = {nullptr, nullptr};
  //@}
};

#endif /* G4TPC_PHG4TPCDISTORTION_H */
