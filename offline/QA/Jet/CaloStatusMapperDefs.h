/// ===========================================================================
/*! \file   CaloStatusMapperDefs.h
 *  \author Derek Anderson
 *  \date   05.22.2024
 *
 *  A Fun4All QA module to plot no. of towers per event
 *  and vs. eta, phi as a function of status.
 */
/// ===========================================================================

#ifndef CLUSTERSTATUSMAPPER_DEFS_H
#define CLUSTERSTATUSMAPPER_DEFS_H

// calobase libraries
#include <calobase/TowerInfov2.h>

// root libraries
#include <TH1.h>
#include <TH2.h>

// c++ utilities
#include <map>
#include <string>
#include <utility>
#include <vector>



// ============================================================================
//! Miscellaneous definitions for CaloStatusMapper
// ============================================================================
/*! A namespace to collect various constants, helper
 *  methods, etc. used in the CaloStatusMapper
 *  module.
 */
namespace CaloStatusMapperDefs
{

  // convenience types
  typedef std::pair<std::string, int> NodeDef;



  // ==========================================================================
  //! Enumeration of calorimeters
  // ==========================================================================
  /*! This enumerates the different calorimeters one can
   *  map in the module.
   */
  enum Calo
  {
    EMCal,  ///!< EMCal geometry
    HCal,   ///!< I/OHCal geometry
    NONE  ///!< Unspecified geometry
  };



  // ==========================================================================
  //! Possible status codes
  // ==========================================================================
  /*! This enumerates possible status codes of towers.
   */
  enum Stat
  {
    Good,
    Hot,
    BadTime,
    BadChi,
    NotInstr,
    NoCalib,
    Unknown
  };



  // ==========================================================================
  //! Maps status codes onto labels
  // ==========================================================================
  inline std::map<Stat, std::string> const& StatLabels()
  {
    static std::map<Stat, std::string> mapStatLabels = {
      {Stat::Good,     "Good"},
      {Stat::Hot,      "Hot"},
      {Stat::BadTime,  "BadTime"},
      {Stat::BadChi,   "BadChi"},
      {Stat::NotInstr, "NotInstr"},
      {Stat::NoCalib,  "NoCalib"},
      {Stat::Unknown,  "Unknown"}
    };
    return mapStatLabels;
  }



  // ==========================================================================
  //! Helper struct to define an axis of a histogram
  // ==========================================================================
  /*! This is a lightweight struct to consolidate information
   *  necessary to define a histogram axis.
   */
  struct AxisDef
  {

    // members
    std::string label {""};  ///! axis label
    std::size_t nBins {10};  ///! no. of bins in an axis
    float       start {0.};  ///! low edge of 1st bin
    float       stop  {1.};  ///! high edge of last bin

  };  // end AxisDef



  // ==========================================================================
  //! Helper struct to define histograms
  // ==========================================================================
  /*! This is a lightweight struct to organize the axis definitions
   *  for status, eta, phi axes and make histograms using them.
   */ 
  template <std::size_t H, std::size_t F, std::size_t S>
  struct HistDef
  {

    // axis definitions
    AxisDef stat {"Status", S, -0.5, S - 0.5};
    AxisDef eta  {"i_{#eta}", H, -0.5, H - 0.5};
    AxisDef phi  {"i_{#phi}", F, -0.5, F - 0.5};

    //! make 1 1d status plot
    TH1D* MakeStatus1D(const std::string& name) const
    {
      const std::string title = ";" + stat.label;
      return new TH1D(name.data(), title.data(), stat.nBins, stat.start, stat.stop);
    }

    //! make a 1d eta plot
    TH1D* MakeEta1D(const std::string& name) const
    {
      const std::string title = ";" + eta.label;
      return new TH1D(name.data(), title.data(), eta.nBins, eta.start, eta.stop);
    }

    //! make a 1d phi plot
    TH1D* MakePhi1D(const std::string& name) const
    {
      const std::string title = ";" + phi.label;
      return new TH1D(name.data(), title.data(), phi.nBins, phi.start, phi.stop);
    }

    //! make a 2d eta-phi plot
    TH2D* MakePhiEta2D(const std::string& name) const
    {
      const std::string title = ";" + eta.label + ";" + phi.label;
      return new TH2D(name.data(), title.data(), eta.nBins, eta.start, eta.stop, phi.nBins, phi.start, phi.stop);
    }

  };  // end HistDef

  // -------------------------------------------------------------------------
  //! Maps for specific calorimeters
  // -------------------------------------------------------------------------
  typedef HistDef<96, 256, 7> EMCalHistDef;
  typedef HistDef<24, 64, 7> HCalHistDef;



  // ==========================================================================
  //! Returns enum corresponding to given tower status
  // ==========================================================================
  /*! This helper methods returns the associated code of
   *  the provided tower.
   */ 
  Stat GetTowerStatus(TowerInfo* tower)
  {

    Stat status = Stat::Unknown;
    if (tower -> get_isHot())
    {
      status = Stat::Hot;
    }
    else if (tower -> get_isBadTime())
    {
      status = Stat::BadTime;
    }
    else if (tower -> get_isBadChi2())
    {
      status = Stat::BadChi;
    }
    else if (tower -> get_isNotInstr())
    {
      status = Stat::NotInstr;
    }
    else if (tower -> get_isNoCalib())
    {
      status = Stat::NoCalib;
    }
    else
    {
      status = Stat::Good;
    }
    return status;

  }  // end 'GetTowerStatus(TowerInfo*)'



  // ==========================================================================
  //! Make QA-compliant histogram name
  // ==========================================================================
  /*! This helper method takes in a base name (e.g. some variable you want
   *  to histogram like "JetEne") and produces a histogram name compliant
   *  w/ the rest of the jet QA.
   *
   *  The format should always be:
   *      h_<module name>_<trigger tag>_<jet tag>_<base name> + <tag>
   *
   *  FIXME this should get moved into JetQADefs.h
   */
  inline std::string MakeQAHistName(
    const std::string& base,
    const std::string& module,
    const std::string& tag = "")
  {

    // set name to base
    std::string name = base;

    // inject module names, tags, etc.
    name.insert(0, "h_" + module + "_");
    if (!tag.empty())
    {
      name.append("_" + tag);
    }
    std::transform(
      name.begin(),
      name.end(),
      name.begin(),
      ::tolower
    );
    return name;

  }  // end 'MakeQAHistNames(std::string& x 3)'

}  // end CaloStatusMapperDefs namespace

#endif

// end ========================================================================

