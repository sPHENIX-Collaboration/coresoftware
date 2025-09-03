/// ===========================================================================
/*! \file    BeamBackgroundFilter.h
 *  \authors Derek Anderson
 *  \date    11.01.2024
 *
 *  Part of the BeamBackgroundFilterAndQA module, this
 *  is an abstract template for any filters run as
 *  part of the module.
 */
/// ===========================================================================

#ifndef BASEBEAMBACKGROUNDFILTER_H
#define BASEBEAMBACKGROUNDFILTER_H

// f4a libraries
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

// root libraries
#include <TH1.h>

// c++ utilities
#include <map>
#include <string>

// forward declarations
class PHCompositeNode;

// ============================================================================
//! Base beam background filter
// ============================================================================
/*! Base class for filters to be applied in the BeamBackgroundFilterAndQA
 *  module. Defines all the machinery common between filters.
 */
class BaseBeamBackgroundFilter
{
 public:
  // ------------------------------------------------------------------------
  //! Apply filter
  // ------------------------------------------------------------------------
  /*! Applies the filter. Should return true if filter finds beam
   *  background, and false if not.
   */
  virtual bool ApplyFilter(PHCompositeNode* /*topNode*/) { return false; }

  // ------------------------------------------------------------------------
  //! Build associated histograms
  // ------------------------------------------------------------------------
  /*! Collects all definitions of histograms.
   */
  virtual void BuildHistograms(const std::string& /*module*/, const std::string& /*tag*/) { return; }

  ///! register histograms
  inline void RegisterHistograms(Fun4AllHistoManager* manager)
  {
    for (auto& hist : m_hists)
    {
      manager->registerHisto(hist.second);
    }
    return;
  }

  ///! Set filter name
  void SetName(const std::string& name) { m_name = name; }

  ///! Get filter name
  const std::string& GetName() { return m_name; }

  ///! default ctor/dtor
  BaseBeamBackgroundFilter(const std::string& name)
    : m_name(name) {};
  virtual ~BaseBeamBackgroundFilter() = default;

 protected:
  // ------------------------------------------------------------------------
  //! Grab relevant input nodes
  // ------------------------------------------------------------------------
  /*! Collects all calls to `findNode`.
   */
  virtual void GrabNodes(PHCompositeNode* /*topNode*/) { return; }

  // ------------------------------------------------------------------------
  //! Histograms
  // ------------------------------------------------------------------------
  /*! All QA histograms for a given filter should be defined in this
   *  map, e.g.
   *
   *  m_hists["hNStreakPhi"] = new TH2D("hNStreakPhi", "", 64, 0., 64., 10, 0., 10.);
   */
  std::map<std::string, TH1*> m_hists;

  ///! filter name
  std::string m_name;

};  // end BaseBeamBackgroundFilter

#endif

// end ========================================================================
