/// ===========================================================================
/*! \file    NullFilter.h
 *  \authors Derek Anderson
 *  \date    11.20.2024
 *
 *  Part of the BeamBackgroundFilterAndQA module, this
 *  filter does nothing, but provides a template for
 *  others to copy and fill in.
 */
/// ===========================================================================

#ifndef NULLFILTER_H
#define NULLFILTER_H

// module components
#include "BaseBeamBackgroundFilter.h"

// c++ utilities
#include <string>

// forward declarations
class PHCompositeNode;

// ============================================================================
//! Null, template filter
// ============================================================================
/*! A beam background filter which does nothing, but
 *  provides a template for other filters.
 */
class NullFilter : public BaseBeamBackgroundFilter
{
 public:
  // ========================================================================
  //! User options for filter
  // ========================================================================
  struct Config
  {
    int verbosity = 0;
    bool debug = true;
    //... additional options go here ...//
  };

  // ctor/dtor
  NullFilter(const std::string& name = "Null");
  NullFilter(const Config& config, const std::string& name = "Null");
  virtual ~NullFilter() = default;

  // inherited methods
  bool ApplyFilter(PHCompositeNode* topNode) override;
  void BuildHistograms(const std::string& module, const std::string& tag = "") override;

 private:
  // inherited methods
  void GrabNodes(PHCompositeNode* /*topNode*/) override;

  // filter-specific methods
  //... go here ...//

  ///! input node(s)
  //... go here ...//

  ///! configuration
  Config m_config;

};  // end NullFilter

#endif

// end ========================================================================
