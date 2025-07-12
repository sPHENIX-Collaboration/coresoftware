// ----------------------------------------------------------------------------
/*! \file    JetQADefs.h
 *  \authors Derek Anderson
 *  \date    06.24.2024
 *
 *  A namespace to hold definitions of various
 *  enums, helper methods, etc. relevant to
 *  all jet QA modules.
 */
// ----------------------------------------------------------------------------

#ifndef JETQADEFS_H
#define JETQADEFS_H

#include <calotrigger/TriggerAnalyzer.h>

#include <phool/phool.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// ----------------------------------------------------------------------------
//! Namespace to hold misc. definitions for the Jet QA
// ----------------------------------------------------------------------------
namespace JetQADefs
{

  // enums -------------------------------------------------------------------

  // -------------------------------------------------------------------------
  //! Trigger flags
  // -------------------------------------------------------------------------
  /*! Enumerates possible triggers for easy reference in QA modules
   *  and Fun4all macros. Note that these numbers ARE NOT necessarily 
   *  matched to the trigger bit numbers
   */
  enum GL1
  {
    Clock           = 0,  /*!< Clock */
    ZDCS            = 1,  /*!< ZDC South */
    ZDCN            = 2,  /*!< ZDC North */
    ZDCNS           = 3,  /*!< ZDC Coincidence */
    HCalSingle      = 4,  /*!< HCAL Singles */
    HCalCoin        = 5,  /*!< HCAL Coincidence */
    MBDS            = 8,  /*!< MBD S >= 1 */
    MBDN            = 9,  /*!< MBD N >= 1 */
    MBDNS1          = 10, /*!< MBD N&S >= 1 */
    MBDNS2          = 11, /*!< MBD N&S >= 2 */
    MBDNSVtx10      = 12, /*!< MBD N&S >= 1, vtx < 10 cm */
    MBDNSVtx30      = 13, /*!< MBD N&S >= 1, vtx < 30 cm */
    MBDNSVtx60      = 14, /*!< MBD N&S >= 1, vtx < 60 cm */
    MBDNSHCalSingle = 15, /*!< HCAL Singles + MBD NS >= 1 */
    MBDNSJet1       = 16, /*!< Jet 6 GeV + MBD NS >= 1 */
    MBDNSJet2       = 17, /*!< Jet 8 GeV + MBD NS >= 1 */
    MBDNSJet3       = 18, /*!< Jet 10 GeV + MBD NS >= 1 */
    MBDNSJet4       = 19, /*!< Jet 12 GeV + MBD NS >= 1 */
    Jet1            = 20, /*!< Jet 6 GeV (no MBD coincidence) */
    Jet2            = 21, /*!< Jet 8 GeV (no MBD coincidence) */
    Jet3            = 22, /*!< Jet 10 GeV (no MBD coincidence) */
    Jet4            = 23, /*!< Jet 12 GeV (no MBD coincidence) */
    MBDNSPhoton1    = 24, /*!< Photon 2 GeV+ MBD NS >= 1 */
    MBDNSPhoton2    = 25, /*!< Photon 3 GeV + MBD NS >= 1 */
    MBDNSPhoton3    = 26, /*!< Photon 4 GeV + MBD NS >= 1 */
    MBDNSPhoton4    = 27, /*!< Photon 5 GeV + MBD NS >= 1 */
    Photon1         = 28, /*!< Photon 2 GeV (no MBD coincidence) */
    Photon2         = 29, /*!< Photon 3 GeV (no MBD coincidence) */
    Photon3         = 30, /*!< Photon 4 GeV (no MBD coincidence) */
    Photon4         = 31, /*!< Photon 5 GeV (no MBD coincidence) */
    MBDNS2Vtx10     = 32, /*!< MBD N&S >= 2, vtx < 10 cm */
    MBDNS2Vtx30     = 33, /*!< MBD N&S >= 2, vtx < 30 cm */
    MBDNS2Vtx60     = 34, /*!< MBD N&S >= 2, vtx < 60 cm */
    MBDNS2Vtx150    = 35, /*!< MBD N&S >= 2, vtx < 150 cm */
    MBDNS2Photon6Vtx10  = 36, /*!< Photon 6 GeV + MBD NS >= 2, vtx < 10 cm */
    MBDNS2Photon8Vtx10  = 37, /*!< Photon 8 GeV + MBD NS >= 2, vtx < 10 cm */
    MBDNS2Photon10Vtx10 = 38, /*!< Photon 10 GeV + MBD NS >= 2, vtx < 10 cm */
    MBDNS2Photon12Vtx10 = 39, /*!< Photon 12 GeV + MBD NS >= 2, vtx < 10 cm */
    MBDNS2Photon6Vtx150 = 40, /*!< Photon 6 GeV + MBD NS >= 2, vtx < 150 cm */
    MBDNS2Photon8Vtx150 = 41, /*!< Photon 8 GeV + MBD NS >= 2, vtx < 150 cm */
    MBDNS2Photon10Vtx150 = 42, /*!< Photon 10 GeV + MBD NS >= 2, vtx < 150 cm */
    MBDNS2Photon12Vtx150 = 43, /*!< Photon 12 GeV + MBD NS >= 2, vtx < 150 cm */
    Inclusive = 44  /*!< inclusive (i.e. NO trigger selection) */
  };

  // constants ----------------------------------------------------------------

  // --------------------------------------------------------------------------
  //! Map of trigger flag onto name
  // --------------------------------------------------------------------------
  /*! This helper method provides a map from the triggers
   *  enumerated by `JetQADefs::GL1` to their official
   *  names. Used by `JetQADefs::DidTriggerFire` to
   *  identify if a trigger fired.
   */
  inline std::map<uint32_t, std::string> MapTriggerFlagToName()
  {
    static const std::map<uint32_t, std::string> mapFlagToName = 
    {
      {Clock           , "Clock"},
      {ZDCS            , "ZDC South"},
      {ZDCN            , "ZDC North"},
      {ZDCNS           , "ZDC Coincidence"},
      {HCalSingle      , "HCAL Singles"},
      {HCalCoin        , "HCAL Coincidence"},
      {MBDS            , "MBD S >= 1"},
      {MBDN            , "MBD N >= 1"},
      {MBDNS1          , "MBD N&S >= 1"},
      {MBDNS2          , "MBD N&S >= 2"},
      {MBDNSVtx10      , "MBD N&S >= 1, vtx < 10 cm"},
      {MBDNSVtx30      , "MBD N&S >= 1, vtx < 30 cm"},
      {MBDNSVtx60      , "MBD N&S >= 1, vtx < 60 cm"},
      {MBDNS2Vtx10     , "MBD N&S >= 2, vtx < 10 cm"},
      {MBDNS2Vtx30     , "MBD N&S >= 2, vtx < 30 cm"},
      {MBDNS2Vtx60     , "MBD N&S >= 2, vtx < 60 cm"},
      {MBDNS2Vtx150    , "MBD N&S >= 2, vtx < 150 cm"},
      {MBDNSHCalSingle , "HCAL Singles + MBD NS >= 1"},
      {MBDNSJet1       , "Jet 6 GeV + MBD NS >= 1"},
      {MBDNSJet2       , "Jet 8 GeV + MBD NS >= 1"},
      {MBDNSJet3       , "Jet 10 GeV + MBD NS >= 1"},
      {MBDNSJet4       , "Jet 12 GeV + MBD NS >= 1"},
      {Jet1            , "Jet 6 GeV"},
      {Jet2            , "Jet 8 GeV"},
      {Jet3            , "Jet 10 GeV"},
      {Jet4            , "Jet 12 GeV"},
      {MBDNSPhoton1    , "Photon 2 GeV+ MBD NS >= 1"},
      {MBDNSPhoton2    , "Photon 3 GeV + MBD NS >= 1"},
      {MBDNSPhoton3    , "Photon 4 GeV + MBD NS >= 1"},
      {MBDNSPhoton4    , "Photon 5 GeV + MBD NS >= 1"},
      {MBDNS2Photon6Vtx10   , "Photon 6 GeV + MBD NS >= 2, vtx < 10 cm"},
      {MBDNS2Photon8Vtx10   , "Photon 8 GeV + MBD NS >= 2, vtx < 10 cm"},
      {MBDNS2Photon10Vtx10  , "Photon 10 GeV + MBD NS >= 2, vtx < 10 cm"},
      {MBDNS2Photon12Vtx10  , "Photon 12 GeV + MBD NS >= 2, vtx < 10 cm"},
      {MBDNS2Photon6Vtx150   , "Photon 6 GeV + MBD NS >= 2, vtx < 150 cm"},
      {MBDNS2Photon8Vtx150   , "Photon 8 GeV + MBD NS >= 2, vtx < 150 cm"},
      {MBDNS2Photon10Vtx150  , "Photon 10 GeV + MBD NS >= 2, vtx < 150 cm"},
      {MBDNS2Photon12Vtx150  , "Photon 12 GeV + MBD NS >= 2, vtx < 150 cm"},
      {Photon1         , "Photon 2 GeV"},
      {Photon2         , "Photon 3 GeV"},
      {Photon3         , "Photon 4 GeV"},
      {Photon4         , "Photon 5 GeV"},
      {Inclusive       , "Inclusive"}
    };
    return mapFlagToName;
  }

  // methods ------------------------------------------------------------------

  // --------------------------------------------------------------------------
  //! Check if a particular trigger fired
  // --------------------------------------------------------------------------
  /*! Helper method to check if a particular trigger based on
   *  it's flag (see `JetQADefs::MapTriggerFlagToName()`).
   *  Returns value of `TriggerAnalyzer::didTriggerFire(
   *  std::string)`.
   *
   *  If unknown trigger flag is provided, method will
   *  cause Fun4All to abort.
   */ 
  inline bool DidTriggerFire(const uint32_t trg, TriggerAnalyzer* analyzer)
  {

    if (MapTriggerFlagToName().count(trg) == 0)
    {
      std::cerr << PHWHERE << ": PANIC: unknown trigger flag (" << trg << ") provided! Aborting" << std::endl;
      exit(1);
    }
    return analyzer->didTriggerFire( MapTriggerFlagToName()[trg] );

  }  // end 'DidTriggerFire(uint32_t, TriggerAnalyzer*)'

  // ==========================================================================
  //! Make a single QA-compliant histogram name
  // ==========================================================================
  /*! This helper method takes in a base name (e.g. some variable you want
   *  to histogram like "JetEne") and produces a histogram name compliant
   *  w/ the rest of the jet QA.
   *
   *  The format should always be:
   *      h_<module name>_<trigger tag>_<jet tag>_<base name> + <tag>
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

  }  // end 'MakeQAHistName(std::string& x 3)'

  // ==========================================================================
  //! Make multiple QA-compliant histogram names
  // ==========================================================================
  /*! This helper method takes in a list of base names (e.g.
   *  some variable you want to histogram like "JetEne") and
   *  produces a list of histogram names compliant w/ the
   *  rest of the jet QA.
   *
   *  The format should always be:
   *    h_<module name>_<trigger tag>_<jet tag>_<base name> + <tag>
   */
  inline std::vector<std::string> MakeQAHistNames(
      const std::vector<std::string>& bases,
      const std::string& module,
      const std::string& tag = "")
  {
    // copy base names to list of hist names
    std::vector<std::string> names = bases;

    // inject module names, tags, etc.
    for (auto& name : names)
    {
      name.insert(0, "h_" + module + "_");
      if (!tag.empty())
      {
        name.append("_" + tag);
      }
      std::transform(
          name.begin(),
          name.end(),
          name.begin(),
          ::tolower);
    }
    return names;

  }  // end 'MakeQAHistNames(std::vector<std::string>&, std::string& x 2)'

}  // namespace JetQADefs

#endif

// end ------------------------------------------------------------------------
