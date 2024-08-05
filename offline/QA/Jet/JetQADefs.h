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

#include <ffarawobjects/Gl1Packet.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include <cassert>
#include <cmath>
#include <iostream>

// ----------------------------------------------------------------------------
//! Namespace to hold misc. definitions for the Jet QA
// ----------------------------------------------------------------------------
namespace JetQADefs
{

  // enums -------------------------------------------------------------------

  // -------------------------------------------------------------------------
  //! Trigger indices in GL1 packets
  // -------------------------------------------------------------------------
  /*! Gives indices of triggers as they are provided in the
   *  GL1 packet. The list is:
   *  | Index | Trigger                       |
   *  | ----- | ----------------------------- |
   *  | 0     | Clock                         |
   *  | 1     | ZDC South                     |
   *  | 2     | ZDC North                     |
   *  | 3     | ZDC Coincidence               |
   *  | 4     | Random                        |
   *  | 5     | HCAL Singles                  |
   *  | 8     | MBD S >= 1                    |
   *  | 9     | MBD N >= 1                    |
   *  | 10    | MBD N&S >= 1                  |
   *  | 11    | MBD N&S >= 2                  |
   *  | 12    | MBD N&S, vtx < 10 cm          |
   *  | 13    | MBD N&S, vtx < 30 cm          |
   *  | 14    | MBD N&S, vtx < 60 cm          |
   *  | 15    | HCAL Singles + MBD NS >= 1    |
   *  | 16    | Jet 1 + MBD NS >= 1           |
   *  | 17    | Jet 2 + MBD NS >= 1           |
   *  | 18    | Jet 3 + MBD NS >= 1           |
   *  | 19    | Jet 4 + MBD NS >= 1           |
   *  | 20    | Jet 1 (no MBD coincidence)    |
   *  | 21    | Jet 2 (no MBD coincidence)    |
   *  | 22    | Jet 3 (no MBD coincidence)    |
   *  | 23    | Jet 4 (no MBD coincidence)    |
   *  | 24    | Photon 1 + MBD NS >= 1        |
   *  | 25    | Photon 2 + MBD NS >= 1        |
   *  | 26    | Photon 3 + MBD NS >= 1        |
   *  | 27    | Photon 4 + MBD NS >= 1        |
   *  | 28    | Photon 1 (no MBD coincidence) |
   *  | 29    | Photon 2 (no MBD coincidence) |
   *  | 30    | Photon 3 (no MBD coincidence) |
   *  | 31    | Photon 4 (no MBD coincidence) |
   */
  enum GL1
  {
    Clock = 0,
    ZDCS = 1,
    ZDCN = 2,
    ZDCNS = 3,
    Random = 4,
    HCalSingle = 5,
    MBDS = 8,
    MBDN = 9,
    MBDNS1 = 10,
    MBDNS2 = 11,
    MBDNSVtx10 = 12,
    MBDNSVtx30 = 13,
    MBDNSVtx60 = 14,
    MBDNSHCalSingle = 15,
    MBDNSJet1 = 16,
    MBDNSJet2 = 17,
    MBDNSJet3 = 18,
    MBDNSJet4 = 19,
    Jet1 = 20,
    Jet2 = 21,
    Jet3 = 22,
    Jet4 = 23,
    MBDNSPhoton1 = 24,
    MBDNSPhoton2 = 25,
    MBDNSPhoton3 = 26,
    MBDNSPhoton4 = 27,
    Photon1 = 28,
    Photon2 = 29,
    Photon3 = 30,
    Photon4 = 31
  };

  // constants ----------------------------------------------------------------

  // --------------------------------------------------------------------------
  //! Max no. of trigger indices
  // --------------------------------------------------------------------------
  inline uint32_t NMaxTrgIndex()
  {
    static const uint32_t nMaxTrgIndex = 32;
    return nMaxTrgIndex;
  }

  // methods ------------------------------------------------------------------

  // --------------------------------------------------------------------------
  //! Check if a particular trigger fired
  // --------------------------------------------------------------------------
  inline bool DidTriggerFire(const uint32_t trg, PHCompositeNode* topNode)
  {
    // grab GL1 packet from node tree
    Gl1Packet* packet = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    if (!packet)
    {
      std::cerr << PHWHERE << ": PANIC: not able to grab GL1 packet! aborting!" << std::endl;
      exit(1);
    }

    // grab trigger bits
    boost::dynamic_bitset<> triggers(NMaxTrgIndex(), packet->getTriggerInput());

    // loop through bits and check if specified one is set
    bool didTrgFire = false;
    for (uint32_t iTrg = 0; iTrg < NMaxTrgIndex(); iTrg++)
    {
      // only consider specified trigger
      if (iTrg != trg) continue;

      // if trigger bit set, break and return true
      if (triggers.test(iTrg))
      {
        didTrgFire = true;
        break;
      }

    }  // end index loop
    return didTrgFire;

  }  // end 'DidTriggerFire(uint32_t, PHCompositeNode*)'

}  // namespace JetQADefs

#endif

// end ------------------------------------------------------------------------
