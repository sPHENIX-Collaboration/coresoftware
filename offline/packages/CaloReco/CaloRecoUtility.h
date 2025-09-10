// $Id: $

/*!
 * \file CaloRecoUtility.h
 * \brief
 * \author Justin Frantz <frantz@ohio.edu>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef CALORECO_CALORECOUTILITY_H
#define CALORECO_CALORECOUTILITY_H

class RawCluster;
class BEmcRec;

/*!
 * \brief CaloRecoUtility

 */

class CaloRecoUtility
{
 public:
  ~CaloRecoUtility();
  CaloRecoUtility();
  // cppcheck: deleting copy ctor and = operator to prevent accidental use
  // if they are used at some point they need to be properly implemented
  // the default does not work for allocated memory
  CaloRecoUtility(const CaloRecoUtility& cru) = delete;
  CaloRecoUtility& operator=(CaloRecoUtility const&) = delete;

  //! corrects cluster Z (implicitly also eta) for updated z vertex
  // assuming
  static void ShowerDepthCorrZVertex(RawCluster* clus, float vz);
  void ProbCorrsZVertex(RawCluster* clus, float vz);
  void LoadProfile();

 private:
  bool _profLoaded {false};
  BEmcRec* _bemc {nullptr};
};

#endif
