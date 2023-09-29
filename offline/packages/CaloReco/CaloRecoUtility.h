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
  CaloRecoUtility(CaloRecoUtility& cru);
  CaloRecoUtility& operator=(CaloRecoUtility const&);

  //! corrects cluster Z (implicitly also eta) for updated z vertex
  // assuming
  void ShowerDepthCorrZVertex(RawCluster* clus, float zvtx);
  void ProbCorrsZVertex(RawCluster* clus, float zvtx);
  void LoadProfile();

 private:
  bool _profLoaded;
  BEmcRec* _bemc = nullptr;
};

#endif
