#ifndef __PHG4SVTXMOMENTUMRECAL_H__
#define __PHG4SVTXMOMENTUMRECAL_H__

//===========================================================
/// \file PHG4SvtxMomentumRecal.h
/// \brief Projects into calorimeters and fills track cal fields
/// \author Mike McCumber
//===========================================================

#include "PHG4HoughTransform.h"
#include  <trackbase_historic/SvtxTrack.h>

// PHENIX includes
#include <fun4all/SubsysReco.h>

#include <TF1.h>

// std includes
#include <vector>
#include <string>

// forward declarations
class PHCompositeNode;

/// \class PHG4SvtxMomentumRecal
///
/// \brief Takes a rescal factor as a function of reco pt
///
class PHG4SvtxMomentumRecal : public SubsysReco
{

 public:
 
  PHG4SvtxMomentumRecal(const std::string &name = "PHG4SvtxMomentumRecal",
			TF1* correction = NULL);
  virtual ~PHG4SvtxMomentumRecal() {}
		
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  
 private:

  TF1* _corr;
};

#endif // __PHG4SVTXMOMENTUMRECAL_H__
