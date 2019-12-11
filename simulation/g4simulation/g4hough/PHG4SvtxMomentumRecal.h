#ifndef G4HOUGH_PHG4SVTXMOMENTUMRECAL_H
#define G4HOUGH_PHG4SVTXMOMENTUMRECAL_H

//===========================================================
/// \file PHG4SvtxMomentumRecal.h
/// \brief Projects into calorimeters and fills track cal fields
/// \author Mike McCumber
//===========================================================

#include <fun4all/SubsysReco.h>

// std includes
#include <cstddef>              // for NULL
#include <string>

// forward declarations
class PHCompositeNode;
class TF1;

/// \class PHG4SvtxMomentumRecal
///
/// \brief Takes a rescal factor as a function of reco pt
///
class PHG4SvtxMomentumRecal : public SubsysReco
{

 public:
 
  PHG4SvtxMomentumRecal(const std::string &name = "PHG4SvtxMomentumRecal",
			TF1* correction = nullptr);
  virtual ~PHG4SvtxMomentumRecal() {}
		
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  
 private:

  TF1* _corr;
};

#endif // __PHG4SVTXMOMENTUMRECAL_H__
