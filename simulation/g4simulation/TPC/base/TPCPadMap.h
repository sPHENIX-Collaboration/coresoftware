#ifndef  __TPCPADMAP_H__
#define __TPCPADMAP_H__

#include <vector>
#include <map>
#include <utility>

#include "TPCDataTypes.h"
#include "TPCConstants.h"

using namespace TPCDataTypes;
using namespace TPCConstants;

class TPCPadMap {
 public:
  TPCPadMap();
  virtual ~TPCPadMap() {};
  Int_t GetNumberOfModules();
  Int_t GetReducedModule(Int_t mod);
  Int_t GetNumberOfPads(Module_t mod);
  Int_t GetNumberOfRows(Module_t mod);
  Int_t GetNumberOfCols(Module_t mod, Pad_t prow);

  Int_t GetRow(Module_t mod, Pad_t pad);
  Int_t GetCol(Module_t mod, Pad_t pad);

  Int_t GetSection(Module_t mod);
  Int_t GetSector(Module_t mod);
  Pad_t GetPad(Module_t mod, Pad_t prow, Pad_t pcol);
  Pad_t GetPad(Module_t mod, Float_t rad, Float_t phi);

  PairOfFloats_t GetXY(Module_t mod, Pad_t pad);
  PairOfFloats_t GetRP(Module_t mod, Pad_t pad);
  PairOfFloats_t GetDXDY(Module_t mod, Pad_t pad);

  PairOfFloats_t GetXY(Module_t mod, Pad_t prow, Pad_t pcol);
  PairOfFloats_t GetRP(Module_t mod, Pad_t prow, Pad_t pcol);
  PairOfFloats_t GetDXDY(Module_t mod, Pad_t prow, Pad_t pcol);

  void DumpModule(Module_t mod);

  ModuleRange_t FindModuleRP(Float_t rad, Float_t phi, Float_t rms, Float_t ele);
  PadQuotaRange_t GetPadQuotasRP(Float_t ele, Module_t mod, Float_t rad, Float_t phi, Float_t rms );
  TimeQuotaRange_t GetTimeQuotas(Float_t ele, Float_t mu, Float_t rms0, Float_t rms1);

 protected:
  void GeneratePadsV1();

  Pad_t fNPadRows[kNModulesPerPlate]; // stores number of rows per module
  PadRange_t fNPadCols[kNModulesPerPlate]; // stores number of cols per module-row
  MapPadXY_t fXY[kNModulesPerPlate]; // center of pad in global cartesian coordinates
  MapPadXY_t fRP[kNModulesPerPlate]; // center of pad in global polar coordinates
  MapPadXY_t fDXDY[kNModulesPerPlate]; // half width and height of pad in global cartesian coordinates
};

#endif /* __TPCPADMAP_H__ */
