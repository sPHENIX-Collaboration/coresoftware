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
  int GetNumberOfModules();
  int GetReducedModule(int mod);
  int GetNumberOfPads(Module_t mod);
  int GetNumberOfRows(Module_t mod);
  int GetNumberOfCols(Module_t mod, Pad_t prow);

  int GetRow(Module_t mod, Pad_t pad);
  int GetCol(Module_t mod, Pad_t pad);

  int GetSection(Module_t mod);
  int GetSector(Module_t mod);
  Pad_t GetPad(Module_t mod, Pad_t prow, Pad_t pcol);
  Pad_t GetPad(Module_t mod, float rad, float phi);

  PairOfFloats_t GetXY(Module_t mod, Pad_t pad);
  PairOfFloats_t GetRP(Module_t mod, Pad_t pad);
  PairOfFloats_t GetDXDY(Module_t mod, Pad_t pad);

  PairOfFloats_t GetXY(Module_t mod, Pad_t prow, Pad_t pcol);
  PairOfFloats_t GetRP(Module_t mod, Pad_t prow, Pad_t pcol);
  PairOfFloats_t GetDXDY(Module_t mod, Pad_t prow, Pad_t pcol);

  void DumpModule(Module_t mod);

  ModuleRange_t FindModuleRP(float rad, float phi, float rms, float ele);
  PadQuotaRange_t GetPadQuotasRP(float ele, Module_t mod, float rad, float phi, float rms );
  TimeQuotaRange_t GetTimeQuotas(float ele, float mu, float rms0, float rms1);

 protected:
  void GeneratePadsV1();

  Pad_t fNPadRows[kNModulesPerPlate]; // stores number of rows per module
  PadRange_t fNPadCols[kNModulesPerPlate]; // stores number of cols per module-row
  MapPadXY_t fXY[kNModulesPerPlate]; // center of pad in global cartesian coordinates
  MapPadXY_t fRP[kNModulesPerPlate]; // center of pad in global polar coordinates
  MapPadXY_t fDXDY[kNModulesPerPlate]; // half width and height of pad in global cartesian coordinates
};

#endif /* __TPCPADMAP_H__ */
