#ifndef MBD_MBDCALIB_H
#define MBD_MBDCALIB_H

#include "MbdDefs.h"

#include <fun4all/Fun4AllBase.h>

#include <phool/recoConsts.h>

#include <array>
#include <vector>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>

#include <mbd/MbdGeom.h>

class TTree;
class CDBInterface;

class MbdCalib : public Fun4AllBase
{
 public:
  MbdCalib();

  // MbdCalib(MbdCalib &other) = delete;
  // void operator=(const MbdCalib &) = delete;

  virtual ~MbdCalib() {}

  float get_qgain(const int ipmt) const { return _qfit_mpv[ipmt]; }
  float get_tq0(const int ipmt) const { return _tqfit_t0mean[ipmt]; }
  int get_sampmax(const int ifeech) const { return _sampmax[ifeech]; }
  std::vector<float> get_shape(const int ifeech) const { return _shape_y[ifeech]; }
  std::vector<float> get_sherr(const int ifeech) const { return _sherr_yerr[ifeech]; }

  int Download_Gains(const std::string& dbfile);
  int Download_TQT0(const std::string& dbfile);
  int Download_TTT0(const std::string& dbfile);
  int Download_Slew(const std::string& dbfile);
  int Download_SampMax(const std::string& dbfile);
  int Download_Shapes(const std::string& dbfile);
  int Download_All();

  int Write_CDB_SampMax(const std::string& dbfile);
  int Write_CDB_Shapes(const std::string& dbfile);
  int Write_CDB_All();

  // void Dump_to_file(const std::string& what = "ALL");

  void Reset();
  // void Print(Option_t* option) const;

 private:
  CDBInterface* _cdb{nullptr};
  recoConsts* _rc{nullptr};

  std::unique_ptr<MbdGeom> _mbdgeom{nullptr};

  int _status{0};
  // int          _run_number {0};
  // uint64_t     _timestamp {0};
  std::string _dbfilename;

  // Assumes Landau fit
  std::array<float, MbdDefs::BBC_N_PMT> _qfit_integ{};
  std::array<float, MbdDefs::BBC_N_PMT> _qfit_mpv{};
  std::array<float, MbdDefs::BBC_N_PMT> _qfit_sigma{};
  std::array<float, MbdDefs::BBC_N_PMT> _qfit_integerr{};
  std::array<float, MbdDefs::BBC_N_PMT> _qfit_mpverr{};
  std::array<float, MbdDefs::BBC_N_PMT> _qfit_sigmaerr{};
  std::array<float, MbdDefs::BBC_N_PMT> _qfit_chi2ndf{};

  // T0 offsets, time channels
  std::array<float, MbdDefs::BBC_N_PMT> _ttfit_t0mean{};
  std::array<float, MbdDefs::BBC_N_PMT> _ttfit_t0meanerr{};
  std::array<float, MbdDefs::BBC_N_PMT> _ttfit_t0sigma{};
  std::array<float, MbdDefs::BBC_N_PMT> _ttfit_t0sigmaerr{};

  // T0 offsets, charge channels
  std::array<float, MbdDefs::BBC_N_PMT> _tqfit_t0mean{};
  std::array<float, MbdDefs::BBC_N_PMT> _tqfit_t0meanerr{};
  std::array<float, MbdDefs::BBC_N_PMT> _tqfit_t0sigma{};
  std::array<float, MbdDefs::BBC_N_PMT> _tqfit_t0sigmaerr{};

  // Slew Correction
  std::array<int, MbdDefs::BBC_N_FEECH>   _slew_npts{};      // num points in template
  std::array<float, MbdDefs::BBC_N_FEECH> _slew_minrange{};  // in template units (samples)
  std::array<float, MbdDefs::BBC_N_FEECH> _slew_maxrange{};  // in template units (samples)
  std::array<std::vector<float>, MbdDefs::BBC_N_FEECH> _slew_y{};

  // Peak of waveform
  std::array<int, MbdDefs::BBC_N_FEECH> _sampmax{};

  // Waveform Template
  int do_templatefit{0};
  std::array<int, MbdDefs::BBC_N_FEECH>   _shape_npts{};     // num points in template
  std::array<float, MbdDefs::BBC_N_FEECH> _shape_minrange{}; // in template units (samples)
  std::array<float, MbdDefs::BBC_N_FEECH> _shape_maxrange{}; // in template units (samples)
  std::array<std::vector<float>, MbdDefs::BBC_N_FEECH> _shape_y{};

  std::array<int, MbdDefs::BBC_N_FEECH> _sherr_npts{};     // num points in template
  std::array<float, MbdDefs::BBC_N_FEECH> _sherr_minrange{}; // in template units (samples)
  std::array<float, MbdDefs::BBC_N_FEECH> _sherr_maxrange{}; // in template units (samples)
  std::array<std::vector<float>, MbdDefs::BBC_N_FEECH> _sherr_yerr{};
};

#endif  // MBD_MBDCALIB_H
