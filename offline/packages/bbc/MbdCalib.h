#ifndef __MBDCALIB_H__
#define __MBDCALIB_H__

#include "MbdDefs.h"
#include <phool/recoConsts.h>

#include <cmath>
#include <cstdint>
#include <cstring>

class TTree;
class CDBInterface;

class MbdCalib
{
public:
  MbdCalib();

  //MbdCalib(MbdCalib &other) = delete;
  //void operator=(const MbdCalib &) = delete;

  virtual ~MbdCalib() {}

  float get_qgain(const int ipmt) const { return _qfit_mpv[ipmt]; }
  float get_tq0(const int ipmt) const { return _tqfit_t0mean[ipmt]; }
  int   get_sampmax(const int ifeech) const { return _sampmax[ifeech]; }

  int  Download_Gains(const std::string& dbfile);
  int  Download_TQT0(const std::string& dbfile);
  int  Download_SampMax(const std::string& dbfile);

  int  Download_All();

  int  StoreInDatabase();

  //void Dump_to_file(const std::string& what = "ALL");

  void Reset();
  //void Print(Option_t* option) const;

private:

  CDBInterface *_cdb {nullptr};
  recoConsts   *_rc {nullptr};

  int          _status {0};
  //int          _run_number {0};
  //uint64_t     _timestamp {0};
  std::string  _dbfilename;

  // Assumes Landau fit
  float _qfit_integ[MbdDefs::BBC_N_PMT];
  float _qfit_mpv[MbdDefs::BBC_N_PMT];
  float _qfit_sigma[MbdDefs::BBC_N_PMT];
  float _qfit_integerr[MbdDefs::BBC_N_PMT];
  float _qfit_mpverr[MbdDefs::BBC_N_PMT];
  float _qfit_sigmaerr[MbdDefs::BBC_N_PMT];
  float _qfit_chi2ndf[MbdDefs::BBC_N_PMT];

  // T0 offsets, charge channels
  float _tqfit_t0mean[MbdDefs::BBC_N_PMT];
  float _tqfit_t0meanerr[MbdDefs::BBC_N_PMT];
  float _tqfit_t0sigma[MbdDefs::BBC_N_PMT];
  float _tqfit_t0sigmaerr[MbdDefs::BBC_N_PMT];

  // Peak of waveform
  int   _sampmax[MbdDefs::BBC_N_FEECH];
};


#endif	// __MBDCALIB_H__

