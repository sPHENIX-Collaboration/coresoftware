#ifndef __BBCCALIB_H__
#define __BBCCALIB_H__

#include "BbcDefs.h"
#include <phool/recoConsts.h>

#include <cmath>
#include <cstdint>
#include <cstring>

class TTree;
class CDBInterface;

class BbcCalib
{
public:
  BbcCalib();

  //BbcCalib(BbcCalib &other) = delete;
  //void operator=(const BbcCalib &) = delete;

  virtual ~BbcCalib() {}

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
  float _qfit_integ[BbcDefs::BBC_N_PMT];
  float _qfit_mpv[BbcDefs::BBC_N_PMT];
  float _qfit_sigma[BbcDefs::BBC_N_PMT];
  float _qfit_integerr[BbcDefs::BBC_N_PMT];
  float _qfit_mpverr[BbcDefs::BBC_N_PMT];
  float _qfit_sigmaerr[BbcDefs::BBC_N_PMT];
  float _qfit_chi2ndf[BbcDefs::BBC_N_PMT];

  // T0 offsets, charge channels
  float _tqfit_t0mean[BbcDefs::BBC_N_PMT];
  float _tqfit_t0meanerr[BbcDefs::BBC_N_PMT];
  float _tqfit_t0sigma[BbcDefs::BBC_N_PMT];
  float _tqfit_t0sigmaerr[BbcDefs::BBC_N_PMT];

  // Peak of waveform
  int   _sampmax[BbcDefs::BBC_N_FEECH];
};


#endif	// __BBCCALIB_H__

