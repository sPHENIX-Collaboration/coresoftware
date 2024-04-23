#ifndef MBD_MBDCALIB_H
#define MBD_MBDCALIB_H

#include "MbdDefs.h"
#include <mbd/MbdGeom.h>

#ifndef ONLINE
#include <fun4all/Fun4AllBase.h>
#include <phool/recoConsts.h>
#endif


#include <array>
#include <vector>
#include <string>
#include <memory>

class TTree;
class CDBInterface;

class MbdCalib 
{
 public:
  MbdCalib();

  // MbdCalib(MbdCalib &other) = delete;
  // void operator=(const MbdCalib &) = delete;

  virtual ~MbdCalib() {}

  float get_qgain(const int ipmt) const { return _qfit_mpv[ipmt]; }
  float get_tq0(const int ipmt) const { return _tqfit_t0mean[ipmt]; }
  float get_tt0(const int ipmt) const { return _ttfit_t0mean[ipmt]; }
  int get_sampmax(const int ifeech) const { return _sampmax[ifeech]; }
  float get_tcorr(const int ifeech, const int tdc) const {
    if (tdc<0)
    {
      std::cout << "bad tdc " << ifeech << " " << tdc << std::endl;
      return _tcorr_y_interp[ifeech][0];
    }
    if (tdc>=_tcorr_maxrange[ifeech])
    {
      std::cout << "bad tdc " << ifeech << " " << tdc << " " << _tcorr_maxrange[ifeech] << std::endl;
      return _tcorr_y_interp[ifeech][_tcorr_maxrange[ifeech]-1];
    }
    //std::cout << "aaa " << _tcorr_y_interp[ifeech].size() << std::endl;
    return _tcorr_y_interp[ifeech][tdc];
  }
  float get_scorr(const int ifeech, const int adc) const {
    if ( _scorr_y_interp[ifeech].size() == 0 ) return 0.; // return 0 if calib doesn't exist
    if (adc<0)
    {
      std::cout << "bad adc " << ifeech << " " << adc << std::endl;
      return _scorr_y_interp[ifeech][0];
    }
    if (adc>=_scorr_maxrange[ifeech])
    {
      std::cout << "bad adc " << ifeech << " " << adc << std::endl;
      return _scorr_y_interp[ifeech][_scorr_maxrange[ifeech]-1];
    }
    return _scorr_y_interp[ifeech][adc];
  }

  std::vector<float> get_shape(const int ifeech) const { return _shape_y[ifeech]; }
  std::vector<float> get_sherr(const int ifeech) const { return _sherr_yerr[ifeech]; }

  void set_sampmax(const int ifeech, const int val) { _sampmax[ifeech] = val; }

  int Download_Gains(const std::string& dbfile);
  int Download_TQT0(const std::string& dbfile);
  int Download_TTT0(const std::string& dbfile);
  int Download_SampMax(const std::string& dbfile);
  int Download_Shapes(const std::string& dbfile);
  int Download_TimeCorr(const std::string& dbfile);
  int Download_SlewCorr(const std::string& dbfile);
  int Download_All();

#ifndef ONLINE
  int Write_CDB_SampMax(const std::string& dbfile);
  int Write_CDB_TTT0(const std::string& dbfile);
  int Write_CDB_TQT0(const std::string& dbfile);
  int Write_CDB_Shapes(const std::string& dbfile);
  int Write_CDB_TimeCorr(const std::string& dbfile);
  int Write_CDB_SlewCorr(const std::string& dbfile);
  int Write_CDB_Gains(const std::string& dbfile);
  int Write_CDB_All();
#endif

  int Write_SampMax(const std::string& dbfile);
  int Write_TQT0(const std::string& dbfile);
  int Write_TTT0(const std::string& dbfile);

  void Reset_TQT0();
  void Reset_TTT0();
  void Reset_Gains();

  void Update_TQT0(const float dz); // update with new z-vertex
  void Update_TTT0(const float dz);

  // void Dump_to_file(const std::string& what = "ALL");

  void Reset();
  // void Print(Option_t* option) const;

  int Verbosity() { return _verbose; }
  void Verbosity(const int v) { _verbose = v; }

 private:
#ifndef ONLINE
  CDBInterface* _cdb{nullptr};
  recoConsts* _rc{nullptr};
#endif

  std::unique_ptr<MbdGeom> _mbdgeom{nullptr};

  int _status{0};
  int _verbose{0};
  // int          _run_number {0};
  // uint64_t     _timestamp {0};
  std::string _dbfilename;

  // Assumes Landau fit
  std::array<float, MbdDefs::MBD_N_PMT> _qfit_integ{};
  std::array<float, MbdDefs::MBD_N_PMT> _qfit_mpv{};
  std::array<float, MbdDefs::MBD_N_PMT> _qfit_sigma{};
  std::array<float, MbdDefs::MBD_N_PMT> _qfit_integerr{};
  std::array<float, MbdDefs::MBD_N_PMT> _qfit_mpverr{};
  std::array<float, MbdDefs::MBD_N_PMT> _qfit_sigmaerr{};
  std::array<float, MbdDefs::MBD_N_PMT> _qfit_chi2ndf{};

  // T0 offsets, time channels
  std::array<float, MbdDefs::MBD_N_PMT> _ttfit_t0mean{};
  std::array<float, MbdDefs::MBD_N_PMT> _ttfit_t0meanerr{};
  std::array<float, MbdDefs::MBD_N_PMT> _ttfit_t0sigma{};
  std::array<float, MbdDefs::MBD_N_PMT> _ttfit_t0sigmaerr{};

  // T0 offsets, charge channels
  std::array<float, MbdDefs::MBD_N_PMT> _tqfit_t0mean{};
  std::array<float, MbdDefs::MBD_N_PMT> _tqfit_t0meanerr{};
  std::array<float, MbdDefs::MBD_N_PMT> _tqfit_t0sigma{};
  std::array<float, MbdDefs::MBD_N_PMT> _tqfit_t0sigmaerr{};

  // Peak of waveform
  std::array<int, MbdDefs::MBD_N_FEECH> _sampmax{};

  // Waveform Template
  int do_templatefit{0};
  std::array<int, MbdDefs::MBD_N_FEECH>   _shape_npts{};      // num points in template
  std::array<float, MbdDefs::MBD_N_FEECH> _shape_minrange{};  // in template units (samples)
  std::array<float, MbdDefs::MBD_N_FEECH> _shape_maxrange{};  // in template units (samples)
  std::array<std::vector<float>, MbdDefs::MBD_N_FEECH> _shape_y{};

  std::array<int, MbdDefs::MBD_N_FEECH>   _sherr_npts{};      // num points in template
  std::array<float, MbdDefs::MBD_N_FEECH> _sherr_minrange{};  // in template units (samples)
  std::array<float, MbdDefs::MBD_N_FEECH> _sherr_maxrange{};  // in template units (samples)
  std::array<std::vector<float>, MbdDefs::MBD_N_FEECH> _sherr_yerr{};

  // Fine Timing Corrections
  std::array<int, MbdDefs::MBD_N_FEECH>   _tcorr_npts{};      // num points in timing corr
  std::array<float, MbdDefs::MBD_N_FEECH> _tcorr_minrange{};  // in units (delta-TDC)
  std::array<float, MbdDefs::MBD_N_FEECH> _tcorr_maxrange{};  // in units (detta-TDC)
  std::array<std::vector<float>, MbdDefs::MBD_N_FEECH> _tcorr_y{};
  std::array<std::vector<float>, MbdDefs::MBD_N_FEECH> _tcorr_y_interp{}; // interpolated tcorr

  // Slew Correction
  std::array<int, MbdDefs::MBD_N_FEECH>   _scorr_npts{};      // num points in slew corr
  std::array<float, MbdDefs::MBD_N_FEECH> _scorr_minrange{};  // in units (delta-TDC)
  std::array<float, MbdDefs::MBD_N_FEECH> _scorr_maxrange{};  // in units (detta-TDC)
  std::array<std::vector<float>, MbdDefs::MBD_N_FEECH> _scorr_y{};
  std::array<std::vector<float>, MbdDefs::MBD_N_FEECH> _scorr_y_interp{}; // interpolated scorr

};

#endif  // MBD_MBDCALIB_H
