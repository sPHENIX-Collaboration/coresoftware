#include "Rossegger.h"

#include <TVector3.h>

#include <cmath>   // for NAN, abs
#include <string>  // for string

class AnalyticFieldModel;
class ChargeMapReader;
class TH3;
class TTree;

template <class T>
class MultiArray;

class AnnularFieldSim
{
 private:
  //units:
  const float cm = 1;        //centimeters -- if you change this, check that all the loading functions are properly agnostic.
  const float m = 100 * cm;  //meters.
  const float mm = cm / 10;
  const float um = mm / 1e3;

  const float C = 1;  //Coulombs
  const float nC = C / 1e9;
  const float fC = C / 1e15;

  const float s = 1;  //seconds
  const float us = s / 1e6;
  const float ns = s / 1e9;

  const float V = 1;  //volts

  const float Tesla = V * s / m / m;  //Tesla=Vs/m^2
  const float kGauss = Tesla / 10;    //kGauss

  const float eps0 = 8.854e-12 * (C / V) / m;    //Farads(=Coulombs/Volts) per meter
  const float epsinv = 1 / eps0;                 //Vcm/C
  const float k_perm = 1 / (4 * 3.1416 * eps0);  //implied units of V*cm/C because we're doing unitful work here.

 public:
  enum BoundsCase
  {
    InBounds,
    OnHighEdge,
    OnLowEdge,
    OutOfBounds
  };  //note that 'OnLowEdge' is qualitatively different from 'OnHighEdge'.  Low means there is a non-zero distance between the point and the edge of the bin.  High applies even if that distance is exactly zero.
  enum LookupCase
  {
    Full3D,
    HybridRes,
    PhiSlice,
    Analytic,
    NoLookup
  };
  //Full3D = uses (nr x nphi x nz)^2 lookup table
  //Hybrid = uses (nr x nphi x nz) x (nr_local x nphi_local x nz_local) + (nr_low x nphi_low x nz_low)^2 set of tables
  //PhiSlice = uses (nr x 1 x nz) x (nr x nphi x nz) lookup table exploiting phi symmetry.
  //Analytic = doesn't use lookup tables -- no memory footprint, uses analytic E field at center of each bin.
  //    Note that this is not the same as analytic propagation, which checks the analytic field integrals in each step.
  //NoLookup = Don't build any structures -- effectively ignores any calculated spacecharge field
  enum ChargeCase
  {
    FromFile,
    AnalyticSpacecharge,
    NoSpacecharge
  };  //load from file, load from AnalyticFieldModel, or set to zero.
  //note that if we set to Zero, we skip the lookup step.

  //debug items
  //
  int debug_printActionEveryN;
  int debug_npercent;
  int debug_printCounter;
  TVector3 debug_distortionScale;

  AnalyticFieldModel *aliceModel = nullptr;

  //the other half of the detector:
  AnnularFieldSim *twin = nullptr;
  bool hasTwin = false;

  //constants of motion, dimensions, etc:
  //
  TVector3 zero_vector;  //a shorthand way to return a vectorial zero.
  //static constexpr float k=8.987e13;//=1/(4*pi*eps0) in N*cm^2/C^2 in a vacuum. N*cm^2/C units, so that we supply space charge in coulomb units.
  //static constexpr float k_perm=8.987e11;//=1/(4*pi*eps0) in (V*cm)/C in a vacuum. so that we supply space charge in Coulombs, distance in cm, and fields in V/cm

  //gas constants:
  double vdrift = NAN;  //gas drift speed in cm/s
  double langevin_T1 = NAN;
  double langevin_T2 = NAN;       //gas tensor drift terms.
  double omegatau_nominal = NAN;  //nominal omegatau value, derived from vdrift and field strengths.
  //double vprime; //first derivative of drift velocity at specific E
  //double vprime2; //second derivative of drift velocity at specific E

  //field constants:
  std::string fieldstring;
  std::string Bfieldname;
  std::string Efieldname;
  //  char fieldstring[300],Bfieldname[100],Efieldname[100];
  std::string chargesourcename;
  char chargestring[300] = {0};  //, chargefilename[100];
  float Enominal = NAN;          //magnitude of the nominal field on which drift speed is based, in V/cm.
  float Bnominal;                //magnitude of the nominal magnetic field on which drift speed is based, in Tesla.

  //physical dimensions
  float phispan;     //angular span of the area in the phi direction, since TVector3 is too smart.
  float rmin, rmax;  //inner and outer radii of the annulus
  float zmin, zmax;  //lower and upper edges of the coordinate system in z (not fully implemented yet)
  //float phimin, phimax;//not implemented at all yet.
  TVector3 dim;       //dimensions of simulated region, in cm
  Rossegger *green;   //stand-alone class to compute greens functions.
  float green_shift;  //how far to offset our position in z when querying our green's functions.

  //variables related to the whole-volume tiling:
  //
  int nr, nphi, nz;       //number of fundamental bins (f-bins) in each direction = dimensions of 3D array covering entire volume
  TVector3 step;          //size of an f-bin in each direction
  LookupCase lookupCase;  //which lookup system to instantiate and use.
  ChargeCase chargeCase;  //which charge model to use
  int truncation_length;  //distance in cells (full 3D metric in units of bins)

  //variables related to the region of interest:
  //
  int rmin_roi, phimin_roi, zmin_roi;  //lower edge of our region of interest, measured in f-bins
  int rmax_roi, phimax_roi, zmax_roi;  //excluded upper edge of our region of interest, measured in f-bins
  int nr_roi, nphi_roi, nz_roi;        //dimensions of our roi in f-bins

  //variables related to the high-res behavior:
  //
  int nr_high = -1;
  int nphi_high = -1;
  int nz_high = -1;  //dimensions, in f-bins of neighborhood of a f-bin in which we calculate the field in full resolution

  //variables related to the low-res behavior:
  //
  int r_spacing = -1;
  int phi_spacing = -1;
  int z_spacing = -1;  //number of f-bins, in each direction, to gang together to make a single low-resolution bin (l-bin)
  int nr_low = -1;
  int nphi_low = -1;
  int nz_low = -1;  //dimensions, in l-bins, of the entire volume
  int rmin_roi_low = -1;
  int phimin_roi_low = -1;
  int zmin_roi_low = -1;  //lowest l-bin that is at least partly in our region of interest
  int rmax_roi_low = -1;
  int phimax_roi_low = -1;
  int zmax_roi_low = -1;  //excluded upper edge l-bin of our region of interest
  int nr_roi_low = -1;
  int nphi_roi_low = -1;
  int nz_roi_low = -1;  //dimensions of our roi in l-bins

  //3- and 6-dimensional arrays to handle bin and bin-to-bin data
  //
  MultiArray<TVector3> *Efield;             //total electric field in each f-bin in the roi for given configuration of charge AND external field.
  MultiArray<TVector3> *Epartial_highres;   //electric field in each f-bin in the roi from charge in a given f-bin or summed bin in the high res region.
  MultiArray<TVector3> *Epartial_lowres;    //electric field in each l-bin in the roi from charge in a given l-bin anywhere in the volume.
  MultiArray<TVector3> *Epartial;           //electric field for the old brute-force model.
  MultiArray<TVector3> *Epartial_phislice;  //electric field in a 2D phi-slice from the full 3D region.
  MultiArray<TVector3> *Eexternal;          //externally applied electric field in each f-bin in the roi
  MultiArray<TVector3> *Bfield;             //magnetic field in each f-bin in the roi

  ChargeMapReader *q;            // //class to read and report charge.
                                 //  MultiArray<double> *q;                    //space charge in each f-bin in the whole volume
  MultiArray<double> *q_local;   //temporary holder of space charge in each f-bin and summed bin of the high-res region.
  MultiArray<double> *q_lowres;  //space charge in each l-bin. = sums over sets of f-bins.

 public:
  //constructors with history for backwards compatibility
  AnnularFieldSim(float rmin, float rmax, float dz, int r, int phi, int z, float vdr);  //abbr. constructor with roi=full region
  AnnularFieldSim(float rin, float rout, float dz,
                  int r, int roi_r0, int roi_r1,
                  int phi, int roi_phi0, int roi_phi1,
                  int z, int roi_z0, int roi_z1,
                  float vdr, LookupCase in_lookupCase = PhiSlice);
  AnnularFieldSim(float in_innerRadius, float in_outerRadius, float in_outerZ,
                  int r, int roi_r0, int roi_r1, int in_rLowSpacing, int in_rHighSize,
                  int phi, int roi_phi0, int roi_phi1, int in_phiLowSpacing, int in_phiHighSize,
                  int z, int roi_z0, int roi_z1, int in_zLowSpacing, int in_zHighSize,
                  float vdr, LookupCase in_lookupCase);
  AnnularFieldSim(float in_innerRadius, float in_outerRadius, float in_outerZ,
                  int r, int roi_r0, int roi_r1, int in_rLowSpacing, int in_rHighSize,
                  int phi, int roi_phi0, int roi_phi1, int in_phiLowSpacing, int in_phiHighSize,
                  int z, int roi_z0, int roi_z1, int in_zLowSpacing, int in_zHighSize,
                  float vdr, LookupCase in_lookupCase, ChargeCase in_chargeCase);
  //! delete copy ctor and assignment opertor (cppcheck)
  explicit AnnularFieldSim(const AnnularFieldSim &) = delete;
  AnnularFieldSim &operator=(const AnnularFieldSim &) = delete;

  //debug functions:
  void UpdateEveryN(int n)
  {
    debug_npercent = n;
    debug_printActionEveryN = 0;
    return;
  };
  bool debugFlag()
  {
    if (debug_printActionEveryN > 0 && debug_printCounter++ >= debug_printActionEveryN)
    {
      debug_printCounter = 0;
      return true;
    }
    return false;
  };
  void SetDistortionScaleRPZ(float a, float b, float c)
  {
    debug_distortionScale.SetXYZ(a, b, c);
    return;
  };
  void SetTruncationDistance(int x)
  {
    truncation_length = x;
    return;
  }

  //getters for internal states:
  const char *GetLookupString();
  const char *GetGasString();
  const char *GetFieldString();
  const char *GetChargeString() { return chargestring; };
  float GetNominalB() { return Bnominal; };
  float GetNominalE() { return Enominal; };
  float GetChargeAt(TVector3 pos);
  TVector3 GetFieldAt(TVector3 pos);
  TVector3 GetBFieldAt(TVector3 pos);
  TVector3 GetFieldStep() { return step; };
  int GetFieldStepsR() { return nr_roi; };
  int GetFieldStepsPhi() { return nphi_roi; };
  int GetFieldStepsZ() { return nz_roi; };
  TVector3 GetInnerEdge() { return TVector3(rmin, 0, zmin); };
  TVector3 GetOuterEdge() { return TVector3(rmax, 0, zmax); };

  //file-writing functions for complex mapping questions:
  void GenerateDistortionMaps(const char *filebase, int r_subsamples = 1, int p_subsamples = 1, int z_subsamples = 1, int z_substeps = 1, bool andCartesian = false);
  void GenerateSeparateDistortionMaps(const char *filebase, int r_subsamples = 1, int p_subsamples = 1, int z_subsamples = 1, int z_substeps = 1, bool andCartesian = false);
  void PlotFieldSlices(const char *filebase, TVector3 pos, char which = 'E');

  void load_spacecharge(const std::string &filename, const std::string &histname, float zoffset = 0, float chargescale = 1, float cmscale = 1, bool isChargeDensity = true);
  void load_spacecharge(TH3 *hist, float zoffset, float chargescale, float cmscale, bool isChargeDensity, const char *inputchargestring = "");

  void load_and_resample_spacecharge(int new_nphi, int new_nr, int new_nz, const std::string &filename, const std::string &histname, float zoffset, float chargescale, float cmscale, bool isChargeDensity);

  void load_and_resample_spacecharge(int new_nphi, int new_nr, int new_nz, TH3 *hist, float zoffset, float chargescale, float cmscale, bool isChargeDensity);
  void save_spacecharge(const std::string &filename);
  void load_analytic_spacecharge(float scalefactor);
  void add_testcharge(float r, float phi, float z, float coulombs);

  void setNominalB(float x)
  {
    Bnominal = x;
    UpdateOmegaTau();
    return;
  };
  void setNominalE(float x)
  {
    Enominal = x;
    UpdateOmegaTau();
    return;
  };
  void setFlatFields(float B, float E);
  void loadEfield(const std::string &filename, const std::string &treename, int zsign = 1);
  void loadBfield(const std::string &filename, const std::string &treename);
  void load3dBfield(const std::string &filename, const std::string &treename, int zsign = 1, float scale = 1.0);

  void loadField(MultiArray<TVector3> **field, TTree *source, float *rptr, float *phiptr, float *zptr, float *frptr, float *fphiptr, float *fzptr, float fieldunit, int zsign);

  void load_rossegger(double epsilon = 1E-4)
  {
    green = new Rossegger(rmin, rmax, zmax, epsilon);
    return;
  };
  void borrow_rossegger(Rossegger *ross, float zshift)
  {
    green = ross;
    green_shift = zshift;
    return;
  };  //get an already-existing rossegger table instead of loading it ourselves.
  void borrow_epartial_from(AnnularFieldSim *sim, float zshift)
  {
    Epartial_phislice = sim->Epartial_phislice;
    green_shift = zshift;
    return;
  };  //get an already-existing rossegger table instead of loading it ourselves.
  void set_twin(AnnularFieldSim *sim)
  {
    twin = sim;
    hasTwin = true;
    return;
  };  //define a twin to handle the negative-z drifting.  If asked to drift something out of range in z, if the twin flag is set we will ask the twin to do the drifting.  Note that the twin does not get linked back in to this side.  It is only the follower.

  TVector3 calc_unit_field(TVector3 at, TVector3 from);
  TVector3 analyticFieldIntegral(float zdest, TVector3 start) { return analyticFieldIntegral(zdest, start, Efield); };

  TVector3 analyticFieldIntegral(float zdest, TVector3 start, MultiArray<TVector3> *field);
  TVector3 interpolatedFieldIntegral(float zdest, TVector3 start) { return interpolatedFieldIntegral(zdest, start, Efield); };
  TVector3 interpolatedFieldIntegral(float zdest, TVector3 start, MultiArray<TVector3> *field);
  double FilterPhiPos(double phi);         //puts phi in 0<phi<2pi
  int FilterPhiIndex(int phi, int range);  //puts phi in bin range 0<phi<range.  defaults to using nphi for range.

  TVector3 GetCellCenter(int r, int phi, int z);
  TVector3 GetRoiCellCenter(int r, int phi, int z);
  TVector3 GetGroupCellCenter(int r0, int r1, int phi0, int phi1, int z0, int z1);
  TVector3 GetWeightedCellCenter(int r, int phi, int z);
  TVector3 fieldIntegral(float zdest, TVector3 start, MultiArray<TVector3> *field);
  void populate_fieldmap();
  //now handled by setting 'analytic' lookup:  void populate_analytic_fieldmap();
  void populate_lookup();
  void populate_full3d_lookup();
  void populate_highres_lookup();
  void populate_lowres_lookup();
  void populate_phislice_lookup();

  void load_phislice_lookup(const char *sourcefile);
  void save_phislice_lookup(const char *destfile);

  TVector3 sum_field_at(int r, int phi, int z);
  TVector3 sum_full3d_field_at(int r, int phi, int z);
  TVector3 sum_local_field_at(int r, int phi, int z);
  TVector3 sum_nonlocal_field_at(int r, int phi, int z);
  TVector3 sum_phislice_field_at(int r, int phi, int z);
  TVector3 swimToInAnalyticSteps(float zdest, TVector3 start, int steps, int *goodToStep);
  TVector3 swimToInSteps(float zdest, TVector3 start, int steps, bool interpolate, int *goodToStep);
  TVector3 swimTo(float zdest, TVector3 start, bool interpolate = true, bool useAnalytic = false);
  TVector3 GetStepDistortion(float zdest, TVector3 start, bool interpolate = true, bool useAnalytic = false);
  TVector3 GetTotalDistortion(float zdest, TVector3 start, int nsteps, bool interpolate = true, int *goodToStep = 0);

 private:
  BoundsCase GetRindexAndCheckBounds(float pos, int *r);
  BoundsCase GetPhiIndexAndCheckBounds(float pos, int *phi);
  BoundsCase GetZindexAndCheckBounds(float pos, int *z);
  int GetRindex(float pos);
  int GetPhiIndex(float pos);
  int GetZindex(float pos);

  void UpdateOmegaTau()
  {
    omegatau_nominal = -Bnominal * vdrift / abs(Enominal);
    return;
  };  //various constants to match internal representation to the familiar formula.  Adding in these factors suggests I should switch to a unitful calculation throughout...
};
