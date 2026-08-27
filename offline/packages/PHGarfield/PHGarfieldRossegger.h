// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHGARFIELD_PHGARFIELDROSSEGGER_H
#define PHGARFIELD_PHGARFIELDROSSEGGER_H

#include <fun4all/SubsysReco.h>

#include <array>
#include <string>
#include <utility>
#include <vector>

class PHCompositeNode;

class PHGarfieldRossegger : public SubsysReco
{
 public:
  explicit PHGarfieldRossegger(const std::string& name = "PHGarfieldRossegger");
  ~PHGarfieldRossegger() override = default;

  int Init(PHCompositeNode*) override;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  void setGeometryCm(double inner_radius, double outer_radius, double half_length);
  void setSourceRadiusCm(double min_radius, double max_radius);
  void setDensity(double reference_density_nC_per_m3, double k_eff, double radial_power_alpha);
  void setPhiModulation(double m1_amplitude, double m1_phase, double m12_amplitude, double m12_phase);
  void setSourceGrid(unsigned int nr, unsigned int nphi, unsigned int nz);
  void setObservationGrid(unsigned int nr, unsigned int nphi, unsigned int nz);
  void setModeTruncation(unsigned int m_phi_max, unsigned int n_radial_modes, unsigned int n_longitudinal_modes);
  void setAutoAxisymmetric(bool value) { m_autoAxisymmetric = value; }
  void setTpcSide(unsigned int side) { m_tpcSide = side; }

  void setUseRealTpcSourceGeometry(bool value) { m_useRealTpcSourceGeometry = value; }
  void setPadPlacementFile(const std::string& filename) { m_padPlacementFile = filename; }
  void setGainMapFile(const std::string& filename) { m_gainMapFile = filename; }
  void setGainHistograms(const std::string& side0_histogram, const std::string& side1_histogram);
  void setDivideChargeByGain(bool value) { m_divideChargeByGain = value; }
  void setNormalizeGainWeightedTotal(bool value) { m_normalizeGainWeightedTotal = value; }

  void setUseDensityMap(bool value) { m_useDensityMap = value; }
  void setDensityMapFile(const std::string& filename, const std::string& side0_histogram, const std::string& side1_histogram);
  void setNormalizeDensityMap(bool value) { m_normalizeDensityMap = value; }

  void setUseFrameChargeModel(bool value) { m_useFrameChargeModel = value; }
  void setFrameReferencePhi(double value) { m_frameReferencePhi = value; }
  void setFrameGeometryFile(const std::string& filename) { m_frameGeometryFile = filename; m_framePolygonsLoaded = false; }
  enum class FrameChargeWeighting
  {
    EqualChargePerPiece,
    ProportionalToArea
  };
  void setFrameBoundaryPotential(double value) { m_frameBoundaryPotential = value; }
  void setFrameChargeWeighting(FrameChargeWeighting mode) { m_frameChargeWeighting = mode; }
  void setFrameEzScale(double value) { m_frameEzScale = value; }

  // Split observation radial bins across condor jobs. Histograms keep full
  // binning, so PART files can be merged with hadd.
  void setRadialJob(unsigned int job_index, unsigned int n_jobs);

  void setGarfieldOutputFile(const std::string& filename) { m_garfieldOutputFile = filename; }
  void setField3DOutputFile(const std::string& filename) { m_field3DOutputFile = filename; }
  void setPHGarfieldField3DOutputFiles(const std::string& side0_filename, const std::string& side1_filename);
  void setWriteField3D(bool value) { m_writeField3D = value; }
  void setWritePHGarfieldField3D(bool value) { m_writePHGarfieldField3D = value; }
  void setVerifyOutput(bool value) { m_verifyOutput = value; }

 private:
  struct SourceGrid
  {
    std::vector<double> r_edges_m;
    std::vector<double> r_centers_m;
    std::vector<int> module_index;
    std::vector<int> layer_index;
    std::vector<bool> is_antenna;
  };

  struct PadInterval
  {
    double start_edge{0.0};
    double angular_width{0.0};
    int direction{1};
  };

  using PadGeometry = std::array<std::array<std::array<PadInterval, 12>, 2>, 3>;
  using GainMap = std::array<std::array<std::array<double, 48>, 12>, 2>;

  struct Point2D
  {
    double x_cm{0.0};
    double y_cm{0.0};
  };

  struct FramePolygon
  {
    std::vector<Point2D> inner;
    std::vector<Point2D> outer;
  };

  using FramePolygons = std::array<FramePolygon, 3>;

  struct FrameBoundaryPattern
  {
    std::vector<double> geometry_fraction;
    std::vector<double> weight;
    std::vector<double> boundary_potential;
  };

  int calculate();
  bool validateConfig() const;
  unsigned int effectivePhiMax() const;
  std::pair<unsigned int, unsigned int> radialRange(unsigned int nr_obs) const;

  double radialBoundaryFunction(double kval, unsigned int mode_m) const;
  double radialBasis(double kval, unsigned int mode_m, double radius_m) const;
  double radialBasisDerivative(double kval, unsigned int mode_m, double radius_m) const;
  std::vector<double> findRadialRoots(unsigned int mode_m, unsigned int n_roots) const;
  std::vector<double> legendreRootsAndWeights(unsigned int n_points, std::vector<double>& weights) const;
  SourceGrid buildSourceGrid() const;
  PadGeometry parsePadPlacementFile() const;
  GainMap readGainMap() const;
  std::vector<double> makeChargeDensity(const SourceGrid& source_grid,
                                        const std::vector<double>& phi_source_edges,
                                        const std::vector<double>& z_source_edges_m,
                                        unsigned int side) const;
  void loadFramePolygons() const;
  static bool pointInPolygon(const std::vector<Point2D>& polygon, double x_cm, double y_cm);
  static double polygonAreaCm2(const std::vector<Point2D>& polygon);
  bool pointInFrameGeometry(const FramePolygon& frame, double r_cm, double phi_rel) const;
  double frameAreaM2(const FramePolygon& frame) const;
  FrameBoundaryPattern makeFrameBoundaryPattern(const SourceGrid& source_grid,
                                                const std::vector<double>& phi_source_edges,
                                                unsigned int side) const;

  void writeGarfieldRootFile(const std::vector<double>& r_edges_m,
                             const std::vector<double>& phi_edges,
                             const std::vector<double>& z_edges_m,
                             const std::vector<double>& er,
                             const std::vector<double>& ephi,
                             const std::vector<double>& ez,
                             unsigned int r_begin,
                             unsigned int r_end) const;

  void writeField3DRootFile(const std::vector<double>& r_source_edges_m,
                            const std::vector<double>& r_obs_edges_m,
                            const std::vector<double>& phi_source_edges,
                            const std::vector<double>& phi_obs_edges,
                            const std::vector<double>& z_source_edges_m,
                            const std::vector<double>& z_obs_edges_m,
                            const std::vector<double>& rho,
                            const std::vector<double>& potential,
                            const std::vector<double>& er,
                            const std::vector<double>& ephi,
                            const std::vector<double>& ez,
                            unsigned int r_begin,
                            unsigned int r_end) const;

  void writeFrameBoundaryField3DRootFile(const std::vector<double>& r_source_edges_m,
                                         const std::vector<double>& r_obs_edges_m,
                                         const std::vector<double>& phi_source_edges,
                                         const std::vector<double>& phi_obs_edges,
                                         const std::vector<double>& z_obs_edges_m,
                                         const FrameBoundaryPattern& frame_pattern,
                                         const std::vector<double>& potential,
                                         const std::vector<double>& er,
                                         const std::vector<double>& ephi,
                                         const std::vector<double>& ez,
                                         unsigned int r_begin,
                                         unsigned int r_end) const;

  void writePHGarfieldField3DRootFile(const std::vector<double>& r_obs_edges_m,
                                      const std::vector<double>& phi_obs_edges,
                                      const std::vector<double>& z_obs_edges_m,
                                      const std::vector<double>& potential,
				      const std::vector<double>& er,
                                      const std::vector<double>& ephi,
                                      const std::vector<double>& ez,
                                      unsigned int side,
                                      unsigned int r_begin,
                                      unsigned int r_end) const;

  bool verifyOutput() const;

  double m_aCm{21.6};
  double m_bCm{76.4};
  double m_lCm{102.325};
  double m_sourceRMinCm{22.8};
  double m_sourceRMaxCm{75.43};

  double m_rhoReferenceNCPerM3{20.0};
  double m_kEff{1.0};
  double m_radialPowerAlpha{1.4};
  double m_m1Amplitude{0.0};
  double m_m1Phase{0.0};
  double m_m12Amplitude{0.0};
  double m_m12Phase{0.0};

  unsigned int m_tpcSide{0};
  unsigned int m_nrSource{18};
  unsigned int m_nphiSource{24};
  unsigned int m_nzSource{20};
  unsigned int m_nrObs{24};
  unsigned int m_nphiObs{24};
  unsigned int m_nzObs{28};
  unsigned int m_mPhiMax{4};
  unsigned int m_nRadialModes{18};
  unsigned int m_nLongitudinalModes{22};
  bool m_autoAxisymmetric{true};
  bool m_useRealTpcSourceGeometry{false};
  bool m_useFrameChargeModel{false};
  double m_frameReferencePhi{0.0};
  std::string m_frameGeometryFile{"tpc_frame_geometry.csv"};
  mutable FramePolygons m_framePolygons;
  mutable bool m_framePolygonsLoaded{false};
  double m_frameBoundaryPotential{1.0};
  FrameChargeWeighting m_frameChargeWeighting{FrameChargeWeighting::ProportionalToArea};
  double m_frameEzScale{1.0};
  std::string m_padPlacementFile{"input/TPC_pad_placement.txt"};
  std::string m_gainMapFile{"input/layer_gain_79513_Mariia_side01.root"};
  std::array<std::string, 2> m_gainHistograms{{"hGainMap_side0_South", "hGainMap_side1_North"}};
  bool m_divideChargeByGain{true};
  bool m_normalizeGainWeightedTotal{true};

  bool m_useDensityMap{false};
  bool m_normalizeDensityMap{true};
  std::string m_densityMapFile;
  std::array<std::string, 2> m_densityMapHistograms{{"h_ibf_final_rphi_side0", "h_ibf_final_rphi_side1"}};

  unsigned int m_jobIndex{0};
  unsigned int m_nJobs{1};

  std::string m_garfieldOutputFile{"sphenix_rossegger_garfield_field_1p4.root"};
  std::string m_field3DOutputFile{"sphenix_rossegger_fields_3d.root"};
  std::array<std::string, 2> m_phGarfieldField3DOutputFiles{{"sphenix_3d_ibf_field_side0_South_v1.root", "sphenix_3d_ibf_field_side1_North_v1.root"}};
  bool m_writeField3D{false};
  bool m_writePHGarfieldField3D{false};
  bool m_verifyOutput{true};
  bool m_done{false};
};

#endif
