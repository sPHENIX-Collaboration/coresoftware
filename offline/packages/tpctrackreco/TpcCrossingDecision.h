#ifndef TPCTRACKRECO_TPCCROSSINGDECISION_H
#define TPCTRACKRECO_TPCCROSSINGDECISION_H


#include <phool/PHObject.h>

#include <iostream>
#include <limits>
#include <string>
#include <vector>

enum class TpcCrossingStatus : unsigned char
{
  Unknown = 0,
  NoInttCrossings = 1,
  NoAllowedCrossing = 2,
  GarfieldInvalid = 3,
  OutsideTpcVolume = 4,
  WrongTpcSide = 5,
  FitFailed = 6,
  SelectedByContainment = 7,
  SelectedByVertex = 8,
  AmbiguousWithoutVertex = 9,
  AmbiguousWithVertex = 10,
  VertexIncompatible = 11,
  NoValidCrossing = 12,
  SelectedByVertexAmbiguous = 13,
  SelectedByVertexLoose = 14,
  SelectedByContainmentAmbiguous = 15
};

enum class TpcCrossingCandidateStage : unsigned char
{
  AvailableFromIntt = 0,
  PassedTimeWindow = 1,
  GarfieldValid = 2,
  InsideTpc = 3,
  CorrectSide = 4,
  FitValid = 5,
  HasVertex = 6,
  VertexCompatible = 7,
  FinalRanking = 8,
  Selected = 9
};

enum TpcCrossingCandidateQABits : unsigned int
{
  FromInttCrossing = 1U << 0,
  PassedTimeWindow = 1U << 1,
  MinEndpointGarfieldOK = 1U << 2,
  MaxEndpointGarfieldOK = 1U << 3,
  EndpointsInsideTPC = 1U << 4,
  EndpointsCorrectSide = 1U << 5,
  FitSuccessful = 1U << 6,
  HasSiliconVertex = 1U << 7,
  PassesVertexDz = 1U << 8,
  IsBestVertexCandidate = 1U << 9,
  IsSelected = 1U << 10
};

class TpcCrossingCandidate
{
 public:
  static float nan() { return std::numeric_limits<float>::quiet_NaN(); }

  short crossing {0};
  bool is_available_from_intt {false};
  bool passes_time_window {false};
  bool was_tested {false};
  bool min_endpoint_garfield_valid {false};
  bool max_endpoint_garfield_valid {false};
  bool endpoints_inside_tpc {false};
  bool endpoints_on_correct_side {false};
  bool fit_valid {false};
  bool tpc_valid {false};
  bool has_silicon_vertex {false};
  bool vertex_compatible {false};
  bool is_selected {false};
  unsigned char confidence_tier {std::numeric_limits<unsigned char>::max()};
  float confidence_score {nan()};
  unsigned char rejection_status {static_cast<unsigned char>(TpcCrossingStatus::Unknown)};

  float min_tbin_time_crossing0_ns {nan()};
  float max_tbin_time_crossing0_ns {nan()};
  float candidate_min_time_ns {nan()};
  float candidate_max_time_ns {nan()};
  float max_lookup_time_ns {nan()};
  float min_time_margin_ns {nan()};
  float max_time_margin_ns {nan()};

  unsigned int min_time_layer {0};
  unsigned int min_time_side {0};
  unsigned int min_time_pad {0};
  unsigned int min_time_tbin {0};
  float min_time_x {nan()};
  float min_time_y {nan()};
  float min_time_z {nan()};
  float min_time_r {nan()};
  float min_time_phi {nan()};

  unsigned int max_time_layer {0};
  unsigned int max_time_side {0};
  unsigned int max_time_pad {0};
  unsigned int max_time_tbin {0};
  float max_time_x {nan()};
  float max_time_y {nan()};
  float max_time_z {nan()};
  float max_time_r {nan()};
  float max_time_phi {nan()};

  unsigned int inner_layer {0};
  unsigned int inner_tbin {0};
  float inner_x {nan()};
  float inner_y {nan()};
  float inner_z {nan()};

  unsigned int outer_layer {0};
  unsigned int outer_tbin {0};
  float outer_x {nan()};
  float outer_y {nan()};
  float outer_z {nan()};

  float min_time_distance_to_central_membrane {nan()};
  float max_time_distance_to_central_membrane {nan()};
  float min_time_distance_to_padplane {nan()};
  float max_time_distance_to_padplane {nan()};
  bool min_time_inside_tpc {false};
  bool max_time_inside_tpc {false};
  bool min_time_correct_side {false};
  bool max_time_correct_side {false};

  std::vector<unsigned int> fit_point_layer;
  std::vector<unsigned int> fit_point_pad;
  std::vector<unsigned int> fit_point_tbin;
  std::vector<float> fit_point_x;
  std::vector<float> fit_point_y;
  std::vector<float> fit_point_z;
  std::vector<float> fit_point_r;
  std::vector<float> fit_point_s;

  unsigned int n_fit_points {0};
  float z_vs_s_slope {nan()};
  float z_vs_s_intercept {nan()};
  float z_fit_chi2 {nan()};
  int z_fit_ndf {-1};
  float s_at_transverse_pca {nan()};
  float minimum_transverse_radius {nan()};
  float tpc_z_at_transverse_pca {nan()};
  float tpc_z_at_r0 {nan()};

  unsigned int n_vertices_at_crossing {0};
  std::vector<unsigned int> candidate_vertex_ids;
  std::vector<float> candidate_vertex_z;
  std::vector<float> candidate_vertex_sigma_z;
  std::vector<unsigned int> candidate_vertex_ntracks;
  std::vector<float> candidate_vertex_delta_z;
  std::vector<float> candidate_vertex_pull_z;

  unsigned int closest_vertex_id {std::numeric_limits<unsigned int>::max()};
  float closest_vertex_z {nan()};
  float closest_vertex_sigma_z {nan()};
  float closest_vertex_delta_z {nan()};
  float closest_vertex_abs_delta_z {nan()};
  float closest_vertex_pull_z {nan()};

  int candidate_rank_by_abs_delta_z {-1};
  int candidate_rank_by_abs_collision_z {-1};
  unsigned int candidate_qa_bits {0};
  unsigned char first_failed_stage {static_cast<unsigned char>(TpcCrossingCandidateStage::AvailableFromIntt)};

 private:
  ClassDefNV(TpcCrossingCandidate, 3)
};

class TpcCrossingDecision : public PHObject
{
 public:
  TpcCrossingDecision() = default;
  ~TpcCrossingDecision() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "TpcCrossingDecision base class" << std::endl;
  }
  void Reset() override {}
  int isValid() const override { return 0; }

  virtual unsigned int get_assembled_track_id() const { return 0; }
  virtual void set_assembled_track_id(unsigned int) {}

  virtual short get_selected_crossing() const { return 0; }
  virtual void set_selected_crossing(short) {}

  virtual unsigned int get_silicon_vertex_id() const { return 0; }
  virtual void set_silicon_vertex_id(unsigned int) {}

  virtual float get_tpc_z0() const { return 0.0F; }
  virtual void set_tpc_z0(float) {}

  virtual float get_silicon_vertex_z() const { return 0.0F; }
  virtual void set_silicon_vertex_z(float) {}

  virtual float get_delta_z() const { return 0.0F; }
  virtual void set_delta_z(float) {}

  virtual float get_best_abs_delta_z() const { return 0.0F; }
  virtual void set_best_abs_delta_z(float) {}

  virtual float get_second_best_abs_delta_z() const { return 0.0F; }
  virtual void set_second_best_abs_delta_z(float) {}

  virtual unsigned char get_selected_tier() const { return std::numeric_limits<unsigned char>::max(); }
  virtual void set_selected_tier(unsigned char) {}

  virtual float get_selected_score() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void set_selected_score(float) {}

  virtual unsigned short get_number_of_available_crossings() const { return 0; }
  virtual void set_number_of_available_crossings(unsigned short) {}

  virtual unsigned short get_number_of_allowed_crossings() const { return 0; }
  virtual void set_number_of_allowed_crossings(unsigned short) {}

  virtual unsigned short get_number_of_tested_crossings() const { return 0; }
  virtual void set_number_of_tested_crossings(unsigned short) {}

  virtual unsigned short get_number_of_tpc_valid_crossings() const { return 0; }
  virtual void set_number_of_tpc_valid_crossings(unsigned short) {}

  virtual unsigned short get_number_of_vertex_compatible_crossings() const { return 0; }
  virtual void set_number_of_vertex_compatible_crossings(unsigned short) {}

  virtual unsigned int get_number_of_candidates() const { return 0; }
  virtual const TpcCrossingCandidate* get_candidate(unsigned int) const { return nullptr; }
  virtual void add_candidate(const TpcCrossingCandidate&) {}
  virtual void clear_candidates() {}

  virtual unsigned char get_status() const { return static_cast<unsigned char>(TpcCrossingStatus::Unknown); }
  virtual void set_status(unsigned char) {}
  virtual void set_status(TpcCrossingStatus status) { set_status(static_cast<unsigned char>(status)); }

 private:
  ClassDefOverride(TpcCrossingDecision, 0)
};

#endif  // TPCTRACKRECO_TPCCROSSINGDECISION_H
