#ifndef EVENTPLANEINFO_EVENTPLANERECOV2_H
#define EVENTPLANEINFO_EVENTPLANERECOV2_H

#include <fun4all/SubsysReco.h>


#include <string>
#include <array>
#include <memory>

class CDBTTree;
class PHCompositeNode;

class EventPlaneRecov2 : public SubsysReco
{
 public:

  explicit EventPlaneRecov2(const std::string &name = "EventPlaneRecov2");
  ~EventPlaneRecov2() override = default;

  // Explicitly disable copying and moving
  EventPlaneRecov2(const EventPlaneRecov2&) = delete;
  EventPlaneRecov2& operator=(const EventPlaneRecov2&) = delete;
  EventPlaneRecov2(EventPlaneRecov2&&) = delete;
  EventPlaneRecov2& operator=(EventPlaneRecov2&&) = delete;

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode *topNode) override;


  void set_inputNode(const std::string &inputNode)
  {
    m_inputNode = inputNode;
  }

  void set_directURL_EventPlaneCalib(const std::string &directURL_EventPlaneCalib)
  {
    m_directURL_EventPlaneCalib = directURL_EventPlaneCalib;
  }

  void set_doAbortNoEventPlaneCalib(bool status = true)
  {
    m_doAbortNoEventPlaneCalib = status;
  }

  void set_sepd_min_channel_charge(double sepd_min_channel_charge)
  {
    m_sepd_min_channel_charge = sepd_min_channel_charge;
  }

 private:

 static int CreateNodes(PHCompositeNode *topNode);

 std::array<std::array<double, 2>, 2> calculate_flattening_matrix(double xx, double yy, double xy, int n, int cent_bin, const std::string& det_label);
 void LoadCalib();

 void print_correction_data();
 void print_QVectors();

 int process_centrality(PHCompositeNode *topNode);
 int process_sEPD(PHCompositeNode *topNode);
 void correct_QVecs();

 int FillNode(PHCompositeNode *topNode);

 std::string m_directURL_EventPlaneCalib;
 bool m_doAbortNoEventPlaneCalib{false};
 bool m_doNotCalib{false};
 bool m_doNotCalibEvent{false};

 double m_cent{0.0};
 double m_globalEvent{0};
 double m_sepd_min_channel_charge{0.2};

 std::string m_calibName{"SEPD_EventPlaneCalib"};
 std::string m_inputNode{"TOWERINFO_CALIB_SEPD"};

  CDBTTree *m_cdbttree {nullptr};

 enum class Subdetector
 {
   S,
   N,
   NS
 };

 struct QVec
 {
    double x{0.0};
    double y{0.0};
 };

 struct CorrectionData
 {
    // Averages of Qx, Qy, Qx^2, Qy^2, Qxy
    QVec avg_Q{};
    double avg_Q_xx{0.0};
    double avg_Q_yy{0.0};
    double avg_Q_xy{0.0};

    // Correction matrix
    std::array<std::array<double, 2>, 2> X_matrix{};
 };

 static constexpr size_t m_cent_bins {80};
 static constexpr std::array<int, 3> m_harmonics = {2, 3, 4};

 // Holds all correction data
 // key: [Harmonic][Cent][Subdetector]
 // Harmonics {2,3,4} -> 3 elements
 // Subdetectors {S,N,NS} -> 3 elements
 std::array<std::array<std::array<CorrectionData, 3>, m_cent_bins>, m_harmonics.size()> m_correction_data;

 // sEPD Q Vectors
 // key: [Harmonic][Subdetector]
 // Subdetectors {S,N,NS} -> 3 elements
 std::array<std::array<QVec, 3>, m_harmonics.size()> m_Q_raw{};
 std::array<std::array<QVec, 3>, m_harmonics.size()> m_Q_recentered{};
 std::array<std::array<QVec, 3>, m_harmonics.size()> m_Q_flat{};
};
#endif
