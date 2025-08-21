#ifndef FLOWAFTERBURNER_AFTERBURNERALGO_H
#define FLOWAFTERBURNER_AFTERBURNERALGO_H

#include <string>
#include <iostream>

namespace CLHEP
{
  class HepRandomEngine;
}

class AfterburnerAlgo
{
 public:
    //! flowAfterburner algorithms
    //! minbias_algorithm: standard flowAfterburner algorithm
    //! minbias_v2_algorithm: flowAfterburner algorithm with v2 only
    //! custom_algorithm: user defined flowAfterburner algorithm
    enum flowAfterburnerAlgorithm
    {
        minbias_algorithm,
        minbias_v2_algorithm,
        custom_algorithm
    }; // flowAfterburner algorithms
    
    AfterburnerAlgo() = default;
    explicit AfterburnerAlgo(flowAfterburnerAlgorithm algorithm = minbias_algorithm);
    
    ~AfterburnerAlgo() = default;

    void print( std::ostream &os = std::cout) const; // debugging output

    // set by an event
    void set_impact_parameter(double b)
    {
      m_impact_parameter = b;
    }


    // to scale the calculated vn values by a constant factor
    // only 1 harmonic ( set to 0 to disable ) or all harmonics
    void set_single_scale_N( const unsigned int n, const float scale ); 
    void set_scale_all( const float scale ); // set scale for all harmonics
    
    void enable_fluctuations(bool enable = true) { _do_fluctuations = enable; }

    void calc_flow(double eta, double pt, CLHEP::HepRandomEngine* engine = nullptr);
    void flucatate( CLHEP::HepRandomEngine* engine, float &v1, float &v2, float &v3, float &v4, float &v5, float &v6) const; // implements event-by-event fluctuations flow

    // getter
    float get_vn(unsigned int n) const;

    static std::string getAlgoName(flowAfterburnerAlgorithm algo);
    static flowAfterburnerAlgorithm getAlgoFromName(const std::string &name);

 private:

    flowAfterburnerAlgorithm m_algorithm = minbias_algorithm; // flowAfterburner algorithm
    double m_impact_parameter = 0.0; // impact parameter in fm
    float m_vn[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // v1 to v6
    float m_vn_scalefactors[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}; // scale factors for v1 to v6
    bool _do_fluctuations = false; // enable or disable event-by-event fluctuations flow

    static float calc_v2(double b, double eta, double pt);


};

#endif // FLOWAFTERBURNER_AFTERBURNERALGO_H