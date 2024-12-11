// Biquad filter module for executing various filtering procedures

#include <math.h>
#include "Biquad.h"

Biquad::Biquad() { // Define Biquad constructor function and its initial variables
  type = bq_type_lowpass; // Initialize type of filter
    a0 = 1.0; 
    a1 = a2 = b1 = b2 = 0.0;
    Fc = 0.50; // Initialize normalized cutoff frequency (Fc = cutoff/sample rate)
    Q = 0.707; // Initialize frequency range of filter 
    peakGain = 0.0; // Initialize gain of the filter at filtering point
    z1 = z2 = 0.0;
}

Biquad::Biquad(int ttype, double tFc, double tQ, double tpeakGainDB) { // Define Biquad constructor function with tunable parameters
  setBiquad(ttype, tFc, tQ, tpeakGainDB); // Inout parameters into setBiquad class
  z1 = z2 = 0.0;
}

Biquad::~Biquad() { // Define Biquad deconstructor function
}

void Biquad::setType(int ttype) { // Define setType function
type = ttype; // Set type to take value given in setBiquad
calcBiquad();
}

void Biquad::setQ(double tQ) { // Define setQ function
    Q = tQ;
    calcBiquad();
}

void Biquad::setFc(double tFc) { // Define setFc function
    Fc = tFc;
    calcBiquad();
}

void Biquad::setPeakGain(double tpeakGainDB) { // Define peakGainDB function
    peakGain = tpeakGainDB;
    calcBiquad();
}
    
void Biquad::setBiquad(int ttype, double tFc, double tQ, double tpeakGainDB) { // Define set Biquad function
    type = ttype;
    Q = tQ;
    Fc = tFc;
    setPeakGain(tpeakGainDB);
}

void Biquad::calcBiquad(void) { // Define calcBiquad function to set parameters of filter
    double norm;
    double V = pow(10, fabs(peakGain) / 20.0);
    double K = tan(M_PI * Fc);
    switch (this->type) { // Choose the 'case type' matching the type given to setBiquad
    case bq_type_lowpass: // Low pass filter Biquad. Filters higher frequencies
            norm = 1 / (1 + K / Q + K * K);
            a0 = K * K * norm;
            a1 = 2 * a0;
            a2 = a0;
            b1 = 2 * (K * K - 1) * norm;
            b2 = (1 - K / Q + K * K) * norm;
            break;
            
    case bq_type_highpass: // High pass filter Biquad. Filters lower frequencies
            norm = 1 / (1 + K / Q + K * K);
            a0 = 1 * norm;
            a1 = -2 * a0;
            a2 = a0;
            b1 = 2 * (K * K - 1) * norm;
            b2 = (1 - K / Q + K * K) * norm;
            break;
            
    case bq_type_bandpass: // Band pass filter. Allows frequencies in a target range
            norm = 1 / (1 + K / Q + K * K);
            a0 = K / Q * norm;
            a1 = 0;
            a2 = -a0;
            b1 = 2 * (K * K - 1) * norm;
            b2 = (1 - K / Q + K * K) * norm;
            break;
            
    case bq_type_notch: // Notch filter. Filters frequencies in a target range
            norm = 1 / (1 + K / Q + K * K);
            a0 = (1 + K * K) * norm;
            a1 = 2 * (K * K - 1) * norm;
            a2 = a0;
            b1 = a1;
            b2 = (1 - K / Q + K * K) * norm;
            break;
            
    case bq_type_peak: // Peak filter. Narrow band pass filter that allows a small range of frequencies
            if (peakGain >= 0) {    // boost
                norm = 1 / (1 + 1/Q * K + K * K);
                a0 = (1 + V/Q * K + K * K) * norm;
                a1 = 2 * (K * K - 1) * norm;
                a2 = (1 - V/Q * K + K * K) * norm;
                b1 = a1;
                b2 = (1 - 1/Q * K + K * K) * norm;
            }
            else {    // cut
                norm = 1 / (1 + V/Q * K + K * K);
                a0 = (1 + 1/Q * K + K * K) * norm;
                a1 = 2 * (K * K - 1) * norm;
                a2 = (1 - 1/Q * K + K * K) * norm;
                b1 = a1;
                b2 = (1 - V/Q * K + K * K) * norm;
            }
            break;
    case bq_type_lowshelf: // Low shelf filter. Boosts or attentuates lower frequency ranges
            if (peakGain >= 0) {    // boost
                norm = 1 / (1 + sqrt(2) * K + K * K);
                a0 = (1 + sqrt(2*V) * K + V * K * K) * norm;
                a1 = 2 * (V * K * K - 1) * norm;
                a2 = (1 - sqrt(2*V) * K + V * K * K) * norm;
                b1 = 2 * (K * K - 1) * norm;
                b2 = (1 - sqrt(2) * K + K * K) * norm;
            }
            else {    // cut
                norm = 1 / (1 + sqrt(2*V) * K + V * K * K);
                a0 = (1 + sqrt(2) * K + K * K) * norm;
                a1 = 2 * (K * K - 1) * norm;
                a2 = (1 - sqrt(2) * K + K * K) * norm;
                b1 = 2 * (V * K * K - 1) * norm;
                b2 = (1 - sqrt(2*V) * K + V * K * K) * norm;
            }
            break;
    case bq_type_highshelf: // High shelf filter. Boosts or attentuates higher frequency ranges 
            if (peakGain >= 0) {    // boost
                norm = 1 / (1 + sqrt(2) * K + K * K);
                a0 = (V + sqrt(2*V) * K + K * K) * norm;
                a1 = 2 * (K * K - V) * norm;
                a2 = (V - sqrt(2*V) * K + K * K) * norm;
                b1 = 2 * (K * K - 1) * norm;
                b2 = (1 - sqrt(2) * K + K * K) * norm;
            }
            else {    // cut
                norm = 1 / (V + sqrt(2*V) * K + K * K);
                a0 = (1 + sqrt(2) * K + K * K) * norm;
                a1 = 2 * (K * K - 1) * norm;
                a2 = (1 - sqrt(2) * K + K * K) * norm;
                b1 = 2 * (K * K - V) * norm;
                b2 = (V - sqrt(2*V) * K + K * K) * norm;
            }
            break;
    }
    
    return;
}

