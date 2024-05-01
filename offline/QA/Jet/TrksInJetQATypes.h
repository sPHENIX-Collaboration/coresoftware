// ----------------------------------------------------------------------------
// 'TrksInJetQATypes.h'
// Derek Anderson
// 04.30.2024
//
// Conveninent type definitions for the TrksInJetQA module.
// ----------------------------------------------------------------------------

#ifndef TRKSINJETSQATYPES_H
#define TRKSINJETSQATYPES_H

// c++ utilities
#include <string>
#include <vector>
#include <utility>
// root libraries
#include <TH1.h>
#include <TH2.h>
// particle flow libraries
#include <particleflowreco/ParticleFlowElement.h>
#include <particleflowreco/ParticleFlowElementContainer.h>



// type definitions -----------------------------------------------------------

typedef std::pair<float, float>                 BinRange;
typedef std::pair<uint32_t, BinRange>           BinDef;
typedef std::tuple<std::string, BinDef>         HistDef1D;
typedef std::tuple<std::string, BinDef, BinDef> HistDef2D;
typedef std::vector<HistDef1D>                  VecHistDef1D;
typedef std::vector<HistDef2D>                  VecHistDef2D;
typedef std::vector<std::vector<TH1D*>>         VecHist1D;
typedef std::vector<std::vector<TH2D*>>         VecHist2D;
typedef std::vector<std::string>                VecHistTypes;
typedef ParticleFlowElement                     PFObject;
typedef ParticleFlowElementContainer            PFObjectStore;

#endif

// end ------------------------------------------------------------------------
