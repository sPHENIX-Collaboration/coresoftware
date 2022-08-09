// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
///
/// @file Selector.cc
/// @brief Implementation  of Selector wrappers
///
#include "HepMC3/Selector.h"

namespace HepMC3 {
const SelectorWrapper<int> StandardSelector::STATUS      = SelectorWrapper<int>([](ConstGenParticlePtr p)->int{return p->status();});
const SelectorWrapper<int> StandardSelector::PDG_ID      = SelectorWrapper<int>([](ConstGenParticlePtr p)->int{return p->pdg_id();});
const SelectorWrapper<double> StandardSelector::PT       = SelectorWrapper<double>([](ConstGenParticlePtr p)->double{return p->momentum().pt();});
const SelectorWrapper<double> StandardSelector::ENERGY   = SelectorWrapper<double>([](ConstGenParticlePtr p)->double{return p->momentum().e();});
const SelectorWrapper<double> StandardSelector::RAPIDITY = SelectorWrapper<double>([](ConstGenParticlePtr p)->double{return p->momentum().rap();});
const SelectorWrapper<double> StandardSelector::ETA      = SelectorWrapper<double>([](ConstGenParticlePtr p)->double{return p->momentum().eta();});
const SelectorWrapper<double> StandardSelector::PHI      = SelectorWrapper<double>([](ConstGenParticlePtr p)->double{return p->momentum().phi();});
const SelectorWrapper<double> StandardSelector::ET       = SelectorWrapper<double>([](ConstGenParticlePtr p)->double{return p->momentum().e() * (p->momentum().pt() / p->momentum().p3mod());});
const SelectorWrapper<double> StandardSelector::MASS     = SelectorWrapper<double>([](ConstGenParticlePtr p)->double{return p->momentum().m();});

ConstSelectorPtr abs(const Selector &input)
{
    return input.abs();
}

AttributeFeature Selector::ATTRIBUTE(const std::string &name) {return AttributeFeature(name);}

} // namespace HepMC3
