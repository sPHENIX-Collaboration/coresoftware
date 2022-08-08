// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
///
/// @file Feature.h
/// @brief Defines Feature interface for selecting Particles according to extracted Features.
///

#ifndef HEPMC3_FEATURE_H
#define HEPMC3_FEATURE_H

#include <functional>
#include <memory>
#include <limits>
#include "HepMC3/GenParticle.h"
#include "HepMC3/Filter.h"


namespace HepMC3 {

//////////////////////////////////////////////////////////////////////

/**
 *  @brief GenericFeature defines the Feature interface
 *  GenericFeature is not intended to be used directly.  The
 *  derived Feature class and its specialisations should be used.
 *
 *  A Feature wraps a function object that can extract a
 *  generic Feature_type from a ConstGenParticlePtr.  Usually the
 *  Feature_type would be something like int (e.g. status) or
 *  double (e.g. pT), but it could in principle be any attribute of a
 *  particle so long as there are well defined <, <=, >, >=, == and
 *  != operators for that attribute, as well as an abs function.
 *
 *  Once a Feature is defined, you can obtain Filters that select
 *  Particles according to that Feature by e.g.
 *  Feature<int> status([](ConstGenParticlePtr p)->int{return p->status();});
 *  bool is_stable = (status == 1)(p);
 *  Filter is_beam = (status == 4);
 *  bool beam = is_beam(p);
 *
 *  An abs function is also defined, so abs(Feature) works as you'd
 *  expect, e.g.
 *  Feature<double> rapidity([](ConstGenParticlePtr p)->double{return p->momentum().rap();});
 *  Filter rapCut = abs(rapidity) < 2.5;
 *
 *  Please also see the Selector interface, which defines an
 *  abstract interface to Feature that is free of the template params
 *  and also includes some standard Features such as
 *
 *  Selector::STATUS;
 *  Selector::PDG_ID;
 *  Selector::PT;
 *  Selector::RAPIDITY;
 */
template<typename Feature_type>
class GenericFeature {
public:
    /// @brief evaluator type
    using Evaluator_type = std::function<Feature_type(ConstGenParticlePtr)>;
    /// @brief shared pointer for evaluator type
    using EvaluatorPtr   =  std::shared_ptr<Evaluator_type>;

    /// @brief access the underlying feature value
    Feature_type operator()(ConstGenParticlePtr input)const {
        return (*m_internal)(input);
    }

    /// @brief greater than operator
    /// @return Filter function
    Filter operator > (Feature_type value) const {
        EvaluatorPtr functor = m_internal;
        return [value, functor](ConstGenParticlePtr input)->bool{return (*functor)(input) >  value;};
    }
    /// @brief less than operator
    /// @return Filter function
    Filter operator < (Feature_type value) const {
        EvaluatorPtr functor = m_internal;
        return [value, functor](ConstGenParticlePtr input)->bool{return (*functor)(input) <  value;};
    }

    /// @brief greater than or equals operator
    /// @return Filter function
    Filter operator >= (Feature_type value) const {
        EvaluatorPtr functor = m_internal;
        return [value, functor](ConstGenParticlePtr input)->bool{return (*functor)(input) >=  value;};
    }

    /// @brief less than or equals operator
    /// @return Filter function
    Filter operator <= (Feature_type value) const {
        EvaluatorPtr functor = m_internal;
        return [value, functor](ConstGenParticlePtr input)->bool{return (*functor)(input) <=  value;};
    }

    /// @brief equality operator
    /// @return Filter function
    virtual Filter operator == (Feature_type value) const {
        EvaluatorPtr functor = m_internal;
        return [value, functor](ConstGenParticlePtr input)->bool{return (*functor)(input) == value;};
    }

    /// @brief inequality operator
    /// @return Filter function
    virtual Filter operator != (Feature_type value) const {
        EvaluatorPtr functor = m_internal;
        return [value, functor](ConstGenParticlePtr input)->bool{return (*functor)(input) != value;};
    }

protected:
    /// Hide the constructor so no one can use GenericFeature directly
    GenericFeature(Evaluator_type functor):m_internal(std::make_shared<Evaluator_type>(functor)) {}

    /// Hide the copy constructor
    GenericFeature(const GenericFeature &copy) : m_internal(copy.m_internal) {}

    /// @brief internal copy of func for evaluation
    /// on the heap so will persist in resulting Filters even if
    /// parent Feature object was destroyed
    EvaluatorPtr m_internal;
};

//////////////////////////////////////////////////////////////////////

/** @brief Expose GenericFeature interface to derived Feature class
 *
 *  This will get used for generic class types that aren't integral or
 *  floating point types.
 *
 *  A Feature wraps a function object that can extract a
 *  generic Feature_type from a ConstGenParticlePtr.  Usually the
 *  Feature_type would be something like int (e.g. status) or
 *  double (e.g. pT), but it could in principle be any attribute of a
 *  particle so long as there are well defined <, <=, >, >=, == and
 *  != operators for that attribute, as well as an abs function.
 *
 *  Once a Feature is defined, you can obtain Filters that select
 *  Particles according to that Feature by e.g.
 *  Feature<int> status([](ConstGenParticlePtr p)->int{return p->status();});
 *  bool is_stable = (status == 1)(p);
 *  Filter is_beam = (status == 4);
 *  bool beam = is_beam(p);
 *
 *  An abs function is also defined, so abs(Feature) works as you'd
 *  expect, e.g.
 *  Feature<double> rapidity([](ConstGenParticlePtr p)->double{return p->momentum().rap();});
 *  Filter rapCut = abs(rapidity) < 2.5;
 *
 *  Please also see the Selector interface, which defines an
 *  abstract interface to Feature that is free of the template params
 *  and also includes some standard Features such as
 *
 *  Selector::STATUS;
 *  Selector::PDG_ID;
 *  Selector::PT;
 *  Selector::RAPIDITY;

 */
template<typename Feature_type, typename Dummy = void>
class Feature : public GenericFeature<Feature_type> {
public:
    using typename GenericFeature<Feature_type>::Evaluator_type;
    using typename GenericFeature<Feature_type>::EvaluatorPtr;
    using GenericFeature<Feature_type>::m_internal;

    using GenericFeature<Feature_type>::operator ();
    using GenericFeature<Feature_type>::operator >;
    using GenericFeature<Feature_type>::operator >=;
    using GenericFeature<Feature_type>::operator <;
    using GenericFeature<Feature_type>::operator <=;
    using GenericFeature<Feature_type>::operator ==;
    using GenericFeature<Feature_type>::operator !=;

    /// @brief  Feature
    Feature(Evaluator_type functor) : GenericFeature<Feature_type>(functor) {}
    /// @brief  Copy
    Feature(const Feature &copy) : GenericFeature<Feature_type>(copy) {}

    /// @brief  Abs function
    Feature<Feature_type> abs() const {
        EvaluatorPtr functor = m_internal;
        Evaluator_type absfunctor = [functor](ConstGenParticlePtr p)->Feature_type{return ::abs((*functor)(p));};
        return Feature<Feature_type>(absfunctor);
    }
};

//////////////////////////////////////////////////////////////////////

/** @brief Specialisation of Feature for integral types
 *
 *  It is a valid operator to compare an int to a float, but the
 *  generic version of these operators in the base class will
 *  first cast input float to an int, then compare that.  In some cases
 *  the comparison will be incorrect because of rounding the float.
 *  e.g. int x=5; float y=5.5; bool result = x<y; would be wrong
 *  because y first gets converted to int 5.
 *
 *  To solve this, we provide specialised comparison operators for
 *  integral type and double.  Note that the opposite specialisation
 *  in which the Feature_type is floating_point is not necessary
 */
template<typename Feature_type>
class Feature<Feature_type, typename std::enable_if<std::is_integral<Feature_type>::value, void>::type> : public GenericFeature<Feature_type> {
public:
    using GenericFeature<Feature_type>::operator ();
    using GenericFeature<Feature_type>::operator >;
    using GenericFeature<Feature_type>::operator >=;
    using GenericFeature<Feature_type>::operator <;
    using GenericFeature<Feature_type>::operator <=;
    using GenericFeature<Feature_type>::operator ==;
    using GenericFeature<Feature_type>::operator !=;

    using typename GenericFeature<Feature_type>::Evaluator_type;
    using typename GenericFeature<Feature_type>::EvaluatorPtr;

    using GenericFeature<Feature_type>::m_internal;

    /// @brief  Feature
    Feature(Evaluator_type functor) : GenericFeature<Feature_type>(functor) {}
    /// @brief  Feature
    Feature(const Feature &copy) : GenericFeature<Feature_type>(copy) {}

    /// @brief abs function
    Feature<Feature_type> abs() const {
        EvaluatorPtr functor = m_internal;
        Evaluator_type absfunctor = [functor](ConstGenParticlePtr p)->Feature_type{return ::abs((*functor)(p));};
        return Feature<Feature_type>(absfunctor);
    }

    /// @brief greater operator
    Filter operator > (double value) const {
        EvaluatorPtr functor = m_internal;
        return [value, functor](ConstGenParticlePtr input)->bool{return (*functor)(input) >  value;};
    }

    /// @brief less operator
    Filter operator < (double value) const {
        EvaluatorPtr functor = m_internal;
        return [value, functor](ConstGenParticlePtr input)->bool{return (*functor)(input) <  value;};
    }

    /// @brief equal operator
    Filter operator == (double value) const {
        EvaluatorPtr functor = m_internal;
        return [value, functor](ConstGenParticlePtr input)->bool{
            Feature_type local = (*functor)(input);
            return std::abs(local - value) <  std::numeric_limits<double>::epsilon();
        };
    }

    /// @brief greater or equal operator
    Filter operator >= (double value) const { return !( (*this) < value );}

    /// @brief less or equal operator
    Filter operator <= (double value) const { return !( (*this) > value );}

    /// @brief not equal operator
    Filter operator != (double value) const {
        return !( (*this) == value );
    }
};

//////////////////////////////////////////////////////////////////////

/** @brief specialisation of Feature for floating point type
 *
 *  Test of equality of floating point types is not safe.  Here we
 *  provide a "reasonable" definition of equality based on the
 *  floating point precision.
 */

template<typename Feature_type>
class Feature<Feature_type, typename std::enable_if<std::is_floating_point<Feature_type>::value, void>::type> : public GenericFeature<Feature_type> {
public:
    using typename GenericFeature<Feature_type>::Evaluator_type;
    using typename GenericFeature<Feature_type>::EvaluatorPtr;

    using GenericFeature<Feature_type>::operator ();
    using GenericFeature<Feature_type>::operator >;
    using GenericFeature<Feature_type>::operator >=;
    using GenericFeature<Feature_type>::operator <;
    using GenericFeature<Feature_type>::operator <=;

    using GenericFeature<Feature_type>::m_internal;

    /// @brief Feature
    Feature(Evaluator_type functor) : GenericFeature<Feature_type>(functor) {}
    /// @brief Copy
    Feature(const Feature &copy) : GenericFeature<Feature_type>(copy) {}

    /// @brief abs function
    Feature<Feature_type> abs() const {
        EvaluatorPtr functor = m_internal;
        Evaluator_type absfunctor = [functor](ConstGenParticlePtr p)->Feature_type{return std::abs((*functor)(p));};
        return Feature<Feature_type>(absfunctor);
    }

    Filter operator == (Feature_type value) const override {
        EvaluatorPtr functor = m_internal;
        return [value, functor](ConstGenParticlePtr input)->bool{
            Feature_type local = (*functor)(input);
            return std::less_equal<Feature_type>{}(fabs(local - value) ,  std::numeric_limits<Feature_type>::epsilon());
        };
    }

    Filter operator != (Feature_type value) const override {
        return !( (*this) == value );
    }
};

//////////////////////////////////////////////////////////////////////

/**
 *  @brief Obtain the absolute value of a Feature.
 *  This works as you'd expect.  If foo is a valid Feature, then
 *  abs(foo) returns a new Feature that corresponds to the absolute
 *  value of the foo feature.  You can construct a Filter from that in
 *  the usual way with e.g. Filter f = abs(foo) > 10.;
 */
template<typename Feature_type>
Feature<Feature_type> abs(const Feature<Feature_type> &input) {
    return input.abs();
}

}

#endif
