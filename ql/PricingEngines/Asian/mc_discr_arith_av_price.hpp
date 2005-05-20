/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2003, 2004 Ferdinando Ametrano
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/reference/license.html>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file mc_discr_arith_av_price.hpp
    \brief Monte Carlo engine for discrete arithmetic average price Asian
*/

#ifndef quantlib_mc_discrete_arithmetic_average_price_asian_engine_hpp
#define quantlib_mc_discrete_arithmetic_average_price_asian_engine_hpp

#include <ql/PricingEngines/Asian/mc_discr_geom_av_price.hpp>
#include <ql/PricingEngines/Asian/analytic_discr_geom_av_price.hpp>

namespace QuantLib {

    //!  Monte Carlo pricing engine for discrete arithmetic average price Asian
    /*!  Monte Carlo pricing engine for discrete arithmetic average price
         Asian options. It can use MCDiscreteGeometricAPEngine (Monte Carlo
         discrete arithmetic average price engine) and
         AnalyticDiscreteGeometricAveragePriceAsianEngine (analytic discrete
         arithmetic average price engine) for control variation.

         \ingroup asianengines

         \test the correctness of the returned value is tested by
               reproducing results available in literature.
    */
    template <class RNG = PseudoRandom, class S = Statistics>
    class MCDiscreteArithmeticAPEngine :
                                public MCDiscreteAveragingAsianEngine<RNG,S> {
      public:
        typedef typename MCDiscreteAveragingAsianEngine<RNG,S>::path_generator_type
            path_generator_type;
        typedef typename MCDiscreteAveragingAsianEngine<RNG,S>::path_pricer_type
            path_pricer_type;
        typedef typename MCDiscreteAveragingAsianEngine<RNG,S>::stats_type
            stats_type;
        // constructor
        MCDiscreteArithmeticAPEngine(Size maxTimeStepPerYear,
                                     bool brownianBridge,
                                     bool antitheticVariate,
                                     bool controlVariate,
                                     Size requiredSamples,
                                     Real requiredTolerance,
                                     Size maxSamples,
                                     BigNatural);
      protected:
        boost::shared_ptr<path_pricer_type> pathPricer() const;
        boost::shared_ptr<path_pricer_type> controlPathPricer() const;
        boost::shared_ptr<PricingEngine> controlPricingEngine() const {
            return boost::shared_ptr<PricingEngine>(new
                AnalyticDiscreteGeometricAveragePriceAsianEngine());
        }
    };


    class ArithmeticAPOPathPricer : public PathPricer<Path> {
      public:
        ArithmeticAPOPathPricer(Option::Type type,
                                Real underlying,
                                Real strike,
                                DiscountFactor discount,
                                Real runningSum = 0.0,
                                Size pastFixings = 0);
        Real operator()(const Path& path) const  {
            #ifndef QL_DISABLE_DEPRECATED
            Size n = path.size();
            #else
            Size n = path.length() - 1;
            #endif
            QL_REQUIRE(n>0, "the path cannot be empty");

            Real price = underlying_, averagePrice = runningSum_;
            Size fixings = n + pastFixings_;
            // not sure the if is correct
            if (path.timeGrid().mandatoryTimes()[0]==0.0) {
                averagePrice = price;
                fixings = n + pastFixings_ + 1;
            }
            #ifndef QL_DISABLE_DEPRECATED
            for (Size i=0; i<n; i++) {
                price *= std::exp(path[i]);
                averagePrice += price;
            }
            #else
            for (Size i=0; i<n; i++) {
                price = path.value(i+1);
                averagePrice += price;
            }
            #endif
            averagePrice = averagePrice/fixings;
            return discount_ * payoff_(averagePrice);
        }
      private:
        Real underlying_;
        PlainVanillaPayoff payoff_;
        DiscountFactor discount_;
        Real runningSum_;
        Size pastFixings_;
    };


    // inline definitions

    template <class RNG, class S>
    inline
    MCDiscreteArithmeticAPEngine<RNG,S>::MCDiscreteArithmeticAPEngine(
                                                    Size maxTimeStepPerYear,
                                                    bool brownianBridge,
                                                    bool antitheticVariate,
                                                    bool controlVariate,
                                                    Size requiredSamples,
                                                    Real requiredTolerance,
                                                    Size maxSamples,
                                                    BigNatural seed)
    : MCDiscreteAveragingAsianEngine<RNG,S>(maxTimeStepPerYear,
                                            brownianBridge,
                                            antitheticVariate,
                                            controlVariate,
                                            requiredSamples,
                                            requiredTolerance,
                                            maxSamples,
                                            seed) {}

    template <class RNG, class S>
    inline
    boost::shared_ptr<QL_TYPENAME MCDiscreteArithmeticAPEngine<RNG,S>::path_pricer_type>
        MCDiscreteArithmeticAPEngine<RNG,S>::pathPricer() const {

        boost::shared_ptr<PlainVanillaPayoff> payoff =
            boost::dynamic_pointer_cast<PlainVanillaPayoff>(
                this->arguments_.payoff);
        QL_REQUIRE(payoff, "non-plain payoff given");

        boost::shared_ptr<EuropeanExercise> exercise =
            boost::dynamic_pointer_cast<EuropeanExercise>(
                this->arguments_.exercise);
        QL_REQUIRE(exercise, "wrong exercise given");

        boost::shared_ptr<BlackScholesProcess> process =
            boost::dynamic_pointer_cast<BlackScholesProcess>(
                                          this->arguments_.stochasticProcess);
        QL_REQUIRE(process, "Black-Scholes process required");

        return boost::shared_ptr<QL_TYPENAME
            MCDiscreteArithmeticAPEngine<RNG,S>::path_pricer_type>(
            new ArithmeticAPOPathPricer(
              payoff->optionType(),
              process->stateVariable()->value(),
              payoff->strike(),
              process->riskFreeRate()->discount(this->timeGrid().back())));
    }

    template <class RNG, class S>
    inline
    boost::shared_ptr<QL_TYPENAME MCDiscreteArithmeticAPEngine<RNG,S>::path_pricer_type>
        MCDiscreteArithmeticAPEngine<RNG,S>::controlPathPricer() const {

        boost::shared_ptr<PlainVanillaPayoff> payoff =
            boost::dynamic_pointer_cast<PlainVanillaPayoff>(
                this->arguments_.payoff);
        QL_REQUIRE(payoff, "non-plain payoff given");

        boost::shared_ptr<EuropeanExercise> exercise =
            boost::dynamic_pointer_cast<EuropeanExercise>(
                this->arguments_.exercise);
        QL_REQUIRE(exercise, "wrong exercise given");

        boost::shared_ptr<BlackScholesProcess> process =
            boost::dynamic_pointer_cast<BlackScholesProcess>(
                                          this->arguments_.stochasticProcess);
        QL_REQUIRE(process, "Black-Scholes process required");

        // for seasoned option the geometric strike might be rescaled
        // to obtain an equivalent arithmetic strike.
        // Any change applied here MUST be applied to the analytic engine too
        return boost::shared_ptr<QL_TYPENAME
            MCDiscreteArithmeticAPEngine<RNG,S>::path_pricer_type>(
            new GeometricAPOPathPricer(
              payoff->optionType(),
              process->stateVariable()->value(),
              payoff->strike(),
              process->riskFreeRate()->discount(this->timeGrid().back())));
    }

    template <class RNG = PseudoRandom, class S = Statistics>
    class MakeMCDiscreteArithmeticAPEngine {
      public:
        MakeMCDiscreteArithmeticAPEngine();
        // named parameters
        MakeMCDiscreteArithmeticAPEngine& withStepsPerYear(Size maxSteps);
        MakeMCDiscreteArithmeticAPEngine& withBrownianBridge(bool b = true);
        MakeMCDiscreteArithmeticAPEngine& withSamples(Size samples);
        MakeMCDiscreteArithmeticAPEngine& withTolerance(Real tolerance);
        MakeMCDiscreteArithmeticAPEngine& withMaxSamples(Size samples);
        MakeMCDiscreteArithmeticAPEngine& withSeed(BigNatural seed);
        MakeMCDiscreteArithmeticAPEngine& withAntitheticVariate(bool b = true);
        MakeMCDiscreteArithmeticAPEngine& withControlVariate(bool b = true);
        // conversion to pricing engine
        operator boost::shared_ptr<PricingEngine>() const;
      private:
        bool antithetic_, controlVariate_;
        Size steps_, samples_, maxSamples_;
        Real tolerance_;
        bool brownianBridge_;
        BigNatural seed_;
    };

    template <class RNG, class S>
    inline MakeMCDiscreteArithmeticAPEngine<RNG,S>::MakeMCDiscreteArithmeticAPEngine()
    : antithetic_(false), controlVariate_(false),
      steps_(Null<Size>()), samples_(Null<Size>()), maxSamples_(Null<Size>()),
      tolerance_(Null<Real>()), brownianBridge_(true), seed_(0) {}

    template <class RNG, class S>
    inline MakeMCDiscreteArithmeticAPEngine<RNG,S>&
    MakeMCDiscreteArithmeticAPEngine<RNG,S>::withStepsPerYear(Size maxSteps) {
        steps_ = maxSteps;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCDiscreteArithmeticAPEngine<RNG,S>&
    MakeMCDiscreteArithmeticAPEngine<RNG,S>::withSamples(Size samples) {
        QL_REQUIRE(tolerance_ == Null<Real>(),
                   "tolerance already set");
        samples_ = samples;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCDiscreteArithmeticAPEngine<RNG,S>&
    MakeMCDiscreteArithmeticAPEngine<RNG,S>::withTolerance(Real tolerance) {
        QL_REQUIRE(samples_ == Null<Size>(),
                   "number of samples already set");
        QL_REQUIRE(RNG::allowsErrorEstimate,
                   "chosen random generator policy "
                   "does not allow an error estimate");
        tolerance_ = tolerance;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCDiscreteArithmeticAPEngine<RNG,S>&
    MakeMCDiscreteArithmeticAPEngine<RNG,S>::withMaxSamples(Size samples) {
        maxSamples_ = samples;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCDiscreteArithmeticAPEngine<RNG,S>&
    MakeMCDiscreteArithmeticAPEngine<RNG,S>::withSeed(BigNatural seed) {
        seed_ = seed;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCDiscreteArithmeticAPEngine<RNG,S>&
    MakeMCDiscreteArithmeticAPEngine<RNG,S>::withBrownianBridge(bool b) {
        brownianBridge_ = b;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCDiscreteArithmeticAPEngine<RNG,S>&
    MakeMCDiscreteArithmeticAPEngine<RNG,S>::withAntitheticVariate(bool b) {
        antithetic_ = b;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCDiscreteArithmeticAPEngine<RNG,S>&
    MakeMCDiscreteArithmeticAPEngine<RNG,S>::withControlVariate(bool b) {
        controlVariate_ = b;
        return *this;
    }

    template <class RNG, class S>
    inline
    MakeMCDiscreteArithmeticAPEngine<RNG,S>::operator boost::shared_ptr<PricingEngine>()
                                                                      const {
        QL_REQUIRE(steps_ != Null<Size>(),
                   "max number of steps per year not given");
        return boost::shared_ptr<PricingEngine>(new
            MCDiscreteArithmeticAPEngine<RNG,S>(steps_,
                                                brownianBridge_,
                                                antithetic_, controlVariate_,
                                                samples_, tolerance_,
                                                maxSamples_,
                                                seed_));
    }



}


#endif
