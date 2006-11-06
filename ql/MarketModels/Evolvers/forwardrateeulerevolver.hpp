/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006 Ferdinando Ametrano
 Copyright (C) 2006 Mark Joshi

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


#ifndef quantlib_forward_rate_euler_evolver_hpp
#define quantlib_forward_rate_euler_evolver_hpp

#include <ql/MarketModels/marketmodelevolver.hpp>
#include <ql/MarketModels/marketmodel.hpp>
#include <ql/MarketModels/evolutiondescription.hpp>
#include <ql/MarketModels/browniangenerator.hpp>
#include <ql/MarketModels/driftcalculator.hpp>

namespace QuantLib {

    //! Euler
    class ForwardRateEulerEvolver : public MarketModelEvolver {
      public:
        ForwardRateEulerEvolver(const boost::shared_ptr<MarketModel>&,
                                const BrownianGeneratorFactory&,
                                const std::vector<Size>& numeraires);
        //! \name MarketModelEvolver interface
        //@{
        const std::vector<Size>& numeraires() const;
        Real startNewPath();
        Real advanceStep();
        Size currentStep() const;
        const CurveState& currentState() const;
        //@}
      private:
        // inputs
        boost::shared_ptr<MarketModel> marketModel_;
        std::vector<Size> numeraires_;
        boost::shared_ptr<BrownianGenerator> generator_;
        // fixed variables
        std::vector<std::vector<Real> > fixedDrifts_;
         // working variables
        Size n_, F_;
        CurveState curveState_;
        Size currentStep_;
        std::vector<Rate> forwards_, displacements_, logForwards_, initialLogForwards_;
        std::vector<Real> drifts1_, initialDrifts_;
        Array brownians_, correlatedBrownians_;
        std::vector<Size> alive_;
        // helper classes
        std::vector<DriftCalculator> calculators_;
    };

}

#endif