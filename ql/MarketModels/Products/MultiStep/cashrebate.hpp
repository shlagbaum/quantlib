/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
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


#ifndef quantlib_marketmodel_rebate_hpp
#define quantlib_marketmodel_rebate_hpp

#include <ql/MarketModels/marketmodelproduct.hpp>

namespace QuantLib {

    class MarketModelCashRebate : public MarketModelMultiProduct {
      public:
        MarketModelCashRebate(const EvolutionDescription& evolution,
                              const std::vector<Time>& paymentTimes,
                              const Matrix& amounts,
                              Size numberOfProducts);
        //! \name MarketModelMultiProduct interface
        //@{
        EvolutionDescription suggestedEvolution() const;
        std::vector<Time> possibleCashFlowTimes() const;
        Size numberOfProducts() const;
        Size maxNumberOfCashFlowsPerProductPerStep() const;
        void reset(); 
        bool nextTimeStep(
                      const CurveState& currentState, 
                      std::vector<Size>& numberCashFlowsThisStep,
                      std::vector<std::vector<CashFlow> >& cashFlowsGenerated);
        //@}
      private:
        EvolutionDescription evolution_;
        std::vector<Time> paymentTimes_;
        Matrix amounts_;
        Size numberOfProducts_;
        // things that vary in a path
        Size currentIndex_;
    };

}

#endif