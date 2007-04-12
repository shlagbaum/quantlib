/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006, 2007 Ferdinando Ametrano
 Copyright (C) 2006 Mark Joshi
 Copyright (C) 2007 StatPro Italia srl

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

#include <ql/marketmodels/models/expcorrflatvol.hpp>
#include <ql/math/pseudosqrt.hpp>
#include <ql/math/interpolations/linearinterpolation.hpp>
#include <ql/marketmodels/models/correlations.hpp>

namespace QuantLib {

    ExpCorrFlatVol::ExpCorrFlatVol(
            const std::vector<Volatility>& volatilities,
            const Matrix& correlations,
            const EvolutionDescription& evolution,
            Size numberOfFactors,
            const std::vector<Rate>& initialRates,
            const std::vector<Spread>& displacements)
    : numberOfFactors_(numberOfFactors),
      numberOfRates_(initialRates.size()),
      numberOfSteps_(evolution.evolutionTimes().size()),
      initialRates_(initialRates),
      displacements_(displacements),
      evolution_(evolution),
      pseudoRoots_(numberOfSteps_, Matrix(numberOfRates_, numberOfFactors_))
    {
        const std::vector<Time>& rateTimes = evolution.rateTimes();
        QL_REQUIRE(numberOfRates_==rateTimes.size()-1,
                   "mismatch between number of rates (" << numberOfRates_ <<
                   ") and rate times");
        QL_REQUIRE(numberOfRates_==displacements.size(),
                   "mismatch between number of rates (" << numberOfRates_ <<
                   ") and displacements (" << displacements.size() << ")");
        QL_REQUIRE(numberOfRates_==volatilities.size(),
                   "mismatch between number of rates (" << numberOfRates_ <<
                   ") and volatilities (" << volatilities.size() << ")");
        QL_REQUIRE(numberOfRates_<=numberOfFactors_*numberOfSteps_,
                   "number of rates (" << numberOfRates_ <<
                   ") greater than number of factors (" << numberOfFactors_
                   << ") times number of steps (" << numberOfSteps_ << ")");

        std::vector<Volatility> stdDev(numberOfRates_);

        Time effStartTime;
        const Matrix& effectiveStopTime = evolution.effectiveStopTime();
        Matrix covariance(numberOfRates_, numberOfRates_);
        for (Size k=0; k<numberOfSteps_; ++k) {
            Size i;
            for (i=0; i<numberOfRates_; ++i) {
                effStartTime = (k>0 ? effectiveStopTime[k-1][i] : 0.0);
                stdDev[i] = volatilities[i] *
                    std::sqrt(effectiveStopTime[k][i]-effStartTime);
            }

            for (i=0; i<numberOfRates_; ++i) {
                for (Size j=i; j<numberOfRates_; ++j) {
                     covariance[i][j] =  covariance[j][i] =
                         stdDev[j] * correlations[i][j] * stdDev[i];
                 }
            }

            pseudoRoots_[k] = rankReducedSqrt(covariance,
                                              numberOfFactors, 1.0,
                                              SalvagingAlgorithm::None);

            QL_ENSURE(pseudoRoots_[k].rows()==numberOfRates_,
                      "step " << k
                      << " flat vol wrong number of rows: "
                      << pseudoRoots_[k].rows()
                      << " instead of " << numberOfRates_);
            QL_ENSURE(pseudoRoots_[k].columns()==numberOfFactors,
                      "step " << k
                      << " flat vol wrong number of columns: "
                      << pseudoRoots_[k].columns()
                      << " instead of " << numberOfFactors_);
        }
    }


    ExpCorrFlatVolFactory::ExpCorrFlatVolFactory(
                                Real longTermCorrelation,
                                Real beta,
                                const std::vector<Time>& times,
                                const std::vector<Volatility>& vols,
                                const Handle<YieldTermStructure>& yieldCurve,
                                Spread displacement)
    : longTermCorrelation_(longTermCorrelation), beta_(beta),
      times_(times), vols_(vols), yieldCurve_(yieldCurve),
      displacement_(displacement) {
        volatility_ = LinearInterpolation(times_.begin(), times_.end(),
                                          vols_.begin());
        volatility_.update();
        registerWith(yieldCurve_);
    }

    boost::shared_ptr<MarketModel>
    ExpCorrFlatVolFactory::create(const EvolutionDescription& evolution,
                                  Size numberOfFactors) const {
        const std::vector<Time>& rateTimes = evolution.rateTimes();
        Size numberOfRates = rateTimes.size()-1;

        std::vector<Rate> initialRates(numberOfRates);
        for (Size i=0; i<numberOfRates; ++i)
            initialRates[i] = yieldCurve_->forwardRate(rateTimes[i],
                                                       rateTimes[i+1],
                                                       Simple);

        std::vector<Volatility> displacedVolatilities(numberOfRates);
        for (Size i=0; i<numberOfRates; ++i) {
            Volatility vol = // to be changes
                volatility_(rateTimes[i]);
            displacedVolatilities[i] =
                initialRates[i]*vol/(initialRates[i]+displacement_);
        }

        std::vector<Spread> displacements(numberOfRates, displacement_);

        Matrix correlations = exponentialCorrelations(evolution.rateTimes(),
                                                      longTermCorrelation_,
                                                      beta_);
        return boost::shared_ptr<MarketModel>(new
            ExpCorrFlatVol(displacedVolatilities,
                           correlations,
                           evolution,
                           numberOfFactors,
                           initialRates,
                           displacements));
    }

    void ExpCorrFlatVolFactory::update() {
        notifyObservers();
    }

}

