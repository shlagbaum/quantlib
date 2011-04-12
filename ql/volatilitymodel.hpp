/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006 Joseph Wang, 2010 Liquidnet Holdings

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file volatilitymodel.hpp
    \brief Volatility term structures
*/

#ifndef quantlib_volatility_model_hpp
#define quantlib_volatility_model_hpp

#include <ql/types.hpp>
#include <ql/timeseries.hpp>


namespace QuantLib {

    template <class T, class Time = Date>
    class LocalVolatilityEstimator {
    public:
    	typedef Time TimeUnits;
    	typedef TimeSeriesBase<T, Time> TimeSeriesIn;
    	typedef TimeSeriesBase<Volatility, Time> TimeSeriesOut;
        virtual ~LocalVolatilityEstimator() {}
        virtual TimeSeriesOut calculate(const TimeSeriesIn &quoteSeries) { return TimeSeriesOut(); }
        virtual TimeSeries<Volatility> calculate(const TimeSeries<T> &quoteSeries) = 0;
    };

    template <class Time = Date>
    class VolatilityCompositor {
      public:
    	typedef Time TimeUnits;
    	typedef TimeSeriesBase<Volatility, Time> TimeSeriesIn;
    	typedef TimeSeriesBase<Volatility, Time> TimeSeriesOut;
        virtual ~VolatilityCompositor() {}
        virtual TimeSeriesOut calculate(const TimeSeriesIn& volatilitySeries) { return TimeSeriesOut(); }
        virtual void calibrate(const TimeSeriesIn& volatilitySeries) {}
        virtual TimeSeries<Volatility> calculate(const TimeSeries<Volatility>& volatilitySeries) = 0;
        virtual void calibrate(const TimeSeries<Volatility>& volatilitySeries) = 0;
    };
}

#endif
