/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006 Joseph Wang

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
    class LocalVolatilityEstimatorTime {
    public:
        virtual ~LocalVolatilityEstimatorTime() {}
        virtual TimeSeriesBase<Volatility, Time>
        calculate(const TimeSeriesBase<T, Time> &quoteSeries) = 0;
    };

    template <class T>
	class LocalVolatilityEstimator : public LocalVolatilityEstimatorTime<T, Date> {
    public:
        virtual TimeSeries<Volatility>
        calculate(const TimeSeries<T> &quoteSeries) = 0;
        virtual TimeSeriesBase<Volatility, Date>
		calculate(const TimeSeriesBase<T, Date> &quoteSeries) {
			return calculate(static_cast<const TimeSeries<T> &>(quoteSeries));
		}
    };

    template <class Time>
	class VolatilityCompositorTime {
      public:
		typedef Time time;
		typedef TimeSeriesBase<Volatility, Time> time_series;

        virtual ~VolatilityCompositorTime() {}
        virtual time_series calculate(const time_series& volatilitySeries) = 0;
        virtual void calibrate(const time_series& volatilitySeries) = 0;
    };

	class VolatilityCompositor : public VolatilityCompositorTime<Date> {
	public:
		typedef Date time;
		typedef TimeSeries<Volatility> time_series;
		typedef TimeSeriesBase<Volatility, Date> time_series_base;

        virtual time_series calculate(const time_series& volatilitySeries) = 0;
        virtual void calibrate(const time_series& volatilitySeries) = 0;
        virtual time_series_base calculate(const time_series_base& volatilitySeries) {
			return calculate(static_cast<const time_series &>(volatilitySeries));
		}
		virtual void calibrate(const time_series_base & volatilitySeries) {
			calibrate(static_cast<const time_series &>(volatilitySeries));
		}
	};

}


#endif
