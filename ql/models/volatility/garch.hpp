/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006 Joseph Wang
 Copyright (C) 2010 Liquidnet Holdings, Inc.

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

/*! \file garch.hpp
    \brief GARCH volatility model
*/

#ifndef quantlib_garch_volatility_model_hpp
#define quantlib_garch_volatility_model_hpp

#include <ql/volatilitymodel.hpp>
#include <ql/math/array.hpp>
#include <ql/math/matrix.hpp>
#include <ql/math/optimization/constraint.hpp>
#include <ql/math/optimization/problem.hpp>
#include <ql/math/optimization/endcriteria.hpp>
#include <ql/errors.hpp>
#include <vector>

namespace QuantLib {

    //! GARCH volatility model
    /*! 
    */

	typedef boost::shared_ptr<Problem> ProblemPtr;
	
	namespace Garch {
	
		//! calibrates GARCH for r^2
		ProblemPtr calibrate_r2 (Natural mode,      // 0 - 1st type initial guess, 1 - 2nd type initial guess, 2 - best of two, 3 - double optimization;
			const std::vector<Volatility> &r2,  // vector of r^2
			Real meanr2,                        // <r^2>
			Real &alpha,                        // results
			Real &beta, 
			Real &omega);
	
		//! the same with user defined optimization method and end criteria
		ProblemPtr calibrate_r2 (Natural mode,
			const std::vector<Volatility> &r2, 
			Real meanr2,
			OptimizationMethod &method, 
			const EndCriteria &endCriteria,
			Real &alpha, 
			Real &beta, 
			Real &omega);
	
		//! the same with user defined initial guess
		ProblemPtr calibrate_r2 (const std::vector<Volatility> &r2,
			Real meanr2,
			OptimizationMethod &method,
			const EndCriteria &endCriteria,
			const Array &initialGuess, 
			Real &alpha, 
			Real &beta, 
			Real &omega);
	
		//! the same with user defined constraints and initial guess
		ProblemPtr calibrate_r2 (const std::vector<Volatility> &r2,
			Real meanr2,
			OptimizationMethod &method,
			Constraint &constraints,
			const EndCriteria &endCriteria,
			const Array &initialGuess, 
			Real &alpha, 
			Real &beta, 
			Real &omega);
	
		//! a helper for calculation of r^2 and <r^2>
		template<typename Iterator>
		Real to_r2 (Iterator begin, Iterator end, std::vector<Volatility> &r2) {
			Real u2(0.0), meanr2(0.0), w(1.0);
			for (; begin != end; ++begin) {
				u2 = *begin; u2 *= u2;
				meanr2 = (1.0 - w) * meanr2 + w * u2;
				r2.push_back(u2);
				w /= (w + 1.0);
			}
			return meanr2;
		}
	}; // namespace Garch
	
	template <class Time>
	class Garch11T : public VolatilityCompositor<Time> {
	public:
		typedef typename VolatilityCompositor<Time>::TimeUnits TimeUnits;
		typedef typename VolatilityCompositor<Time>::TimeSeriesIn TimeSeriesIn;
		typedef typename VolatilityCompositor<Time>::TimeSeriesOut TimeSeriesOut;
		typedef typename TimeSeriesIn::const_iterator ts_const_iterator;
		typedef typename TimeSeriesIn::const_value_iterator ts_value_iterator;
	
		Garch11T(Real a = 0.0, Real b = 0.0, Real vl = 0.0) :
			alpha_(a), beta_(b), gamma_ (1 - a - b), vl_(vl), logLikelihood_(0), mode_(2) {
		}
	
		Garch11T(const TimeSeriesIn& qs, Natural mode = 2) : alpha_(0), beta_(0), vl_(0), logLikelihood_(0), mode_(mode) {
			calibrate(qs);
		};
	
		Real alpha() const { return alpha_; }
	
		Real beta() const { return beta_; }
	
		Real omega() const { return vl_ * gamma_; }
	
		Real ltVol() const { return vl_; }
	
		Real logLikelihood() const { return logLikelihood_; }
	
		Natural mode() const { return mode_; }
	
		void setMode(Natural mode) { mode_ = mode; }
	
		TimeSeriesOut calculate(const TimeSeriesIn &quoteSeries) {
			return calculate(quoteSeries, alpha_, beta_, gamma_* vl_);
		}
	
		static TimeSeriesOut
		calculate(const TimeSeriesIn &quoteSeries, Real alpha, Real beta, Real omega) {
			TimeSeriesOut retval;
			ts_const_iterator cur = quoteSeries.cbegin();
			Real u = cur->second;
			Real sigma2 = u*u;
			retval[cur->first] = std::sqrt(sigma2);
			++cur;
			while (cur != quoteSeries.end()) {
				sigma2 = omega + alpha * u * u + beta * sigma2;
				retval[cur->first] = std::sqrt(sigma2);
				u = cur->second;
				++cur;
			}
			return retval;
		}
	
		Real forecast (Real r, Real sigma2) const {
			return gamma_* vl_ + alpha_ * r * r + beta_ * sigma2;
		}
	
		void calibrate(const TimeSeriesIn &quoteSeries) {
			calibrate (quoteSeries.cbegin_values(), quoteSeries.cend_values());
		}
	
		void calibrate(const TimeSeriesIn &quoteSeries,
			OptimizationMethod &method, 
			const EndCriteria &endCriteria) {
			calibrate (quoteSeries.cbegin_values(), quoteSeries.cend_values(), method, endCriteria);
		}
	
		void calibrate(const TimeSeriesIn &quoteSeries,
			OptimizationMethod &method, 
			const EndCriteria &endCriteria,
			const Array &initialGuess) {
			calibrate (quoteSeries.cbegin_values(), quoteSeries.cend_values(), method, endCriteria, initialGuess);
		}
	
		template<typename Iterator>
		void calibrate (Iterator begin, Iterator end) {
			std::vector<Volatility> r2;
			Real meanr2 = Garch::to_r2(begin, end, r2);
			ProblemPtr p = Garch::calibrate_r2 (mode_, r2, meanr2, alpha_, beta_, vl_);
			gamma_ = 1 - alpha_ - beta_;
			vl_ /= gamma_;
			logLikelihood_ = p ? -p->functionValue() : -costFunction(begin, end);
		}
	
		template<typename Iterator, class Optimizer, class EndCriteria>
		void calibrate (Iterator begin, Iterator end, 
			Optimizer &method,
			EndCriteria endCriteria) {
			std::vector<Volatility> r2;
			Real meanr2 = Garch::to_r2(begin, end, r2);
			ProblemPtr p = Garch::calibrate_r2 (mode_, r2, meanr2, method, endCriteria, alpha_, beta_, vl_);
			gamma_ = 1 - alpha_ - beta_;
			vl_ /= gamma_;
			logLikelihood_ = p ? -p->functionValue() : -costFunction(begin, end);
		}
	
		template<typename Iterator, class Optimizer, class EndCriteria>
		void calibrate (Iterator begin, Iterator end, 
			Optimizer &method,
			EndCriteria endCriteria, 
			const Array &initialGuess) {
			std::vector<Volatility> r2;
			Garch::to_r2(begin, end, r2);
			ProblemPtr p = Garch::calibrate_r2 (r2, method, endCriteria, initialGuess, alpha_, beta_, vl_);
			gamma_ = 1 - alpha_ - beta_;
			vl_ /= gamma_;
			logLikelihood_ = p ? -p->functionValue() : -costFunction(begin, end);
		}
	
		template<class Iterator>
		Real costFunction(Iterator begin, Iterator end, 
			Real alpha = -1, Real beta = -1, Real omega = -1) const {
				if (alpha < 0) alpha = alpha_;
				if (beta < 0) beta = beta_;
				if (omega < 0) omega = vl_ * gamma_;
	
				Real retval(0.0);
				Real u2(0.0), sigma2(0.0);
				Size N = 0;
				for (; begin != end; ++begin, ++N) {
					sigma2 = omega + alpha * u2 + beta * sigma2;
					u2 = *begin; u2 *= u2;
					retval += std::log(sigma2) + u2 / sigma2;
				}
				return N > 0 ? retval / (2*N) : 0.0;
		}

        TimeSeries<Volatility> calculate(const TimeSeries<Volatility>& volatilitySeries) { return TimeSeries<Volatility>(); }
        void calibrate(const TimeSeries<Volatility>& volatilitySeries) {}

	
		Real costFunction(const TimeSeriesIn &quoteSeries, 
			Real alpha = -1, Real beta = -1, Real omega = -1) const {
				return costFunction(quoteSeries.cbegin_values(), quoteSeries.cend_values(), alpha, beta, omega);
		}
	
	private:
		Real alpha_, beta_, gamma_, vl_;
		Real logLikelihood_;
		Natural mode_;
	};
	
	typedef Garch11T<Date> Garch11;

	// GARCH 2D

	class Garch11Diag {
	public:
		Garch11Diag (Size dim) : omega_(dim, 0.0), alpha_(dim, dim, 0.0), beta_(dim, dim, 0.0) {
		}

		Garch11Diag (const Array &omega, const Matrix &alpha, const Matrix &beta) :
			omega_(omega), alpha_(alpha), beta_(beta) {
		}

		const Disposable<Array> forecast (const Array &r, const Array &sigma) const;
		const Array & omega() const { return omega_; }
		const Matrix & alpha() const { return alpha_; }
		const Matrix & beta() const { return beta_; }

		Array & omega() { return omega_; }
		Matrix & alpha() { return alpha_; }
		Matrix & beta() { return beta_; }

	private:
		Array omega_;
		Matrix alpha_, beta_;
	};

	namespace Garch {
		void initGuess (const std::vector<Volatility> &r1, const std::vector<Volatility> &r2,
				Array &params);
		void initGuess (const std::vector<Volatility> &r1, const std::vector<Volatility> &r2,
				Real omega1, Real alpha1, Real beta1, Array &params);
		Real calibrate (const std::vector<Volatility> &r1, const std::vector<Volatility> &r2, Array &params);
		Real calibrate (const std::vector<Volatility> &r1, const std::vector<Volatility> &r2,
				Array &params, OptimizationMethod &method, EndCriteria &endCriteria);
		Disposable<Matrix> forecast (const Array &r, const Matrix &cov, const Garch11Diag &varModel,
				Real omegaCov, Real alphaCov, Real betaCov);
		void params2Model (const Array &params, Garch11Diag &varModel);
		void model2Params (const Garch11Diag &varModel, Array &params);
		Disposable<Matrix> calcCovMtx (const Array &g, Real beta);
	}

	//! GARCH(1,1) 2D volatility model
	/*!
	*/
	template <class Time = Date>
	class Garch2_11 {
	public:
		  typedef Garch11T< Time > Garch11;
		  typedef typename Garch11::TimeSeriesIn TimeSeriesIn;
		  typedef typename TimeSeriesIn::const_iterator ts_const_iterator;
		  typedef typename TimeSeriesIn::value_iterator ts_value_iterator;

		  Garch2_11() : varModel_(2), covModel_()  {
		  }

		  Garch2_11(const Garch11Diag & varModel, const Garch11 & covarModel) :
			  varModel_(varModel), covModel_(covarModel)  {
		  }

		  const Garch11Diag & diagModel() const { return varModel_; }

		  const Garch11 & covarModel() const { return covModel_; }

		  Disposable<Matrix> ltVol() const {
			  Matrix I(2, 2, 0.0);
			  I[0][0] = I[1][1] = 1.0;
			  Array g = inverse(I - varModel_.alpha() - varModel_.beta())*varModel_.omega();
			  Real ltCov = covModel_.ltVol();
			  Real beta = ltCov / g[0];
			  return Garch::calcCovMtx (g, beta);
		  }

		  Real calibrate (const TimeSeriesIn &r1, const TimeSeriesIn &r2) {
			  return calibrate (r1.begin_values(), r1.end_values(), r2.begin_values(), r2.end_values());
		  }

		  template <typename Iterator1, typename Iterator2>
		  Real calibrate (Iterator1 begin1, Iterator1 end1, Iterator2 begin2, Iterator2 end2) {
			  std::vector<Volatility> vr1(std::distance(begin1, end1));
			  std::vector<Volatility> vr2(std::distance(begin2, end2));
			  std::copy (begin1, end1, vr1.begin());
			  std::copy (begin2, end2, vr2.begin());
			  return calibrate(vr1, vr2);
		  }

		  template <typename Iterator1, typename Iterator2>
		  Real calibrate (Iterator1 begin1, Iterator1 end1, Iterator2 begin2, Iterator2 end2,
				  OptimizationMethod &method, EndCriteria &endCriteria) {
			  std::vector<Volatility> vr1(std::distance(begin1, end1));
			  std::vector<Volatility> vr2(std::distance(begin2, end2));
			  std::copy (begin1, end1, vr1.begin());
			  std::copy (begin2, end2, vr2.begin());
			  return calibrate(vr1, vr2, method, endCriteria);
		  }

		  Real calibrate (const std::vector<Volatility> &r1, const std::vector<Volatility> &r2) {
				Array params(11, 0.0);
				Garch::initGuess (r1, r2, params);
				Real res = Garch::calibrate (r1, r2, params);
				Garch::params2Model (params, varModel_);
				covModel_ = Garch11(params[9], params[10], params[8] / (1-params[9]-params[10]));
				return res;
		  }

		  Real calibrate (const std::vector<Volatility> &r1, const std::vector<Volatility> &r2,
				  OptimizationMethod &method, EndCriteria &endCriteria) {
				Array params(11, 0.0);
				Garch::initGuess (r1, r2, params);
				Real res = Garch::calibrate (r1, r2, params, method, endCriteria);
				Garch::params2Model (params, varModel_);
				covModel_ = Garch11(params[9], params[10], params[8] / (1-params[9]-params[10]));
				return res;
		  }

		  template <typename Iterator1, typename Iterator2>
		  Real calibrate (Iterator1 begin1, Iterator1 end1, Iterator2 begin2, Iterator2 end2,
				  const Garch11 &model1) {
			  std::vector<Volatility> vr1(std::distance(begin1, end1));
			  std::vector<Volatility> vr2(std::distance(begin2, end2));
			  std::copy (begin1, end1, vr1.begin());
			  std::copy (begin2, end2, vr2.begin());
			  return calibrate(vr1, vr2, model1);
		  }

		  template <typename Iterator1, typename Iterator2>
		  Real calibrate (Iterator1 begin1, Iterator1 end1, Iterator2 begin2, Iterator2 end2,
				  const Garch11 &model1, OptimizationMethod &method, EndCriteria &endCriteria) {
			  std::vector<Volatility> vr1(std::distance(begin1, end1));
			  std::vector<Volatility> vr2(std::distance(begin2, end2));
			  std::copy (begin1, end1, vr1.begin());
			  std::copy (begin2, end2, vr2.begin());
			  return calibrate(vr1, vr2, model1, method, endCriteria);
		  }

		  Real calibrate (const std::vector<Volatility> &r1, const std::vector<Volatility> &r2,
				  const Garch11 &model1) {
				Array params(11, 0.0);
				Garch::initGuess (r1, r2, model1.omega(), model1.alpha(), model1.beta(), params);
				Real res = Garch::calibrate (r1, r2, params);
				Garch::params2Model (params, varModel_);
				covModel_ = Garch11(params[9], params[10], params[8] / (1-params[9]-params[10]));
				return res;
		  }

		  Real calibrate (const std::vector<Volatility> &r1, const std::vector<Volatility> &r2,
				  const Garch11 &model1, OptimizationMethod &method, EndCriteria &endCriteria) {
				Array params(11, 0.0);
				Garch::initGuess (r1, r2, model1.omega(), model1.alpha(), model1.beta(), params);
				Real res = Garch::calibrate (r1, r2, params, method, endCriteria);
				Garch::params2Model (params, varModel_);
				covModel_ = Garch11(params[9], params[10], params[8] / (1-params[9]-params[10]));
				return res;
		  }

		  Disposable<Matrix> forecast (const Array &r, const Matrix &cov) const {
			  return Garch::forecast(r, cov, varModel_, covModel_.omega(), covModel_.alpha(), covModel_.beta());
		  }

	private:
		Garch11Diag varModel_;
		Garch11 covModel_;
	};

} // namespace QuantLib


#endif
