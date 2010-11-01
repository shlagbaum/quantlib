/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006 Joseph Wang
 Copyright (C) 2010 Liquidnet Holdings

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

#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <ql/models/volatility/garch.hpp>
#include <ql/math/array.hpp>
#include <ql/math/optimization/costfunction.hpp>
#include <ql/math/optimization/problem.hpp>
#include <ql/math/optimization/simplex.hpp>
#include <ql/math/optimization/leastsquare.hpp>
#include <ql/math/linearleastsquaresregression.hpp>
#include <ql/math/functional.hpp>
#include <ql/math/solvers1d/brent.hpp>
#include <ql/experimental/math/autocovariance.hpp>
#include <ql/errors.hpp>
#include <ql/qldefines.hpp>

namespace QuantLib { namespace Garch {

	const Real tol_level = 1.0e-8;

	class Garch11CostFunction : public CostFunction {
	public:
		Garch11CostFunction (const std::vector<Volatility> &);
		virtual Real value(const Array& x) const;
		virtual Disposable<Array> values(const Array& x) const;
		virtual void gradient(Array& grad, const Array& x) const;
		virtual Real valueAndGradient(Array& grad, const Array& x) const;
	private:
		const std::vector<Volatility> &r2_;
	};

	class Garch11Constraint : public Constraint {
	  private:
		class Impl : public Constraint::Impl {
			Real gammaLower_, gammaUpper_;
		  public:
			  Impl (Real gammaLower, Real gammaUpper) : gammaLower_(gammaLower), gammaUpper_(gammaUpper) {
			  }
			bool test(const Array &x) const {
		      QL_REQUIRE(x.size() >= 3, "size of parameters vector < 3");
		      return x[0] > 0 && x[1] >= 0 && x[2] >= 0 && x[1] + x[2] < gammaUpper_ && x[1] + x[2] >= gammaLower_;
			}
		};
	  public:
		  Garch11Constraint(Real gammaLower, Real gammaUpper) : Constraint(boost::shared_ptr<Constraint::Impl>(
												   new Garch11Constraint::Impl(gammaLower, gammaUpper))) {}
	};

  Garch11CostFunction::Garch11CostFunction (const std::vector<Volatility> &r2) : r2_(r2) {

  }

  Real Garch11CostFunction::value(const Array& x) const {
  	  Real retval(0.0);
  	  Real sigma2 = 0;
  	  Real u2 = 0;
  	  BOOST_FOREACH (Volatility r2, r2_) {
  		  sigma2 = x[0] + x[1] * u2 + x[2] * sigma2;
  		  u2 = r2;
  		  retval += std::log(sigma2) + u2 / sigma2;
  	  }
  	  return retval / (2.0*r2_.size());
  }

  Disposable<Array> Garch11CostFunction::values(const Array& x) const {
      Array retval (r2_.size());
  	  Real sigma2 = 0;
  	  Real u2 = 0;
  	  Size i = 0;
  	  BOOST_FOREACH (Volatility r2, r2_) {
  		  sigma2 = x[0] + x[1] * u2 + x[2] * sigma2;
  		  u2 = r2;
  		  retval[i++] = (std::log(sigma2) + u2 / sigma2)/(2.0*r2_.size());
  	  }
  	  return retval;
  }

  void Garch11CostFunction::gradient(Array& grad, const Array& x) const {
  	  std::fill (grad.begin(), grad.end(), 0.0);
  	  Real sigma2 = 0;
  	  Real u2 = 0;
  	  Real sigma2prev = sigma2;
  	  Real u2prev = u2;
  	  Real norm = 2.0 * r2_.size();
  	  BOOST_FOREACH (Volatility r2, r2_) {
  	    sigma2 = x[0] + x[1] * u2 + x[2] * sigma2;
  	    u2 = r2;
  	    Real w = (sigma2 - u2) / (sigma2*sigma2);
  	    grad[0] += w;
  	    grad[1] += u2prev * w;
  	    grad[2] += sigma2prev * w;
  	    u2prev = u2;
  	    sigma2prev = sigma2;
  	  }
  	  std::transform(grad.begin(), grad.end(), grad.begin(),
  			  std::bind2nd(std::divides<Real>(), norm));
  }

  Real Garch11CostFunction::valueAndGradient(Array& grad, const Array& x) const {
  	  std::fill (grad.begin(), grad.end(), 0.0);
  	  Real retval(0.0);
  	  Real sigma2 = 0;
  	  Real u2 = 0;
  	  Real sigma2prev = sigma2;
  	  Real u2prev = u2;
  	  Real norm = 2.0 * r2_.size();
  	  BOOST_FOREACH (Volatility r2, r2_) {
  	    sigma2 = x[0] + x[1] * u2 + x[2] * sigma2;
  	    u2 = r2;
  	    retval += std::log(sigma2) + u2 / sigma2;
  	    Real w = (sigma2 - u2) / (sigma2*sigma2);
  	    grad[0] += w;
  	    grad[1] += u2prev * w;
  	    grad[2] += sigma2prev * w;
  	    u2prev = u2;
  	    sigma2prev = sigma2;
  	  }
  	  std::transform(grad.begin(), grad.end(), grad.begin(),
  			  std::bind2nd(std::divides<Real>(), norm));
  	  return retval / norm;
  }

  class FitAcfProblem : public LeastSquareProblem {
  public:
  	FitAcfProblem(Real A2, const Array &acf, const std::vector<size_t> &idx);
  	virtual Size size();
  	virtual void targetAndValue(const Array& x, Array& target, Array& fct2fit);
  	virtual void targetValueAndGradient(const Array& x, Matrix& grad_fct2fit, Array& target, Array& fct2fit);
  private:
  	Real A2_;
  	Array acf_;
  	std::vector<size_t> idx_;
  };

  class FitAcfConstraint : public Constraint {
    private:
  	class Impl : public Constraint::Impl {
  		Real gammaLower_, gammaUpper_;
  	  public:
  		Impl(Real gammaLower, Real gammaUpper) : gammaLower_(gammaLower), gammaUpper_(gammaUpper) {

  		}
  		bool test(const Array &x) const {
  	      QL_REQUIRE(x.size() >= 2, "size of parameters vector < 2");
  	      return x[0] >= gammaLower_ && x[0] < gammaUpper_ && x[1] >= 0 && x[1] <= x[0];
  		}
  	};
    public:
  	  FitAcfConstraint(Real gammaLower, Real gammaUpper) :
  		  Constraint(boost::shared_ptr<Constraint::Impl>(
  			new FitAcfConstraint::Impl(gammaLower, gammaUpper))) {}
  };

  Real fGamma (Real gamma, Real A, Real B) {
  	Real beta = gamma * (1 - A) - B;
  	return 3*A*(1 - gamma*gamma) - 1 + 3*gamma*gamma + 2*beta*beta - 4*beta*gamma;
  }

  // Initial guess based on fitting ACF - initial guess for fitting acf is
  // a moment matching estimates for mean(r2), acf(0), and acf(1).
  Real initialGuess1 (const Array &acf, Real meanr2, Real &alpha, Real &beta, Real &omega) {
    Real A21 = acf[1];
    Real A4 = acf[0] + meanr2*meanr2;

    Real A = meanr2*meanr2/A4; // 1/sigma^2
    Real B = A21 / A4; // rho(1)

    Real gammaLower = A <= 1./3. - tol_level ? std::sqrt((1 - 3*A)/(3 - 3*A)) + tol_level : tol_level;
	Garch11Constraint constraints(gammaLower, 1.0 - tol_level);

    //Real fGamma0 = fGamma (0, B4, B21);
    //Real fGammaLower = fGamma (gammaLower, B4, B21);
    //Real fGamma1 = fGamma (1.0, B4, B21);
    Real gamma = gammaLower + (1 - gammaLower) * 0.5;
    beta = std::min(gamma, std::max(gamma * (1 - A) - B, 0.0));
    alpha = gamma - beta;
    omega = meanr2 * (1 - gamma);

    if (std::fabs(A-0.5) < QL_EPSILON) {
    	gamma = std::max(gammaLower, -(1+4*B*B)/(4*B));
    	beta = std::min(gamma, std::max(gamma * (1 - A) - B, 0.0));
        alpha = gamma - beta;
        omega = meanr2 * (1 - gamma);
    } else
    if (A > 1.0 - QL_EPSILON) {
    	gamma = std::max(gammaLower, -(1+B*B)/(2*B));
    	beta = std::min(gamma, std::max(gamma * (1 - A) - B, 0.0));
        alpha = gamma - beta;
        omega = meanr2 * (1 - gamma);
    } else {
		Real D = (3*A-1)*(2*B*B+(1-A)*(2*A-1));
		if (D >= 0) {
			Real d = std::sqrt(D);
			Real b = (B - d)/(2*A-1);
			Real g = 0;
			if (b >= tol_level && b <= 1.0 - tol_level) {
				g = (b + B) / (1 - A);
			}
			if (g < gammaLower) {
				b = (B + d)/(2*A-1);
				if (b >= tol_level && b <= 1.0 - tol_level) {
					g = (b + B) / (1 - A);
				}
			}
			if (g >= gammaLower) {
				gamma = g;
		    	beta = std::min(gamma, std::max(gamma * (1 - A) - B, 0.0));
		        alpha = gamma - beta;
		        omega = meanr2 * (1 - gamma);
			}
		}
    }

    std::vector<size_t> idx;
    size_t nCov = acf.size() - 1;
    for (size_t i = 0; i <= nCov; ++i) {
  	  if (i < 2 || i > 1 && acf[i] > 0 && acf[i-1] > 0 && acf[i-1] > acf[i]) {
  		  idx.push_back(i);
  	  }
    }

    Array x(2);
    x[0] = gamma;
    x[1] = beta;

    try {
		FitAcfConstraint c(gammaLower, 1.0 - tol_level);
		NonLinearLeastSquare nnls(c);
		nnls.setInitialValue(x);
		FitAcfProblem pr(meanr2, acf, idx);
		x = nnls.perform(pr);
		Array guess(3);
		guess[0] = meanr2 * (1 - x[0]);
		guess[1] = x[0] - x[1];
		guess[2] = x[1];
		if (constraints.test(guess)) {
			omega = guess[0];
			alpha = guess[1];
			beta = guess[2];
		}
    } catch (const std::exception &) {
    	// failed -- returning initial values
    }
    return gammaLower;
  }

  // Initial guess based on fitting ACF - initial guess for fitting acf is
  // an estimate of gamma = alpfa+beta based on a property: acf(i+1) = gamma*acf(i) for i > 1.
  Real initialGuess2 (const Array &acf, Real meanr2, Real &alpha, Real &beta, Real &omega) {
    Real A21 = acf[1];
    Real A4 = acf[0] + meanr2*meanr2;
    Real A = meanr2*meanr2/A4; // 1/sigma^2
    Real B = A21 / A4; // rho(1)
    Real gammaLower = A <= 1./3. - tol_level ? std::sqrt((1 - 3*A)/(3 - 3*A)) + tol_level : tol_level;
	Garch11Constraint constraints(gammaLower, 1.0 - tol_level);

    // ACF
    Real gamma = 0;
    size_t nn = 0;
    std::vector<size_t> idx;
    size_t nCov = acf.size() - 1;
    for (size_t i = 0; i <= nCov; ++i) {
  	  if (i < 2) idx.push_back(i);
  	  if (i > 1 && acf[i] > 0 && acf[i-1] > 0 && acf[i-1] > acf[i]) {
  		  gamma += acf[i]/acf[i-1];
  		  nn++;
  		  idx.push_back(i);
  	  }
    }
    if (nn > 0)
    	gamma /= nn;
    if (gamma < gammaLower) gamma = gammaLower;
    beta = std::min(gamma, std::max(gamma * (1 - A) - B, 0.0));
    omega = meanr2 * (1 - gamma);

    Array x(2);
    x[0] = gamma;
    x[1] = beta;

    try {
		FitAcfConstraint c(gammaLower, 1 - tol_level);
		NonLinearLeastSquare nnls(c);
		nnls.setInitialValue(x);
		FitAcfProblem pr(meanr2, acf, idx);
		x = nnls.perform(pr);
		Array guess(3);
		guess[0] = meanr2 * (1 - x[0]);
		guess[1] = x[0] - x[1];
		guess[2] = x[1];
		if (constraints.test(guess)) {
			omega = guess[0];
			alpha = guess[1];
			beta = guess[2];
		}
    } catch (const std::exception &) {
    	// failed -- returning initial values
    }
    return gammaLower;
  }

   ProblemPtr calibrate_r2 (Natural mode, const std::vector<Volatility> &r2, Real meanr2,
		  Real &alpha, Real &beta, Real &omega) {
	EndCriteria endCriteria(10000, 500, tol_level, tol_level, tol_level);
	Simplex method(0.001);
    return calibrate_r2 (mode, r2, meanr2, method, endCriteria, alpha, beta, omega);
  }

ProblemPtr calibrate_r2 (const std::vector<Volatility> &r2,
	  OptimizationMethod &method,
	  Constraint &constraints,
	  const EndCriteria &endCriteria,
	  const Array &initGuess, Real &alpha, Real &beta, Real &omega);

ProblemPtr calibrate_r2 (Natural mode, const std::vector<Volatility> &r2, Real meanr2,
		  OptimizationMethod &method, const EndCriteria &endCriteria,
		  Real &alpha, Real &beta, Real &omega) {
	Real dataSize = Real(r2.size());
	alpha = 0.0;
	beta = 0.0;
	omega = 0.0;
	QL_REQUIRE (dataSize >= 4, "Data series is too short to fit GARCH model");
	QL_REQUIRE (meanr2 > 0, "Data series is constant");
	omega = meanr2 * dataSize / (dataSize - 1);
	// ACF
	Size maxLag = (Size)std::sqrt(dataSize);
	Array acf(maxLag+1);
	std::vector<Volatility> tmp(r2.size());
	std::transform (r2.begin(), r2.end(), tmp.begin(), std::bind2nd(std::minus<double>(), meanr2));
	autocovariances (tmp.begin(), tmp.end(), acf.begin(), maxLag);
	QL_REQUIRE (acf[0] > 0, "Data series is constant");

	Garch11CostFunction cost (r2);

	// 2 initial guesses based on fitting ACF
	Real gammaLower = 0.0;
	Array opt1(3);
	Real fCost1 = QL_MAX_REAL;
	if (0 <= mode && mode != 1) {
		gammaLower = initialGuess1 (acf, meanr2, opt1[1], opt1[2], opt1[0]);
		fCost1 = cost.value(opt1);
	}

	Array opt2(3);
	Real fCost2 = QL_MAX_REAL;
	if (1 <= mode) {
		gammaLower = initialGuess2 (acf, meanr2, opt2[1], opt2[2], opt2[0]);
		fCost2 = cost.value(opt2);
	}

	Garch11Constraint constraints(gammaLower, 1.0 - tol_level);

	ProblemPtr ret;
	if (mode <= 2) {
		try {
			ret = calibrate_r2 (r2, method, constraints, endCriteria, fCost1 <= fCost2 ? opt1 : opt2, alpha, beta, omega);
		} catch (const std::exception &) {
			if (fCost1 <= fCost2) {
				alpha = opt1[1];
				beta = opt1[2];
				omega = opt1[0];
			} else {
				alpha = opt2[1];
				beta = opt2[2];
				omega = opt2[0];
			}
		}
	} else {
		ProblemPtr ret1, ret2;
		try {
			ret1 = calibrate_r2 (r2, method, constraints, endCriteria, opt1, alpha, beta, omega);
			opt1[1] = alpha;
			opt1[2] = beta;
			opt1[0] = omega;
			double fCost = QL_MAX_REAL;
			if (constraints.test(opt1) && (fCost = cost.value(opt1)) < fCost1)
				fCost1 = fCost;
		} catch (const std::exception &) {
			fCost1 = QL_MAX_REAL;
		}

		try {
			ret2 = calibrate_r2 (r2, method, constraints, endCriteria, opt2, alpha, beta, omega);
			opt2[1] = alpha;
			opt2[2] = beta;
			opt2[0] = omega;
			double fCost = QL_MAX_REAL;
			if (constraints.test(opt2) && (fCost = cost.value(opt2)) < fCost2)
				fCost2 = fCost;
		} catch (const std::exception &) {
			fCost2 = QL_MAX_REAL;
		}

		if (fCost1 <= fCost2) {
			alpha = opt1[1];
			beta = opt1[2];
			omega = opt1[0];
			ret = ret1;
		} else {
			alpha = opt2[1];
			beta = opt2[2];
			omega = opt2[0];
			ret = ret2;
		}
	}
	return ret;
  }

  ProblemPtr calibrate_r2 (const std::vector<Volatility> &r2,
		  OptimizationMethod &method,
		  Constraint &constraints,
		  const EndCriteria &endCriteria,
		  const Array &initGuess, Real &alpha, Real &beta, Real &omega) {
	  Garch11CostFunction cost(r2);
	  ProblemPtr problem = ProblemPtr(new Problem(cost, constraints, initGuess));
	  EndCriteria::Type ret = method.minimize(*problem, endCriteria);
	  const Array &optimum = problem->currentValue();
	  alpha = optimum[1];
	  beta = optimum[2];
	  omega = optimum[0];
	  return problem;
  }

  ProblemPtr calibrate_r2 (const std::vector<Volatility> &r2, Real meanr2,
		  OptimizationMethod &method,
		  Constraint &constraints,
		  const EndCriteria &endCriteria,
		  const Array &initGuess, Real &alpha, Real &beta, Real &omega) {
		std::vector<Volatility> tmp(r2.size());
		std::transform (r2.begin(), r2.end(), tmp.begin(), std::bind2nd(std::minus<double>(), meanr2));
		return calibrate_r2 (tmp, method, constraints, endCriteria, initGuess, alpha, beta, omega);
  }

  ProblemPtr calibrate_r2 (const std::vector<Volatility> &r2, Real meanr2,
		  OptimizationMethod &method,
		  const EndCriteria &endCriteria,
		  const Array &initGuess, Real &alpha, Real &beta, Real &omega) {
		Garch11Constraint constraints(0.0, 1.0 - tol_level);
		std::vector<Volatility> tmp(r2.size());
		std::transform (r2.begin(), r2.end(), tmp.begin(), std::bind2nd(std::minus<double>(), meanr2));
		return calibrate_r2 (tmp, method, constraints, endCriteria, initGuess, alpha, beta, omega);
  }

  // class FitAcfProblem : public LeastSquareProblem {
  FitAcfProblem::FitAcfProblem(Real A2, const Array &acf, const std::vector<size_t> &idx) :
  	A2_(A2), acf_(acf), idx_(idx) {
  }

  Size FitAcfProblem::size() { return idx_.size(); }

  void FitAcfProblem::targetAndValue(const Array& x, Array& target, Array& fct2fit) {
  	Real A4 = acf_[0] + A2_*A2_;
  	Real gamma = x[0];
  	Real beta = x[1];
  	target[0] = A2_*A2_/A4;
  	fct2fit[0] = (1 - 3*gamma*gamma - 2*beta*beta + 4*beta*gamma) / (3*(1 - gamma*gamma));
  	target[1] = acf_[1] / A4;
  	fct2fit[1] = gamma * (1 - fct2fit[0]) - beta;
  	for (size_t i = 2; i < idx_.size(); ++i) {
  		target[i] = acf_[idx_[i]] / A4;
  		fct2fit[i] = std::pow(gamma, (int)idx_[i]-1)* fct2fit[1];
  	}
  }

  void FitAcfProblem::targetValueAndGradient(const Array& x, Matrix& grad_fct2fit, Array& target, Array& fct2fit) {
  	Real A4 = acf_[0] + A2_*A2_;
  	Real gamma = x[0];
  	Real beta = x[1];
  	target[0] = A2_*A2_/A4;
  	Real w1 = (1 - 3*gamma*gamma - 2*beta*beta + 4*beta*gamma);
  	Real w2 = (1 - gamma*gamma);
  	fct2fit[0] = w1 / (3*w2);
  	grad_fct2fit[0][0] = (2.0/3.0) * ((2*beta-3*gamma)*w2 + 2*w1*gamma) / (w2*w2);
  	grad_fct2fit[0][1] = (4.0/3.0) * (gamma - beta) / w2;
  	target[1] = acf_[1] / A4;
  	fct2fit[1] = gamma * (1 - fct2fit[0]) - beta;
  	grad_fct2fit[1][0] = (1 - fct2fit[0]) - gamma * grad_fct2fit[0][0];
  	grad_fct2fit[1][1] = -gamma * grad_fct2fit[0][1] - 1;
  	for (size_t i = 2; i < idx_.size(); ++i) {
  		target[i] = acf_[idx_[i]] / A4;
  		w1 = std::pow(gamma, (int)idx_[i]-1);
  		fct2fit[i] = w1 * fct2fit[1];
  		grad_fct2fit[i][0] = (idx_[i]-1) * (w1/gamma)*fct2fit[1] + w1*grad_fct2fit[1][0];
  		grad_fct2fit[i][1] = w1 * grad_fct2fit[1][1];
  	}
  }

// GARCH(1, 1) 2D
const Real tol_level2 = 1.0e-6;

class Garch2_11CostFunction : public CostFunction {
public:
	Garch2_11CostFunction (const std::vector<Volatility> &r1, const std::vector<Volatility> &r2,
			Garch11Diag &varModel);
	virtual Real value(const Array& x) const;
	virtual Disposable<Array> values(const Array& x) const;
	virtual void gradient(Array& grad, const Array& x) const;
	virtual Real valueAndGradient(Array& grad, const Array& x) const;
private:
	friend class Garch2_11Constraint;
	const std::vector<Volatility> &r1_;
	const std::vector<Volatility> &r2_;
	Garch11Diag &varModel_;
	Size size_;
};

class Garch2_11Constraint : public Constraint {
  private:
	class Impl : public Constraint::Impl {
		Garch11Diag &varModel_;
	  public:
		  Impl (Garch11Diag &varModel) : varModel_(varModel) {
		  }
		bool test(const Array &x) const {
	      params2Model (x, varModel_);
	      const Array &omega = varModel_.omega();
	      const Matrix &alpha = varModel_.alpha();
	      const Matrix &beta = varModel_.beta();
	      Real gamma11 = alpha[0][0] + beta[0][0];
	      Real gamma21 = alpha[1][0] + beta[1][0];
	      Real gamma22 = alpha[1][1] + beta[1][1];
	      Real gamma12 = x[9] + x[10];
	      return
			  omega[0] > 0 && omega[1] > 0
			  // && alpha[0][0] >= 0 && alpha[0][0] < 1
			  // && alpha[1][1] >= 0 && alpha[1][1] < 1
			  // && beta[0][0] >= 0 && beta[0][0] < 1
			  // && beta[1][1] >= 0 && beta[1][1] < 1
			  // && x[9] >= 0 && x[9] < 1
			  // && x[10] >= 0 && x[10] < 1
			  && gamma11 >= 0 && gamma11 < 1
			  && gamma22 >= 0 && gamma22 < 1
			  && gamma21 > -omega[1]*(1.0-gamma11)/omega[0] && gamma21 < 1
			  && gamma12 >= 0 && gamma12 < 1;
		}
	};
  public:
	  Garch2_11Constraint(Garch11Diag &varModel) : Constraint(boost::shared_ptr<Constraint::Impl>(
							new Garch2_11Constraint::Impl(varModel))) {}
};

void initGuess (const std::vector<Volatility> &r1, const std::vector<Volatility> &r2, Array &params) {
	Size n = std::min(r1.size(), r2.size());
	std::vector<Volatility> vr1(n), vr2(n);
	std::vector<boost::function1<Volatility, Volatility> > f;
    f.push_back(identity<Volatility>());

    LinearLeastSquaresRegression<> lr(r1, r2, f);
    Volatility bta = lr.coefficients()[0];

	Real mean1 = 0;
	Real mean2 = 0;
	Real w = 1.0, u1 = 0, u2 = 0;
	std::vector<Volatility>::const_iterator itr1 = r1.begin(), itr2 = r2.begin();
	std::vector<Volatility>::iterator itvr1 = vr1.begin(), itvr2 = vr2.begin();
	for (Size i = 0; i < n; ++i, ++itr1, ++itr2, ++itvr1, ++itvr2) {
		u1 = *itr1;
		u2 = *itr2 - bta*u1;
		u1 *= u1;
		u2 *= u2;
		mean1 = (1.0 - w) * mean1 + w * u1;
		mean2 = (1.0 - w) * mean2 + w * u2;
		*itvr1 = u1;
		*itvr2 = u2;
		w /= (w + 1.0);
	}

	Garch11Diag varModel(2);
	Array &omega = varModel.omega();
	Matrix &alpha = varModel.alpha();
	Matrix &beta = varModel.beta();

	try {
		calibrate_r2 (2, vr1, mean1, alpha[0][0], beta[0][0], omega[0]);
		calibrate_r2 (2, vr2, mean2, alpha[1][1], beta[1][1], omega[1]);
	} catch (const std::exception &) {
		// TODO: what to do?
	}
	model2Params (varModel, params);
	params[8] = bta * omega[0];
	params[9] = alpha[0][0];
	params[10] = beta[0][0];
}

void initGuess (const std::vector<Volatility> &r1, const std::vector<Volatility> &r2,
		Real omega1, Real alpha1, Real beta1, Array &params) {
	Size n = std::min(r1.size(), r2.size());
	std::vector<Volatility> vr2(n);
	std::vector<boost::function1<Volatility, Volatility> > f;
    f.push_back(identity<Volatility>());

    LinearLeastSquaresRegression<> lr(r1, r2, f);
    Volatility bta = lr.coefficients()[0];

	Real mean2 = 0;
	Real w = 1.0, u2 = 0;
	std::vector<Volatility>::const_iterator itr1 = r1.begin(), itr2 = r2.begin();
	std::vector<Volatility>::iterator itvr2 = vr2.begin();
	for (Size i = 0; i < n; ++i, ++itr1, ++itr2, ++itvr2) {
		u2 = *itr2 - bta*(*itr1);
		u2 *= u2;
		mean2 = (1.0 - w) * mean2 + w * u2;
		*itvr2 = u2;
		w /= (w + 1.0);
	}

	Garch11Diag varModel(2);
	Array &omega = varModel.omega();
	Matrix &alpha = varModel.alpha();
	Matrix &beta = varModel.beta();

	omega[0] = omega1;
	alpha[0][0] = alpha1;
	beta[0][0] = beta1;

	try {
		calibrate_r2 (2, vr2, mean2, alpha[1][1], beta[1][1], omega[1]);
	} catch (const std::exception &) {
	}
	model2Params (varModel, params);
	params[8] = bta * omega[0];
	params[9] = alpha[0][0];
	params[10] = beta[0][0];
}

ProblemPtr calibrate (const std::vector<Volatility> &r1, const std::vector<Volatility> &r2, Array &params) {
	EndCriteria endCriteria(10000, 500, tol_level2, tol_level2, tol_level2);
	Simplex method(0.001);

	return calibrate (r1, r2, params, method, endCriteria);
}

ProblemPtr calibrate (const std::vector<Volatility> &r1, const std::vector<Volatility> &r2,
		Array &params, OptimizationMethod &method, EndCriteria &endCriteria) {
	ProblemPtr res;
	Garch11Diag varModel(2);

	try {
		Garch2_11CostFunction cost (r1, r2, varModel);
		Garch2_11Constraint constraints(varModel);
		res = ProblemPtr( new Problem(cost, constraints, params));
		method.minimize(*res, endCriteria);
		const Array &optimum = res->currentValue();
		std::copy (optimum.begin(), optimum.end(), params.begin());
	} catch (const std::exception &) {
	}
	return res;
}

void params2Model (const Array &params, Garch11Diag &varModel) {
	Array &omega = varModel.omega();
	Matrix &alpha = varModel.alpha();
	Matrix &beta = varModel.beta();

	omega[0] = params[0];
	omega[1] = params[1];

	alpha[0][0] = params[2];
	alpha[0][1] = 0;
	alpha[1][0] = params[3];
	alpha[1][1] = params[4];

	beta[0][0] = params[5];
	beta[0][1] = 0;
	beta[1][0] = params[6];
	beta[1][1] = params[7];
}

void model2Params (const Garch11Diag &varModel, Array &params) {
	const Array &omega = varModel.omega();
	const Matrix &alpha = varModel.alpha();
	const Matrix &beta = varModel.beta();

	params[0] = omega[0];
	params[1] = omega[1];

	params[2] = alpha[0][0];
	params[3] = alpha[1][0];
	params[4] = alpha[1][1];

	params[5] = beta[0][0];
	params[6] = beta[1][0];
	params[7] = beta[1][1];
}

Disposable<Matrix> calcCovMtx (const Array &g, Real beta) {
	  Matrix res(2, 2);
	  res[0][0] = g[0];
	  res[1][1] = g[1] + beta*beta*g[0];
	  res[1][0] = res[0][1] = beta*g[0];
	  return res;
}

Garch2_11CostFunction::Garch2_11CostFunction (const std::vector<Volatility> &r1, const std::vector<Volatility> &r2,
		Garch11Diag &varModel) : r1_(r1), r2_(r2), varModel_(varModel),
		size_(std::min(r1.size(), r2.size())) {

}

Disposable<Matrix> forecast (const Array &r, const Matrix &cov, const Garch11Diag &varModel,
		Real omegaCov, Real alphaCov, Real betaCov, Real maxBeta)
{
	  Real beta = maxBeta > 0 ? cov[0][0]*maxBeta > std::fabs(cov[1][0]) ? cov[1][0] / cov[0][0] : 0 : cov[0][0] > 0 ? cov[1][0] / cov[0][0] : 0;
	  Array g(2);
	  g[0] = cov[0][0];
	  g[1] = cov[1][1] - beta*cov[1][0];
	  Array b2(2);
	  b2[0] = r[0]*r[0];
	  b2[1] = r[1] - beta * r[0];
	  b2[1] *= b2[1];
	  //std::cerr << cov[0][0] << "," << cov[1][0] << ", beta = " << beta << ", g = " << g[0] << "," << g[1] << ", b2 = " << b2[0] << "," << b2[1] << ",";
	  g = varModel.forecast (b2, g);
	  if (g[1] < QL_EPSILON)
		  g[1] = QL_EPSILON;
	  Real sigma12 = omegaCov + alphaCov*r[0]*r[1] + betaCov*cov[1][0];
	  Real sigma11 = varModel.omega()[0] + varModel.alpha()[0][0]*r[0]*r[0] + varModel.beta()[0][0]*cov[0][0];
	  beta = maxBeta > 0 ? sigma11*maxBeta > std::fabs(sigma12) ? sigma12 / sigma11 : 0 : sigma11 > 0 ? sigma12 / sigma11 : 0;
	  //std::cerr << sigma11 << "," << sigma12 << ", beta = " << beta << ", g = " << g[0] << "," << g[1] << std::endl;
	  return calcCovMtx (g, beta);
}

Real Garch2_11CostFunction::value(const Array& x) const {
	Real retval(0.0);
	params2Model (x, varModel_);
	Array g(2, 0.0), b2(2, 0.0);
	Real sigma21 = 0;
	Real omega21 = x[8];
	Real alpha21 = x[9];
	Real beta21 = x[10];
	Real q21 = 0, g11 = 0, g22 = 0, r1, r2;

	std::vector<Volatility>::const_iterator itr1 = r1_.begin(), itr2 = r2_.begin();
	for (Size t = 0; t < size_; ++t, ++itr1, ++itr2) {
		r1 = *itr1;
		r2 = *itr2;
		sigma21 = omega21 + alpha21 * r1*r2 + beta21 * sigma21;
		b2[0] = r1;
		b2[1] = r2 - q21 * r1;
		b2[0] *= b2[0];
		b2[1] *= b2[1];
		g = varModel_.forecast (b2, g);

		g11 = g[0];
		q21 = sigma21 / g11;

		if (g[1] > QL_EPSILON) {
			g22 = g[1];
			//rho = std::sqrt (g11*q21*q21/(g22 + g11*q21*q21));
			//if (q21 < 0) rho = -rho;
		} else {
			g22 = QL_EPSILON;
		}

		retval += std::log (g11*g22) + b2[0] / g11 + b2[1] / g22;
	}

	return 0.5 * retval / size_;
}

Disposable<Array> Garch2_11CostFunction::values(const Array& x) const {
    Array retval (size_);
	params2Model (x, varModel_);
	Array g(2, 0.0), b2(2, 0.0);
	Real sigma21 = 0;
	Real omega21 = x[8];
	Real alpha21 = x[9];
	Real beta21 = x[10];
	Real q21 = 0, g11 = 0, g22 = 0, rho = 0;

	for (Size t = 0; t < size_; ++t) {
		sigma21 = omega21 + alpha21 * r1_[t]*r2_[t] + beta21 * sigma21;
		b2[0] = r1_[t];
		b2[1] = r2_[t] - q21 * r1_[t];
		b2[0] *= b2[0];
		b2[1] *= b2[1];
		g = varModel_.forecast (b2, g);

		g11 = g[0];
		q21 = sigma21 / g11;

		if (g[1] > QL_EPSILON) {
			g22 = g[1];
			rho = std::sqrt (g11*q21*q21/(g22 + g11*q21*q21));
			if (q21 < 0) rho = -rho;
		} else {
			g22 = QL_EPSILON;
		}

		retval[t] = std::log (g11) + std::log (g22) + b2[0] / g11 + b2[1] / g22;
	}

	return retval;
}

void Garch2_11CostFunction::gradient(Array& grad, const Array& x) const {
	std::fill (grad.begin(), grad.end(), 0);
}

Real Garch2_11CostFunction::valueAndGradient(Array& grad, const Array& x) const {
	gradient(grad, x);
	return value (x);
}

} // namespace Garch 

const Disposable<Array> Garch11Diag::forecast (const Array &r2, const Array &sigma2) const {
	Array res = omega_ + alpha_*r2 + beta_*sigma2;
	return res;
}

} // namespace QuantLib


