/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2009 Andreas Gaida
 Copyright (C) 2009 Ralph Schreyer
 Copyright (C) 2009 Klaus Spanderen

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

/*! \file fdmbackwardsolver.hpp
*/

#ifndef quantlib_fdm_backward_solver_hpp
#define quantlib_fdm_backward_solver_hpp

#include <ql/experimental/finitedifferences/fdmdirichletboundary.hpp>

namespace QuantLib {

    class FdmLinearOpComposite;
    class FdmStepConditionComposite;

    class FdmBackwardSolver {
      public:
        enum FdmSchemeType { Hundsdorfer, Douglas, 
                             CraigSneyd, ModifiedCraigSneyd, 
                             ImplicitEuler, ExplicitEuler };
        typedef FdmLinearOp::array_type array_type;
        
        FdmBackwardSolver(
          const boost::shared_ptr<FdmLinearOpComposite>& map,
          const FdmBoundaryConditionSet& bcSet,
          const boost::shared_ptr<FdmStepConditionComposite> condition,
          FdmSchemeType schemeType = Hundsdorfer,
          Real theta = 0.5+std::sqrt(3.0)/6,
          Real mu    = 0.5);

        void rollback(array_type& a, 
                      Time from, Time to,
                      Size steps, Size dampingSteps);

      protected:
        const boost::shared_ptr<FdmLinearOpComposite> map_;
        const FdmBoundaryConditionSet bcSet_;
        const boost::shared_ptr<FdmStepConditionComposite> condition_;
        const FdmSchemeType schemeType_;
        const Real theta_, mu_;
    };
}

#endif
