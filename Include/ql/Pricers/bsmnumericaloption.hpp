
/*
 * Copyright (C) 2000-2001 QuantLib Group
 *
 * This file is part of QuantLib.
 * QuantLib is a C++ open source library for financial quantitative
 * analysts and developers --- http://quantlib.sourceforge.net/
 *
 * QuantLib is free software and you are allowed to use, copy, modify, merge,
 * publish, distribute, and/or sell copies of it under the conditions stated
 * in the QuantLib License.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the license for more details.
 *
 * You should have received a copy of the license along with this file;
 * if not, contact ferdinando@ametrano.net
 * The license is also available at http://quantlib.sourceforge.net/LICENSE.TXT
 *
 * The members of the QuantLib Group are listed in the Authors.txt file, also
 * available at http://quantlib.sourceforge.net/Authors.txt
*/

/*! \file bsmnumericaloption.hpp
    \brief common code for numerical option evaluation

    $Source$
    $Log$
    Revision 1.6  2001/05/24 13:57:51  nando
    smoothing #include xx.hpp and cutting old Log messages

    Revision 1.5  2001/05/24 12:52:02  nando
    smoothing #include xx.hpp

    Revision 1.4  2001/05/23 19:30:27  nando
    smoothing #include xx.hpp

    Revision 1.3  2001/05/22 13:23:04  marmar
    Method controlVariateCorrection added

    Revision 1.2  2001/04/26 16:04:52  marmar
    underlying_ not mutable anymore

    Revision 1.1  2001/04/09 14:05:48  nando
    all the *.hpp moved below the Include/ql level

    Revision 1.3  2001/04/06 18:46:20  nando
    changed Authors, Contributors, Licence and copyright header

*/

#ifndef BSM_numerical_option_pricer_h
#define BSM_numerical_option_pricer_h

#include "ql/Pricers/bsmoption.hpp"
#include "ql/FiniteDifferences/bsmoperator.hpp"

namespace QuantLib {

    namespace Pricers {

        class BSMNumericalOption : public BSMOption {
          public:
                BSMNumericalOption(Type type, double underlying, double strike,
                    Rate dividendYield, Rate riskFreeRate, Time residualTime,
                    double volatility, int gridPoints);
                // accessors
                virtual void calculate() const = 0;
                double value() const;
                double delta() const;
                double gamma() const;
                double theta() const;
                Array getGrid() const{return grid_;}

          protected:
            // methods
            virtual void setGridLimits(double center, double timeDelay) const;
            virtual void initializeGrid() const;
            virtual void initializeInitialCondition() const;
            virtual void initializeOperator() const;
            // input data
            int gridPoints_;
            // results
            mutable double delta_, gamma_, theta_;

            mutable Array grid_;
            mutable FiniteDifferences::BSMOperator finiteDifferenceOperator_;
            mutable Array initialPrices_;
            // temporaries
            mutable double sMin_, center_, sMax_;
          private:
            // temporaries
            mutable double gridLogSpacing_;
            int safeGridPoints(int gridPoints, Time residualTime);
        };

        //! This is a safety check to be sure we have enough grid points.
        #define QL_NUM_OPT_MIN_GRID_POINTS            100
        //! This is a safety check to be sure we have enough grid points.
        #define QL_NUM_OPT_GRID_POINTS_PER_YEAR        50

            // The following is a safety check to be sure we have enough grid
            // points.
        inline int BSMNumericalOption::safeGridPoints(int gridPoints,
                                                        Time residualTime){
            return QL_MAX(gridPoints,
              residualTime>1.0 ? (int)(QL_NUM_OPT_MIN_GRID_POINTS +
              (residualTime-1.0)*QL_NUM_OPT_GRID_POINTS_PER_YEAR) :
              QL_NUM_OPT_MIN_GRID_POINTS);
        }

    }

}


#endif
