/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006 Chiara Fornarola

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

/*! \file euriborswapfixtgm.hpp
    \brief %euriborswapfixtgm index
*/

#ifndef quantlib_euriborswapfixifr_hpp
#define quantlib_euriborswapfixifr_hpp

#include <ql/Indexes/swapindex.hpp>
#include <ql/Indexes/euribor.hpp>
#include <ql/Calendars/target.hpp>
#include <ql/DayCounters/thirty360.hpp>
#include <ql/Currencies/europe.hpp>

namespace QuantLib {

  //EuriborSwapFix index published by IFR Markets and distributed by Reuters page TGM42281 and
  //       by Telerate. For more info see http://www.ifrmarkets.com
  //  

    class EuriborSwapFixIFR : public SwapIndex {
      public:
        EuriborSwapFixIFR(Integer years,
                        const Handle<YieldTermStructure>& h =
                                    Handle<YieldTermStructure>())
        : SwapIndex("EuriborSwapFixIFR", // familyName
                    years,
                    2, // settlementDays
                    EURCurrency(),
                    TARGET(), 
                    Annual, // fixedLegFrequency
                    Unadjusted, // fixedLegConvention
                    Thirty360(Thirty360::BondBasis), // fixedLegDaycounter 
                    boost::shared_ptr<Xibor>(new Euribor6M(h))) {}
    };



    //! 1-year %EuriborSwapFixIFR index
    class EuriborSwapFixIFR1Y : public EuriborSwapFixIFR {
      public:
        EuriborSwapFixIFR1Y(const Handle<YieldTermStructure>& h)
        : EuriborSwapFixIFR(1,h) {}
    };

    //! 2-year %EuriborSwapFixIFR index
    class EuriborSwapFixIFR2Y : public EuriborSwapFixIFR {
      public:
        EuriborSwapFixIFR2Y(const Handle<YieldTermStructure>& h)
        : EuriborSwapFixIFR(2,h) {}
    };

    //! 3-year %EuriborSwapFixIFR index
    class EuriborSwapFixIFR3Y : public EuriborSwapFixIFR {
      public:
        EuriborSwapFixIFR3Y(const Handle<YieldTermStructure>& h)
        : EuriborSwapFixIFR(3,h) {}
    };

    //! 4-year %EuriborSwapFixIFR index
    class EuriborSwapFixIFR4Y : public EuriborSwapFixIFR {
      public:
        EuriborSwapFixIFR4Y(const Handle<YieldTermStructure>& h)
        : EuriborSwapFixIFR(4,h) {}
    };

    //! 5-year %EuriborSwapFixIFR index
    class EuriborSwapFixIFR5Y : public EuriborSwapFixIFR {
      public:
        EuriborSwapFixIFR5Y(const Handle<YieldTermStructure>& h)
        : EuriborSwapFixIFR(5,h) {}
    };

    //! 6-year %EuriborSwapFixIFR index
    class EuriborSwapFixIFR6Y : public EuriborSwapFixIFR {
      public:
        EuriborSwapFixIFR6Y(const Handle<YieldTermStructure>& h)
        : EuriborSwapFixIFR(6,h) {}
    };
    
    //! 7-year %EuriborSwapFixIFR index
    class EuriborSwapFixIFR7Y : public EuriborSwapFixIFR {
      public:
        EuriborSwapFixIFR7Y(const Handle<YieldTermStructure>& h)
        : EuriborSwapFixIFR(7,h) {}
    };

    //! 8-year %EuriborSwapFixIFR index
    class EuriborSwapFixIFR8Y : public EuriborSwapFixIFR {
      public:
        EuriborSwapFixIFR8Y(const Handle<YieldTermStructure>& h)
        : EuriborSwapFixIFR(8,h) {}
    };
    
    //! 9-year %EuriborSwapFixIFR index
    class EuriborSwapFixIFR9Y : public EuriborSwapFixIFR {
      public:
        EuriborSwapFixIFR9Y(const Handle<YieldTermStructure>& h)
        : EuriborSwapFixIFR(9,h) {}
    };

    //! 10-year %EuriborSwapFixIFR index
    class EuriborSwapFixIFR10Y : public EuriborSwapFixIFR {
      public:
        EuriborSwapFixIFR10Y(const Handle<YieldTermStructure>& h)
        : EuriborSwapFixIFR(10,h) {}
    };

    //! 12-year %EuriborSwapFixIFR index
    class EuriborSwapFixIFR12Y : public EuriborSwapFixIFR {
      public:
        EuriborSwapFixIFR12Y(const Handle<YieldTermStructure>& h)
        : EuriborSwapFixIFR(12,h) {}
    };

    //! 15-year %EuriborSwapFixIFR index
    class EuriborSwapFixIFR15Y : public EuriborSwapFixIFR {
      public:
        EuriborSwapFixIFR15Y(const Handle<YieldTermStructure>& h)
        : EuriborSwapFixIFR(15,h) {}
    };

    //! 20-year %EuriborSwapFixIFR index
    class EuriborSwapFixIFR20Y : public EuriborSwapFixIFR {
      public:
        EuriborSwapFixIFR20Y(const Handle<YieldTermStructure>& h)
        : EuriborSwapFixIFR(20,h) {}
    };

    //! 25-year %EuriborSwapFixIFR index
    class EuriborSwapFixIFR25Y : public EuriborSwapFixIFR {
      public:
        EuriborSwapFixIFR25Y(const Handle<YieldTermStructure>& h)
        : EuriborSwapFixIFR(25,h) {}
    };

    //! 30-year %EuriborSwapFixIFR index
    class EuriborSwapFixIFR30Y : public EuriborSwapFixIFR {
      public:
        EuriborSwapFixIFR30Y(const Handle<YieldTermStructure>& h)
        : EuriborSwapFixIFR(30,h) {}
    };
   
}


#endif