
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

/*! \file calendar.cpp
    \brief Abstract calendar class

    $Source$
    $Log$
    Revision 1.20  2001/05/24 13:57:51  nando
    smoothing #include xx.hpp and cutting old Log messages

    Revision 1.19  2001/05/08 17:23:47  lballabio
    removed unnecessary if branch (although more convoluted, it had the same effect of the else branch

    Revision 1.18  2001/04/09 14:13:33  nando
    all the *.hpp moved below the Include/ql level

    Revision 1.17  2001/04/06 18:46:21  nando
    changed Authors, Contributors, Licence and copyright header

    Revision 1.16  2001/04/04 11:07:24  nando
    Headers policy part 1:
    Headers should have a .hpp (lowercase) filename extension
    All *.h renamed to *.hpp

    Revision 1.15  2001/03/02 15:43:36  lballabio
    Fixed a bug in advance() with a negative number of days

    Revision 1.14  2001/03/01 11:37:07  lballabio
    Fixed bug in advance(...,Days)

    Revision 1.13  2001/01/24 13:17:46  marmar
    style redefined

    Revision 1.12  2001/01/17 14:37:56  nando
    tabs removed

    Revision 1.11  2000/12/20 18:26:15  enri
    test

    Revision 1.10  2000/12/20 16:42:38  enri
    commit test

    Revision 1.9  2000/12/14 12:32:31  lballabio
    Added CVS tags in Doxygen file documentation blocks

*/

#include "ql/calendar.hpp"

namespace QuantLib {

    Date Calendar::roll(const Date& d , bool modified) const {
        Date d1 = d;
        while (isHoliday(d1))
            d1++;
        if (modified && d1.month() != d.month()) {
            d1 = d;
            while (isHoliday(d1))
                d1--;
        }
        return d1;
    }

    Date Calendar::advance(const Date& d, int n, TimeUnit unit,
      bool modified) const {
        if (n == 0) {
            return roll(d,modified);
        } else if (unit == Days) {
            Date d1 = d;
            if (n > 0) {
                while (n > 0) {
                    d1++;
                    while (isHoliday(d1))
                        d1++;
                    n--;
                }
            } else {
                while (n < 0) {
                    d1--;
                    while(isHoliday(d1))
                        d1--;
                    n++;
                }
            }
            return d1;
        } else {
            Date d1 = d.plus(n,unit);
            return roll(d1,modified);
        }
        QL_DUMMY_RETURN(Date());
    }

}
