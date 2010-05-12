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

/*! \file timeseries.hpp
    \brief Container for historical data
*/

#ifndef quantlib_timeseries_hpp
#define quantlib_timeseries_hpp

#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <ql/time/date.hpp>
#include <ql/utilities/null.hpp>
#include <ql/errors.hpp>
#include <map>
#include <vector>

namespace QuantLib {

    //! Container for historical data
    /*! This class acts as a generic repository for a set of
        historical data.  Any single datum can be accessed through its
        date, while sets of consecutive data can be accessed through
        iterators.

        \pre The <c>Container</c> type must satisfy the requirements
             set by the C++ standard for associative containers.
    */
    template <class T, class Time = Date, class Container = std::map<Time, T> >
    class TimeSeries {
      public:
        typedef Time key_type;
        typedef T value_type;
      private:
        mutable Container values_;
      public:
        /*! Default constructor */
        TimeSeries() {}
        /*! This constructor initializes the history with a set of
            values passed as two sequences, the first containing dates
            and the second containing corresponding values.
        */
        template <class TimeIterator, class ValueIterator>
        TimeSeries(TimeIterator tBegin, TimeIterator tEnd,
                   ValueIterator vBegin) {
            while (tBegin != tEnd)
                values_[*(tBegin++)] = *(vBegin++);
        }
        /*! This constructor initializes the history with a set of
            values. Such values are assigned to a corresponding number
            of consecutive dates starting from <b><i>firstDate</i></b>
            included.
        */
        template <class ValueIterator>
        TimeSeries(const Time& tBegin,
                   ValueIterator begin, ValueIterator end) {
            Time t = tBegin;
            while (begin != end)
                values_[t++] = *(begin++);
        }
        //! \name Inspectors
        //@{
        //! returns the first date for which a historical datum exists
        Time firstDate() const;
        Time earliest() const { return firstDate(); }
        //! returns the last date for which a historical datum exists
        Time lastDate() const;
        Time latest() const { return lastDate(); }
        //! returns the number of historical data including null ones
        Size size() const;
        //! returns whether the series contains any data
        bool empty() const;
        //@}
        //! \name Historical data access
        //@{
        //! returns the (possibly null) datum corresponding to the given date
        T operator[](const Time& t) const {
            if (values_.find(t) != values_.end())
                return values_[t];
            else
                return Null<T>();
        }
        T& operator[](const Time& t) {
            if (values_.find(t) == values_.end())
                values_[t] = Null<T>();
            return values_[t];
        }
        //@}

        //! \name Iterators
        //@{
        typedef typename Container::const_iterator const_iterator;
        typedef typename const_iterator::iterator_category iterator_category;
        typedef typename Container::const_reverse_iterator
                                           const_reverse_iterator;
        const_iterator begin() const;
        const_iterator end() const;
        const_reverse_iterator rbegin() const;
        const_reverse_iterator rend() const;
        //@}

        //! \name Projection iterators
        //@{
        template<class U>
        class projection_iterator_base : public boost::iterator_facade<
			projection_iterator_base<U>, U, boost::bidirectional_traversal_tag> {
        public:
        	typedef U value_type;
        	explicit projection_iterator_base(const_iterator it) : it_(it) {
        	}
        private:
        	friend class boost::iterator_core_access;

        	const_iterator it_;

        	void increment() { ++it_; }

        	void decrement() { --it_; }

        	bool equal(projection_iterator_base<U> const& other) const {
        	        return this->it_ == other.it_;
        	}

        	U &dereference() const { return dereference(U()); }
        	const T &dereference (const T &) const { return this->it_->second; }
        	const Time &dereference (const Time &) const { return this->it_->first; }
        };

        template<class Base>
        class projection_iterator : public boost::iterator_adaptor<
			projection_iterator<Base>, Base, typename Base::value_type,
			boost::bidirectional_traversal_tag> {

        	typedef boost::iterator_adaptor<
				projection_iterator<Base>, Base, typename Base::value_type,
				boost::bidirectional_traversal_tag> iterator_adaptor_;
        public:
        	explicit projection_iterator(const_iterator it) :
        		projection_iterator::iterator_adaptor_(Base(it)) {
        	}
       	};

        typedef projection_iterator<projection_iterator_base<const T> > value_iterator;
        typedef projection_iterator<projection_iterator_base<const Time> > time_iterator;

        value_iterator begin_values() const { return value_iterator(values_.begin()); }
        value_iterator end_values() const { return value_iterator(values_.end()); }

        time_iterator begin_time() const { return time_iterator(values_.begin()); }
        time_iterator end_time() const { return time_iterator(values_.end()); }
        //@}

        //! \name Utilities
        //@{
        const_iterator find(const Time&);
        const_iterator find(const Time&t) const { return values_.find(t); }
        void clear() { values_.clear(); }
        //! returns the dates for which historical data exist
        std::vector<Time> dates() const;
        //! returns the historical data
        std::vector<T> values() const;
        //@}
    };

    // inline definitions

    template <class T, class Time, class C>
    inline Time TimeSeries<T,Time,C>::firstDate() const {
        QL_REQUIRE(!values_.empty(), "empty timeseries");
        return values_.begin()->first;
    }

    template <class T, class Time, class C>
    inline Time TimeSeries<T,Time,C>::lastDate() const {
        QL_REQUIRE(!values_.empty(), "empty timeseries");
        return values_.rbegin()->first;
    }

    template <class T, class Time, class C>
    inline Size TimeSeries<T,Time,C>::size() const {
        return values_.size();
    }

    template <class T, class Time, class C>
    inline bool TimeSeries<T,Time,C>::empty() const {
        return values_.empty();
    }

    template <class T, class Time, class C>
    inline typename TimeSeries<T,Time,C>::const_iterator
    TimeSeries<T,Time,C>::begin() const {
        return values_.begin();
    }

    template <class T, class Time, class C>
    inline typename TimeSeries<T,Time,C>::const_iterator
    TimeSeries<T,Time,C>::end() const {
        return values_.end();
    }

    template <class T, class Time, class C>
    inline typename TimeSeries<T,Time,C>::const_reverse_iterator
    TimeSeries<T,Time,C>::rbegin() const {
        return values_.rbegin();
    }

    template <class T, class Time, class C>
    inline typename TimeSeries<T,Time,C>::const_reverse_iterator
    TimeSeries<T,Time,C>::rend() const {
        return values_.rend();
    }

    template <class T, class Time, class C>
    inline typename TimeSeries<T,Time,C>::const_iterator
    TimeSeries<T,Time,C>::find(const Time& t) {
        const_iterator i = values_.find(t);
        if (i == values_.end()) {
            values_[t] = Null<T>();
            i = values_.find(t);
        }
        return i;
    }

    template <class T, class Time, class C>
    std::vector<Time> TimeSeries<T,Time,C>::dates() const {
        std::vector<Time> v;
        v.reserve(size());
        for (const_iterator i = begin(); i != end(); ++i)
            v.push_back(i->first);
        return v;
    }

    template <class T, class Time, class C>
    std::vector<T> TimeSeries<T,Time,C>::values() const {
        std::vector<T> v;
        v.reserve(size());
        for (const_iterator i = begin(); i != end(); ++i)
            v.push_back(i->second);
        return v;
    }
}

#endif
