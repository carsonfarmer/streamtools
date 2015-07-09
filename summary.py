#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Methods for computing various summary statistics on streams of data."""

# Copyright (C) 2015, Carson Farmer <carsonfarmer@gmail.com>
# Based on code from http://www.johndcook.com/blog/skewness_kurtosis/
# All rights reserved. MIT Licensed.

import math


class RunningRegression(object):
    """Class to compute simple linear regressions on data streams.

    The RunningRegression class is the analog of the RunningStats class. You
    add pairs of (x, y) values by using the Push. At any point along the way
    you can call the slope, intercept, or correlation functions to see the
    current value of these statistics.

    Notes
    -----
    You can combine two RunningRegression objects by using the + and +=
    operators. For example, you might accrue data on several different threads
    in parallel then add their RunningRegression objects together.
    """

    def __init__(self):
        self.x_stats = RunningStats()
        self.y_stats = RunningStats()
        self.S_xy = 0.
        self.n = 0

    def clear(self):
        self.x_stats.clear()
        self.y_stats.clear()
        self.S_xy = 0.
        self.n = 0

    def add(self, x, y):
        self.S_xy += (self.x_stats.mean() - x)*(self.y_stats.mean() - y) * \
            float(self.n)/float(self.n+1)

        self.x_stats.add(x)
        self.y_stats.add(y)
        self.n += 1

    def num_samples(self):
        return self.n

    def slope(self):
        S_xx = self.x_stats.variance()*(self.n - 1.0)
        return self.S_xy / S_xx

    def intercept(self):
        return self.y_stats.mean() - self.slope()*self.x_stats.mean()

    def correlation(self):
        t = self.x_stats.stddev() * self.y_stats.stddev()
        return self.S_xy / ((self.n-1) * t)

    def __add__(self, other):
        a, b = self, other
        c = RunningRegression()

        c.x_stats = a.x_stats + b.x_stats
        c.y_stats = a.y_stats + b.y_stats
        c.n = a.n + b.n

        delta_x = b.x_stats.mean() - a.x_stats.mean()
        delta_y = b.y_stats.mean() - a.y_stats.mean()
        c.S_xy = a.S_xy + b.S_xy + float(a.n*b.n)*delta_x*delta_y/float(c.n)

        return c

    def __iadd__(self, other):
        comb = self + other
        self = comb
        return self


class RunningStats(object):
    """Class to compute summary statistics on data streams.

    This code is an extension of the method of Knuth and Welford for computing
    standard deviation in one pass through the data. It computes skewness and
    kurtosis as well with a similar interface. In addition to only requiring
    one pass through the data, the algorithm is numerically stable and accurate.

    To use the class, create an instance of `RunningStats`, and use the `add`
    method to insert new data values. After accruing your data, you can obtain
    sample statistics by calling one of these methods:

        mean
        variance
        stddev
        skewness
        kurtosis

    Notes
    -----
    You can combine two RunningStats objects by using the + and += operators.
    For example, you might accrue data on several different processes in
    parallel then add their RunningStats objects together to create a single
    object with the state that it would have had if all the data had been
    accumulated by it alone.

    References
    ----------
    [1] Algorithms for calculating variance
        http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    [2] Timothy B. Terriberry. Computing Higher-Order Moments Online.
        https://people.xiph.org/~tterribe/notes/homs.html
    [3] Philippe Pébay. SANDIA REPORT SAND2008-6212 (2008). Formulas for Robust,
        One-Pass Parallel Computation of Co-variances and Arbitrary-Order
        Statistical Moments.
        http://prod.sandia.gov/techlib/access-control.cgi/2008/086212.pdf
    """

    def __init__(self):
        self.clear()

    def clear(self):
        """Reset the collection of sufficient statistics."""
        self.n = 0
        self.M1 = 0.
        self.M2 = 0.
        self.M3 = 0.
        self.M4 = 0.
        self.min_val = float('inf')
        self.max_val = float('-inf')

    def add(self, x):
        """Add an additional data value.

        Parameters
        ----------
        x : numeric
            A new data value from the input data stream.
        """
        n1 = self.n
        self.n += 1
        self.min_val = min(self.min_val, x)
        self.max_val = max(self.max_val, x)
        delta = x - self.M1
        delta_n = delta / self.n
        delta_n2 = delta_n * delta_n
        term1 = delta * delta_n * n1
        self.M1 += delta_n
        self.M4 += term1 * delta_n2 * (self.n*self.n - 3*self.n + 3) + 6 * \
            delta_n2 * self.M2 - 4 * delta_n * self.M3
        self.M3 += term1 * delta_n * (self.n - 2) - 3 * delta_n * self.M2
        self.M2 += term1

    def num_samples(self):
        """The current number of data values."""
        return self.n

    def mean(self):
        """The current mean of the input data values."""
        return self.M1

    def variance(self):
        """The current variance of the input data values."""
        return self.M2/(self.n-1.0)

    def stddev(self):
        """The current standard deviation of the input data values."""
        return math.sqrt(self.variance())

    def skewness(self):
        """The current skewness of the input data values."""
        return math.sqrt(float(self.n))*self.M3/math.pow(self.M2, 1.5)

    def kurtosis(self):
        """The current kurtosis of the input data values."""
        return float(self.n)*self.M4 / (self.M2*self.M2) - 3.0

    def __add__(self, other):
        """Merge this RunningStats object with another.

        Parameters
        ----------
        other : RunningStats
            The RunningStats instance to merge with this one.
        """
        a, b = self, other
        c = RunningStats()

        c.n = a.n + b.n

        delta = b.M1 - a.M1
        delta2 = delta*delta
        delta3 = delta*delta2
        delta4 = delta2*delta2

        c.M1 = (a.n*a.M1 + b.n*b.M1) / c.n
        c.M2 = a.M2 + b.M2 + delta2 * a.n * b.n / c.n
        c.M3 = a.M3 + b.M3 + delta3 * a.n * b.n * (a.n - b.n) / (c.n*c.n)
        c.M3 += 3.0*delta * (a.n*b.M2 - b.n*a.M2) / c.n
        c.M4 = a.M4 + b.M4 + delta4*a.n*b.n * \
            (a.n*a.n - a.n*b.n + b.n*b.n) / (c.n*c.n*c.n)
        c.M4 += 6.0*delta2 * (a.n*a.n*b.M2 + b.n*b.n*a.M2) / \
            (c.n*c.n) + 4.0*delta*(a.n*b.M3 - b.n*a.M3) / c.n

        return c

    def __iadd__(self, other):
        """Merge another RunningStats object into this one.

        Parameters
        ----------
        other : RunningStats
            The RunningStats instance to merge into this one
        """
        comb = self + other
        self = comb
        return self
