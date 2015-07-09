#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Exponentially-weighted moving average.

The class in this module wraps a Exponentially Weighted Moving Average
algorithm (EWMA [1]), based on a talk on "Quantifying Abnormal Behavior" by
VividCortex [2]. They provide an excellent discussion of EWMAs, along with an
implementation in Go, in their github repo [3].

This module contains a generic implementation of the EWMA algorithm. It has a
default warm-up period of 1 and it uses an exponential decay. It supports a
custom age which must be stored, and thus uses more memory. It will report a
value of 0.0 until you have added the required number of samples to pass the
warmup period. It uses some memory to store the number of samples added to it.
As a result it uses a little over twice the memory of a simpler EWMA.

The current implementation assumes an implicit time interval of 1.0 between
every sample added. That is, the passage of time is treated as though it's the
same as the arrival of samples. If you need time-based decay when samples are
not arriving precisely at set intervals, then this package will not support your
needs at present.

Notes
-----
An exponentially weighted moving average is a way to continuously compute a
type of average for a series of numbers, as the numbers arrive. After a value
in the series is added to the average, its weight in the average decreases
exponentially over time. This biases the average towards more recent data.
EWMAs are useful for several reasons, chiefly their inexpensive computational
and memory cost, as well as the fact that they represent the recent central
tendency of the series of values.

The overall algorithm works thus, in pseudocode:
    * Multiply the next number in the series by alpha.
    * Multiply the current value of the average by 1 minus alpha.
    * Add the result of steps 1 and 2, and store it as the new current value of
      the average.
    * Repeat for each number in the series.

There are special-case behaviors for how to initialize the current value, and
these vary between implementations. One approach is to start with the first
value in the series; another is to average the first n or so values in the
series using an arithmetic average, and then begin the incremental updating of
the average. Each method has pros and cons.

Extras
------
The EWMA algorithm requires a decay factor, alpha. The larger the alpha, the
more the average is biased towards recent history. The alpha must be between 0
and 1, and is typically a fairly small number, such as 0.04. Consider a
fixed-size sliding-window moving average (not an exponentially weighted moving
average) that averages over the previous N samples. What is the average age of
each sample? It is N/2. Now suppose that you wish to construct a EWMA whose
samples have the same average age. The formula to compute the alpha required
for this is: alpha = 2/(N+1). Proof is in the book "Production and Operations
Analysis" by Steven Nahmias. So, for example, if you have a time-series with
samples once per second, and you want to get the moving average over the
previous minute, you should use an alpha of 0.032786885 (i.e., 2 / (60 + 1)).

References
----------
[1]: http://en.wikipedia.org/wiki/Moving_average#Exponential_moving_average
[2]: https://vividcortex.com/blog/2013/07/23/a-fast-go-library-for-exponential
     -moving-averages/
[3]: https://github.com/VividCortex/ewma
"""

# Copyright (C) 2013 VividCortex <https://vividcortex.com>
# Copyright (C) 2014, Carson Farmer <carsonfarmer@gmail.com>
# All rights reserved. MIT Licensed.


class EWMA(object):
    """Class to calculate a moving average.

    Computes a moving average over a time-series stream of numbers. The average
    may be over a window or exponentially decaying. If no age is given, this
    function will return a default exponentially weighted implementation
    that consumes minimal memory. The age is related to the decay factor
    alpha by 2.0 / (`age` + 1). It signifies the average age of the samples as
    time goes to infinity.

    Parameters
    ----------
    age : numeric (default=30)
        By default, we average over a one-minute period, which means the
        average age of the metrics in the period is 30 seconds.
    warmup : numeric (default=1)
        For best results, the moving average should not be initialized to
        the samples it sees immediately. The book "Production and Operations
        Analysis" by Steven Nahmias suggests initializing the moving average
        to the mean of the first 10 samples. Until the Ewma has seen
        this many samples, it is not "ready" to be queried for the value of
        the moving average.
    """

    def __init__(self, age=30.0, warmup=1):
        super(EWMA, self).__init__()
        self._decay = 2.0 / (age + 1)
        self._value = 0.0
        # self._count = 0
        # self._warmup = warmup if warmup > 0 else 1

    def add(self, value):
        """Adds a value to the series and updates the moving average."""
        # if self._count < self._warmup:
        #     self._count += 1
        #     self._value += value
        # elif self._count == self._warmup:
        #     self._value = self._value / float(self._warmup)
        #     self._count += 1
        # else:
        self._value = (value*self._decay) + (self._value*(1 - self._decay))

    # @property
    # def count(self):
    #     """The number of samples added to this instance."""
    #     return self._count

    @property
    def value(self):
        """Returns the current value of the moving average."""
        # if self._count <= self._warmup:
        #     return 0.0
        return self._value

    @value.setter
    def value(self, value):
        """Sets the EWMA's value."""
        self._value = value
        # if self._count <= self._warmup:
        #     self._count = self._warmup + 1  # Force warmup period
