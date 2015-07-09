#!/usr/bin/env python

from math import copysign  # , fabs, sqrt
import random
# import sys

# http://www.seancassidy.me/on-accepting-interview-question-answers.html


def calcP2(qp1, q, qm1, d, np1, n, nm1):
    d = float(d)
    n = float(n)
    np1 = float(np1)
    nm1 = float(nm1)

    outer = d / (np1 - nm1)
    inner_left = (n - nm1 + d) * (qp1 - q) / (np1 - n)
    inner_right = (np1 - n - d) * (q - qm1) / (n - nm1)

    return q + outer * (inner_left + inner_right)


class Quantile:
    LEN = 5

    def __init__(self, p):
        """ Constructs a single quantile object """
        self.dn = [0, p/2, p, (1 + p)/2, 1]
        self.npos = [1, 1 + 2*p, 1 + 4*p, 3 + 2*p, 5]
        self.pos = list(range(1, self.LEN + 1))
        self.heights = []
        self.initialized = False
        self.p = p

    def add(self, item):
        """ Adds another datum """
        if len(self.heights) != self.LEN:
            self.heights.append(item)
        else:
            if self.initialized is False:
                self.heights.sort()
                self.initialized = True

            # find cell k
            if item < self.heights[0]:
                self.heights[0] = item
                k = 1
            else:
                for i in range(1, self.LEN):
                    if self.heights[i - 1] <= item and item < self.heights[i]:
                        k = i
                        break
                else:
                    k = 4
                    if self.heights[-1] < item:
                        self.heights[-1] = item

            # increment all positions greater than k
            self.pos = [j if i < k else j + 1 for i, j in enumerate(self.pos)]
            self.npos = [x + y for x, y in zip(self.npos, self.dn)]

            self.__adjust()

    def __adjust(self):
        for i in range(1, self.LEN - 1):
            n = self.pos[i]
            q = self.heights[i]

            d = self.npos[i] - n

            if (d >= 1 and self.pos[i + 1] - n > 1) or \
               (d <= -1 and self.pos[i - 1] - n < -1):
                d = int(copysign(1, d))

                qp1 = self.heights[i + 1]
                qm1 = self.heights[i - 1]
                np1 = self.pos[i + 1]
                nm1 = self.pos[i - 1]
                qn = calcP2(qp1, q, qm1, d, np1, n, nm1)

                if qm1 < qn and qn < qp1:
                    self.heights[i] = qn
                else:
                    # use linear form
                    self.heights[i] = q + d * (self.heights[i + d] - q) / \
                        (self.pos[i + d] - n)

                self.pos[i] = n + d

    def quantile(self):
        if self.initialized:
            return self.heights[2]
        else:
            self.heights.sort()
            l = len(self.heights)
            # make sure we don't overflow on p == 1 or underflow on p == 0
            return self.heights[int(min(max(l - 1, 0), l * self.p))]


class Frugal(object):

    def __init__(self, guess=0, q=0.5, f=lambda step: 1, kind="two"):
        """
        Parameters
        ----------
        func : function
            Suggestions:
            f(step) = 1
            f(step) = 10
            f(step) = step
            f(step) = 2*step [default]
        References
        ----------
        [1] research.neustar.biz/2013/09/16/sketch-of-the-day-frugal-streaming/

        """
        self.m = guess
        self.q = q
        self.f = f
        self.step = 1
        self.sign = 1
        self.kind = kind

    def quantile(self):
        return self.m

    def add(self, item):
        r = random.random()
        if self.kind.startswith("one"):
            if item > self.m and r > 1 - 0.75:
                self.m += 1
            elif item < self.m and r > 0.75:
                self.m -= 1
        else:  # Use Frugal-2U by default
            if item > self.m and r > 1 - self.q:
                # Increment the step size if and only if the estimate keeps
                # moving in the same direction. Step size is incremented by the
                # result of applying the specified step function to the previous
                # step size.
                if self.sign > 0:
                    self.step += self.f(self.step)
                else:
                    self.step -= self.f(self.step)
                    # Increment the estimate by step size if step is positive.
                    # Otherwise, increment the step size by one.
                self.m += self.step if self.step > 0 else 1

                # If the estimate overshot the item in the stream, pull the
                # estimate back and re-adjust the step size.
                if self.m > item:
                    self.step += (item - self.m)
                    self.m = item

                # Mark that the estimate increased this step
                self.sign = 1
            # If the item is less than the stream, follow all of the same steps
            # as above, with signs reversed.
            elif item < self.m and r > self.q:
                if self.sign < 0:
                    self.step += self.f(self.step)
                else:
                    self.step -= self.f(self.step)

                self.m -= self.step if self.step > 0 else 1

                if self.m < item:
                    self.step += (self.m - item)
                    self.m = item

                self.sign = -1

            # Damp down the step size to avoid oscillation.
            if (self.m - item) * self.sign < 0 and self.step > 1:
                self.step = 1


class LiveStats:
    def __init__(self, p=[0.5], kind="frugal"):
        """ Constructs a LiveStream object

        Keyword arguments:

        p -- A list of quantiles to track, by default, [0.5]

        Some ideas for allowing the removal of data values:
        alias-i.com/lingpipe/src/com/aliasi/stats/OnlineNormalEstimator.java

        """
        self.min_val = float('inf')
        self.max_val = float('-inf')
        self.var_m2 = 0.0
        self.kurt_m4 = 0.0
        self.skew_m3 = 0.0
        self.average = 0.0
        self.count = 0
        self.tiles = {}
        self.initialized = False
        for i in p:
            if kind in "p-square":
                self.tiles[i] = Quantile(i)
            else:  # Frugal-2U (or 1U)
                self.tiles[i] = Frugal(i)

    def add(self, item):
        """ Adds another datum """
        delta = item - self.average

        self.min_val = min(self.min_val, item)
        self.max_val = max(self.max_val, item)

        # Average
        self.average = (self.count * self.average + item) / (self.count + 1)
        self.count = self.count + 1

        # Variance (except for the scale)
        self.var_m2 = self.var_m2 + delta * (item - self.average)

        # tiles
        for perc in list(self.tiles.values()):
            perc.add(item)

        # Kurtosis
        self.kurt_m4 = self.kurt_m4 + (item - self.average)**4.0

        # Skewness
        self.skew_m3 = self.skew_m3 + (item - self.average)**3.0

    def quantiles(self):
        """ Returns a list of tuples of the quantile and its location """
        return [(key, val.quantile()) for key, val in self.tiles.items()]

    def maximum(self):
        """ Returns the maximum value given """
        return self.max_val

    def mean(self):
        """ Returns the cumulative moving average of the data """
        return self.average

    def minimum(self):
        """ Returns the minimum value given """
        return self.min_val

    def num(self):
        """ Returns the number of items added so far"""
        return self.count

    def variance(self):
        """ Returns the sample variance of the data given so far"""
        if self.count > 1:
            return self.var_m2 / (self.count - 1)
        else:
            return float('NaN')

    def kurtosis(self):
        """ Returns the sample kurtosis of the data given so far"""
        if self.count > 1:
            return self.kurt_m4 / (self.count * self.variance()**2.0) - 3.0
        else:
            return float('NaN')

    def skewness(self):
        """ Returns the sample skewness of the data given so far"""
        if self.count > 1:
            return self.skew_m3 / (self.count * self.variance()**1.5)
        else:
            return float('NaN')

    def __add__(self, other):
        if not isinstance(type(self), other):
            raise(TypeError("cannot merge '%s' and '%s' obejcts" % (type(self),
                  type(other))))
        if not self.tiles.keys() == other.tiles.keys():
            raise(ValueError("cannot merge '%s' objects with differing "
                  "quantiles" % type(self)))

        a, b = self, other
        c = LiveStats()

        c.count = a.count + b.count
        c.average = (a.count*a.average + c.count*c.average) / c.count

        delta = b.average - a.average
        c.var_m2 = a.var_m2 + b.var_m2 + delta**2.0 * a.count * b.count / \
            c.count

        # Skewness
        c.skew_m3 = a.skew_m3 + b.skew_m3 + delta**3.0 * a.count * b.count * \
            (a.count - b.count) / c.count**2.0
        c.skew_m3 += 3.0*delta * (a.count*b.var_m2 - b.count*a.var_m2) / c.count

        # Kurtosis
        c.kurt_m4 = a.kurt_m4 + b.kurt_m4 + delta**4.0 * a.count * b.count * \
            (a.count*a.count - a.count*b.count + b.count*b.count) / \
            (c.count**3.0)
        c.kurt_m4 += 6.0*delta**2.0 * (a.count**2.0*b.var_m2 +
                                       b.count**2.0*a.var_m2) / \
            (c.count*c.count) + 4.0*delta*(a.count*b.skew_m3 -
                                           b.count*a.skew_m3) / c.count
